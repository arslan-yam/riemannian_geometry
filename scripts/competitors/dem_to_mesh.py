"""A real DEM (GeoTIFF, e.g. swissALTI3D 0.5 m or USGS 3DEP 1 m) to a terrain mesh (OFF), as a regular grid TIN
or as an irregular TIN with the same relief.

  python3 scripts/competitors/dem_to_mesh.py dem.tif out.off [--window X0 Y0 W H] [--step S]
                                             [--irregular N] [--alpha 1.0] [--seed 1] [--json out.json]

  --window    a pixel window of the raster (default: all of it); several tiles can be mosaicked beforehand
  --step      keep every S-th pixel in both directions (the grid spacing becomes S pixels)
  --irregular keep N points instead of the whole grid, with probability growing with the local curvature
              (|discrete Laplacian| ** alpha, plus a floor so that flat areas keep some points), all boundary points of
              the window, and triangulate them by Delaunay in the chart: an irregular TIN of varying density
Heights and chart coordinates are in metres (the pixel size from ModelPixelScaleTag). No-data pixels (GDAL_NODATA, or
values below -1e4) are not allowed inside the window. Run scripts/competitors/prepare_mesh.py on the result, as for
every other terrain. The JSON records the source, window, spacing and the slope / aspect statistics.
"""
import json, os, sys
import numpy as np
import tifffile
from scipy.spatial import Delaunay


def read_dem(path):
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        z = page.asarray().astype(np.float64)
        tags = {t.name: t.value for t in page.tags.values()}
    sx, sy = (tags.get("ModelPixelScaleTag") or (1.0, 1.0, 0.0))[:2]
    nod = tags.get("GDAL_NODATA")
    nodata = float(str(nod).strip("\x00 ")) if nod is not None else None
    return z, float(sx), float(sy), nodata


def grid_faces(h, w):
    idx = np.arange(h * w).reshape(h, w)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel()
    return np.concatenate([np.stack([a, c, b], 1), np.stack([b, c, d], 1)])


def stats(V, F):
    P = V[F]
    e = np.stack([np.linalg.norm(P[:, 1] - P[:, 0], axis=1), np.linalg.norm(P[:, 2] - P[:, 1], axis=1),
                  np.linalg.norm(P[:, 0] - P[:, 2], axis=1)], 1)
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    area2 = np.linalg.norm(n, axis=1)
    slope = np.degrees(np.arccos(np.clip(np.abs(n[:, 2]) / np.maximum(area2, 1e-300), 0, 1)))
    aspect = e.max(1) ** 2 / np.maximum(area2, 1e-300)  # longest edge over height
    return dict(vertices=len(V), faces=len(F), edge_min=float(e.min()), edge_mean=float(e.mean()), edge_max=float(e.max()),
                slope_p50_deg=float(np.percentile(slope, 50)), slope_p99_deg=float(np.percentile(slope, 99)),
                aspect_p50=float(np.percentile(aspect, 50)), aspect_p99=float(np.percentile(aspect, 99)), aspect_max=float(aspect.max()))


def write_off(path, V, F):
    with open(path, "w") as f:
        f.write(f"OFF\n{len(V)} {len(F)} 0\n")
        np.savetxt(f, V, fmt="%.4f")
        np.savetxt(f, np.c_[np.full(len(F), 3), F], fmt="%d")


def main():
    a = sys.argv[1:]
    src, out = a[0], a[1]
    opt = dict(window=None, step=1, irregular=0, alpha=1.0, seed=1, json=None)
    i = 2
    while i < len(a):
        k = a[i]
        if k == "--window": opt["window"] = [int(x) for x in a[i + 1:i + 5]]; i += 5; continue
        if k == "--step": opt["step"] = int(a[i + 1])
        elif k == "--irregular": opt["irregular"] = int(a[i + 1])
        elif k == "--alpha": opt["alpha"] = float(a[i + 1])
        elif k == "--seed": opt["seed"] = int(a[i + 1])
        elif k == "--json": opt["json"] = a[i + 1]
        else: sys.exit(f"unknown option {k}")
        i += 2
    z, sx, sy, nodata = read_dem(src)
    if opt["window"]:
        x0, y0, w, h = opt["window"]; z = z[y0:y0 + h, x0:x0 + w]
    s = opt["step"]; z = z[::s, ::s]
    bad = (z < -1e4) | (np.isclose(z, nodata) if nodata is not None else False)
    if bad.any():
        sys.exit(f"{int(bad.sum())} no-data pixels in the window: choose another window")
    h, w = z.shape
    yy, xx = np.mgrid[0:h, 0:w]
    X, Y = xx * sx * s, (h - 1 - yy) * sy * s  # row 0 is north
    if opt["irregular"]:
        lap = np.zeros_like(z)
        lap[1:-1, 1:-1] = np.abs(z[:-2, 1:-1] + z[2:, 1:-1] + z[1:-1, :-2] + z[1:-1, 2:] - 4 * z[1:-1, 1:-1])
        wgt = lap ** opt["alpha"]
        wgt = wgt + 0.05 * wgt.mean() + 1e-12  # the floor: flat areas keep a sparse sampling
        border = np.zeros_like(z, bool); border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
        rng = np.random.default_rng(opt["seed"])
        free = np.flatnonzero(~border.ravel())
        n_free = max(0, opt["irregular"] - int(border.sum()))
        p = wgt.ravel()[free]; p = p / p.sum()
        pick = rng.choice(free, size=min(n_free, len(free)), replace=False, p=p)
        keep = np.sort(np.r_[np.flatnonzero(border.ravel()), pick])
        V = np.c_[X.ravel()[keep], Y.ravel()[keep], z.ravel()[keep]]
        tri = Delaunay(V[:, :2])
        F = tri.simplices.astype(np.int64)
        # counter-clockwise in the chart
        P = V[F]; cr = (P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1])
        F[cr < 0] = F[cr < 0][:, [0, 2, 1]]
        F = F[np.abs(cr) > 1e-12 * (sx * s) ** 2]  # drop degenerate hull slivers
    else:
        V = np.c_[X.ravel(), Y.ravel(), z.ravel()]
        F = grid_faces(h, w)
    write_off(out, V, F)
    st = stats(V, F)
    info = dict(source=os.path.basename(src), pixel_m=[sx, sy], window=opt["window"], step=s,
                kind="irregular" if opt["irregular"] else "grid", alpha=opt["alpha"] if opt["irregular"] else None, **st)
    print(json.dumps(info))
    if opt["json"]:
        json.dump(info, open(opt["json"], "w"), indent=1)


if __name__ == "__main__":
    main()
