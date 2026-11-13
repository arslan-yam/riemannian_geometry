"""Adversarial variants of a prepared terrain (review item 5), each aimed at one assumption of the analysis:
  steep5, steep10  heights times 5 / 10: steep slopes, a large metric speed sigma_C along the cuts (longer capped
                   tails) and larger sub-meshes for the exact sweeps
  aniso            the chart stretched 10:1 along x (coordinate rescaling): sliver triangles in the chart, boxes of
                   high aspect ratio (a large window ratio kappa_C)
  slivers          the surface resampled on a grid whose columns alternate between 0.05 and 1.95 times the mean
                   spacing: half the triangles are needles with aspect ratio ~40
  cluster          the surface resampled at random points, 90% of them in 20 Gaussian clusters (density contrast
                   ~100:1), plus pairs of near-coincident points (1e-4 of the mean spacing apart); Delaunay in the
                   chart: dense near tables in the clusters, tiny edges next to long ones
  ties             heights rounded to multiples of 5 m on the original grid: plateaus and terraces with many exactly
                   equal distances (the tie cases that broke VTP's queue), and near ties for kNN
Obstacles with many cut runs: scripts/competitors/make_obstacles.py --kind comb.
The output keeps about the input's vertex count (resampled kinds: --n vertices, default the input's); run
scripts/competitors/prepare_mesh.py on it afterwards like on every terrain.
usage: python3 scripts/competitors/make_adversarial.py in.off out.off --kind steep5 [--n 250000] [--seed 1]
writes out.off and out.json (what was done, mesh statistics)
"""
import argparse, json, os, sys
import numpy as np
from matplotlib.tri import Triangulation, LinearTriInterpolator
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rgeo.terrain import read_off, write_off

ap = argparse.ArgumentParser()
ap.add_argument("inp"); ap.add_argument("out")
ap.add_argument("--kind", choices=["steep5", "steep10", "aniso", "slivers", "cluster", "ties"], required=True)
ap.add_argument("--n", type=int, default=0)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

V, F = read_off(a.inp)
V = np.asarray(V, dtype=float); F = np.asarray(F, dtype=np.int64)
n = a.n or len(V)
rng = np.random.default_rng(a.seed)
lo, hi = V[:, :2].min(0), V[:, :2].max(0)
info = {"input": a.inp, "kind": a.kind, "seed": a.seed}


def surface():
    tri = Triangulation(V[:, 0], V[:, 1], F)
    interp = LinearTriInterpolator(tri, V[:, 2])
    def z(x, y):
        out = np.asarray(interp(x, y).filled(np.nan))
        m = np.isnan(out)
        if m.any():  # on the hull boundary the interpolator can miss: nearest source vertex
            from scipy.spatial import cKDTree
            out[m] = V[cKDTree(V[:, :2]).query(np.column_stack([x[m], y[m]]))[1], 2]
        return out
    return z


def grid_faces(W, H):
    r, c = np.meshgrid(np.arange(H - 1), np.arange(W - 1), indexing="ij")
    v = (r * W + c).ravel()
    return np.concatenate([np.column_stack([v, v + 1, v + W]), np.column_stack([v + 1, v + 1 + W, v + W])])


if a.kind in ("steep5", "steep10"):
    f = 5.0 if a.kind == "steep5" else 10.0
    V2, F2 = V.copy(), F
    V2[:, 2] = (V[:, 2] - V[:, 2].min()) * f
    info["height_factor"] = f
elif a.kind == "aniso":
    V2, F2 = V.copy(), F
    V2[:, 0] = (V[:, 0] - lo[0]) * 10.0
    info["x_factor"] = 10.0
elif a.kind == "ties":
    V2, F2 = V.copy(), F
    V2[:, 2] = np.round(V[:, 2] / 5.0) * 5.0
    info["height_step"] = 5.0
    info["distinct_heights"] = int(len(np.unique(V2[:, 2])))
elif a.kind == "slivers":
    z = surface()
    ar = (hi[0] - lo[0]) / (hi[1] - lo[1])
    H = int(round(np.sqrt(n / ar))); W = int(round(n / H))
    W += W % 2 == 0  # odd: pairs of columns plus the last
    s = (hi[0] - lo[0]) / (W - 1)  # mean spacing; columns alternate 0.05 s and 1.95 s apart
    gaps = np.where(np.arange(W - 1) % 2 == 0, 0.05 * s, 1.95 * s)
    xs = lo[0] + np.concatenate([[0], np.cumsum(gaps)]); xs[-1] = min(xs[-1], hi[0])
    ys = np.linspace(lo[1], hi[1], H)
    X, Y = np.meshgrid(xs, ys)
    V2 = np.column_stack([X.ravel(), Y.ravel(), z(X.ravel(), Y.ravel())])
    F2 = grid_faces(W, H)
    info.update(columns=W, rows=H, spacing=float(s))
else:  # cluster
    from scipy.spatial import Delaunay
    z = surface()
    side = hi - lo
    nc = 20
    centres = lo + rng.random((nc, 2)) * side
    sig = 0.02 * side.max() * rng.uniform(0.5, 1.5, nc)
    n_cl, n_bg = int(0.9 * n), n - int(0.9 * n)
    k = rng.integers(0, nc, n_cl)
    P = centres[k] + rng.normal(size=(n_cl, 2)) * sig[k, None]
    P = P[(P >= lo).all(1) & (P <= hi).all(1)]
    B = lo + rng.random((n_bg, 2)) * side
    corners = np.array([[lo[0], lo[1]], [hi[0], lo[1]], [lo[0], hi[1]], [hi[0], hi[1]]])
    # near-coincident pairs: 200 points copied 1e-4 mean spacings away
    spacing = np.sqrt(side.prod() / n)
    idx = rng.choice(len(P), 200, replace=False)
    ang = rng.random(200) * 2 * np.pi
    Q = P[idx] + 1e-4 * spacing * np.column_stack([np.cos(ang), np.sin(ang)])
    XY = np.unique(np.concatenate([corners, P, B, Q]), axis=0)
    tri = Delaunay(XY)
    F2 = tri.simplices.astype(np.int64)
    # drop degenerate faces (zero chart area) that Delaunay can leave on collinear hull points
    p = XY[F2]
    area = 0.5 * ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    F2 = F2[np.abs(area) > 1e-12 * side.prod()]
    F2[area[np.abs(area) > 1e-12 * side.prod()] < 0] = F2[area[np.abs(area) > 1e-12 * side.prod()] < 0][:, [0, 2, 1]]
    used = np.unique(F2)
    remap = -np.ones(len(XY), dtype=np.int64); remap[used] = np.arange(len(used))
    XY, F2 = XY[used], remap[F2]
    V2 = np.column_stack([XY, z(XY[:, 0], XY[:, 1])])
    info.update(clusters=nc, cluster_points=int(len(P)), background_points=int(n_bg), near_coincident_pairs=200)

E = np.concatenate([F2[:, [0, 1]], F2[:, [1, 2]], F2[:, [2, 0]]])
el = np.linalg.norm(V2[E[:, 0]] - V2[E[:, 1]], axis=1)
p = V2[F2]
cr = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
area3 = 0.5 * np.linalg.norm(cr, axis=1)
lmax = np.max(np.stack([np.linalg.norm(p[:, 1] - p[:, 0], axis=1), np.linalg.norm(p[:, 2] - p[:, 1], axis=1), np.linalg.norm(p[:, 0] - p[:, 2], axis=1)]), 0)
aspect = lmax ** 2 / np.maximum(2 * area3, 1e-300)  # longest edge over the height on it (1 for a right isosceles ~)
slope = np.degrees(np.arccos(np.clip(np.abs(cr[:, 2]) / np.maximum(np.linalg.norm(cr, axis=1), 1e-300), 0, 1)))
info.update(vertices=int(len(V2)), faces=int(len(F2)), edge_mean=float(el.mean()), edge_min=float(el.min()), edge_max=float(el.max()),
            aspect_p50=float(np.median(aspect)), aspect_p99=float(np.percentile(aspect, 99)), aspect_max=float(aspect.max()),
            slope_p50_deg=float(np.median(slope)), slope_p99_deg=float(np.percentile(slope, 99)))
write_off(a.out, V2, F2)
json.dump(info, open(os.path.splitext(a.out)[0] + ".json", "w"), indent=1)
print(f"{os.path.basename(a.out)}: {a.kind}, {len(V2)} vertices, {len(F2)} faces, edge {el.min():.3g}-{el.max():.3g} (mean {el.mean():.3g}), "
      f"aspect p50 {np.median(aspect):.1f} p99 {np.percentile(aspect, 99):.1f}, slope p50 {np.median(slope):.0f} p99 {np.percentile(slope, 99):.0f} deg")
