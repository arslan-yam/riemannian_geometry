"""OFF mesh input/output and repair of folded height fields."""
import numpy as np


def read_off(path):
    with open(path) as f:
        head = f.readline().strip()
        if head != "OFF":
            raise ValueError(f"not an OFF file: {path}")
        nv, nf, _ = (int(x) for x in f.readline().split())
        V = np.array([[float(x) for x in f.readline().split()[:3]] for _ in range(nv)])
        F = np.array([[int(x) for x in f.readline().split()[1:4]] for _ in range(nf)], dtype=int)
    return V, F


def write_off(path, V, F):
    with open(path, "w") as f:
        f.write("OFF\n")
        f.write(f"{len(V)} {len(F)} 0\n")
        for v in V:
            f.write(f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for t in F:
            f.write(f"3 {t[0]} {t[1]} {t[2]}\n")


def projected_areas(xy, F):
    P = xy[F]
    return 0.5 * ((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1])
                  - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))


def repair_height_field(V, F, verbose=True):
    #a terrain is only a height field if no face folds over in the (x, y) projection. the public
    #meshes contain the odd folded sliver (BigMount has one), which makes point location ambiguous
    #and the triangulation invalid. unfold it by sliding one vertex of the sliver toward the
    #centroid of its 1-ring, by the smallest fraction that makes every incident face positively
    #oriented with a margin. z is kept, so the surface changes only inside that vertex's 1-ring
    V = np.array(V, dtype=float); F = np.asarray(F, dtype=int)
    xy = V[:, :2]
    area = projected_areas(xy, F)
    if (area < 0).sum() > (area > 0).sum():
        F = F[:, ::-1].copy(); area = -area #whole mesh wound clockwise: flip everything
    margin = 0.1 * np.median(np.abs(area))
    for _ in range(100):
        bad = np.nonzero(area <= 1e-9 * np.median(np.abs(area)))[0] #folded or degenerate only
        if len(bad) == 0:
            break
        f = bad[0]
        best = None
        for v in F[f]:
            inc = np.nonzero((F == v).any(1))[0]
            ring = np.unique(F[inc]); ring = ring[ring != v]
            cen = xy[ring].mean(0)
            lo, hi = 0.0, 1.0
            if projected_areas(np.where(np.arange(len(xy))[:, None] == v, cen, xy), F)[inc].min() < margin:
                continue
            for _ in range(30): #bisection on the slide fraction
                s = 0.5 * (lo + hi)
                xy2 = xy.copy(); xy2[v] = xy[v] + s * (cen - xy[v])
                if projected_areas(xy2, F)[inc].min() >= margin: hi = s
                else: lo = s
            disp = hi * np.linalg.norm(cen - xy[v])
            if best is None or disp < best[0]:
                best = (disp, v, xy[v] + hi * (cen - xy[v]))
        if best is None:
            raise ValueError(f"cannot unfold face {f}")
        disp, v, new = best
        if verbose:
            print(f"repair: face {f} folded, moved vertex {v} by {disp:.1f} (median edge "
                  f"{np.median(np.linalg.norm(xy[F[:, 0]] - xy[F[:, 1]], axis=1)):.1f})")
        xy[v] = new
        area = projected_areas(xy, F)
    V[:, :2] = xy
    return V, F
