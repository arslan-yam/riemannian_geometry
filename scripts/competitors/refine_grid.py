"""A finer regular grid on the same terrain surface, for scaling beyond the largest UP-Oracle terrain (2.5M faces).

The UP-Oracle terrains are regular grids of one region at several resolutions, resampled from a common elevation model
(their finer grids differ from a linear interpolation of the coarser ones by 0.3 m on average, at most 1.3 m). We
resample the finest one the same way: an N x N grid over its bounding box, each height linearly interpolated inside the
source triangle that contains it. The new mesh lies on the source surface and adds no relief; quads are split along the
same diagonal as in the source grids.

usage: python3 scripts/competitors/refine_grid.py src.off out.off --side 1584      # 1584^2 vertices, 5,011,778 faces
"""
import argparse
import numpy as np
from matplotlib.tri import Triangulation, LinearTriInterpolator

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out"); ap.add_argument("--side", type=int, required=True)
a = ap.parse_args()

with open(a.src) as f:
    assert f.readline().strip() == "OFF"
    nv, nf, _ = map(int, f.readline().split())
    V = np.loadtxt(f, max_rows=nv)
    F = np.loadtxt(f, max_rows=nf, dtype=np.int64)[:, 1:4]
tri = Triangulation(V[:, 0], V[:, 1], F)
interp = LinearTriInterpolator(tri, V[:, 2])
lo, hi = V[:, :2].min(0), V[:, :2].max(0)
W = a.side
xs = np.linspace(lo[0], hi[0], W); ys = np.linspace(lo[1], hi[1], W)
X, Y = np.meshgrid(xs, ys)                    # row-major: x fastest, as in the source grids
Z = np.asarray(interp(X, Y).filled(np.nan))
if np.isnan(Z).any():                         # points on the hull boundary that the interpolator misses: nearest source vertex
    from scipy.spatial import cKDTree
    m = np.isnan(Z)
    Z[m] = V[cKDTree(V[:, :2]).query(np.column_stack([X[m], Y[m]]))[1], 2]
r, c = np.meshgrid(np.arange(W - 1), np.arange(W - 1), indexing="ij")
v = (r * W + c).ravel()
faces = np.concatenate([np.column_stack([v, v + 1, v + W]), np.column_stack([v + 1, v + 1 + W, v + W])])
with open(a.out, "w") as f:
    f.write(f"OFF\n{W * W} {len(faces)} 0\n")
    np.savetxt(f, np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]), fmt="%.6f")
    np.savetxt(f, np.column_stack([np.full(len(faces), 3), faces]), fmt="%d")
print(f"{a.out}: {W * W} vertices, {len(faces)} faces, heights {np.nanmin(Z):.1f}..{np.nanmax(Z):.1f}")
