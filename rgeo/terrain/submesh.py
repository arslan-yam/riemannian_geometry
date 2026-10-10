"""Fast sub-mesh extraction, Steiner insertion and early-terminating exact sweeps.

Three facts make exact labels cheap:
* `FaceGrid` answers "which faces meet this box" in time proportional to the answer, so a cut
  only ever touches the part of the terrain its geodesics can visit.
* `insert_points` inserts many Steiner points at once (vectorised 1->3 splits, a small Delaunay
  triangulation for faces that receive several points), without rebuilding a TerrainMesh.
* pygeodesic's MMP only honours its stop-at-targets condition when `max_distance` is finite (the
  default is infinite, so every sweep floods the whole mesh). With `max_distance=0` the sweep
  stops as soon as every target's distance is final: identical distances, 10-200x less work.
"""
import numpy as np
from scipy.spatial import Delaunay


class FaceGrid:
    #uniform grid over the chart; every face is registered in each cell its bounding box meets
    def __init__(self, tm, cell=None):
        self.tm = tm
        P = tm.xy[tm.F]
        self.flo = P.min(1)
        self.fhi = P.max(1)
        e = float(np.median(np.linalg.norm(P[:, 0] - P[:, 1], axis=1)))
        self.cell = float(cell or 2.0 * e)
        self.org = tm.xy.min(0) - 1e-9 * max(1.0, self.cell)
        i0 = np.floor((self.flo - self.org) / self.cell).astype(np.int64)
        i1 = np.floor((self.fhi - self.org) / self.cell).astype(np.int64)
        self.nc = i1.max(0) + 1
        nx = i1[:, 0] - i0[:, 0] + 1
        ny = i1[:, 1] - i0[:, 1] + 1
        cnt = nx * ny
        fid = np.repeat(np.arange(len(tm.F)), cnt)
        k = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        nyr = np.repeat(ny, cnt)
        cx = np.repeat(i0[:, 0], cnt) + k // nyr
        cy = np.repeat(i0[:, 1], cnt) + k % nyr
        cid = cx * self.nc[1] + cy
        order = np.argsort(cid, kind="stable")
        self.cell_faces = fid[order]
        self.cell_ptr = np.searchsorted(cid[order], np.arange(self.nc[0] * self.nc[1] + 1))

    def faces_in_box(self, x0, x1, y0, y1):
        a0 = np.clip(np.floor((np.array([x0, y0]) - self.org) / self.cell).astype(np.int64), 0, self.nc - 1)
        a1 = np.clip(np.floor((np.array([x1, y1]) - self.org) / self.cell).astype(np.int64), 0, self.nc - 1)
        rows = []
        for cx in range(a0[0], a1[0] + 1):
            lo = self.cell_ptr[cx * self.nc[1] + a0[1]]
            hi = self.cell_ptr[cx * self.nc[1] + a1[1] + 1]
            rows.append(self.cell_faces[lo:hi])
        f = np.unique(np.concatenate(rows)) if rows else np.zeros(0, dtype=np.int64)
        keep = (self.flo[f, 0] <= x1) & (self.fhi[f, 0] >= x0) & (self.flo[f, 1] <= y1) & (self.fhi[f, 1] >= y0)
        return f[keep]


def extract(tm, faces):
    #the sub-mesh made of the given global faces: (V, F) in local numbering, the global ids of the
    #local vertices, and the global -> local vertex map (-1 outside)
    Fg = tm.F[faces]
    used = np.unique(Fg)
    remap = np.full(len(tm.V), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    return tm.V[used].copy(), remap[Fg], used, remap


def _bary(P, a, b, c):
    #barycentric coordinates of the points P (k,2) in the triangles (a, b, c) (each (k,2))
    det = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    l1 = ((P[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (P[:, 1] - a[:, 1])) / det
    l2 = ((b[:, 0] - a[:, 0]) * (P[:, 1] - a[:, 1]) - (P[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])) / det
    return np.stack([1.0 - l1 - l2, l1, l2], 1)


def insert_points(tm, V, F, faces, X, along=None, tol=1e-7):
    """Insert the chart points X as Steiner vertices of the sub-mesh (V, F) made of the global faces
    `faces`. Returns (V2, F2, idx, Xused): new mesh, the local vertex index of every point, and the
    positions actually inserted (a point within `tol` (relative) of a face edge is moved along
    `along` -- e.g. the cut direction, so a portal stays on its cut -- or towards the face centroid).
    The surface is unchanged: every new vertex lies on its face's plane."""
    X = np.atleast_2d(np.asarray(X, dtype=float)).copy()
    k = len(X)
    idx = np.full(k, -1, dtype=np.int64)
    if k == 0:
        return V, F, idx, X
    gf = tm.faces_of(X)
    fmap = np.full(len(tm.F), -1, dtype=np.int64)
    fmap[faces] = np.arange(len(faces))
    lf = np.where(gf >= 0, fmap[np.maximum(gf, 0)], -1)
    ok = lf >= 0
    #keep points off the edges of their face
    tri = tm.xy[tm.F[np.maximum(gf, 0)]]
    for _ in range(8):
        B = _bary(X, tri[:, 0], tri[:, 1], tri[:, 2])
        bad = ok & (B.min(1) < tol)
        if not bad.any():
            break
        scale = np.linalg.norm(tri[bad, 1] - tri[bad, 0], axis=1)
        if along is not None:
            d = np.asarray(along, dtype=float)
            #step along `along`, in the direction that increases the smallest barycentric weight
            plus = X[bad] + 1e-5 * scale[:, None] * d
            minus = X[bad] - 1e-5 * scale[:, None] * d
            bp = _bary(plus, tri[bad, 0], tri[bad, 1], tri[bad, 2]).min(1)
            bm = _bary(minus, tri[bad, 0], tri[bad, 1], tri[bad, 2]).min(1)
            X[bad] = np.where((bp >= bm)[:, None], plus, minus)
        else:
            c = tri[bad].mean(1)
            X[bad] = X[bad] + 1e-5 * (c - X[bad])
    #drop duplicates within a face
    z = tm.a[np.maximum(gf, 0)] * X[:, 0] + tm.b[np.maximum(gf, 0)] * X[:, 1] + tm.c[np.maximum(gf, 0)]
    newV = np.column_stack([X, z])
    nV0 = len(V)
    idx[ok] = nV0 + np.arange(ok.sum())
    pts_local = np.nonzero(ok)[0]
    V2 = np.vstack([V, newV[pts_local]])
    order = np.argsort(lf[pts_local], kind="stable")
    pl = pts_local[order]
    lfs = lf[pl]
    uniq, start, cnt = np.unique(lfs, return_index=True, return_counts=True)
    keep_face = np.ones(len(F), dtype=bool)
    keep_face[uniq] = False
    new_faces = []
    one = cnt == 1
    if one.any():
        f1 = uniq[one]; p1 = idx[pl[start[one]]]
        a, b, c = F[f1, 0], F[f1, 1], F[f1, 2]
        new_faces.append(np.stack([np.stack([a, b, p1], 1), np.stack([b, c, p1], 1), np.stack([c, a, p1], 1)], 1).reshape(-1, 3))
    for fi, s0, cn in zip(uniq[~one], start[~one], cnt[~one]):
        corners = F[fi]
        pids = idx[pl[s0:s0 + cn]]
        ids = np.concatenate([corners, pids])
        P2 = V2[ids, :2]
        dl = Delaunay(P2)
        T = ids[dl.simplices]
        Q = V2[T][:, :, :2]
        area = (Q[:, 1, 0] - Q[:, 0, 0]) * (Q[:, 2, 1] - Q[:, 0, 1]) - (Q[:, 2, 0] - Q[:, 0, 0]) * (Q[:, 1, 1] - Q[:, 0, 1])
        T = np.where((area < 0)[:, None], T[:, [0, 2, 1]], T)
        new_faces.append(T[np.abs(area) > 0])
    F2 = np.vstack([F[keep_face]] + new_faces) if new_faces else F
    return V2, F2, idx, X


class ExactSolver:
    #MMP (pygeodesic) on a mesh, with early termination at the targets
    def __init__(self, V, F):
        import pygeodesic.geodesic as pg
        self.geo = pg.PyGeodesicAlgorithmExact(np.ascontiguousarray(V, dtype=np.float64),
                                                np.ascontiguousarray(F, dtype=np.int32))

    def dist(self, src, targets):
        d, _ = self.geo.geodesicDistances(np.array([int(src)], dtype=np.int32),
                                          np.asarray(targets, dtype=np.int32), max_distance=0.0)
        return d
