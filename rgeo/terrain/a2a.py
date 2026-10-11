import numpy as np
from .mesh import TerrainMesh
from .exact import ExactTerrain


# --------------------------------------------------------------------------------------------
# arbitrary-point queries. the index stores labels for mesh vertices; a surface point p inside
# face f is attached to them in one of two ways:
#   mode="snap":   D(p,q) <= min_{a,b} d_loc(p,a) + D(a,b) + d_loc(b,q), an admissible upper
#                  bound. rings=0 takes the three vertices of the face (d_loc = straight 3-d
#                  segment inside the face), rings=1 the one-ring fan with d_loc from an exact
#                  MMP sweep on the fan. first order: the corner at the vertex costs
#                  |p-a|(1-cos theta), tens of metres on these meshes
#   mode="interp": bilinear (barycentric x barycentric) interpolation of the 3 x 3 matrix of
#                  vertex-pair answers D(a, b) over the two faces. the distance field is smooth
#                  inside a face away from the cut locus (hessian ~ 1/d), so the error is
#                  O(edge^2 / d): second order, two-sided like the parabola interpolation
# pairs whose faces are within `near` rings of each other are answered exactly by one MMP sweep
# on the local sub-mesh with both points inserted (interpolation is no good when d ~ edge)
# --------------------------------------------------------------------------------------------
class TerrainA2A:
    def __init__(self, tm, index, mode="interp", rings=0, near=2):
        self.tm = tm
        self.index = index
        self.mode = mode
        self.rings = rings
        self.near = near
        #faces within `near` rings of every face, for the exact local route of close pairs
        vf = tm.vfaces
        self.near_faces = []
        for f in range(len(tm.F)):
            fs = {f}
            for _ in range(near):
                fs |= {g for h in fs for v in tm.F[h] for g in vf[v]}
            self.near_faces.append(fs)

    def bary(self, p, f):
        tri = self.tm.F[int(f)]
        (x0, y0), (x1, y1), (x2, y2) = self.tm.xy[tri]
        det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        l1 = ((p[0] - x0) * (y2 - y0) - (x2 - x0) * (p[1] - y0)) / det
        l2 = ((x1 - x0) * (p[1] - y0) - (p[0] - x0) * (y1 - y0)) / det
        return tri, np.clip(np.array([1.0 - l1 - l2, l1, l2]), 0.0, 1.0)

    def local(self, p, f):
        #vertices near face f and exact surface distances from p (chart point inside f) to them
        tm = self.tm
        f = int(f)
        tri = tm.F[f]
        if self.rings == 0:
            z = tm.a[f] * p[0] + tm.b[f] * p[1] + tm.c[f]
            return tri, np.linalg.norm(tm.V[tri] - np.array([p[0], p[1], z]), axis=1)
        import pygeodesic.geodesic as pg
        faces = sorted({g for v in tri for g in tm.vfaces[v]})
        used = np.unique(tm.F[faces])
        remap = -np.ones(len(tm.V), dtype=int); remap[used] = np.arange(len(used))
        sub = ExactTerrain(TerrainMesh(tm.V[used].copy(), remap[tm.F[faces]], shift=False))
        V2, F2, idx = sub.insert_points(np.asarray(p, dtype=float)[None])
        if idx[0] < 0:
            z = tm.a[f] * p[0] + tm.b[f] * p[1] + tm.c[f]
            return tri, np.linalg.norm(tm.V[tri] - np.array([p[0], p[1], z]), axis=1)
        geo = pg.PyGeodesicAlgorithmExact(V2, F2)
        d, _ = geo.geodesicDistances(np.array([int(idx[0])]), np.arange(len(used)))
        return used, np.asarray(d)

    def exact_near(self, p, fp, q, fq):
        #both points inserted into the sub-mesh around the two faces, one exact sweep
        import pygeodesic.geodesic as pg
        tm = self.tm
        faces = sorted(self.near_faces[int(fp)] | self.near_faces[int(fq)])
        used = np.unique(tm.F[faces])
        remap = -np.ones(len(tm.V), dtype=int); remap[used] = np.arange(len(used))
        sub = ExactTerrain(TerrainMesh(tm.V[used].copy(), remap[tm.F[faces]], shift=False))
        V2, F2, idx = sub.insert_points(np.array([p, q], dtype=float))
        if (idx < 0).any():
            return np.inf
        geo = pg.PyGeodesicAlgorithmExact(V2, F2)
        return float(geo.geodesicDistance(int(idx[0]), int(idx[1]))[0])

    def query(self, p, fp, q, fq):
        ix = self.index
        fp, fq = int(fp), int(fq)
        if fq in self.near_faces[fp]:
            return self.exact_near(p, fp, q, fq)
        if self.mode == "snap":
            Np, dp = self.local(p, fp)
            Nq, dq = self.local(q, fq)
        else:
            Np, wp = self.bary(p, fp)
            Nq, wq = self.bary(q, fq)
        #the vertex-pair oracle on every (a, b): its own leaf table, halo, portals, interpolation.
        #(interpolating the label vectors first and taking one portal minimum would be 9x cheaper
        #but has a jensen gap, first order in the spread of the vertices' argmins: 0.8% vs 0.18%)
        T = np.array([[0.0 if a == b else ix.query(int(a), int(b)) for b in Nq] for a in Np])
        if self.mode == "snap":
            return float((dp[:, None] + T + dq[None, :]).min())
        return float(wp @ T @ wq)


# --------------------------------------------------------------------------------------------
# certified arbitrary-point queries on top of the Certified Cut Index. For points p in face f and
# q in face g (w_f, w_g their weights, a / b the corners of f / g):
#   * |p - q| <= radius (default: the index's near radius R = g/(eps lambda)): exact. One MMP sweep
#     on the sub-mesh of faces meeting bbox(p, q) dilated by UB/(2 lambda) -- by the ellipse lemma
#     every geodesic of length <= UB lies inside -- with both points inserted. LB = UB = d (rounding)
#   * otherwise, from the certified vertex-pair intervals [LB(a,b), UB(a,b)] and d(p,a) <= w_f|pa|
#     (the straight segment inside the face is a path):
#       UB = min_ab w_f|pa| + UB(a,b) + w_g|qb|            (the cost of a real path)
#       LB = max(max_ab LB(a,b) - w_f|pa| - w_g|qb|, lambda |PQ|_3d)   (triangle inequality)
#     estimate = bilinear interpolation of the 3 x 3 vertex estimates, clamped into [LB, UB].
# Guarantee: LB <= d <= UB for every pair; UB <= (1+eps) d + 2(1+eps)(w_f r_p + w_g r_q) with r the
# distance from a point to the nearest corner of its face, and exact answers below the radius. With
# `width` (default 2 eps) a pair whose interval is wider than 1 + width is computed exactly too, so
# every answer has UB <= (1 + width) LB: a hard (1 + 2 eps) guarantee for arbitrary points.
# --------------------------------------------------------------------------------------------
class CertifiedA2A:
    #width: certify-or-refine. if the corner interval is wider than UB/LB = 1 + width, the pair is
    #computed exactly as well (ellipse sub-mesh of the corner UB), so every answer satisfies
    #UB <= (1 + width) LB, i.e. a hard (1 + width) guarantee; None keeps the corner bounds as they are
    def __init__(self, tm, index, query=None, radius=None, width="auto"):
        self.tm = tm
        self.ix = index
        self.q = query if query is not None else index #anything with .bounds(s, t) -> (lb, ub, est)
        self.radius = float(index.Rq if radius is None else radius)
        self.lam = float(index.lam)
        self.exact_ok = bool(getattr(index.oracle, "exact", False))
        self.width = 2.0 * float(index.eps) if width == "auto" else width
        self.refined = 0

    def _lift(self, p, f):
        tm = self.tm
        return np.array([p[0], p[1], tm.a[f] * p[0] + tm.b[f] * p[1] + tm.c[f]])

    def _vertex_bounds(self, p, fp, q, fq):
        tm = self.tm
        P, Q = self._lift(p, fp), self._lift(q, fq)
        A, B = tm.F[fp], tm.F[fq]
        da = tm.w[fp] * np.linalg.norm(tm.V[A] - P, axis=1) #upper bounds of d(p, a)
        db = tm.w[fq] * np.linalg.norm(tm.V[B] - Q, axis=1)
        L = np.empty((3, 3)); U = np.empty((3, 3)); E = np.empty((3, 3))
        for i, a in enumerate(A):
            for j, b in enumerate(B):
                L[i, j], U[i, j], E[i, j] = (0.0, 0.0, 0.0) if a == b else self.q.bounds(int(a), int(b))
        ub = float((da[:, None] + U + db[None, :]).min())
        lb = max(float((L - da[:, None] - db[None, :]).max()), self.lam * float(np.linalg.norm(P - Q)))
        _, wp = _bary1(tm, p, fp); _, wq = _bary1(tm, q, fq)
        est = float(wp @ E @ wq)
        return lb, ub, est

    def _exact(self, p, fp, q, fq, ub):
        from .submesh import extract, insert_points, ExactSolver
        tm = self.tm
        R = ub / (2.0 * self.lam) + 1e-4 * self.ix.edge #+ the inserted points may move by ~1e-5 edges
        lo = np.minimum(p, q) - R; hi = np.maximum(p, q) + R
        faces = self.ix.grid.faces_in_box(lo[0], hi[0], lo[1], hi[1])
        V, F, used, remap = extract(tm, faces)
        X = np.array([p, q], dtype=float)
        V2, F2, idx, Xu = insert_points(tm, V, F, faces, X)
        if (idx < 0).any():
            return None
        d = float(ExactSolver(V2, F2).dist(int(idx[0]), [int(idx[1])])[0])
        #a point within 1e-7 of an edge is moved ~1e-5 edge lengths: account for the move
        slack = sum(float(np.linalg.norm(self._lift(x, f) - self._lift(xu, f))) * tm.w[f]
                    for x, xu, f in ((p, Xu[0], fp), (q, Xu[1], fq)))
        return d, slack

    def bounds(self, p, fp, q, fq):
        p = np.asarray(p, dtype=float); q = np.asarray(q, dtype=float)
        fp, fq = int(fp), int(fq)
        lb, ub, est = self._vertex_bounds(p, fp, q, fq)
        loose = self.width is not None and ub > (1.0 + self.width) * lb
        if self.exact_ok and (np.hypot(*(p - q)) <= self.radius or loose):
            self.refined += int(loose and np.hypot(*(p - q)) > self.radius)
            r = self._exact(p, fp, q, fq, ub)
            if r is not None:
                d, slack = r
                lb2, ub2 = d * (1 - 1e-9) - slack, d * (1 + 1e-9) + slack
                lb, ub = max(lb, lb2), min(ub, ub2)
                est = d
        return lb, ub, min(max(est, lb), ub)

    def query(self, p, fp, q, fq):
        return self.bounds(p, fp, q, fq)[2]


def _bary1(tm, p, f):
    tri = tm.F[int(f)]
    (x0, y0), (x1, y1), (x2, y2) = tm.xy[tri]
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    l1 = ((p[0] - x0) * (y2 - y0) - (x2 - x0) * (p[1] - y0)) / det
    l2 = ((x1 - x0) * (p[1] - y0) - (p[0] - x0) * (y1 - y0)) / det
    return tri, np.clip(np.array([1.0 - l1 - l2, l1, l2]), 0.0, 1.0)
