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
