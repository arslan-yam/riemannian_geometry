import numpy as np
from .mesh import TerrainMesh


# --------------------------------------------------------------------------------------------
# exact labels on a terrain. the refinement pipeline is the general-manifold route; on a
# polyhedral surface the geodesic distance from one point to every vertex is computable
# exactly by the MMP algorithm (pygeodesic, the same code the terrain-oracle papers use as
# ground truth). a portal is not a vertex, so it is inserted as a steiner vertex first (the
# face containing it is split in three, which changes nothing about the surface)
# --------------------------------------------------------------------------------------------
class ExactTerrain:
    def __init__(self, tm):
        if getattr(tm, "weighted", False):
            raise ValueError("no exact solver for weighted terrains: use refined labels (exact=None)")
        self.tm = tm

    def insert_points(self, X):
        #returns (V, F, idx) with the points of X added as vertices idx[i] on the surface.
        #several points in one face are handled by splitting sub-triangles recursively
        tm = self.tm
        X = np.atleast_2d(np.asarray(X, dtype=float))
        faces = tm.faces_of(X)
        V = list(map(tuple, tm.V))
        subs = {} #original face -> list of sub-triangles (vertex index triples)
        idx = np.full(len(X), -1, dtype=int)
        xy = [tuple(p) for p in tm.xy]
        for i, (p, f) in enumerate(zip(X, faces)):
            if f < 0:
                continue
            z = tm.a[f] * p[0] + tm.b[f] * p[1] + tm.c[f]
            V.append((p[0], p[1], z)); xy.append((p[0], p[1]))
            vi = len(V) - 1
            idx[i] = vi
            tris = subs.get(int(f), [tuple(int(v) for v in tm.F[f])])
            #find the sub-triangle containing p (largest minimum barycentric coordinate)
            best, bestw = 0, -np.inf
            for k, (a, b, c) in enumerate(tris):
                (x0, y0), (x1, y1), (x2, y2) = xy[a], xy[b], xy[c]
                det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
                if det == 0:
                    continue
                l1 = ((p[0] - x0) * (y2 - y0) - (x2 - x0) * (p[1] - y0)) / det
                l2 = ((x1 - x0) * (p[1] - y0) - (p[0] - x0) * (y1 - y0)) / det
                w = min(l1, l2, 1.0 - l1 - l2)
                if w > bestw:
                    best, bestw = k, w
            a, b, c = tris.pop(best)
            tris += [(a, b, vi), (b, c, vi), (c, a, vi)]
            subs[int(f)] = tris
        F = []
        for f in range(len(tm.F)):
            F += subs.get(f, [tuple(int(v) for v in tm.F[f])])
        V = np.array(V)
        V[:, :2] += tm.origin #pygeodesic sees the original coordinates, like every other tool
        return V, np.array(F, dtype=int), idx

    def distances_from_points(self, X, targets=None):
        #exact geodesic distance from every point of X to every vertex of the terrain (rows)
        import pygeodesic.geodesic as pg
        V, F, idx = self.insert_points(X)
        geo = pg.PyGeodesicAlgorithmExact(V, F)
        n = len(self.tm.V)
        out = np.full((len(X), n), np.inf)
        tgt = np.arange(n) if targets is None else np.asarray(targets, dtype=int)
        for i, vi in enumerate(idx):
            if vi < 0:
                continue
            d, _ = geo.geodesicDistances(np.array([vi]), tgt)
            out[i, tgt] = d
        return out

    def region_geo(self, box, portals=()):
        #exact solver on the part of the terrain that any geodesic between points of `box`
        #(and portals within its padded window) can visit. a path of length L stays within
        #euclidean distance L of its start (lambda_min = 1 on a terrain), and any two points
        #of the box are at most lambda_max(box) * diagonal apart, so faces within that radius
        #of the box suffice. deep cuts of the separator tree get small sub-meshes, which turns
        #the O(n) sweep per portal into O(region), and the total build from O(n^2) to ~O(n log n)
        import pygeodesic.geodesic as pg
        tm = self.tm
        x0, x1, y0, y1 = box
        portals = np.atleast_2d(np.asarray(portals, dtype=float)) if len(portals) else np.zeros((0, 2))
        if len(portals):
            x0, x1 = min(x0, portals[:, 0].min()), max(x1, portals[:, 0].max())
            y0, y1 = min(y0, portals[:, 1].min()), max(y1, portals[:, 1].max())
        c = tm.xy[tm.F].mean(1)
        inbox = (c[:, 0] >= x0) & (c[:, 0] <= x1) & (c[:, 1] >= y0) & (c[:, 1] <= y1)
        slope = np.sqrt(tm.a ** 2 + tm.b ** 2)
        lam = np.sqrt(1.0 + (slope[inbox].max() if inbox.any() else slope.max()) ** 2)
        R = lam * np.hypot(x1 - x0, y1 - y0)
        near = (tm.xy[:, 0] >= x0 - R) & (tm.xy[:, 0] <= x1 + R) & (tm.xy[:, 1] >= y0 - R) & (tm.xy[:, 1] <= y1 + R)
        faces = np.nonzero(near[tm.F].any(1))[0]
        used = np.unique(tm.F[faces])
        remap = -np.ones(len(tm.V), dtype=int); remap[used] = np.arange(len(used))
        V = tm.V[used].copy()
        F = remap[tm.F[faces]]
        pidx = np.full(len(portals), -1, dtype=int)
        if len(portals):
            #steiner-insert the portals into the sub-mesh (same recursive face split as insert_points).
            #V is in the chart of tm; a sub-mesh built with shift=False keeps that chart
            base = self if len(faces) == len(tm.F) else ExactTerrain(TerrainMesh(V, F, shift=False))
            V2, F2, pidx = base.insert_points(portals)
            V2[:, :2] -= base.tm.origin #back to the chart of tm
            V, F = V2, F2
        Vg = V.copy(); Vg[:, :2] += tm.origin
        return pg.PyGeodesicAlgorithmExact(Vg, F), remap, pidx

    def submesh(self, vertices, rings=3):
        #faces within `rings` vertex rings of the given vertices, reindexed. distances on it are
        #admissible (paths stay on the surface) and exact whenever the true geodesic stays inside
        tm = self.tm
        seed = set(int(v) for v in vertices)
        for _ in range(rings):
            grow = set()
            for v in seed:
                for f in tm.vfaces[v]:
                    grow.update(int(u) for u in tm.F[f])
            seed |= grow
        faces = sorted({f for v in seed for f in tm.vfaces[v]})
        used = np.unique(tm.F[faces])
        remap = -np.ones(len(tm.V), dtype=int); remap[used] = np.arange(len(used))
        V = tm.V[used].copy(); V[:, :2] += tm.origin
        return V, remap[tm.F[faces]], remap

    def local_table(self, members, rings=3, halo=None):
        #exact pairwise distances among `members` on their local sub-mesh. with `halo`, a set of
        #extra vertices, also returns each member's distances to those (same sweeps, more targets)
        import pygeodesic.geodesic as pg
        V, F, remap = self.submesh(members, rings)
        geo = pg.PyGeodesicAlgorithmExact(V, F)
        members = np.asarray(members, dtype=int)
        loc = remap[members]
        extra = np.array(sorted(v for v in (halo or ()) if remap[v] >= 0), dtype=int)
        tgt = np.concatenate([loc, remap[extra]]) if len(extra) else loc
        D = np.full((len(members), len(members)), np.inf)
        H = np.full((len(members), len(extra)), np.inf)
        for i, vi in enumerate(loc):
            d, _ = geo.geodesicDistances(np.array([vi]), tgt)
            D[i] = d[:len(loc)]
            H[i] = d[len(loc):]
        D = np.minimum(D, D.T)
        return (D, extra, H) if halo is not None else D
