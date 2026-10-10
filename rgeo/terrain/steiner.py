"""Steiner base graph for weighted terrains, and labels from it.

On a weighted terrain (face costs w_f) shortest paths refract at edges (Snell's law) and hug cheap
edges, so the mesh graph is the wrong skeleton: it routes ~9% longer on HorseMount with costs in
(1, 5). The Steiner graph of the weighted-region literature places `m` points on every edge and
connects, inside every face, every pair of its boundary points by a straight segment costing
w_f x 3-d length (points on the same edge are chained along it at the cheaper adjacent cost). A path
in this graph is a real surface curve, so graph distances are admissible upper bounds, and they
converge to the weighted geodesic distance as m grows.

`SteinerLabels` plugs the graph into `CertifiedCutIndex` in place of exact MMP labels: a portal is
a source seeded at the boundary points of its face, labels and near tables come from a multi-seed
Dijkstra that stops at its last target (`rgeo/graphs/cdijkstra.py`, C) on the sub-graph the window
lemma allows.
"""
import numpy as np
from scipy.sparse import csr_matrix
from ..oracles.certified import f32_up
from ..graphs.cdijkstra import multi_dijkstra


class SteinerGraph:
    def __init__(self, tm, m=5):
        self.tm, self.m = tm, int(m)
        n = len(tm.V)
        F = tm.F
        #unique edges, with the index of each face's three edges and their orientation
        E = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
        key = np.sort(E, 1)
        uk, inv = np.unique(key[:, 0].astype(np.int64) * n + key[:, 1], return_inverse=True)
        self.edges = np.stack([uk // n, uk % n], 1)
        ne = len(self.edges)
        face_edge = inv.reshape(3, -1).T #(nf, 3): edge id of (v0v1, v1v2, v2v0)
        #face weights per edge: min over the (one or two) adjacent faces
        wmin = np.full(ne, np.inf)
        np.minimum.at(wmin, inv, np.tile(tm.w, 3))
        m = self.m
        t = np.arange(1, m + 1) / (m + 1.0)
        A, B = self.edges[:, 0], self.edges[:, 1]
        Pst = tm.V[A][:, None, :] + t[None, :, None] * (tm.V[B] - tm.V[A])[:, None, :] #(ne, m, 3)
        self.P = np.vstack([tm.V, Pst.reshape(-1, 3)]) #node coordinates (chart x, y, z)
        self.N = len(self.P)
        steiner_id = n + np.arange(ne * m).reshape(ne, m)
        rows, cols, vals = [], [], []
        #chains along every edge: A, s_1, ..., s_m, B at the cheaper adjacent face cost
        chain = np.column_stack([A, steiner_id, B]) #(ne, m+2)
        a, b = chain[:, :-1].ravel(), chain[:, 1:].ravel()
        w = np.repeat(wmin, m + 1) * np.linalg.norm(self.P[a] - self.P[b], axis=1)
        rows += [a]; cols += [b]; vals += [w]
        #inside every face: all pairs of boundary points not lying on a common edge
        nf = len(F)
        bp = np.zeros((nf, 3 + 3 * m), dtype=np.int64) #boundary points of each face
        lab = np.zeros((3 + 3 * m, 3), dtype=bool) #which of the face's edges each point lies on
        bp[:, :3] = F
        lab[0, [0, 2]] = lab[1, [0, 1]] = lab[2, [1, 2]] = True
        for k in range(3):
            eid = face_edge[:, k]
            u0 = F[:, k]
            fwd = self.edges[eid, 0] == u0 #edge stored as (u0 -> u1)? else reversed
            ids = steiner_id[eid]
            ids = np.where(fwd[:, None], ids, ids[:, ::-1])
            bp[:, 3 + k * m:3 + (k + 1) * m] = ids
            lab[3 + k * m:3 + (k + 1) * m, k] = True
        I, J = np.triu_indices(3 + 3 * m, 1)
        keep = ~(lab[I] & lab[J]).any(1)
        I, J = I[keep], J[keep]
        a, b = bp[:, I].ravel(), bp[:, J].ravel()
        w = np.repeat(tm.w, len(I)) * np.linalg.norm(self.P[a] - self.P[b], axis=1)
        rows += [a]; cols += [b]; vals += [w]
        r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(vals)
        self.G = csr_matrix((np.r_[v, v], (np.r_[r, c], np.r_[c, r])), shape=(self.N, self.N))
        self.face_points = bp #(nf, 3+3m) node ids of every face's boundary points
        self.edge_steiner = steiner_id

    def nodes_of_faces(self, faces):
        return np.unique(self.face_points[faces])


class SteinerLabels:
    name = "Steiner graph"
    exact = False

    def __init__(self, tm, grid, m=5, graph=None):
        self.tm, self.grid = tm, grid
        self.sg = graph or SteinerGraph(tm, m)

    def _subgraph(self, faces):
        nodes = self.sg.nodes_of_faces(faces)
        return nodes, self.sg.G[nodes][:, nodes]

    def cut(self, faces, X, along, members):
        #labels from every portal (a point inside its face: seeded at the face's boundary points
        #at w_f x straight distance) to every member, on the sub-graph of the given faces
        tm, sg = self.tm, self.sg
        X = np.atleast_2d(np.asarray(X, dtype=float)).copy()
        nodes, Gs = self._subgraph(faces)
        loc = np.full(sg.N, -1, dtype=np.int64)
        loc[nodes] = np.arange(len(nodes))
        gf = tm.faces_of(X)
        fs = np.maximum(gf, 0)
        X3 = np.column_stack([X, tm.a[fs] * X[:, 0] + tm.b[fs] * X[:, 1] + tm.c[fs]])
        seeds = []
        for j in range(len(X)):
            if gf[j] < 0: #not located: no seeds, the label stays inf and the cut drops the portal
                seeds.append((np.zeros(0, dtype=np.int64), np.zeros(0)))
                continue
            q = sg.face_points[gf[j]]
            ql = loc[q]
            ok = ql >= 0
            seeds.append((ql[ok], tm.w[gf[j]] * np.linalg.norm(sg.P[q[ok]] - X3[j], axis=1)))
        tgt = loc[members]
        if (tgt < 0).any():
            raise RuntimeError("a member is outside the label sub-graph")
        D = multi_dijkstra(Gs, seeds, tgt) #(portals, members), each search stops at its last member
        return f32_up(D.T), X

    def near(self, faces, members, targets):
        sg = self.sg
        nodes, Gs = self._subgraph(faces)
        loc = np.full(sg.N, -1, dtype=np.int64)
        loc[nodes] = np.arange(len(nodes))
        seeds = [(np.array([loc[s]], dtype=np.int64), np.zeros(1)) for s in members]
        res = multi_dijkstra(Gs, seeds, [loc[np.asarray(T, dtype=np.int64)] for T in targets])
        return [f32_up(r) if len(T) else np.zeros(0, np.float32) for r, T in zip(res, targets)]
