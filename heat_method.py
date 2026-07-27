import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.spatial import Delaunay, cKDTree
from graph_build import blocked


class HeatMethod:
    def __init__(self, manifold, points, t_factor=1.0, obstacles=(), max_edge=None):
        self.manifold = manifold
        self.points = np.asarray([np.asarray(p, dtype=float) for p in points])
        self.n = len(self.points)

        tri = Delaunay(self.points)
        keep, grads, areas = [], [], []
        spacing = max_edge if max_edge is not None else 3.0 * self.typical_spacing()

        for f in tri.simplices:
            p0, p1, p2 = self.points[f]
            if max(np.linalg.norm(p1 - p0), np.linalg.norm(p2 - p1), np.linalg.norm(p0 - p2)) > spacing:
                continue
            centre = (p0 + p1 + p2) / 3.0
            if not manifold.in_domain(centre) or blocked(centre, obstacles):
                continue
            if any(blocked(0.5 * (a + b), obstacles) for a, b in ((p0, p1), (p1, p2), (p2, p0))):
                continue
            
            dm = np.array([p1 - p0, p2 - p0])
            det = np.linalg.det(dm)
            if abs(det) < 1e-14:
                continue
            inv = np.linalg.inv(dm)
            g = manifold.metric(centre)
            keep.append(f)
            grads.append(np.array([inv @ [-1.0, -1.0], inv @ [1.0, 0.0], inv @ [0.0, 1.0]]))
            areas.append(0.5 * abs(det) * np.sqrt(max(np.linalg.det(g), 0.0)))

        self.faces = np.array(keep, dtype=int)
        self.grads = np.array(grads)
        self.areas = np.array(areas)
        self.ginv = np.array([np.linalg.inv(manifold.metric(self.points[f].mean(0))) for f in self.faces])
        rows, cols, vals = [], [], []
        mass = np.zeros(self.n)
        
        for t, f in enumerate(self.faces):
            for a in range(3):
                mass[f[a]] += self.areas[t] / 3.0
                for b in range(3):
                    rows.append(f[a])
                    cols.append(f[b])
                    vals.append(self.areas[t] * (self.grads[t][a] @ self.ginv[t] @ self.grads[t][b]))

        self.L = sp.csc_matrix((vals, (rows, cols)), shape=(self.n, self.n))
        self.mass = np.where(mass > 0, mass, mass[mass > 0].mean() if (mass > 0).any() else 1.0)
        self.t = t_factor * self.mean_edge() ** 2
        eps = 1e-10 * self.L.diagonal().max()
        self.heat_solve = spl.factorized((sp.diags(self.mass) + self.t * self.L).tocsc())
        self.poisson_solve = spl.factorized((self.L + eps * sp.identity(self.n)).tocsc())

    def typical_spacing(self):
        if self.n < 2:
            return 1.0
        d, idxs = cKDTree(self.points).query(self.points, k=2)
        return float(np.median(d[:, 1]))

    def mean_edge(self):
        if not len(self.faces):
            return 1.0
        p = self.points[self.faces]
        e = np.concatenate([p[:, 1] - p[:, 0], p[:, 2] - p[:, 1], p[:, 0] - p[:, 2]])
        return float(np.mean(np.linalg.norm(e, axis=1)))

    def distances(self, src):
        b = np.zeros(self.n)
        b[src] = 1.0
        u = self.heat_solve(b)
        gu = np.einsum("tai,ta->ti", self.grads, u[self.faces])
        norm = np.sqrt(np.maximum(np.einsum("ti,tij,tj->t", gu, self.ginv, gu), 0.0))
        X = np.where(norm[:, None] > 0, -np.einsum("tij,tj->ti", self.ginv, gu) / np.maximum(norm, 1e-300)[:, None], 0.0)
        contrib = self.areas[:, None] * np.einsum("tai,ti->ta", self.grads, X)
        div = np.bincount(self.faces.ravel(), weights=contrib.ravel(), minlength=self.n)
        phi = self.poisson_solve(div)
        phi = phi - phi[src]
        
        if phi.sum() < 0:
            phi = -phi
        return np.maximum(phi, 0.0)
