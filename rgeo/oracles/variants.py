"""Variants of the separator oracle kept as recorded negative results (see notebooks/exps.ipynb).

* SeparatorHeatIndex: labels from the heat method instead of Dijkstra (smooth but biased).
* SeparatorEikonalLeafIndex / SeparatorEikonalIndex: leaf tables / labels from fast marching on a
  grid; O(h) absolute error and a grid fine enough to beat the graph is far too slow.
* SeparatorContinuousIndex: allows the path to move along the band between two portals
  (min_ij a_i + D_ij + b_j). Provably identical to the plain index for graph labels, since
  a_j <= a_i + D_ij by the triangle inequality: same errors, twice the query time.
* SeparatorBFSRefinedIndex: balanced BFS frontiers (true vertex separators) instead of coordinate
  cuts; no escape round a cut end, but the frontier is not a line, so no interpolation.
"""
import numpy as np
from ..solvers.heat_method import HeatMethod
from ..solvers.eikonal import solve_eikonal, sample_points
from .separator import SeparatorIndex, SeparatorRefinedIndex


class SeparatorHeatIndex(SeparatorIndex):
    def __init__(self, points, adj, manifold, k=16, leaf_size=8, band_frac=1.1, max_depth=32, t_factor=1.0, obstacles=()):
        self.heat = HeatMethod(manifold, points, t_factor, obstacles)
        self.heat_cache = {}
        super().__init__(points, adj, k, leaf_size, band_frac, max_depth)

    def portal_distances(self, src):
        if src not in self.heat_cache:
            self.heat_cache[src] = self.heat.distances(src)
        return self.heat_cache[src]

    def query_path(self, s, t):
        raise NotImplementedError("heat labels carry no predecessors, use SeparatorIndex for paths")


class SeparatorEikonalLeafIndex(SeparatorIndex):
    #the leaf table is built by fast marching on the leaf bounding box instead of by dijkstra
    #on the leaf subgraph, so same leaf distances are continuum ones and carry no stencil error.
    #a leaf covers a small box, so the grid is fine there for the same resolution.
    #NOTE: same leaf pairs are short, and the fmm error is O(h) in absolute terms, so the
    #relative error is large unless leaf_res is high. accurate from ~121 but very slow
    def __init__(self, points, adj, manifold, obstacles=(), leaf_res=181, pad_frac=0.15,
                 k=16, leaf_size=8, band_frac=1.1, max_depth=32, halo=1):
        self.manifold = manifold
        self.obstacles = list(obstacles)
        self.leaf_res = leaf_res
        self.pad_frac = pad_frac
        super().__init__(points, adj, k, leaf_size, band_frac, max_depth, True, halo)

    def build_leaf_tables(self):
        self.in_leaf = np.full(self.n, -1, dtype=int)
        self.leaf_table = []
        for leaf_id, members in enumerate(self.leaves):
            for i, g in enumerate(members):
                self.in_leaf[g] = i
            region = sorted(self.leaf_region(leaf_id))
            q = self.points[region]
            pad = self.pad_frac * max(q[:, 0].max() - q[:, 0].min(), q[:, 1].max() - q[:, 1].min())
            box = (q[:, 0].min() - pad, q[:, 0].max() + pad, q[:, 1].min() - pad, q[:, 1].max() + pad)
            D = np.full((len(members), len(members)), np.inf)
            for i, s in enumerate(members):
                xs, ys, u = solve_eikonal(self.manifold, box, self.leaf_res, self.leaf_res,
                                          [self.points[s]], self.obstacles)
                D[i] = sample_points(xs, ys, u, self.points[members])
                D[i, i] = 0.0
            self.leaf_table.append(np.minimum(D, D.T)) #symmetrize, each row is its own solve


class SeparatorEikonalIndex(SeparatorEikonalLeafIndex):
    def __init__(self, points, adj, manifold, bounds, res=121, leaf_res=None, pad_frac=0.15,
                 k=16, leaf_size=8, band_frac=1.1, max_depth=32, obstacles=(), halo=1):
        self.bounds = bounds
        self.res = res
        self.eik_cache = {}
        self.coords = np.array([np.asarray(p, dtype=float) for p in points])
        super().__init__(points, adj, manifold, obstacles, leaf_res if leaf_res else res,
                         pad_frac, k, leaf_size, band_frac, max_depth, halo)

    def portal_distances(self, src):
        if src not in self.eik_cache:
            xs, ys, u = solve_eikonal(self.manifold, self.bounds, self.res, self.res, [self.coords[src]], self.obstacles)
            self.eik_cache[src] = sample_points(xs, ys, u, self.coords)
        return self.eik_cache[src]

    def query_path(self, s, t):
        raise NotImplementedError("eikonal labels carry no predecessors, use SeparatorIndex for paths")


class SeparatorContinuousIndex(SeparatorIndex):
    def __init__(self, points, adj, k=16, leaf_size=8, band_frac=1.1, max_depth=32):
        super().__init__(points, adj, k, leaf_size, band_frac, max_depth)
        self.cross = []
        for portals in self.sep_portals:
            D = np.empty((len(portals), len(portals)))
            for i, p in enumerate(portals):
                D[i] = self.portal_distances(int(p))[portals]
            self.cross.append(D)

    def query(self, s, t):
        if s == t:
            return 0.0
        sep = self.lcs(s, t)
        if sep is None:
            return float("inf")
        a = self.label[s][sep]
        b = self.label[t][sep]
        cross = float((a[:, None] + self.cross[sep] + b[None, :]).min())
        return min(cross, self.leaf_bound(s, t))

    def query_path(self, s, t):
        raise NotImplementedError("the cut-crossing walk is not reconstructed, use SeparatorIndex for paths")


class SeparatorBFSRefinedIndex(SeparatorRefinedIndex):
    #a coordinate cut is not a separator of the GRAPH: the line x=sv is crossed somewhere by
    #the intermediate value theorem, but possibly outside the region, so a path can loop round
    #the end of the cut. a balanced BFS frontier is a true vertex separator, so deleting it
    #disconnects the two sides and every path must pass through it. no escape by construction
    def frontier_bfs(self, root, member_set):
        level = {root: 0}
        order = [root]
        dq = [root]
        head = 0
        while head < len(dq):
            v = dq[head]
            head += 1
            lv = level[v] + 1
            for u in self.adj[v]:
                if u in member_set and u not in level:
                    level[u] = lv
                    order.append(u)
                    dq.append(u)
        return level, order

    def build(self, points, depth, bbox, parent=-1):
        if len(points) <= self.leaf_size or depth >= self.max_depth:
            leaf_id = len(self.leaves)
            self.leaves.append(points)
            for s in points:
                self.leaf_of[s] = leaf_id
            return
        member_set = set(points)
        _, order0 = self.frontier_bfs(points[0], member_set)
        level, _ = self.frontier_bfs(order0[-1], member_set) #peripheral root by double bfs
        reached = list(level)
        others = [v for v in points if v not in level]
        maxlvl = max(level.values())

        best_i, best_score = None, None
        for i in range(1, maxlvl):
            below = sum(1 for v in reached if level[v] < i)
            sep = sum(1 for v in reached if level[v] == i)
            above = len(reached) - below - sep + len(others)
            if below == 0 or above == 0:
                continue
            score = (abs(below - above), sep)
            if best_score is None or score < best_score:
                best_score, best_i = score, i
        if best_i is None:
            leaf_id = len(self.leaves)
            self.leaves.append(points)
            for s in points:
                self.leaf_of[s] = leaf_id
            return

        S = [v for v in reached if level[v] == best_i]
        left = [v for v in reached if level[v] < best_i] + S #keep S in a child so every
        right = [v for v in reached if level[v] > best_i] + others #vertex still has a leaf

        portals = S if len(S) <= self.k else [S[i] for i in
                  np.unique(np.linspace(0, len(S) - 1, self.k).round().astype(int))]
        sep_id = len(self.sep_portals)
        self.sep_portals.append(np.array(portals, dtype=int))
        self.sep_meta.append((0, 0.0, depth)) #no axis or split value for a graph frontier
        self.sep_bbox.append(bbox)
        self.sep_points.append(points)
        self.sep_parent.append(parent)
        for s in points:
            self.anc[s].append(sep_id)
        self.build(left, depth + 1, bbox, sep_id)
        self.build(right, depth + 1, bbox, sep_id)
