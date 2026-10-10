"""Separator-tree distance oracles on a geometric graph.

`SeparatorIndex` is the common skeleton: recursive median coordinate cuts, portals chosen from
graph vertices in a band around each cut, labels = graph distances from every region member to
the portals of every cut above it, O(1) lowest common separator, exact leaf tables.
`SeparatorRefinedIndex` keeps the same portals but stores the *refined* (continuum) length of
each graph path instead of its graph length. The cut-portal index that the project is about
(`SeparatorCutIndex`) is in `cut_index.py`; recorded negative results are in `variants.py`.
"""
import heapq
import numpy as np
from scipy.spatial import cKDTree
from ..graphs.dijkstra import dijkstra, get_path
from ..refine import refine_path


class SeparatorIndex:
    def __init__(self, points, adj, k=16, leaf_size=8, band_frac=1.1, max_depth=32, leaf_table=True, halo=1):
        self.points = np.array([np.asarray(p, dtype=float) for p in points])
        self.n = len(points)
        self.adj = adj
        self.k = k #portals per separator
        self.leaf_size = leaf_size #min points in a region
        self.max_depth = max_depth
        self.use_leaf_table = leaf_table #cuts never separate a same-leaf pair, so store them exactly
        self.halo = halo
        self.band = band_frac * self.median_spacing()

        self.sep_points = [] #points in separator region
        self.sep_portals = [] #portals in separator regions
        self.sep_meta = [] #separator data: (axis, split_value, depth)
        self.sep_bbox = [] #data for drawing
        self.sep_parent = [] #parent separator id (-1 for root), for O(1) lcs
        self.leaf_of = np.full(self.n, -1)
        self.anc = [[] for _ in range(self.n)]
        self.label = [dict() for _ in range(self.n)]
        self.leaves = []

        bound_box = (float(self.points[:, 0].min()), float(self.points[:, 0].max()), float(self.points[:, 1].min()), float(self.points[:, 1].max()))
        self.build(list(range(self.n)), 0, bound_box)
        self.rep = [self.anc[s][-1] if self.anc[s] else -1 for s in range(self.n)]
        self.build_sep_lca()
        self.cache = {}
        self.build_labels()
        self.build_leaf_tables()
        self.num_seps = len(self.sep_portals)
        self.max_portals = max((len(p) for p in self.sep_portals), default=0)

    def median_spacing(self):
        if self.n < 2:
            return 1.0
        d, idxs = cKDTree(self.points).query(self.points, k=2)
        return float(np.median(d[:, 1]))

    def subsample(self, candidates, other_axis):
        if len(candidates) <= self.k:
            return list(candidates)
        order = sorted(candidates, key=lambda s: self.points[s, other_axis])
        idx = np.unique(np.linspace(0, len(order) - 1, self.k).round().astype(int))
        return [order[i] for i in idx]

    def build(self, points, depth, bbox, parent=-1):
        x0, x1, y0, y1 = bbox
        if len(points) <= self.leaf_size or depth >= self.max_depth:
            leaf_id = len(self.leaves)
            self.leaves.append(points)
            for s in points:
                self.leaf_of[s] = leaf_id
            return

        coords = self.points[points]
        axis = int(np.argmax(coords.max(0) - coords.min(0)))
        sv = float(np.median(coords[:, axis]))
        left = [s for s in points if self.points[s, axis] < sv]
        right = [s for s in points if self.points[s, axis] >= sv]

        if not left or not right:
            ms = sorted(points, key=lambda s: self.points[s, axis])
            mid = len(ms) // 2
            left = ms[:mid]
            right = ms[mid:]
            sv = float(self.points[ms[mid], axis])

        candidates = [s for s in points if abs(self.points[s, axis] - sv) <= self.band]
        if not candidates:
            candidates = sorted(points, key=lambda s: abs(self.points[s, axis] - sv))[:self.k]
        portals = self.subsample(candidates, 1 - axis)

        sep_id = len(self.sep_portals)
        self.sep_portals.append(np.array(portals, dtype=int))
        self.sep_meta.append((axis, sv, depth))
        self.sep_bbox.append(bbox)
        self.sep_points.append(points)
        self.sep_parent.append(parent)
        for s in points:
            self.anc[s].append(sep_id)

        if axis == 0:
            self.build(left, depth + 1, (x0, sv, y0, y1), sep_id)
            self.build(right, depth + 1, (sv, x1, y0, y1), sep_id)
        else:
            self.build(left, depth + 1, (x0, x1, y0, sv), sep_id)
            self.build(right, depth + 1, (x0, x1, sv, y1), sep_id)

    def build_sep_lca(self):
        m = len(self.sep_parent)
        self.sep_depth = np.zeros(m, dtype=np.int64)
        children = [[] for _ in range(m)]
        for v in range(m):
            p = self.sep_parent[v]
            if p >= 0:
                children[p].append(v)
                self.sep_depth[v] = self.sep_depth[p] + 1

        if m == 0:
            self.sep_euler = []
            self.sep_first = []
            self.sep_dep = []
            self.sep_sp = [[]]
            return

        euler = [0]
        first = np.full(m, -1, dtype=np.int64)
        first[0] = 0
        stack = [(0, iter(children[0]))]
        
        while stack:
            v, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                stack.pop()
                if stack:
                    euler.append(stack[-1][0])
            else:
                if first[nxt] < 0:
                    first[nxt] = len(euler)
                euler.append(nxt)
                stack.append((nxt, iter(children[nxt])))
        euler = np.array(euler, dtype=np.int64)
        dep = self.sep_depth[euler].astype(np.int64)

        e = len(euler)
        j = 1
        sp = np.zeros((e.bit_length(), e), dtype=np.int64)
        sp[0] = np.arange(e)
        while (1 << j) <= e:
            span, half = 1 << j, 1 << (j - 1)
            l = sp[j - 1, :e - span + 1]
            r = sp[j - 1, half:e - span + 1 + half]
            sp[j, :e - span + 1] = np.where(dep[l] <= dep[r], l, r)
            j += 1

        self.sep_euler = euler.tolist()
        self.sep_first = first.tolist()
        self.sep_dep = dep.tolist()
        self.sep_sp = sp.tolist()

    def sep_lca(self, u, v):
        l, r = self.sep_first[u], self.sep_first[v]
        if l > r:
            l, r = r, l
        j = (r - l + 1).bit_length() - 1
        a = self.sep_sp[j][l]
        b = self.sep_sp[j][r - (1 << j) + 1]
        return self.sep_euler[a if self.sep_dep[a] <= self.sep_dep[b] else b]

    def dijkstra(self, src):
        if src not in self.cache:
            self.cache[src] = dijkstra(self.adj, src)
        return self.cache[src]

    def portal_distances(self, src):
        return self.dijkstra(src)[0]

    def build_labels(self):
        for sep_id, points in enumerate(self.sep_points):
            portals = self.sep_portals[sep_id]
            dmat = np.empty((len(points), len(portals)))
            for j, p in enumerate(portals):
                dist = self.portal_distances(int(p))
                dmat[:, j] = dist[points]
            for i, s in enumerate(points):
                self.label[s][sep_id] = dmat[i]

    def bounded_dijkstra(self, src, limit):
        #dijkstra over the WHOLE graph but cut off at a distance we already know is an upper
        #bound on every target of interest. no region assumption, and provably exact for any
        #t with d(src,t) <= limit, because the shortest path never leaves that ball
        dist = {src: 0.0}
        pred = {src: -1}
        pq = [(0.0, src)]
        while pq:
            d, v = heapq.heappop(pq)
            if d > dist[v]:
                continue
            if d > limit:
                break
            for u, w in self.adj[v].items():
                nd = d + w
                if nd <= limit and nd < dist.get(u, np.inf):
                    dist[u] = nd
                    pred[u] = v
                    heapq.heappush(pq, (nd, u))
        return dist, pred

    def local_dijkstra(self, src, allowed):
        dist = {src: 0.0}
        pred = {src: -1}
        pq = [(0.0, src)]
        while pq:
            d, v = heapq.heappop(pq)
            if d > dist[v]:
                continue
            for u, w in self.adj[v].items():
                if u in allowed and d + w < dist.get(u, np.inf):
                    dist[u] = d + w
                    pred[u] = v
                    heapq.heappush(pq, (d + w, u))
        return dist, pred

    def leaf_region(self, leaf_id):
        allowed = set(self.leaves[leaf_id])
        for _ in range(self.halo):
            allowed = allowed | {u for v in allowed for u in self.adj[v]}
        return allowed

    def build_leaf_tables(self):
        self.in_leaf = np.full(self.n, -1, dtype=int)
        self.leaf_table = []
        for leaf_id, members in enumerate(self.leaves):
            for i, g in enumerate(members):
                self.in_leaf[g] = i
            if not self.use_leaf_table:
                self.leaf_table.append(None)
                continue
            allowed = self.leaf_region(leaf_id)
            D = np.full((len(members), len(members)), np.inf)
            for i, s in enumerate(members):
                #pass 1 gives a valid upper bound, pass 2 uses it as a certified cutoff so the
                #result no longer depends on the halo being large enough
                first, _ = self.local_dijkstra(s, allowed)
                reach = [first[g] for g in members if g in first]
                limit = max(reach) if reach else 0.0
                dist, pred = self.bounded_dijkstra(s, limit)
                for j, g in enumerate(members):
                    if g in dist:
                        D[i, j] = dist[g]
            self.leaf_table.append(D)

    def leaf_bound(self, s, t):
        if not self.use_leaf_table:
            return float("inf")
        leaf_id = self.leaf_of[s]
        if leaf_id < 0 or leaf_id != self.leaf_of[t]:
            return float("inf")
        return float(self.leaf_table[leaf_id][self.in_leaf[s], self.in_leaf[t]])

    def lcs(self, s, t):
        u, v = self.rep[s], self.rep[t]
        if u < 0 or v < 0:
            return None
        return self.sep_lca(u, v)

    def portal_bound(self, s, t):
        sep = self.lcs(s, t)
        if sep is None:
            return float("inf")
        dists = self.label[s][sep] + self.label[t][sep]
        return float(dists.min())

    def query(self, s, t):
        if s == t:
            return 0.0
        return min(self.portal_bound(s, t), self.leaf_bound(s, t))

    def query_path(self, s, t):
        if s == t:
            return 0.0, [s]
        sep = self.lcs(s, t)
        if sep is None:
            return float("inf"), None

        if self.leaf_bound(s, t) < self.portal_bound(s, t):
            dist, pred = self.local_dijkstra(s, self.leaf_region(self.leaf_of[s]))
            if t in dist:
                path = [t]
                while path[-1] != s:
                    path.append(pred[path[-1]])
                return float(dist[t]), path[::-1]

        portals = self.sep_portals[sep]
        dists = self.label[s][sep] + self.label[t][sep]
        i = int(dists.argmin())
        if not np.isfinite(dists[i]):
            return float("inf"), None

        p = int(portals[i])
        dist, pred = self.dijkstra(p)
        l = get_path(pred, p, s)
        r = get_path(pred, p, t)

        if l is None or r is None:
            return float(dists[i]), None
        return float(dists[i]), l[::-1] + r[1:]


class SeparatorRefinedIndex(SeparatorIndex):
    def __init__(self, points, adj, manifold, obstacles=(), k=16, leaf_size=8, band_frac=1.1,
                 max_depth=32, leaf_table=True, halo=1, rounds=1, use_relax=True, segments=12):
        self.manifold = manifold
        self.obstacles = list(obstacles)
        self.rounds = rounds
        self.use_relax = use_relax
        self.segments = segments
        super().__init__(points, adj, k, leaf_size, band_frac, max_depth, leaf_table, halo)
        self.refine_labels()

    def build_leaf_tables(self):
        #the leaf table must be refined too, otherwise same leaf pairs are answered from raw
        #graph distances while every portal label is a continuum one, and the min() mixes
        #two different quantities. the paths come free from the local dijkstra predecessors
        self.in_leaf = np.full(self.n, -1, dtype=int)
        self.leaf_table = []
        for leaf_id, members in enumerate(self.leaves):
            for i, g in enumerate(members):
                self.in_leaf[g] = i
            if not self.use_leaf_table:
                self.leaf_table.append(None)
                continue
            allowed = self.leaf_region(leaf_id)
            D = np.full((len(members), len(members)), np.inf)
            for i, s in enumerate(members):
                #pass 1 gives a valid upper bound, pass 2 uses it as a certified cutoff so the
                #result no longer depends on the halo being large enough
                first, _ = self.local_dijkstra(s, allowed)
                reach = [first[g] for g in members if g in first]
                limit = max(reach) if reach else 0.0
                dist, pred = self.bounded_dijkstra(s, limit)
                for j, g in enumerate(members):
                    if g not in dist:
                        continue
                    D[i, j] = dist[g]
                    path = [g]
                    while path[-1] != s:
                        path.append(pred[path[-1]])
                    if len(path) < 3:
                        continue
                    d, coords = refine_path(self.manifold, self.points, path[::-1], self.obstacles,
                                            rounds=self.rounds, use_relax=self.use_relax, segments=self.segments)
                    if d < D[i, j]:
                        D[i, j] = d
            self.leaf_table.append(np.minimum(D, D.T)) #each direction is refined separately

    def refine_labels(self):
        self.refined = 0
        for sep_id, members in enumerate(self.sep_points):
            portals = self.sep_portals[sep_id]
            for j, p in enumerate(portals):
                dist, pred = self.dijkstra(int(p))
                for s in members:
                    if s == int(p) or not np.isfinite(dist[s]):
                        continue
                    path = get_path(pred, int(p), s)
                    if path is None or len(path) < 3:
                        continue
                    d, coords = refine_path(self.manifold, self.points, path, self.obstacles,
                                            rounds=self.rounds, use_relax=self.use_relax, segments=self.segments)
                    if d < self.label[s][sep_id][j]:
                        self.label[s][sep_id][j] = d
                        self.refined += 1
