import numpy as np
from scipy.spatial import cKDTree
from dijkstra import dijkstra, get_path

class SeparatorIndex:
    def __init__(self, points, adj, k=16, leaf_size=8, band_frac=1.1, max_depth=32):
        self.points = np.array([np.asarray(p, dtype=float) for p in points])
        self.n = len(points)
        self.adj = adj
        self.k = k #portals per separator
        self.leaf_size = leaf_size #min points in a region
        self.max_depth = max_depth
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
        self.rep = [self.anc[s][-1] if self.anc[s] else -1 for s in range(self.n)] #deepest separator per point
        self.build_sep_lca()
        self.cache = {}
        self.build_labels()
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
        idx = np.unique(np.linspace(0, len(order) - 1, self.k).round().astype(int)) #even coverage of the cut line
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

    def build_labels(self):
        for sep_id, points in enumerate(self.sep_points):
            portals = self.sep_portals[sep_id]
            dmat = np.empty((len(points), len(portals)))
            for j, p in enumerate(portals):
                dist, pred = self.dijkstra(int(p))
                dmat[:, j] = dist[points]
            for i, s in enumerate(points):
                self.label[s][sep_id] = dmat[i]

    def lcs(self, s, t):
        u, v = self.rep[s], self.rep[t]
        if u < 0 or v < 0:
            return None
        return self.sep_lca(u, v)

    def query(self, s, t):
        if s == t:
            return 0.0
        sep = self.lcs(s, t)
        if sep is None:
            return float("inf")
        dists = self.label[s][sep] + self.label[t][sep]
        return float(dists.min())

    def query_path(self, s, t):
        if s == t:
            return 0.0, [s]
        sep = self.lcs(s, t)
        if sep is None:
            return float("inf"), None

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
