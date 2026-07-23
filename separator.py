import numpy as np
from scipy.spatial import cKDTree
from dijkstra import dijkstra, get_path

class SeparatorIndex:
    def __init__(self, points, adj, k=16, leaf_size=8, band_frac=1.6, max_depth=32):
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
        self.leaf_of = np.full(self.n, -1)
        self.anc = [[] for _ in range(self.n)]
        self.label = [dict() for _ in range(self.n)]
        self.leaves = []
        
        bound_box = (float(self.points[:, 0].min()), float(self.points[:, 0].max()), float(self.points[:, 1].min()), float(self.points[:, 1].max()))
        self.build(list(range(self.n)), 0, bound_box)
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
    
    def build(self, points, depth, bbox):
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
        for s in points:
            self.anc[s].append(sep_id)
        
        if axis == 0:
            self.build(left, depth + 1, (x0, sv, y0, y1))
            self.build(right, depth + 1, (sv, x1, y0, y1))
        else:
            self.build(left, depth + 1, (x0, x1, y0, sv))
            self.build(right, depth + 1, (x0, x1, sv, y1))
            
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
        a = self.anc[s]
        b = self.anc[t]
        l = min(len(a), len(b))
        m = 0
        
        while m < l and a[m] == b[m]:
            m += 1
            
        if m > 0:
            return a[m - 1]
        else:
            return None

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
    
    
        