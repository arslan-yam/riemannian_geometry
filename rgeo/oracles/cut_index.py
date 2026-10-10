"""The Separator Cut Index: portals placed ON the cut lines, continuum labels, interpolation.

See the class comment of `SeparatorCutIndex` and notes/error_bound.md for the analysis. The
module-level worker functions exist so that `multiprocessing` (fork) can call bound methods of
the index object inherited by the children.
"""
import numpy as np
from scipy.spatial import cKDTree
from ..graphs.dijkstra import get_path
from ..obstacles import blocked
from ..refine import refine_polyline, refine_batch, segment_lengths
from .separator import SeparatorIndex

#numpy 2.0 renamed trapz to trapezoid, keep working on both
trapezoid = getattr(np, "trapezoid", None) or np.trapz

_WORKER_INDEX = None


def _portal_worker(task):
    sep_id, j = task
    return sep_id, j, _WORKER_INDEX._label_portal(sep_id, j)


def _leaf_worker(leaf_id):
    return leaf_id, _WORKER_INDEX._leaf_table(leaf_id)


def _exact_worker(sep_id):
    return sep_id, _WORKER_INDEX._label_separator_exact(sep_id)


class SeparatorCutIndex(SeparatorIndex):
    #portals placed exactly ON the cut line instead of at nearby graph vertices, and the
    #crossing point interpolated rather than snapped to a sample.
    #
    #E_portal is 1d minimisation error: f(x) = d(s,x) + d(x,t) along the cut is smooth with
    #a nondegenerate minimum, so a discrete min over k samples costs O(h^2) = O(1/k^2).
    #Fitting a parabola through the argmin and its neighbours reports the vertex instead,
    #which is O(h^3). This needs (a) f smooth, which refined continuum labels give and raw
    #graph labels do not, and (b) f parameterised by one variable, which requires the
    #portals to lie on a line rather than scattered through a 2d band.
    #Labels need no graph vertex: d(s,x) for an arbitrary cut point x is the refined length
    #of the graph path to the nearest node, extended to x.
    def _leaf_table(self, leaf_id):
        #continuum leaf table of one leaf: graph paths between all member pairs from the local
        #dijkstras, then ONE batched refinement for the whole leaf. each unordered pair is
        #refined once (the energy is symmetric, so the reversed polyline gives the same curve)
        members = self.leaves[leaf_id]
        m = len(members)
        D = np.full((m, m), np.inf)
        if m == 0:
            return D
        if self.exact is not None:
            #exact distances among the members and, from the same sweeps, to every vertex
            #within halo_rings hops of the leaf: that is what removes the close-pairs-across-a-
            #cut outliers (see notes). the solver lives on the sub-mesh the geodesics can visit
            near = set(int(g) for g in members)
            for _ in range(self.halo_rings):
                near |= {u for v in near for u in self.adj[v]}
            near -= set(int(g) for g in members)
            extra = np.array(sorted(near), dtype=int)
            P = self.points[np.concatenate([members, extra])] if len(extra) else self.points[members]
            box = (P[:, 0].min(), P[:, 0].max(), P[:, 1].min(), P[:, 1].max())
            geo, remap, _ = self.exact.region_geo(box)
            loc = remap[np.asarray(members, dtype=int)]
            tgt = np.concatenate([loc, remap[extra]]) if len(extra) else loc
            H = np.full((m, len(extra)), np.inf)
            for i, vi in enumerate(loc):
                d, _ = geo.geodesicDistances(np.array([int(vi)]), tgt)
                D[i] = d[:m]
                H[i] = d[m:]
            D = np.minimum(D, D.T)
            return (D, extra, H) if self.halo_rings else D
        allowed = self.leaf_region(leaf_id)
        pairs, lines = [], []
        for i, s in enumerate(members):
            first, _ = self.local_dijkstra(s, allowed)
            reach = [first[g] for g in members if g in first]
            limit = max(reach) if reach else 0.0
            dist, pred = self.bounded_dijkstra(s, limit)
            for j, g in enumerate(members):
                if g not in dist:
                    continue
                D[i, j] = dist[g]
                if j <= i:
                    continue
                path = [g]
                while path[-1] != s:
                    path.append(pred[path[-1]])
                if len(path) < 3:
                    continue
                pairs.append((i, j))
                lines.append(self.points[path[::-1]])
        if lines:
            res = refine_batch(self.manifold, lines, self.obstacles,
                               rounds=self.rounds, use_relax=self.use_relax, segments=self.segments)
            for (i, j), (d, _) in zip(pairs, res):
                if d < D[i, j]:
                    D[i, j] = d
        return np.minimum(D, D.T)

    def build_leaf_tables(self):
        #same tables as SeparatorRefinedIndex.build_leaf_tables, but batched per leaf and
        #leaves fanned out across workers: the leaf tables were a third of the build
        import time
        t0 = time.perf_counter()
        self.timing = getattr(self, "timing", {})
        self.in_leaf = np.full(self.n, -1, dtype=int)
        for members in self.leaves:
            for i, g in enumerate(members):
                self.in_leaf[g] = i
        self.leaf_table = [None] * len(self.leaves)
        self.halo_ids = [None] * self.n #sorted vertex ids with an exact local distance
        self.halo_d = [None] * self.n
        if not self.use_leaf_table:
            return

        def store(leaf_id, res):
            if isinstance(res, tuple):
                D, extra, H = res
                for i, s in enumerate(self.leaves[leaf_id]):
                    self.halo_ids[s] = extra
                    self.halo_d[s] = H[i].astype(self.label_dtype)
            else:
                D = res
            self.leaf_table[leaf_id] = D

        if self.workers > 1 and len(self.leaves) > 1:
            import multiprocessing as mp
            global _WORKER_INDEX
            _WORKER_INDEX = self
            order = sorted(range(len(self.leaves)), key=lambda i: -len(self.leaves[i]))
            with mp.get_context("fork").Pool(self.workers) as pool:
                for leaf_id, res in pool.imap_unordered(_leaf_worker, order, chunksize=1):
                    store(leaf_id, res)
            _WORKER_INDEX = None
        else:
            for leaf_id in range(len(self.leaves)):
                store(leaf_id, self._leaf_table(leaf_id))
        self.timing["leaf_tables"] = time.perf_counter() - t0

    def build_labels(self):
        #graph portal labels (one dijkstra per vertex portal) are only a safety net for the
        #refined labels; with exact labels they never win, so skip them and their memory
        if self.exact is None:
            super().build_labels()

    def portal_bound(self, s, t):
        if self.exact is not None:
            return float("inf")
        return super().portal_bound(s, t)

    def halo_bound(self, s, t):
        ids = self.halo_ids[s]
        if ids is None or len(ids) == 0:
            return float("inf")
        j = int(np.searchsorted(ids, t))
        if j < len(ids) and ids[j] == t:
            return float(self.halo_d[s][j])
        return float("inf")

    def __init__(self, points, adj, manifold, obstacles=(), k=16, leaf_size=8, band_frac=1.1,
                 max_depth=32, halo=1, rounds=2, use_relax=True, segments=20, interpolate=True,
                 alpha=0.0, scale_k=True, k_min=4, k_max=None, pad="diam", k_pad=1, workers=1, batched=True,
                 pad_speed=None, exact=None, exact_rings=3, halo_rings=2, label_dtype=np.float32,
                 edge_portals=0):
        self.manifold = manifold
        self.obstacles = list(obstacles)
        self.rounds = rounds
        self.use_relax = use_relax
        self.segments = segments
        self.interpolate = interpolate
        self.alpha = alpha #portal density ~ c(x)^-alpha, c = tangential metric speed on the cut.
        #alpha=-1 is riemannian arc length, alpha=0 is coordinate uniform, alpha>0 puts portals
        #where the cut is cheap to traverse, which is where geodesics actually tend to cross
        self.scale_k = scale_k #portals proportional to cut length, so density is uniform across depths
        self.k_min = k_min
        self.k_max = k_max if k_max else 2 * k
        #a path may cross the cut LINE outside the region (intermediate value theorem gives
        #x=sv somewhere, not necessarily inside sep_bbox). "diam" extends by half the region
        #diameter, which is provably enough: a crossing D outside costs >= 2D extra, and the
        #optimum is at most diam(region). "domain" spans the whole line.
        #measured on a wall bed: diam gained 1.2% and domain 2.4% for 4.3x and 5.5x the build,
        #and the max error was identical, so escapes are rare. default off, enable if the
        #geometry makes looping round a cut end plausible
        self.pad = pad
        self.k_pad = k_pad #portals per padded tail, on top of the region budget
        self.pad_speed = pad_speed #lambda_min of the metric: makes the "diam" window a proof, not a heuristic
        #exact labels when the manifold has an exact point-to-all solver (terrain.ExactTerrain):
        #then no graph path is refined at all, labels and leaf tables are exact geodesics
        self.exact = exact
        self.exact_rings = exact_rings
        self.halo_rings = halo_rings #exact distances to all vertices this many hops away (exact only)
        self.label_dtype = label_dtype #float32 halves the index; 1e-7 relative is far below the portal error
        #extra portals where mesh edges cross the cut (manifold.edge_crossings), at most this many
        #per cut: on weighted terrains geodesics hug cheap edges and cross a cut where those do
        self.edge_portals = edge_portals
        self.model_reach = 2.0 #model interpolation when an endpoint is within this many 3-portal windows of the cut
        self.model_fit = "parabola" #how the smooth remainder is fitted: "parabola" or "linear"
        self.workers = workers #separators are independent, so their labels build in parallel
        self.batched = batched #relax every member's polyline to a portal in one optimiser call
        super().__init__(points, adj, k, leaf_size, band_frac, max_depth, True, halo)
        self.dom = (self.points[:, 0].min(), self.points[:, 0].max(),
                    self.points[:, 1].min(), self.points[:, 1].max())
        self.tree = cKDTree(self.points)
        self.build_cut_labels()

    def cut_arclength(self, sep_id, samples=200):
        #parameterise the cut by riemannian arc length. the same coordinate spacing can be
        #wildly non uniform in true distance, on the poincare disk by an order of magnitude
        axis, sv, depth = self.sep_meta[sep_id]
        #the separator guarantees that a path inside the REGION crosses the cut, so the cut
        #must be sampled across the whole region. the member bounding box can be strictly
        #smaller when obstacles or in_domain carve the region up, and then a legitimate
        #crossing corridor gets no portals at all
        x0, x1, y0, y1 = self.sep_bbox[sep_id]
        rlo, rhi = (y0, y1) if axis == 0 else (x0, x1)
        dom_lo = self.dom[2] if axis == 0 else self.dom[0]
        dom_hi = self.dom[3] if axis == 0 else self.dom[1]

        #the padding proof (D <= 0.5*diam(region)) bounds a RIEMANNIAN distance, so turning it
        #into a sampling window needs a riemannian arc length, not a coordinate offset: on the
        #poincare disk the metric speed is >=2 everywhere and >10x near the boundary, so adding
        #a riemannian length directly to a coordinate wildly overshoots (or, if clipped to the
        #domain, undershoots) the intended window. probe the metric speed once across the whole
        #domain extent on this axis and invert its cumulative arc length instead
        probe = np.linspace(dom_lo, dom_hi, max(samples, 400))
        pline = np.zeros((len(probe), 2))
        pline[:, axis] = sv
        pline[:, 1 - axis] = probe
        pspeed = np.sqrt(np.maximum(self.manifold.metric_batch(pline)[:, 1 - axis, 1 - axis], 1e-300))
        raw_arc = np.concatenate([[0.0], np.cumsum(0.5 * (pspeed[1:] + pspeed[:-1]) * np.diff(probe))])
        arc_rlo = float(np.interp(rlo, probe, raw_arc))
        arc_rhi = float(np.interp(rhi, probe, raw_arc))

        if self.pad in ("diam", "domain"):
            #riemannian length of the bbox diagonal, via the same quadrature used everywhere
            #else in the codebase for continuum lengths -- NOT np.hypot, which is coordinate-only
            margin = 0.5 * self.manifold.segment_lenght(np.array([x0, y0]), np.array([x1, y1]))
        else:
            margin = 1e-12
        if self.pad == "diam" and self.pad_speed is not None:
            #provable window: a crossing at coordinate offset c from the region costs at least
            #2*c*lambda_min extra, so c <= margin/lambda_min. exact for terrains (lambda_min = 1)
            lo = rlo - margin / self.pad_speed
            hi = rhi + margin / self.pad_speed
        elif self.pad == "diam":
            #invert raw_arc: find the coordinate where the accumulated riemannian arc length
            #reaches arc_rlo-margin / arc_rhi+margin, clipping to the domain edge if the domain
            #itself runs out before the margin does (falls back to the "domain" behaviour there).
            #this uses the speed ALONG the cut as the local lambda_min, which is only a heuristic:
            #a detour off the cut can be cheaper than the cut itself. pad_speed gives the proof
            lo = float(np.interp(arc_rlo - margin, raw_arc, probe))
            hi = float(np.interp(arc_rhi + margin, raw_arc, probe))
        elif self.pad == "domain":
            lo, hi = dom_lo, dom_hi
        else:
            lo, hi = rlo, rhi
        lo = max(lo, dom_lo) #never leave the point cloud
        hi = min(hi, dom_hi)

        dense = np.linspace(lo, hi, samples)
        line = np.zeros((samples, 2))
        line[:, axis] = sv
        line[:, 1 - axis] = dense
        g = self.manifold.metric_batch(line)
        speed = np.sqrt(np.maximum(g[:, 1 - axis, 1 - axis], 1e-300)) #|t|_g, works for any metric
        #a crossing D outside the region needs d(s,t) >= 2D, so only a shrinking fraction of
        #pairs can use it. weight the padded tail by that fraction instead of giving it the
        #same density as the region, which would spend most of the budget where crossings are rare.
        #D itself must be riemannian too, so reuse the probe's arc length rather than coordinate out
        raw = np.interp(dense, probe, raw_arc)
        out = np.maximum(np.maximum(arc_rlo - raw, raw - arc_rhi), 0.0)
        feas = np.clip(1.0 - out / max(margin, 1e-12), 0.0, 1.0)
        w = (speed ** (-self.alpha)) * feas
        arc = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(dense))])
        #the budget must be set by the REGION span only. padding is a soundness margin, not a
        #place crossings concentrate; letting it into the length inflated the longest cut and
        #diluted every interior budget (measured: +35% error on poincare, same max)
        inside = (dense >= rlo) & (dense <= rhi)
        region_speed = np.where(inside, speed, 0.0)
        return dense, arc, region_speed, inside

    def build_cut_labels(self):
        self.cut_pos = [] #along cut coordinate of each virtual portal
        self.cut_seg = [] #free sub interval id of each portal, fitting must never cross a gap
        self.cut_xy = []
        self.cut_label = [dict() for _ in range(self.n)]

        geom = [self.cut_arclength(i) for i in range(len(self.sep_points))]
        #scale_k must use TRUE arc length, not the alpha weighted measure, or the budget
        #allocation across separators would change with alpha
        lengths = [float(trapezoid(sp, d)) if len(d) > 1 else 0.0 for d, _, sp, _ in geom]
        longest = max(lengths) if lengths else 1.0

        for sep_id, members in enumerate(self.sep_points):
            axis, sv, depth = self.sep_meta[sep_id]
            dense, arc, speed, inside = geom[sep_id]
            kk = self.k
            if self.scale_k and longest > 0:
                kk = int(np.clip(round(self.k * lengths[sep_id] / longest), self.k_min, self.k_max))

            #nothing ever crosses the cut inside an obstacle, so portals on a blocked stretch
            #are pure waste. find the free sub intervals and share the budget between them
            line = np.zeros((len(dense), 2))
            line[:, axis] = sv
            line[:, 1 - axis] = dense
            free = np.array([self.manifold.in_domain(q) and not blocked(q, self.obstacles) for q in line])
            runs = []
            i = 0
            while i < len(free):
                if free[i]:
                    j = i
                    while j + 1 < len(free) and free[j + 1]:
                        j += 1
                    if j > i:
                        runs.append((i, j))
                    i = j + 1
                else:
                    i += 1
            if not runs: #cut fully blocked inside the bbox, keep the old placement as fallback
                runs = [(0, len(dense) - 1)]

            #split every free run into its in-region part and its padded tail parts. the kk
            #budget goes to the region parts only; each tail part gets a fixed small allowance,
            #so soundness padding never dilutes the interior (measured: letting it share the
            #budget cost 35-50% accuracy on poincare for no change in the max error)
            pieces = [] #(run id, i0, i1, is_region)
            for r, (i0, i1) in enumerate(runs):
                m = inside[i0:i1 + 1]
                start, cur = i0, bool(m[0])
                for t in range(i0 + 1, i1 + 2):
                    flag = bool(m[t - i0]) if t <= i1 else None
                    if flag != cur:
                        if t - 1 > start:
                            pieces.append((r, start, t - 1, cur))
                        start, cur = t - 1, flag
            region = [pc for pc in pieces if pc[3]]
            tails = [pc for pc in pieces if not pc[3]]
            if not region: #region entirely blocked: fall back to whatever is free
                region, tails = [(r, i0, i1, True) for r, (i0, i1) in enumerate(runs)], []
            weights = [arc[i1] - arc[i0] for _, i0, i1, _ in region]
            total = sum(weights)
            pos, seg = [], []
            for (r, i0, i1, _), w in zip(region, weights):
                nr = max(1, int(round(kk * w / total))) if total > 0 else 1
                a = arc[i0:i1 + 1] - arc[i0]
                d = dense[i0:i1 + 1]
                if nr == 1 or a[-1] <= 0:
                    pos.append(float(np.interp(0.5 * max(a[-1], 0.0), a, d)))
                    seg.append(r)
                else:
                    for v in np.interp(np.linspace(0.0, a[-1], nr), a, d):
                        pos.append(float(v))
                        seg.append(r)
            for r, i0, i1, _ in tails:
                a = arc[i0:i1 + 1] - arc[i0]
                d = dense[i0:i1 + 1]
                if a[-1] <= 0:
                    continue
                #interior points of the tail only; its region-side end is already a portal
                for v in np.interp(np.linspace(0.0, a[-1], self.k_pad + 2)[1:-1], a, d):
                    pos.append(float(v))
                    seg.append(r)
            if self.edge_portals and hasattr(self.manifold, "edge_crossings"):
                cross = self.manifold.edge_crossings(axis, sv, dense[0], dense[-1])
                if len(cross) > self.edge_portals: #keep a uniform subset
                    cross = cross[np.linspace(0, len(cross) - 1, self.edge_portals).round().astype(int)]
                for v in cross:
                    #the free run this crossing falls in (all runs are free on a terrain)
                    j = int(np.searchsorted(dense, v))
                    r_of = [r for r, (i0, i1) in enumerate(runs) if i0 <= j <= i1 + 1]
                    if r_of:
                        pos.append(float(v)); seg.append(r_of[0])
            order = np.argsort(pos) #cut_bound fits a parabola over neighbours, keep them sorted
            pos = np.array(pos)[order]
            seg = [seg[i] for i in order]
            xy = np.zeros((len(pos), 2))
            xy[:, axis] = sv
            xy[:, 1 - axis] = pos
            self.cut_pos.append(pos)
            self.cut_seg.append(np.array(seg))
            self.cut_xy.append(xy)
            for s in members:
                self.cut_label[s][sep_id] = np.full(len(pos), np.inf, dtype=self.label_dtype)

        #geometry needed by the model interpolation: the cut's arc length, and for every member
        #its foot on the cut line and the riemannian length of the drop (one batched walk per cut)
        self.cut_arc = [(d, a) for d, a, _, _ in geom]
        self.foot_r = [dict() for _ in range(self.n)]
        for sep_id, members in enumerate(self.sep_points):
            if not members:
                continue
            axis, sv, depth = self.sep_meta[sep_id]
            P = self.points[np.asarray(members, dtype=int)]
            Q = P.copy(); Q[:, axis] = sv
            r = segment_lengths(self.manifold, P, Q)
            for s, rv in zip(members, r):
                self.foot_r[s][sep_id] = float(rv) if np.isfinite(rv) else float(np.hypot(*(P[0] * 0)))
        #everything above is cheap placement. the refinement below is all of the cost, and the
        #separators are independent, so fan them out across processes when asked to
        import time
        t0 = time.perf_counter()
        if self.exact is not None:
            self._label_exact()
        elif self.workers > 1 and len(self.sep_points) > 1:
            self._label_parallel()
        else:
            for sep_id in range(len(self.sep_points)):
                self._merge_labels(sep_id, self._label_separator(sep_id))
        self.timing = getattr(self, "timing", {})
        self.timing["cut_labels"] = time.perf_counter() - t0
        self.cache = {} #the dijkstra trees were construction scaffolding

    def _label_exact(self):
        #one task per cut: an exact solver on the sub-mesh that its geodesics can visit, with
        #the cut's portals inserted as steiner vertices, then one sweep per portal. tasks are
        #fanned out by fork; the top of the tree (whole-mesh sweeps) goes first
        tasks = list(range(len(self.sep_points)))
        if self.workers > 1 and len(tasks) > 1:
            import multiprocessing as mp
            global _WORKER_INDEX
            _WORKER_INDEX = self
            tasks.sort(key=lambda i: -len(self.sep_points[i]))
            with mp.get_context("fork").Pool(self.workers) as pool:
                for sep_id, cols in pool.imap_unordered(_exact_worker, tasks, chunksize=1):
                    self._merge_labels(sep_id, cols)
            _WORKER_INDEX = None
        else:
            for sep_id in tasks:
                self._merge_labels(sep_id, self._label_separator_exact(sep_id))

    def _label_separator_exact(self, sep_id):
        members = np.asarray(self.sep_points[sep_id], dtype=int)
        xy = self.cut_xy[sep_id]
        if len(members) == 0 or len(xy) == 0:
            return {}
        geo, remap, pidx = self.exact.region_geo(self.sep_bbox[sep_id], xy)
        tgt = remap[members]
        cols = {}
        for j, vi in enumerate(pidx):
            if vi < 0: #portal outside the surface
                continue
            d, _ = geo.geodesicDistances(np.array([int(vi)]), tgt)
            cols[j] = {int(s): float(v) for s, v in zip(members, d)}
        return cols

    def _label_portal(self, sep_id, j):
        #refined continuum length from every member of the region to portal j of its cut.
        #one dijkstra from the portal's anchor, then ALL members' polylines refined in one
        #batched optimiser call. pure function of the index state, so it runs in a forked worker
        members = self.sep_points[sep_id]
        x = self.cut_xy[sep_id][j]
        if not self.manifold.in_domain(x) or blocked(x, self.obstacles):
            return {}
        anchor = int(self.tree.query(x)[1]) #reach the cut point through the nearest node
        dist, pred = self.dijkstra(anchor)
        owners, lines = [], []
        for s in members:
            if not np.isfinite(dist[s]):
                continue
            path = get_path(pred, anchor, s)
            if path is None:
                continue
            line = np.vstack([x, self.points[path]]) #x -> anchor -> ... -> s
            keep = np.concatenate([[True], np.linalg.norm(np.diff(line, axis=0), axis=1) > 1e-12])
            line = line[keep] #x may coincide with its anchor node
            if len(line) < 2:
                continue
            owners.append(s)
            lines.append(line)
        if not lines:
            return {}
        if self.batched:
            res = refine_batch(self.manifold, lines, self.obstacles,
                               rounds=self.rounds, use_relax=self.use_relax, segments=self.segments)
            return {s: d for s, (d, _) in zip(owners, res)}
        out = {}
        for s, line in zip(owners, lines):
            d, c = refine_polyline(self.manifold, line, self.obstacles,
                                   rounds=self.rounds, use_relax=self.use_relax, segments=self.segments)
            out[s] = d
        return out

    def _label_separator(self, sep_id):
        cols = {}
        for j in range(len(self.cut_xy[sep_id])):
            cols[j] = self._label_portal(sep_id, j)
        return cols

    def _merge_portal(self, sep_id, j, col):
        for s, d in col.items():
            self.cut_label[s][sep_id][j] = d

    def _merge_labels(self, sep_id, cols):
        for j, col in cols.items():
            self._merge_portal(sep_id, j, col)

    def _label_parallel(self):
        import multiprocessing as mp
        #fork so the children inherit the whole index without pickling it; only small result
        #dicts come back. tasks are (separator, portal) pairs, not separators: the root
        #separator alone holds most of the refinement work, so per-separator tasks cannot
        #balance. largest first so the tail of the schedule is short
        global _WORKER_INDEX
        _WORKER_INDEX = self
        tasks = [(i, j) for i in range(len(self.sep_points)) for j in range(len(self.cut_xy[i]))]
        tasks.sort(key=lambda t: -len(self.sep_points[t[0]]))
        with mp.get_context("fork").Pool(self.workers) as pool:
            for sep_id, j, col in pool.imap_unordered(_portal_worker, tasks, chunksize=1):
                self._merge_portal(sep_id, j, col)
        _WORKER_INDEX = None

    def cut_bound(self, s, t):
        sep = self.lcs(s, t)
        if sep is None:
            return float("inf")
        return self.portal_min(sep, self.cut_label[s][sep] + self.cut_label[t][sep], s, t)

    def foot(self, sep, p):
        #(cut coordinate of the foot of p on the cut line, riemannian length of the drop).
        #p is a node id (precomputed at build) or a chart point (measured now)
        axis, sv, _ = self.sep_meta[sep]
        if isinstance(p, tuple): #(cut coordinate, drop) already known, e.g. interpolated for a face point
            return p
        if isinstance(p, (int, np.integer)):
            r = self.foot_r[p].get(sep)
            if r is not None:
                return float(self.points[p][1 - axis]), r
            p = self.points[p]
        q = np.array(p, dtype=float); q[axis] = sv
        return float(q[1 - axis]), float(self.manifold.segment_lenght(np.asarray(p, dtype=float), q))

    def portal_min(self, sep, f, p=None, q=None):
        #minimum of the sampled crossing cost f over the portals of cut sep, with the optional
        #interpolation. also used for arbitrary-point queries (terrain.TerrainA2A).
        #
        #interpolate=True fits a parabola through the argmin sample and its neighbours.
        #interpolate="model" subtracts what makes f sharp first: an endpoint at distance r from
        #the cut sees d(p, x) ~ sqrt(r^2 + (x - x_p)^2), whose curvature 1/r is huge when r ~ h,
        #and that is exactly where the parabola fails (measured: every error above 1% on
        #HorseMount had an endpoint within a few portal gaps of the cut, at any d). the
        #remainder g(x) = f(x) - m_p(x) - m_q(x) is smooth, gets the parabola, and the model
        #sum is minimised in 1-d by newton steps. reduces to the parabola when r >> h
        i = int(f.argmin())
        best = float(f[i])
        if not self.interpolate or not (0 < i < len(f) - 1):
            return best
        g = self.cut_seg[sep]
        if not (g[i - 1] == g[i] == g[i + 1]): #f jumps across an obstacle gap, no parabola there
            return best
        x = self.cut_pos[sep][i - 1:i + 2]
        y = f[i - 1:i + 2].astype(np.float64) #labels may be float32; the fit wants full precision
        if not np.isfinite(y).all():
            return best
        if self.interpolate == "model" and p is not None and q is not None:
            #the model matters only when an endpoint is within a few portal gaps of the cut;
            #otherwise both hyperbolas are straight lines over the window and the parabola is
            #the same fit for a third of the cost
            xp, rp = self.foot(sep, p)
            xq, rq = self.foot(sep, q)
            arc = self._arc(sep, x)
            if min(rp, rq) < self.model_reach * (arc[2] - arc[0]):
                return self._model_min(sep, x, y, (xp, rp), (xq, rq), best, arc)
        d1 = (y[1] - y[0]) / (x[1] - x[0]) #newton form of the parabola through the three
        d2 = (y[2] - y[1]) / (x[2] - x[1]) #samples, polyfit does the same thing 5x slower
        a = (d2 - d1) / (x[2] - x[0])
        if a <= 0:
            return best
        xv = 0.5 * (x[0] + x[1]) - d1 / (2.0 * a)
        if x[0] <= xv <= x[2]:
            fv = float(y[0] + d1 * (xv - x[0]) + a * (xv - x[0]) * (xv - x[1]))
            #for a resolved parabola the vertex sits at most a quarter of the smaller neighbour
            #gap below the argmin sample (exact in the symmetric case). a deeper vertex means
            #the three samples do not share one basin, and the fit is not trusted past that
            gain_cap = 0.25 * min(y[0] - y[1], y[2] - y[1])
            best = min(best, max(fv, y[1] - gain_cap))
        return best

    def _model_min(self, sep, x, y, p, q, best, arc):
        (xp, rp), (xq, rq) = p, q
        #the cut coordinate is a chart coordinate; the model needs the same units as the labels,
        #so scale offsets along the cut by the local riemannian speed of the cut
        sp = float((arc[2] - arc[0]) / (x[2] - x[0])) if x[2] > x[0] else 1.0
        mp = lambda u: np.sqrt(rp * rp + (sp * (u - xp)) ** 2)
        mq = lambda u: np.sqrt(rq * rq + (sp * (u - xq)) ** 2)
        g = y - mp(x) - mq(x)
        d1 = (g[1] - g[0]) / (x[1] - x[0]); d2 = (g[2] - g[1]) / (x[2] - x[1])
        if self.model_fit == "linear":
            #least-squares line through the three residuals: a parabola through three noisy
            #residuals (steep terrain) can bend the wrong way and overshoot the minimum
            xm, gm = x.mean(), g.mean()
            b = float(((x - xm) * (g - gm)).sum() / ((x - xm) ** 2).sum()); a = 0.0; c = gm - b * xm
        else:
            a = (d2 - d1) / (x[2] - x[0]); b = d1 - a * (x[0] + x[1]); c = g[0] - a * x[0] * x[0] - b * x[0]
        #minimise m(u) = mp + mq + a u^2 + b u + c on [x0, x2]. m is a sum of two hyperbolas
        #and a parabola: nearly linear away from the feet, so plain newton overshoots by the
        #whole window. bracketed newton on the derivative (bisection whenever the step leaves
        #the bracket) is robust and converges in ~10 steps
        def deriv(u):
            tp = sp * (u - xp); tq = sp * (u - xq)
            Mp = np.sqrt(rp * rp + tp * tp); Mq = np.sqrt(rq * rq + tq * tq)
            d = sp * tp / Mp + sp * tq / Mq + 2 * a * u + b
            dd = sp * sp * rp * rp / Mp ** 3 + sp * sp * rq * rq / Mq ** 3 + 2 * a
            return d, dd
        lo, hi = float(x[0]), float(x[2])
        dlo, _ = deriv(lo); dhi, _ = deriv(hi)
        if dlo >= 0:
            u = lo
        elif dhi <= 0:
            u = hi
        else:
            u = float(x[1])
            for _ in range(40):
                d, dd = deriv(u)
                if d < 0:
                    lo = u
                else:
                    hi = u
                un = u - d / dd if dd > 0 else 0.5 * (lo + hi)
                if not (lo < un < hi):
                    un = 0.5 * (lo + hi)
                if abs(un - u) < 1e-7 * max(1.0, abs(u)):
                    u = un; break
                u = un
        val = float(mp(u) + mq(u) + a * u * u + b * u + c)
        #never below the plain triangle-inequality floor of the samples' own model
        return min(best, max(val, float(mp(u) + mq(u))))

    def _arc(self, sep, x):
        #riemannian arc length along the cut at the given cut coordinates
        dense, arc = self.cut_arc[sep]
        return np.interp(x, dense, arc)

    def query(self, s, t):
        if s == t:
            return 0.0
        return min(self.cut_bound(s, t), self.portal_bound(s, t), self.leaf_bound(s, t), self.halo_bound(s, t))

    def query_path(self, s, t):
        raise NotImplementedError("the interpolated crossing point is not a graph vertex")
