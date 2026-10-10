"""Certified Separator Cut Index: a (1+eps) terrain distance oracle with per-query certificates.

The construction follows the theorem in notes/guarantee.md step by step:

* separator tree of axis-parallel cut lines through the chart, each placed strictly between vertex
  coordinates (no vertex ever lies on a cut); recursion stops when the members' bounding box has
  chart diagonal <= R_near, so every same-leaf pair is a near pair;
* window lemma: the shortest path between two members on opposite sides crosses the cut line
  within chart offset mu/lambda_min of the members' box, mu = an upper bound of half the largest
  member-pair distance (from one graph Dijkstra out of a hub member);
* portals on the cut, spaced by Riemannian arc length: gap <= g inside the box, and in the tails
  gap <= max(g, 2 eps lambda_min c) with c the chart offset of the interval's near end; every run
  of the line inside the mesh has portals at both ends;
* labels ell(u, p) = exact geodesic distance (MMP sweep from p, stopped at the members) computed on
  the sub-mesh that provably contains every member-portal geodesic (ellipse bound);
* near tables: exact distances to every vertex within chart distance R_near = g/(eps lambda_min).

Query: near pairs exactly; otherwise UB = min_j ell(s,p_j)+ell(t,p_j) (the length of an actual
path), LB = max(chord, min_j (f_j+f_{j+1})/2 - A_j) <= d(s,t), and an interpolated estimate clamped
to [LB, UB]. Theorem: (1-eps) d <= LB <= d <= UB <= (1+eps) d for every pair.
"""
import math
import time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
from ..terrain.submesh import FaceGrid, extract, insert_points, ExactSolver


def mesh_graph_csr(tm):
    #the mesh graph as a symmetric CSR matrix; an edge costs its 3-d length times the cheaper of
    #its adjacent face weights (walking on the cheap side), so every graph path is a surface curve
    #whose cost is at least the true weighted distance: graph distances are upper bounds
    n = len(tm.V)
    E = np.vstack([tm.F[:, [0, 1]], tm.F[:, [1, 2]], tm.F[:, [2, 0]]])
    fe = np.tile(np.arange(len(tm.F)), 3)
    E = np.sort(E, 1)
    key = E[:, 0].astype(np.int64) * n + E[:, 1]
    order = np.argsort(key, kind="stable")
    key, E, fe = key[order], E[order], fe[order]
    first = np.r_[True, key[1:] != key[:-1]]
    grp = np.cumsum(first) - 1
    wmin = np.full(int(first.sum()), np.inf)
    np.minimum.at(wmin, grp, tm.w[fe])
    Eu = E[first]
    L = np.linalg.norm(tm.V[Eu[:, 0]] - tm.V[Eu[:, 1]], axis=1) * wmin
    G = csr_matrix((np.r_[L, L], (np.r_[Eu[:, 0], Eu[:, 1]], np.r_[Eu[:, 1], Eu[:, 0]])), shape=(n, n))
    return G, Eu


def f32_up(d):
    #float32 copy rounded UP, so a stored distance is never below the exact one (labels and near
    #distances stay upper bounds bit-exactly; lower bounds subtract the matching relative margin)
    d = np.asarray(d, dtype=np.float64)
    x = d.astype(np.float32)
    return np.where(x.astype(np.float64) < d, np.nextafter(x, np.float32(np.inf)), x).astype(np.float32)


F32_REL = 2.0 ** -22 #two float32 ulps, relative


class ExactLabels:
    #exact labels on an unweighted terrain: MMP sweeps (pygeodesic) with early termination
    name = "exact MMP"
    exact = True

    def __init__(self, tm, grid):
        if getattr(tm, "weighted", False):
            raise ValueError("MMP is exact only on unweighted terrains; use SteinerLabels")
        self.tm, self.grid = tm, grid

    def cut(self, faces, X, along, members):
        V, F, used, remap = extract(self.tm, faces)
        V2, F2, pidx, Xu = insert_points(self.tm, V, F, faces, X, along=along)
        tgt = remap[members]
        if (tgt < 0).any():
            raise RuntimeError("a member is outside the label sub-mesh")
        sol = ExactSolver(V2, F2)
        L = np.full((len(members), len(X)), np.inf, dtype=np.float32)
        for j, p in enumerate(pidx):
            if p >= 0:
                L[:, j] = f32_up(sol.dist(p, tgt))
        return L, Xu

    def near(self, faces, members, targets):
        V, F, used, remap = extract(self.tm, faces)
        sol = ExactSolver(V, F)
        out = []
        for s, T in zip(members, targets):
            out.append(f32_up(sol.dist(remap[s], remap[T])) if len(T) else np.zeros(0, np.float32))
        return out


class CertifiedCutIndex:
    def __init__(self, tm, eps=0.1, gap=None, gap_factor=0.5, oracle=None, interpolate="model",
                 hub_slack=4.0, verbose=False):
        tt = time.perf_counter
        self.timing = {}
        t0 = tt()
        self.tm = tm
        self.n = len(tm.V)
        self.xy = tm.xy
        self.V3 = tm.V
        self.lam = float(np.min(tm.w)) #lambda_min of the metric (face weights; 1 when unweighted)
        P = tm.xy[tm.F]
        self.edge = float(np.median(np.linalg.norm(P[:, 0] - P[:, 1], axis=1)))
        self.g = float(gap) if gap else float(gap_factor) * self.edge
        self.eps = float(eps)
        self.Rn = self.g / (self.eps * self.lam)
        #near tables and the near test use R_near (1 + 1e-9): every pair inside one leaf (chart
        #diagonal <= R_near) is a near pair even after rounding
        self.Rq = self.Rn * (1.0 + 1e-9)
        self.interpolate = interpolate
        self.hub_slack = hub_slack
        self.verbose = verbose
        self.grid = FaceGrid(tm)
        self.G, self.E = mesh_graph_csr(tm)
        self.oracle = oracle or ExactLabels(tm, self.grid)
        self.kd = cKDTree(tm.xy)
        self.ucoord = [np.unique(tm.xy[:, 0]), np.unique(tm.xy[:, 1])]
        self.timing["setup"] = tt() - t0
        t0 = tt(); self._tree(); self.timing["tree"] = tt() - t0
        t0 = tt(); self._cut_geometry(); self.timing["cut_geometry"] = tt() - t0
        t0 = tt(); self._labels(); self.timing["labels"] = tt() - t0
        t0 = tt(); self._near(); self.timing["near_tables"] = tt() - t0
        t0 = tt(); self._feet(); self.timing["feet"] = tt() - t0

    # ------------------------------------------------------------------ separator tree
    def _split_value(self, c, axis):
        #a value strictly between the two middle member coordinates and different from EVERY vertex
        #coordinate of the mesh on that axis, so no vertex lies on the cut line
        cs = np.sort(c)
        k = len(cs) // 2
        lo, hi = cs[k - 1], cs[k]
        if lo == hi: #ties at the median: move to the nearest distinct split
            d = np.nonzero(cs != cs[k])[0]
            if len(d) == 0:
                return None
            j = d[np.argmin(np.abs(d - k))]
            if j < k:
                lo, hi = cs[j], cs[j + 1]
            else:
                lo, hi = cs[j - 1], cs[j]
        U = self.ucoord[axis]
        inner = U[(U > lo) & (U < hi)]
        pts = np.concatenate([[lo], inner, [hi]])
        gaps = np.diff(pts)
        i = int(np.argmax(gaps))
        return 0.5 * (pts[i] + pts[i + 1])

    def _tree(self):
        n = self.n
        maxd = 64
        self.anc = np.full((n, maxd), -1, dtype=np.int32)
        self.row = np.full((n, maxd), -1, dtype=np.int32)
        self.depth_of = np.zeros(n, dtype=np.int32)
        self.cut_axis, self.cut_v, self.cut_box, self.cut_members, self.cut_depth = [], [], [], [], []
        self.leaf_of = np.full(n, -1, dtype=np.int64)
        self.leaves = []
        stack = [(np.arange(n), 0)]
        while stack:
            M, depth = stack.pop()
            Pm = self.xy[M]
            lo, hi = Pm.min(0), Pm.max(0)
            v = None
            if len(M) > 1 and np.hypot(*(hi - lo)) > self.Rn and depth < maxd:
                a = int(np.argmax(hi - lo))
                v = self._split_value(Pm[:, a], a)
                if v is None:
                    a = 1 - a
                    v = self._split_value(Pm[:, a], a)
            if v is None:
                self.leaf_of[M] = len(self.leaves)
                self.leaves.append(M)
                self.depth_of[M] = depth
                continue
            cid = len(self.cut_axis)
            self.cut_axis.append(a); self.cut_v.append(float(v)); self.cut_box.append((lo, hi))
            self.cut_members.append(M); self.cut_depth.append(depth)
            self.anc[M, depth] = cid
            self.row[M, depth] = np.arange(len(M), dtype=np.int32)
            stack.append((M[Pm[:, a] > v], depth + 1))
            stack.append((M[Pm[:, a] < v], depth + 1))
        md = int(self.depth_of.max()) + 1
        self.anc = np.ascontiguousarray(self.anc[:, :md])
        self.row = np.ascontiguousarray(self.row[:, :md])

    # ------------------------------------------------------------------ cut lines
    def _line_pieces(self, a, v, lo, hi):
        #exact geometry of the cut line {x_a = v}, u in [lo, hi]: the pieces inside each face,
        #their lifted speed, grouped into runs (maximal intervals inside the mesh)
        tm = self.tm
        b = 1 - a
        box = (v, v, lo, hi) if a == 0 else (lo, hi, v, v)
        faces = self.grid.faces_in_box(*box)
        if len(faces) == 0:
            return []
        Q = self.xy[tm.F[faces]] #(m, 3, 2)
        da = Q[:, :, a] - v
        us = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            cross = da[:, i] * da[:, j] < 0
            t = np.where(cross, da[:, i] / np.where(cross, da[:, i] - da[:, j], 1.0), np.nan)
            us.append(np.where(cross, Q[:, i, b] + t * (Q[:, j, b] - Q[:, i, b]), np.nan))
        U = np.stack(us, 1)
        cnt = np.isfinite(U).sum(1)
        keep = cnt == 2 #no vertex on the line, so a crossed face is crossed through two edges
        faces, U = faces[keep], U[keep]
        u0 = np.nanmin(U, 1); u1 = np.nanmax(U, 1)
        u0c = np.maximum(u0, lo); u1c = np.minimum(u1, hi)
        ok = u1c > u0c
        faces, u0c, u1c = faces[ok], u0c[ok], u1c[ok]
        if len(faces) == 0:
            return []
        grad = tm.b[faces] if a == 0 else tm.a[faces]
        speed = tm.w[faces] * np.sqrt(1.0 + grad * grad)
        o = np.argsort(u0c)
        faces, u0c, u1c, speed = faces[o], u0c[o], u1c[o], speed[o]
        runs = []
        tol = 1e-9 * max(1.0, self.edge)
        start = 0
        for i in range(1, len(faces) + 1):
            if i == len(faces) or u0c[i] > u1c[i - 1] + tol:
                ub = np.concatenate([u0c[start:i], [u1c[i - 1]]])
                ab = np.concatenate([[0.0], np.cumsum((u1c[start:i] - u0c[start:i]) * speed[start:i])])
                runs.append((ub, ab, float(speed[start:i].max())))
                start = i
        return runs

    def _place(self, ub, ab, r_lo, r_hi):
        #portal coordinates on one run (breakpoints ub, cumulative arc ab), following the gap rule
        U0, U1 = ub[0], ub[-1]
        arc = lambda u: np.interp(u, ub, ab)
        inv = lambda s: np.interp(s, ab, ub)
        pos = [U0, U1]
        g, k2 = self.g, 2.0 * self.eps * self.lam
        rl, rh = max(U0, r_lo), min(U1, r_hi)
        if rh > rl: #in-box part: uniform in arc length, gap <= g
            A0, A1 = arc(rl), arc(rh)
            m = max(1, int(np.ceil((A1 - A0) / g)))
            pos += list(inv(np.linspace(A0, A1, m + 1)))
        if U0 < r_lo: #left tail, walking away from the box
            cur = min(U1, r_lo)
            pos.append(cur)
            while cur > U0:
                c = r_lo - cur
                nxt = arc(cur) - max(g, k2 * c)
                if nxt <= 0.0:
                    break
                cur = float(inv(nxt))
                pos.append(cur)
        if U1 > r_hi: #right tail
            cur = max(U0, r_hi)
            pos.append(cur)
            while cur < U1:
                c = cur - r_hi
                nxt = arc(cur) + max(g, k2 * c)
                if nxt >= ab[-1]:
                    break
                cur = float(inv(nxt))
                pos.append(cur)
        pos = np.unique(np.clip(np.array(pos), U0, U1))
        keep = np.r_[True, np.diff(pos) > 1e-9 * max(1.0, self.edge)]
        return pos[keep]

    def _cut_geometry(self):
        lo_dom, hi_dom = self.xy.min(0), self.xy.max(0)
        self.cut_pos, self.cut_arc, self.cut_run, self.cut_gap, self.cut_xy = [], [], [], [], []
        self.cut_runs, self.cut_mu, self.cut_hub = [], [], []
        for cid, M in enumerate(self.cut_members):
            a, v = self.cut_axis[cid], self.cut_v[cid]
            b = 1 - a
            lo, hi = self.cut_box[cid]
            r_lo, r_hi = lo[b], hi[b]
            #hub upper bounds: d(s,t) <= d(s,h) + d(h,t) <= 2 mu
            ctr = 0.5 * (lo + hi)
            hub = int(M[np.argmin(np.linalg.norm(self.xy[M] - ctr, axis=1))])
            lim = self.hub_slack * np.hypot(*(hi - lo)) * 2.0 + 10 * self.edge
            dh = dijkstra(self.G, indices=hub, limit=lim)
            if not np.isfinite(dh[M]).all():
                dh = dijkstra(self.G, indices=hub)
            mu = float(dh[M].max())
            off = mu / self.lam
            wlo, whi = max(r_lo - off, lo_dom[b]), min(r_hi + off, hi_dom[b])
            runs = self._line_pieces(a, v, wlo, whi)
            pos, arcs, rid = [], [], []
            for r, (ub, ab, smax) in enumerate(runs):
                p = self._place(ub, ab, r_lo, r_hi)
                pos.append(p); arcs.append(np.interp(p, ub, ab)); rid.append(np.full(len(p), r))
            pos = np.concatenate(pos) if pos else np.zeros(0)
            arcs = np.concatenate(arcs) if arcs else np.zeros(0)
            rid = np.concatenate(rid) if rid else np.zeros(0, dtype=int)
            gap = np.where(rid[1:] == rid[:-1], np.diff(arcs), np.inf) if len(pos) > 1 else np.zeros(0)
            X = np.zeros((len(pos), 2)); X[:, a] = v; X[:, b] = pos
            self.cut_pos.append(pos); self.cut_arc.append(arcs); self.cut_run.append(rid)
            self.cut_gap.append(gap); self.cut_xy.append(X); self.cut_runs.append(runs)
            self.cut_mu.append(mu); self.cut_hub.append(hub)

    # ------------------------------------------------------------------ labels
    def _labels(self):
        tm = self.tm
        self.L = []
        self.label_submesh = []
        self.cut_end = [None] * len(self.cut_members)
        self.cut_lb_ok = np.ones(len(self.cut_members), dtype=bool)
        for cid, M in enumerate(self.cut_members):
            X = self.cut_xy[cid]
            a = self.cut_axis[cid]
            hub, mu = self.cut_hub[cid], self.cut_mu[cid]
            lo, hi = self.cut_box[cid]
            if len(X) == 0:
                self.L.append(None); self._finalize_cut(cid, np.zeros((len(M), 0), np.float32), X)
                self.label_submesh.append(0); continue
            #upper bound of every portal-member distance through the hub
            lim = self.hub_slack * np.hypot(*(hi - lo)) * 2.0 + 10 * self.edge + 2 * mu
            dh = dijkstra(self.G, indices=hub, limit=lim)
            Ux = tm.segment_lengths_batch(X, np.repeat(self.xy[hub][None], len(X), 0))
            _, vx = self.kd.query(X)
            Uv = tm.segment_lengths_batch(X, self.xy[vx]) + dh[vx]
            Ux = np.minimum(Ux, Uv)
            if not np.isfinite(Ux).all():
                dh = dijkstra(self.G, indices=hub)
                Ux = np.minimum(Ux, tm.segment_lengths_batch(X, self.xy[vx]) + dh[vx])
            R = (float(np.max(Ux)) + mu) / (2.0 * self.lam)
            Pall = np.vstack([self.xy[M], X])
            blo, bhi = Pall.min(0) - R, Pall.max(0) + R
            faces = self.grid.faces_in_box(blo[0], bhi[0], blo[1], bhi[1])
            along = np.zeros(2); along[1 - a] = 1.0
            Lc, Xu = self.oracle.cut(faces, X, along, M)
            self.L.append(None)
            self._finalize_cut(cid, Lc, Xu)
            self.label_submesh.append(len(faces))
            if self.verbose and cid % 200 == 0:
                print(f"  cut {cid}/{len(self.cut_members)}: {len(M)} members, {len(X)} portals, sub-mesh {len(faces)} faces", flush=True)

    def _finalize_cut(self, cid, Lc, Xu):
        #arc positions, gaps and run ends from the portals as inserted (a portal within 1e-7 of a
        #face edge is moved along the cut by ~1e-5 edge lengths, so the gaps of the placement are
        #recomputed). a portal without a finite label is dropped and its neighbours' gap merges.
        #if a run's first (last) portal is not at the run's end, the stretch in between needs its
        #own lower-bound term f_j - 2 * (arc of the stretch); a run left without portals switches
        #the cut's interval bound off (LB = chord), the certificate stays sound either way
        a = self.cut_axis[cid]
        b = 1 - a
        runs = self.cut_runs[cid]
        ok = np.isfinite(Lc).all(0) if Lc.shape[0] else np.ones(Lc.shape[1], dtype=bool)
        Lc, Xu = Lc[:, ok], Xu[ok]
        pos = Xu[:, b] if len(Xu) else np.zeros(0)
        order = np.argsort(pos, kind="stable")
        Lc, Xu, pos = Lc[:, order], Xu[order], pos[order]
        starts = np.array([ub[0] for ub, _, _ in runs]) if runs else np.zeros(0)
        rid = np.clip(np.searchsorted(starts, pos, side="right") - 1, 0, max(len(runs) - 1, 0))
        arcs = np.zeros(len(pos))
        end_i, end_g = [], []
        lb_ok = True
        for r, (ub, ab, _) in enumerate(runs):
            m = np.nonzero(rid == r)[0]
            if len(m) == 0:
                lb_ok = False
                continue
            arcs[m] = np.interp(pos[m], ub, ab)
            lead, trail = arcs[m[0]] - ab[0], ab[-1] - arcs[m[-1]]
            if lead > 0:
                end_i.append(m[0]); end_g.append(2.0 * lead)
            if trail > 0:
                end_i.append(m[-1]); end_g.append(2.0 * trail)
        gap = np.where(rid[1:] == rid[:-1], np.diff(arcs), np.inf) if len(pos) > 1 else np.zeros(0)
        self.L[cid] = np.ascontiguousarray(Lc)
        self.cut_xy[cid] = Xu
        self.cut_pos[cid] = pos
        self.cut_arc[cid] = arcs
        self.cut_run[cid] = rid
        self.cut_gap[cid] = gap
        self.cut_end[cid] = (np.array(end_i, dtype=np.int64), np.array(end_g, dtype=np.float64))
        self.cut_lb_ok[cid] = lb_ok

    # ------------------------------------------------------------------ near tables
    def _near(self):
        n = self.n
        ids_list = [None] * n
        d_list = [None] * n
        for M in self.leaves:
            T = self.kd.query_ball_point(self.xy[M], self.Rq)
            T = [np.array(sorted(set(t) - {int(s)}), dtype=np.int64) for s, t in zip(M, T)]
            allt = np.unique(np.concatenate([M] + [t for t in T if len(t)])) if len(M) else M
            #graph upper bounds of the near distances -> dilation of the sub-mesh
            D = dijkstra(self.G, indices=M, limit=self.hub_slack * self.Rq + 10 * self.edge)
            Umax = 0.0
            for i, t in enumerate(T):
                if len(t):
                    u = D[i, t]
                    if not np.isfinite(u).all():
                        u = dijkstra(self.G, indices=int(M[i]))[t]
                    Umax = max(Umax, float(u.max()))
            R = Umax / (2.0 * self.lam)
            Pall = self.xy[allt]
            blo, bhi = Pall.min(0) - R, Pall.max(0) + R
            faces = self.grid.faces_in_box(blo[0], bhi[0], blo[1], bhi[1])
            res = self.oracle.near(faces, M, T)
            for s, t, d in zip(M, T, res):
                ids_list[s] = t.astype(np.int32)
                d_list[s] = d
        self.near_ptr = np.zeros(n + 1, dtype=np.int64)
        self.near_ptr[1:] = np.cumsum([len(x) for x in ids_list])
        self.near_ids = np.concatenate(ids_list).astype(np.int32)
        self.near_d = np.concatenate(d_list).astype(np.float32)

    def near_lookup(self, s, t):
        a, b = self.near_ptr[s], self.near_ptr[s + 1]
        j = a + int(np.searchsorted(self.near_ids[a:b], t))
        if j < b and self.near_ids[j] == t:
            return float(self.near_d[j])
        return None

    # ------------------------------------------------------------------ feet (model interpolation)
    def _feet(self):
        tm = self.tm
        self.foot_r, self.foot_arc, self.foot_run = [], [], []
        for cid, M in enumerate(self.cut_members):
            a, v = self.cut_axis[cid], self.cut_v[cid]
            b = 1 - a
            P = self.xy[M]
            Q = P.copy(); Q[:, a] = v
            r = tm.segment_lengths_batch(P, Q)
            r = np.where(np.isfinite(r), r, self.lam * np.abs(P[:, a] - v))
            fa = np.full(len(M), np.nan); fr = np.full(len(M), -1, dtype=np.int16)
            for k, (ub, ab, _) in enumerate(self.cut_runs[cid]):
                inside = (P[:, b] >= ub[0]) & (P[:, b] <= ub[-1])
                fa[inside] = np.interp(P[inside, b], ub, ab); fr[inside] = k
            self.foot_r.append(r.astype(np.float32)); self.foot_arc.append(fa.astype(np.float32))
            self.foot_run.append(fr)

    # ------------------------------------------------------------------ queries
    def _cut_of(self, s, t):
        diff = np.nonzero(self.anc[s] != self.anc[t])[0]
        k = int(diff[0]) if len(diff) else min(int(self.depth_of[s]), int(self.depth_of[t]))
        return int(self.anc[s, k - 1]), int(self.row[s, k - 1]), int(self.row[t, k - 1])

    def bounds(self, s, t):
        #(LB, UB, estimate) with LB <= d(s,t) <= UB guaranteed
        if s == t:
            return 0.0, 0.0, 0.0
        dv = self.xy[s] - self.xy[t]
        if dv @ dv <= self.Rq * self.Rq:
            d = self.near_lookup(s, t)
            if d is not None:
                lo = d * (1.0 - F32_REL) if self.oracle.exact else self.lam * float(np.linalg.norm(self.V3[s] - self.V3[t]))
                return lo, d, d
        C, rs, rt = self._cut_of(s, t)
        LC = self.L[C]
        f = np.add(LC[rs], LC[rt], dtype=np.float64)
        i = int(np.argmin(f))
        ub = float(f[i])
        lb = 0.0
        if self.oracle.exact and self.cut_lb_ok[C]: #the interval bound needs exact labels
            gaps = self.cut_gap[C]
            lb = math.inf
            if len(f) > 1: #consecutive portals of the same run; run changes have gap = inf and are skipped
                mid = np.where(np.isfinite(gaps), 0.5 * (f[:-1] + f[1:]) - gaps, np.inf)
                lb = float(np.min(mid))
            ei, eg = self.cut_end[C]
            if len(ei): #stretches between a run's end and its first / last portal
                lb = min(lb, float(np.min(f[ei] - eg)))
            if not math.isfinite(lb):
                lb = 0.0
        dd = self.V3[s] - self.V3[t]
        chord = self.lam * math.sqrt(float(dd @ dd))
        lb = max(lb * (1.0 - 2 * F32_REL), chord)
        est = ub
        if self.interpolate:
            est = self._interp(C, f, i, rs, rt)
        return lb, ub, min(max(est, lb), ub)

    def _interp(self, C, f, i, rs, rt):
        run = self.cut_run[C]
        if not (0 < i < len(f) - 1) or not (run[i - 1] == run[i] == run[i + 1]):
            return float(f[i])
        x = self.cut_arc[C][i - 1:i + 2].astype(np.float64)
        y = f[i - 1:i + 2]
        if self.interpolate == "model":
            fr = self.foot_run[C]
            if fr[rs] == run[i] and fr[rt] == run[i]:
                xs, r1 = float(self.foot_arc[C][rs]), float(self.foot_r[C][rs])
                xt, r2 = float(self.foot_arc[C][rt]), float(self.foot_r[C][rt])
                if min(r1, r2) < 2.0 * (x[2] - x[0]):
                    return _model_min(x, y, xs, r1, xt, r2)
        d1 = (y[1] - y[0]) / (x[1] - x[0]); d2 = (y[2] - y[1]) / (x[2] - x[1])
        a2 = (d2 - d1) / (x[2] - x[0])
        if a2 <= 0:
            return float(y[1])
        xv = 0.5 * (x[0] + x[1]) - d1 / (2.0 * a2)
        if not (x[0] <= xv <= x[2]):
            return float(y[1])
        fv = y[0] + d1 * (xv - x[0]) + a2 * (xv - x[0]) * (xv - x[1])
        cap = 0.25 * min(y[0] - y[1], y[2] - y[1])
        return float(max(fv, y[1] - cap))

    def query(self, s, t):
        return self.bounds(s, t)[2]

    def query_ub(self, s, t):
        return self.bounds(s, t)[1]

    def nbytes(self):
        b = sum(L.nbytes for L in self.L)
        b += self.near_ids.nbytes + self.near_d.nbytes + self.near_ptr.nbytes
        b += sum(x.nbytes + y.nbytes + z.nbytes + w.nbytes for x, y, z, w in zip(self.cut_arc, self.cut_run, self.cut_gap, self.cut_xy))
        b += sum(i.nbytes + g.nbytes for i, g in self.cut_end) + self.cut_lb_ok.nbytes
        b += sum(x.nbytes + y.nbytes + z.nbytes for x, y, z in zip(self.foot_r, self.foot_arc, self.foot_run))
        b += self.anc.nbytes + self.row.nbytes
        return b

    def stats(self):
        return dict(cuts=len(self.cut_members), leaves=len(self.leaves), portals=int(sum(len(p) for p in self.cut_pos)),
                    labels=int(sum(L.size for L in self.L)), near_entries=int(len(self.near_ids)),
                    depth=int(self.anc.shape[1]), R_near=self.Rn, gap=self.g)


def _model_min(x, y, xs, r1, xt, r2):
    #f(u) ~ sqrt(r1^2 + (u-xs)^2) + sqrt(r2^2 + (u-xt)^2) + smooth remainder (parabola through the
    #residuals); minimised on [x0, x2] by bracketed newton. heuristic estimate, clamped by the caller
    m = lambda u, xc, r: np.sqrt(r * r + (u - xc) ** 2)
    g = y - m(x, xs, r1) - m(x, xt, r2)
    d1 = (g[1] - g[0]) / (x[1] - x[0]); d2 = (g[2] - g[1]) / (x[2] - x[1])
    a = (d2 - d1) / (x[2] - x[0]); b = d1 - a * (x[0] + x[1]); c = g[0] - a * x[0] * x[0] - b * x[0]
    def der(u):
        M1, M2 = m(u, xs, r1), m(u, xt, r2)
        return (u - xs) / M1 + (u - xt) / M2 + 2 * a * u + b, r1 * r1 / M1 ** 3 + r2 * r2 / M2 ** 3 + 2 * a
    lo, hi = float(x[0]), float(x[2])
    if der(lo)[0] >= 0:
        u = lo
    elif der(hi)[0] <= 0:
        u = hi
    else:
        u = float(x[1])
        for _ in range(40):
            d, dd = der(u)
            if d < 0:
                lo = u
            else:
                hi = u
            un = u - d / dd if dd > 0 else 0.5 * (lo + hi)
            if not (lo < un < hi):
                un = 0.5 * (lo + hi)
            if abs(un - u) < 1e-9 * max(1.0, abs(u)):
                u = un
                break
            u = un
    val = m(u, xs, r1) + m(u, xt, r2) + a * u * u + b * u + c
    return float(min(float(y[1]), max(val, m(u, xs, r1) + m(u, xt, r2))))
