"""Exact two-point shortest-path queries among obstacles: hyperbolic (H^2) versus flat discs.

Checks notes/hyperbolic.md section 7: under bounded corner density and bounded visibility, exact distance labels on the
tangent visibility graph of the obstacle corners have O(R log m) = O(log^2 m) entries per corner in H^2 (radial cuts have
length R), against about sqrt(m) in a flat disc. Any two free points are then answered exactly through the corners they see.

  python3 bench/hyperbolic_exact.py case '<json spec>'   # one case, prints a json line (used by `run`)
  python3 bench/hyperbolic_exact.py run [--quick]        # the sweep -> results/hyperbolic_exact.json (resumable)
  python3 bench/hyperbolic_exact.py report               # -> results/hyperbolic_exact.png, notes/hyperbolic_exact.md

Obstacles are regular geodesic k-gons of circumradius r_o. Their centres are hard-core Poisson: dart throwing at target
intensity lam per unit area in the disc of radius R - r_o, accepted if no accepted centre is within 2 r_o (1 + gap). The
domain is the disc minus the obstacles. The disc is convex, so shortest paths bend only at obstacle corners and never use
the rim.

Geometry. Points are stored in polar coordinates (r, theta). To look around a point u the isometry that moves u to the
origin is applied and Klein coordinates are used: geodesics are straight there, so visibility around the origin is the
same Euclidean computation in both geometries. Distances use the cancellation-free polar formulas of hyperbolic_h2h.py.
"""
import argparse, json, os, subprocess, sys, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, maximum_bipartite_matching, connected_components
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from bench.hyperbolic_h2h import Flat, Hyp, pairs_within, footprint_mb, TWO_PI  # noqa: E402

OUT_JSON = os.path.join(ROOT, "results", "hyperbolic_exact.json")
OUT_PNG = os.path.join(ROOT, "results", "hyperbolic_exact.png")
OUT_MD = os.path.join(ROOT, "notes", "hyperbolic_exact.md")
GEOM = {"flat": Flat, "hyp": Hyp}
EPS_ANG = 1e-12
EPS_CROSS = 1e-15


def wrap(a):  # to (-pi, pi]
    return np.pi - np.mod(np.pi - a, TWO_PI)


def cross(ax, ay, bx, by):
    return ax * by - ay * bx


# ---------------------------------------------------------------- geometry
def local_coords(geom, ru, tu, r, t):
    """Coordinates of points (r, t) after the isometry moving (ru, tu) to the origin: Klein (hyperbolic) or Euclidean.
    The isometry also rotates so that the direction away from the global origin is +x; visibility does not care."""
    d = t - tu
    if geom == "flat":
        return r * np.cos(d) - ru, r * np.sin(d)
    ch, sh = np.cosh(r), np.sinh(r)
    x1, x2 = sh * np.cos(d), sh * np.sin(d)
    x0 = ch * np.cosh(ru) - x1 * np.sinh(ru)
    y1 = x1 * np.cosh(ru) - ch * np.sinh(ru)
    return y1 / x0, x2 / x0


def global_xy(geom, r, t):
    """Global chart in which geodesics are straight: Klein (hyperbolic) or the plane (flat)."""
    s = np.tanh(r) if geom == "hyp" else r
    return s * np.cos(t), s * np.sin(t)


def within_local(geom, x, y, rho):
    """|local point| < rho in the metric of the geometry (Klein: d(0, x) = artanh |x|)."""
    q = x * x + y * y
    lim = np.tanh(rho) if geom == "hyp" else rho
    return q < lim * lim


def offset_point(geom, rc, tc, rho, psi):
    """The point at distance rho from (rc, tc) in direction psi, measured from the direction away from the origin."""
    if geom == "flat":
        cx, cy = rc * np.cos(tc), rc * np.sin(tc)
        x = cx + rho * np.cos(tc + psi)
        y = cy + rho * np.sin(tc + psi)
        return np.hypot(x, y), np.arctan2(y, x)
    ch, sh = np.cosh(rho), np.sinh(rho)
    x1 = sh * np.cos(psi) * np.cosh(rc) + ch * np.sinh(rc)
    x2 = sh * np.sin(psi) + 0.0 * rc
    return np.arcsinh(np.hypot(x1, x2)), tc + np.arctan2(x2, x1)


# ---------------------------------------------------------------- obstacles
class Domain:
    """Disc of radius R minus n_obs disjoint regular geodesic k-gons. corners: (n_obs, k) global ids, ccw order."""

    def __init__(self, geom, R, lam, r_o, k, gap, rng, lam_eff=None):
        g = GEOM[geom]
        self.geom, self.g, self.R, self.r_o, self.k = geom, g, R, r_o, k
        Rc = R - r_o
        N = rng.poisson(lam * g.area(Rc))
        rc, tc = g.sample(Rc, N, rng)
        dmin = 2.0 * r_o * (1.0 + gap)
        I, J, _ = pairs_within(g, rc, tc, dmin)
        nb = [[] for _ in range(N)]
        for a, b in zip(I.tolist(), J.tolist()):
            nb[a].append(b); nb[b].append(a)
        removed = np.zeros(N, dtype=bool)
        acc = []
        for i in range(N):  # samples come in random order
            if not removed[i]:
                acc.append(i)
                removed[nb[i]] = True
        acc = np.array(acc, dtype=np.int64)
        if lam_eff is not None:  # thin the (near-jammed) hard-core sample to a set intensity
            want = int(round(lam_eff * g.area(Rc)))
            if want > len(acc):
                raise RuntimeError(f"lam_eff {lam_eff} above the hard-core sample's {len(acc) / g.area(Rc):.3f}")
            acc = np.sort(rng.choice(acc, want, replace=False))
        self.rc, self.tc = rc[acc], tc[acc]
        n = len(acc)
        self.n_obs = n
        self.lam_target, self.lam_real = lam, n / g.area(Rc)
        phi = TWO_PI * rng.random(n)
        psi = phi[:, None] + TWO_PI * np.arange(k)[None, :] / k
        r, t = offset_point(geom, self.rc[:, None], self.tc[:, None], r_o, psi)
        self.r, self.t = r.ravel(), t.ravel()          # corner j of obstacle o has id o*k + j
        self.m = n * k
        self.corners = np.arange(self.m).reshape(n, k)
        self.obs_of = np.repeat(np.arange(n), k)
        ids = np.arange(self.m)
        self.prev = (ids // k) * k + (ids % k - 1) % k
        self.next = (ids // k) * k + (ids % k + 1) % k
        self.ctree = cKDTree(g.kd_xy(self.rc, self.tc)) if n else None
        gx, gy = global_xy(geom, self.r, self.t)
        self.gx, self.gy = gx.reshape(n, k), gy.reshape(n, k)

    def dist(self, r1, t1, r2, t2):
        return self.g.dist(r1, t1, r2, t2)

    def obstacles_near(self, ro, to, rho):
        """Obstacles whose centre is within rho + r_o of (ro, to): every obstacle meeting the rho-ball."""
        if self.n_obs == 0:
            return np.zeros(0, dtype=np.int64)
        xy, B = self.g.kd_query(np.array([ro]), np.array([to]), rho + self.r_o)
        return np.array(self.ctree.query_ball_point(xy[0], B[0] * 1.0001 + 1e-12), dtype=np.int64)

    def rho_all(self, ro):
        return ro + self.R + self.r_o

    def inside_any(self, r, t):
        """Is (r, t) inside some obstacle (strictly)? Checked in the global straight chart."""
        x, y = global_xy(self.geom, np.array([r]), np.array([t]))
        for o in self.obstacles_near(r, t, 0.0):
            px, py = self.gx[o], self.gy[o]
            qx, qy = np.roll(px, -1), np.roll(py, -1)
            if np.all(cross(qx - px, qy - py, x[0] - px, y[0] - py) > 0):
                return True
        return False

    def sample_free(self, n, rng):
        out = []
        while len(out) < n:
            r, t = self.g.sample(self.R, 1, rng)
            if not self.inside_any(r[0], t[0]):
                out.append((float(r[0]), float(t[0])))
        return np.array(out)

    # ------------------------------------------------------------ visibility around one point
    def view(self, ro, to, own=-1, rho0=1.0):
        """Corners visible from (ro, to) and, if own >= 0 (the point is that corner), which of them are tangent.

        Returns (vis ids, tangent mask over vis, last search radius). Around the point (moved to the origin, straight
        chart) every other obstacle j is convex and misses the origin: it spans an open angular interval (a_j, b_j) with
        tangent corners T1, T2. A corner v of another obstacle is hidden by j iff its direction is strictly inside
        (a_j, b_j) and it lies strictly beyond the chord T1 T2; it is hidden by its own obstacle iff it is beyond that
        obstacle's chord (on the far chain). From a corner, directions strictly inside the own obstacle's cone are
        blocked, and of the own obstacle only the two neighbours are visible.

        Sweep: U is the set of directions not yet blocked by an obstacle lying wholly within radius rho_prev. At radius
        rho only obstacles whose interval meets U matter (a blocker of a segment in direction theta has theta in its
        interval), and only corners at distance in [rho_prev, rho) with direction in U are new. rho doubles until U is
        empty or the ball covers the disc. This is exact; it is fast when sightlines are short."""
        k = self.k
        own_obs = -1 if own < 0 else own // k
        U = [(0.0, TWO_PI)]
        if own >= 0:
            nx, ny = local_coords(self.geom, ro, to, self.r[self.next[own]], self.t[self.next[own]])
            px, py = local_coords(self.geom, ro, to, self.r[self.prev[own]], self.t[self.prev[own]])
            c0 = float(np.arctan2(ny, nx)) % TWO_PI
            U = subtract_intervals(U, [(c0, c0 + float(np.mod(np.arctan2(py, px) - c0, TWO_PI)))])
        vis_all, tang_all = [], []
        rho_prev, rho = 0.0, rho0
        while True:
            full = rho + self.r_o >= self.rho_all(ro)
            obs = self.obstacles_near(ro, to, self.rho_all(ro) if full else rho)
            if own_obs >= 0:
                obs = obs[obs != own_obs]
            C = self.corners[obs]
            X, Y = local_coords(self.geom, ro, to, self.r[C], self.t[C])
            phi = np.arctan2(Y.mean(1), X.mean(1)) if len(obs) else np.zeros(0)
            rel = wrap(np.arctan2(Y, X) - phi[:, None])
            ar = np.arange(len(obs))
            i1, i2 = rel.argmin(1), rel.argmax(1)
            a, b = rel[ar, i1], rel[ar, i2]
            lo = np.mod(phi + a, TWO_PI); wd = b - a
            sel = overlaps(lo, wd, U)
            obs, C, X, Y, phi, a, b, lo, wd = obs[sel], C[sel], X[sel], Y[sel], phi[sel], a[sel], b[sel], lo[sel], wd[sel]
            i1, i2 = i1[sel], i2[sel]
            ar = np.arange(len(obs))
            T1x, T1y, T2x, T2y = X[ar, i1], Y[ar, i1], X[ar, i2], Y[ar, i2]
            ex, ey = T2x - T1x, T2y - T1y
            s0 = cross(ex, ey, -T1x, -T1y)
            vx, vy = X.ravel(), Y.ravel()
            vob = np.repeat(ar, k)
            band = ~within_local(self.geom, vx, vy, rho_prev)
            if not full:
                band &= within_local(self.geom, vx, vy, rho)
            band &= in_intervals(np.mod(np.arctan2(vy, vx), TWO_PI), U)
            band &= cross(ex[vob], ey[vob], vx - T1x[vob], vy - T1y[vob]) * s0[vob] >= -EPS_CROSS   # near chain
            if own >= 0:
                band &= ~((cross(nx, ny, vx, vy) > EPS_CROSS) & (cross(vx, vy, px, py) > EPS_CROSS))
            cand = np.nonzero(band)[0]
            if len(cand):
                rl = wrap(np.arctan2(vy[cand], vx[cand])[:, None] - phi[None, :])
                ins = (rl > a[None, :] + EPS_ANG) & (rl < b[None, :] - EPS_ANG)
                sv = cross(ex[None, :], ey[None, :], vx[cand][:, None] - T1x[None, :], vy[cand][:, None] - T1y[None, :])
                ins &= sv * s0[None, :] < -EPS_CROSS
                ins[np.arange(len(cand)), vob[cand]] = False
                cand = cand[~ins.any(1)]
                vis = C.ravel()[cand]
                tang = np.ones(len(vis), dtype=bool)
                if own >= 0:
                    wx, wy = vx[cand], vy[cand]
                    tang &= cross(wx, wy, nx, ny) * cross(wx, wy, px, py) >= -EPS_CROSS
                    ax_, ay_ = local_coords(self.geom, ro, to, self.r[self.prev[vis]], self.t[self.prev[vis]])
                    bx_, by_ = local_coords(self.geom, ro, to, self.r[self.next[vis]], self.t[self.next[vis]])
                    tang &= cross(wx, wy, ax_ - wx, ay_ - wy) * cross(wx, wy, bx_ - wx, by_ - wy) >= -EPS_CROSS
                vis_all.append(vis); tang_all.append(tang)
            if full:
                break
            inside = within_local(self.geom, X, Y, rho).all(1)
            U = subtract_intervals(U, [(l, l + w) for l, w in zip(lo[inside].tolist(), wd[inside].tolist())])
            if not U:
                break
            rho_prev, rho = rho, 2.0 * rho
        vis = np.concatenate(vis_all) if vis_all else np.zeros(0, dtype=np.int64)
        tang = np.concatenate(tang_all) if tang_all else np.zeros(0, dtype=bool)
        if own >= 0:
            vis = np.concatenate([vis, [self.prev[own], self.next[own]]])
            tang = np.concatenate([tang, [True, True]])
        return vis, tang, rho

    def segment_free(self, r1, t1, r2, t2):
        """Is the geodesic segment between two free points free of obstacle interiors?"""
        d = float(self.dist(r1, t1, r2, t2))
        obs = self.obstacles_near(r1, t1, d)
        if len(obs) == 0:
            return True
        C = self.corners[obs]
        X, Y = local_coords(self.geom, r1, t1, self.r[C], self.t[C])
        qx, qy = local_coords(self.geom, r1, t1, np.array([r2]), np.array([t2]))
        return not clip_hits(np.zeros(1), np.zeros(1), qx, qy, X, Y).any()



def split_circ(ivs):
    """Circular intervals (lo in [0, 2pi), hi = lo + width) -> sorted merged non-wrapping segments in [0, 2pi]."""
    seg = []
    for l, h in ivs:
        if h > TWO_PI:
            seg.append((l, TWO_PI)); seg.append((0.0, h - TWO_PI))
        else:
            seg.append((l, h))
    seg.sort()
    out = []
    for l, h in seg:
        if out and l <= out[-1][1] + 1e-12:
            out[-1] = (out[-1][0], max(out[-1][1], h))
        else:
            out.append((l, h))
    return out


def subtract_intervals(U, ivs):
    """U (disjoint segments in [0, 2pi]) minus the union of circular intervals; gaps under 1e-12 are dropped."""
    B = split_circ(ivs)
    out = []
    for l, h in U:
        cur = l
        for bl, bh in B:
            if bh <= cur or bl >= h:
                continue
            if bl > cur:
                out.append((cur, bl))
            cur = max(cur, bh)
            if cur >= h:
                break
        if cur < h:
            out.append((cur, h))
    return [(l, h) for l, h in out if h - l > 1e-12]


def overlaps(lo, wd, U):
    """Which circular intervals [lo, lo + wd] meet the set U (disjoint segments)?"""
    hit = np.zeros(len(lo), dtype=bool)
    for l, h in U:
        hit |= (np.mod(l - lo, TWO_PI) < wd) | (np.mod(lo - l, TWO_PI) < h - l)
    return hit


def in_intervals(th, U):
    hit = np.zeros(len(th), dtype=bool)
    for l, h in U:
        hit |= (th >= l) & (th <= h)
    return hit

def covers_circle(ivs):
    """Do the angular intervals (lo in [0, 2pi), hi = lo + width, width < 2pi) cover the circle up to 1e-12 gaps?"""
    if not ivs:
        return False
    seg = []
    for l, h in ivs:
        if h > TWO_PI:
            seg.append((l, TWO_PI)); seg.append((0.0, h - TWO_PI))
        else:
            seg.append((l, h))
    seg.sort()
    reach = 0.0
    for l, h in seg:
        if l > reach + 1e-12:
            return False
        reach = max(reach, h)
    return reach >= TWO_PI - 1e-12


def clip_hits(ax, ay, bx, by, PX, PY, tol=1e-9):
    """For segments a_i b_i (arrays of length s) and convex ccw polygons (n_p, k): (s, n_p) mask of segments that pass
    through a polygon's interior over a parameter length > tol (Cyrus-Beck clipping). Used as an independent check."""
    QX, QY = np.roll(PX, -1, axis=1), np.roll(PY, -1, axis=1)
    ex, ey = QX - PX, QY - PY                                    # (n_p, k)
    dx, dy = (bx - ax)[:, None, None], (by - ay)[:, None, None]
    # inside half-plane: cross(e, x - P) > 0 ; along the segment x = a + s d
    num = cross(ex[None], ey[None], ax[:, None, None] - PX[None], ay[:, None, None] - PY[None])   # value at s=0
    den = cross(ex[None], ey[None], dx, dy)                                                        # slope in s
    with np.errstate(divide="ignore", invalid="ignore"):
        s = -num / den
    lo = np.where(den > 0, s, -np.inf).max(2)
    hi = np.where(den < 0, s, np.inf).min(2)
    par = (den == 0) & (num <= 0)
    lo = np.maximum(lo, 0.0); hi = np.minimum(hi, 1.0)
    return (hi - lo > tol) & ~par.any(2)


# ---------------------------------------------------------------- tangent visibility graph
def build_graph(dom, rng, n_check=0):
    """Tangent visibility graph on the corners. Returns (I, J, W, stats). Pairs are collected from both ends and kept
    only if both ends agree (visibility and tangency are symmetric; a disagreement is a numerical tie and is counted)."""
    k = dom.k
    pairs, nvis, rhos = {}, np.zeros(dom.m, dtype=np.int64), np.zeros(dom.m)
    for u in range(dom.m):
        vis, tang, rho = dom.view(dom.r[u], dom.t[u], own=u)
        nvis[u] = len(vis); rhos[u] = rho
        for v in vis[tang].tolist():
            key = (u, v) if u < v else (v, u)
            pairs[key] = pairs.get(key, 0) + 1
    both = [p for p, c in pairs.items() if c == 2]
    one = sum(1 for c in pairs.values() if c == 1)
    I = np.array([p[0] for p in both], dtype=np.int64)
    J = np.array([p[1] for p in both], dtype=np.int64)
    W = dom.dist(dom.r[I], dom.t[I], dom.r[J], dom.t[J])
    st = dict(edges=len(I), asym_pairs=one, mean_deg=2.0 * len(I) / max(dom.m, 1),
              max_deg=int(np.bincount(np.concatenate([I, J]), minlength=dom.m).max()) if len(I) else 0,
              vis_mean=float(nvis.mean()), vis_max=int(nvis.max()), max_edge=float(W.max()) if len(W) else 0.0,
              rho_mean=float(rhos.mean()), rho_max=float(rhos.max()))
    if n_check:  # independent check of visibility: Cyrus-Beck clipping in the global chart, random corner pairs
        u = rng.integers(0, dom.m, n_check)
        v = rng.integers(0, dom.m, n_check)
        keep = (dom.obs_of[u] != dom.obs_of[v])
        u, v = u[keep], v[keep]
        mism = 0
        visset = {}
        for a, b in zip(u.tolist(), v.tolist()):
            if a not in visset:
                visset[a] = set(dom.view(dom.r[a], dom.t[a], own=a)[0].tolist())
            ax, ay = global_xy(dom.geom, dom.r[a], dom.t[a])
            bx, by = global_xy(dom.geom, dom.r[b], dom.t[b])
            hit = clip_hits(np.array([ax]), np.array([ay]), np.array([bx]), np.array([by]), dom.gx, dom.gy)[0]
            free = not hit.any()
            mism += free != (b in visset[a])
        st.update(vis_checked=len(u), vis_mismatch=int(mism))
    return I, J, W, st


# ---------------------------------------------------------------- cut-tree labels
def min_vertex_cover(A_lr, nl, nr):
    """Minimum vertex cover of a bipartite graph (CSR nl x nr) by Koenig's theorem. Returns (left idx, right idx)."""
    match = maximum_bipartite_matching(A_lr, perm_type="column")    # match[i] = right partner of left i, or -1
    rmatch = -np.ones(nr, dtype=np.int64)
    ml = np.nonzero(match >= 0)[0]
    rmatch[match[ml]] = ml
    visL = np.zeros(nl, dtype=bool); visR = np.zeros(nr, dtype=bool)
    stack = list(np.nonzero(match < 0)[0])
    visL[stack] = True
    indptr, indices = A_lr.indptr, A_lr.indices
    while stack:   # alternating paths from unmatched left vertices
        i = stack.pop()
        for j in indices[indptr[i]:indptr[i + 1]]:
            if not visR[j]:
                visR[j] = True
                i2 = rmatch[j]
                if i2 >= 0 and not visL[i2]:
                    visL[i2] = True
                    stack.append(i2)
    return np.nonzero(~visL)[0], np.nonzero(visR)[0]


class CutLabels:
    """Exact labels from a recursive cut tree. Each piece P has a separator S(P), a minimum vertex cover of the graph
    edges between its two sides, and stores d_P(s, v) for s in S(P), v in P, inside the subgraph induced by P. The
    children are the two sides of P minus S(P). For u, w: d(u, w) = min over common pieces P, s in S(P) of
    d_P(s, u) + d_P(s, w); the deepest piece holding a whole shortest u-w path has its separator on that path.

    mode 'radial': split the disc by the diameter that halves the vertex count, then each wedge by the radius at the
    median angle (every cut is a radius, length R). mode 'kd': split by the median line of the longer side of the
    bounding box in the straight chart (Klein for H^2; lines there are geodesics)."""

    def __init__(self, n, I, J, W, ang, gx, gy, mode, leaf=24):
        A = csr_matrix((np.concatenate([W, W]), (np.concatenate([I, J]), np.concatenate([J, I]))), shape=(n, n))
        self.n = n
        self.chain = [[] for _ in range(n)]       # per vertex: list of (piece id, position in piece)
        self.D, self.S, self.depth, self.size = [], [], [], []
        stack = [(np.arange(n), None, 0)]
        while stack:
            mem, wedge, dep = stack.pop()
            if len(mem) == 0:
                continue
            pid = len(self.D)
            sub = A[mem][:, mem]
            if len(mem) <= leaf:
                Sloc = np.arange(len(mem)); side = None
            else:
                if mode == "radial":
                    side, wl, wr = radial_split(ang[mem], wedge)
                else:
                    side = kd_split(gx[mem], gy[mem]); wl = wr = None
                L = np.nonzero(side == 0)[0]; Rr = np.nonzero(side == 1)[0]
                A_lr = sub[L][:, Rr]
                A_lr.data[:] = 1.0
                cl, cr = min_vertex_cover(A_lr.tocsr(), len(L), len(Rr))
                Sloc = np.concatenate([L[cl], Rr[cr]])
            Dp = dijkstra(sub, directed=False, indices=Sloc) if len(Sloc) else np.zeros((0, len(mem)))
            self.D.append(Dp); self.S.append(mem[Sloc]); self.depth.append(dep); self.size.append(len(mem))
            for pos, v in enumerate(mem.tolist()):
                self.chain[v].append((pid, pos))
            if side is not None:
                rest = np.ones(len(mem), dtype=bool); rest[Sloc] = False
                stack.append((mem[rest & (side == 0)], wl, dep + 1))
                stack.append((mem[rest & (side == 1)], wr, dep + 1))
        self.label_len = np.array([sum(len(self.S[p]) for p, _ in ch) for ch in self.chain])
        self.entries = int(self.label_len.sum())
        self.height = max(self.depth) + 1

    def query(self, u, w):
        best = np.inf
        for (p, a), (q, b) in zip(self.chain[u], self.chain[w]):
            if p != q:
                break
            Dp = self.D[p]
            if len(Dp):
                best = min(best, float((Dp[:, a] + Dp[:, b]).min()))
        return best

    def query_points(self, Vp, dp, Vq, dq):
        """min over u in Vp, w in Vq of dp[u] + d(u, w) + dq[w], through the pieces (cost O((|Vp| + |Vq|) label))."""
        bp, bq = {}, {}
        for V, dd, out in ((Vp, dp, bp), (Vq, dq, bq)):
            for v, x in zip(V.tolist(), dd.tolist()):
                for p, a in self.chain[v]:
                    Dp = self.D[p]
                    if len(Dp):
                        col = Dp[:, a] + x
                        out[p] = np.minimum(out[p], col) if p in out else col
        best = np.inf
        for p, col in bp.items():
            if p in bq:
                best = min(best, float((col + bq[p]).min()))
        return best

    def sep_by_depth(self):
        d = np.array(self.depth); s = np.array([len(x) for x in self.S]); z = np.array(self.size)
        inner = z > 0
        out = []
        for k in range(self.height):
            sel = (d == k) & inner
            if sel.any():
                out.append([k, int(sel.sum()), float(s[sel].mean()), int(s[sel].max()), float(z[sel].mean())])
        return out


def radial_split(ang, wedge):
    """Split by a ray from the centre. wedge None: by a whole diameter (two rays) halving the count."""
    n = len(ang)
    a = np.mod(ang, TWO_PI)
    if wedge is None:
        srt = np.sort(a)
        cand = srt[np.linspace(0, n - 1, min(n, 512)).astype(int)]
        ext = np.concatenate([srt, srt + TWO_PI])
        cnt = np.searchsorted(ext, cand + np.pi) - np.searchsorted(ext, cand)
        th = cand[np.argmin(np.abs(cnt - n / 2))]
        side = (np.mod(a - th, TWO_PI) >= np.pi).astype(np.int8)
        return side, (th, th + np.pi), (th + np.pi, th + TWO_PI)
    lo, hi = wedge
    rel = np.mod(a - lo, TWO_PI)
    o = np.argsort(rel, kind="stable")
    side = np.ones(n, dtype=np.int8); side[o[: n // 2]] = 0
    cut = 0.5 * (rel[o[n // 2 - 1]] + rel[o[n // 2]])
    return side, (lo, lo + cut), (lo + cut, hi)


def kd_split(x, y):
    n = len(x)
    c = x if (x.max() - x.min()) >= (y.max() - y.min()) else y
    o = np.argsort(c, kind="stable")
    side = np.ones(n, dtype=np.int8); side[o[: n // 2]] = 0
    return side


# ---------------------------------------------------------------- one case
def run_case(spec):
    rng = np.random.default_rng(spec["seed"])
    out = dict(spec=spec, status="ok")
    t0 = time.perf_counter()
    R = spec["R"]
    if spec["geom"] == "flat" and spec.get("match_area_R"):   # flat disc with the area of the hyperbolic one
        R = float(np.sqrt(Hyp.area(spec["match_area_R"]) / np.pi))
    dom = Domain(spec["geom"], R, spec["lam"], spec["r_o"], spec["k"], spec.get("gap", 0.1), rng, spec.get("lam_eff"))
    out.update(R=R, n_obs=dom.n_obs, m=dom.m, lam_real=dom.lam_real, area=float(dom.g.area(R)))
    I, J, W, st = build_graph(dom, rng, n_check=spec.get("n_check", 0))
    out.update(st); out["graph_s"] = time.perf_counter() - t0
    n = dom.m
    nc, lab = connected_components(csr_matrix((np.ones(len(I)), (I, J)), shape=(n, n)), directed=False)
    out["components"] = int(nc)
    ang = dom.t
    gx, gy = global_xy(dom.geom, dom.r, dom.t)
    labels = {}
    for mode in spec.get("modes", ["radial", "kd"]):
        if mode == "kd" or (mode == "radial"):
            t1 = time.perf_counter()
            L = CutLabels(n, I, J, W, ang, gx, gy, mode)
            labels[mode] = L
            out[mode] = dict(build_s=time.perf_counter() - t1, label_mean=float(L.label_len.mean()),
                             label_max=int(L.label_len.max()), entries=L.entries, height=L.height,
                             sep_by_depth=L.sep_by_depth())
    if "h2h" in spec.get("modes", []):
        from rgeo.oracles.h2h import H2HIndex
        t1 = time.perf_counter()
        adj = [dict() for _ in range(n)]
        for a, b, w in zip(I.tolist(), J.tolist(), W.tolist()):
            adj[a][b] = w; adj[b][a] = w
        big = np.argmax(np.bincount(lab))   # H2H needs one component
        keep = np.nonzero(lab == big)[0]
        new = -np.ones(n, dtype=np.int64); new[keep] = np.arange(len(keep))
        adj2 = [{int(new[b]): w for b, w in adj[a].items()} for a in keep.tolist()]
        H = H2HIndex(adj2)
        out["h2h"] = dict(build_s=time.perf_counter() - t1, label_mean=float((H.depth + 1).mean()),
                          width=int(H.width), height=int(H.h), n=len(keep))
    # validation: corner pairs against Dijkstra on G; free points against Dijkstra on G + p + q
    A = csr_matrix((np.concatenate([W, W]), (np.concatenate([I, J]), np.concatenate([J, I]))), shape=(n, n))
    nv = spec.get("n_val", 20)
    src = rng.integers(0, n, nv)
    Dref = dijkstra(A, directed=False, indices=src)
    err = {m_: 0.0 for m_ in labels}
    for i, u in enumerate(src.tolist()):
        for w in rng.integers(0, n, 10).tolist():
            ref = Dref[i, w]
            for m_, L in labels.items():
                got = L.query(u, w)
                e = 0.0 if (np.isinf(ref) and np.isinf(got)) else abs(got - ref) / max(ref, 1e-12)
                err[m_] = max(err[m_], e)
    out["val_corner_relerr"] = err
    npts = spec.get("n_points", 30)
    P = dom.sample_free(2 * npts, rng)
    vp_sizes, vp_rho, perr, lb_ok = [], [], {m_: 0.0 for m_ in labels}, True
    tq = {m_: 0.0 for m_ in labels}
    for i in range(npts):
        (r1, t1_), (r2, t2_) = P[2 * i], P[2 * i + 1]
        tv = time.perf_counter()
        Vp, _, rp = dom.view(r1, t1_); Vq, _, _ = dom.view(r2, t2_)
        dp = dom.dist(r1, t1_, dom.r[Vp], dom.t[Vp]); dq = dom.dist(r2, t2_, dom.r[Vq], dom.t[Vq])
        direct = float(dom.dist(r1, t1_, r2, t2_))
        free = dom.segment_free(r1, t1_, r2, t2_)
        tview = time.perf_counter() - tv
        vp_sizes.append(len(Vp)); vp_rho.append(rp)
        # reference: Dijkstra on G plus the two points
        rows = np.concatenate([I, J, np.full(len(Vp), n), Vp, np.full(len(Vq), n + 1), Vq])
        cols = np.concatenate([J, I, Vp, np.full(len(Vp), n), Vq, np.full(len(Vq), n + 1)])
        ww = np.concatenate([W, W, dp, dp, dq, dq])
        if free:
            rows = np.concatenate([rows, [n, n + 1]]); cols = np.concatenate([cols, [n + 1, n]])
            ww = np.concatenate([ww, [direct, direct]])
        ref = float(dijkstra(csr_matrix((ww, (rows, cols)), shape=(n + 2, n + 2)), directed=False, indices=n)[n + 1])
        lb_ok &= ref >= direct * (1 - 1e-12)
        for m_, L in labels.items():
            tl = time.perf_counter()
            got = min(direct if free else np.inf, L.query_points(Vp, dp, Vq, dq))
            tq[m_] += tview + time.perf_counter() - tl
            e = 0.0 if (np.isinf(ref) and np.isinf(got)) else abs(got - ref) / max(ref, 1e-12)
            perr[m_] = max(perr[m_], e)
    out.update(val_point_relerr=perr, point_vis_mean=float(np.mean(vp_sizes)), point_vis_max=int(np.max(vp_sizes)),
               point_rho_max=float(np.max(vp_rho)), lower_bound_ok=bool(lb_ok),
               point_query_ms={m_: 1e3 * v / npts for m_, v in tq.items()})
    out["total_s"] = time.perf_counter() - t0
    out["peak_mb"] = footprint_mb()[1]
    return out


# ---------------------------------------------------------------- sweep
BASE = dict(r_o=0.15, k=3, gap=0.02, lam=40.0, modes=["radial", "kd"], n_val=20, n_points=30)


def plan(quick=False):
    """Main series: near-jammed obstacles (lam_eff None), hyperbolic R and a flat disc of the same area.
    Density series: thinned to lam_eff (mean obstacle width 0.248, so lam_eff * width crosses 1 near lam_eff 4)."""
    specs = []
    Rs = [3.0, 4.0, 5.0] if quick else [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]
    seeds = [1] if quick else [1, 2, 3]
    for R in Rs:
        for s in seeds:
            if R >= 6.0 and s > 1:
                continue
            specs.append(dict(BASE, geom="hyp", R=R, seed=s, series="main"))
            specs.append(dict(BASE, geom="flat", R=0.0, match_area_R=R, seed=s, series="main"))
    lams = [2.0, 4.0] if quick else [2.0, 3.0, 4.0, 5.0]
    for le in lams:
        for R in ([3.0, 4.0] if quick else [3.0, 4.0, 5.0, 5.5]):
            specs.append(dict(BASE, geom="hyp", R=R, seed=1, lam_eff=le, series="density"))
            specs.append(dict(BASE, geom="flat", R=0.0, match_area_R=R, seed=1, lam_eff=le, series="density"))
    specs[0]["n_check"] = 300
    return specs


def key(spec):
    return json.dumps({k: spec[k] for k in sorted(spec) if k not in ("n_check",)}, sort_keys=True)


def load_results():
    if os.path.exists(OUT_JSON):
        return json.load(open(OUT_JSON))
    return {}


def cmd_run(a):
    res = load_results()
    for spec in plan(a.quick):
        kk = key(spec)
        if kk in res and res[kk].get("status") == "ok":
            continue
        t0 = time.time()
        try:
            p = subprocess.run([sys.executable, os.path.abspath(__file__), "case", json.dumps(spec)], capture_output=True,
                               text=True, timeout=a.timeout)
            out = json.loads(p.stdout.strip().splitlines()[-1]) if p.returncode == 0 else \
                dict(spec=spec, status="error", err=p.stderr[-2000:])
        except subprocess.TimeoutExpired:
            out = dict(spec=spec, status="timeout")
        res[kk] = out
        os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
        json.dump(res, open(OUT_JSON, "w"), indent=1)
        print(f"{spec['series']:7s} {spec['geom']:4s} R={spec.get('match_area_R') or spec['R']} lam_eff={spec.get('lam_eff')} "
              f"seed={spec['seed']}: {out['status']} m={out.get('m')} radial={out.get('radial', {}).get('label_mean')} "
              f"({time.time() - t0:.0f} s)", flush=True)


# ---------------------------------------------------------------- report
def ok_rows(res, series, geom, lam_eff=None):
    rows = [o for o in res.values() if o.get("status") == "ok" and o["spec"]["series"] == series
            and o["spec"]["geom"] == geom and o["spec"].get("lam_eff") == lam_eff]
    out = {}
    for o in rows:
        Rk = o["spec"].get("match_area_R") or o["spec"]["R"]
        out.setdefault(Rk, []).append(o)
    return dict(sorted(out.items()))


def agg(lst, f):
    v = np.array([f(o) for o in lst], dtype=float)
    return float(v.mean()), (float(v.std(ddof=1)) if len(v) > 1 else 0.0)


def fit2(x, y):
    """y = a + b x: relative RMS error."""
    A = np.column_stack([np.ones_like(x), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.sqrt(np.mean(((A @ c) - y) ** 2)) / np.mean(y)), c


def slope(x, y):
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])


def cmd_report(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = load_results()
    H, F = ok_rows(res, "main", "hyp"), ok_rows(res, "main", "flat")
    lines = ["| R | hyperbolic m | label radial | label kd | top sep. / cut length | flat m | flat label radial | flat label kd |"
             " flat top sep. / cut length | flat / hyp (kd) |", "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for R in H:
        if R not in F:
            continue
        h, f = H[R], F[R]
        hm, fm = agg(h, lambda o: o["m"]), agg(f, lambda o: o["m"])
        hr, hk = agg(h, lambda o: o["radial"]["label_mean"]), agg(h, lambda o: o["kd"]["label_mean"])
        fr, fk = agg(f, lambda o: o["radial"]["label_mean"]), agg(f, lambda o: o["kd"]["label_mean"])
        hs = agg(h, lambda o: o["kd"]["sep_by_depth"][0][2] / (2 * o["R"]))
        fs = agg(f, lambda o: o["kd"]["sep_by_depth"][0][2] / (2 * o["R"]))
        lines.append(f"| {R:g} | {hm[0]:,.0f} | {hr[0]:.0f} ± {hr[1]:.0f} | {hk[0]:.0f} ± {hk[1]:.0f} | {hs[0]:.1f} | "
                     f"{fm[0]:,.0f} | {fr[0]:.0f} ± {fr[1]:.0f} | {fk[0]:.0f} ± {fk[1]:.0f} | {fs[0]:.1f} | {fk[0] / hk[0]:.2f} |")
    main_tab = "\n".join(lines)
    # fits over the main series (means per R)
    fit_lines = ["| series | label | fit a + b log^2 m | fit a + b sqrt m | free exponent in m |", "|---|---|---:|---:|---:|"]
    for name, D in (("hyperbolic", H), ("flat", F)):
        if len(D) < 3:
            continue
        m = np.array([agg(D[R], lambda o: o["m"])[0] for R in D])
        for mode in ("radial", "kd"):
            y = np.array([agg(D[R], lambda o: o[mode]["label_mean"])[0] for R in D])
            e1, _ = fit2(np.log(m) ** 2, y); e2, _ = fit2(np.sqrt(m), y)
            fit_lines.append(f"| {name} | {mode} | {100 * e1:.1f}% | {100 * e2:.1f}% | {slope(m, y):.2f} |")
    fit_tab = "\n".join(fit_lines)
    # visibility in the main series
    vis_lines = ["| geometry | R | lam_real | mean degree of G | visible corners (mean / max) | longest edge | free point sees (mean / max) |"
                 " exact (max rel. err.) | visibility check |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for name, D in (("hyperbolic", H), ("flat", F)):
        for R, lst in D.items():
            o = lst[0]
            err = max(max(o["val_corner_relerr"].values()), max(o["val_point_relerr"].values()))
            chk = f"{o['vis_mismatch']}/{o['vis_checked']}" if "vis_checked" in o else ""
            vis_lines.append(f"| {name} | {R:g} | {o['lam_real']:.2f} | {o['mean_deg']:.1f} | {o['vis_mean']:.1f} / {o['vis_max']} | "
                             f"{o['max_edge']:.2f} | {o['point_vis_mean']:.1f} / {o['point_vis_max']} | {err:.1e} | {chk} |")
    vis_tab = "\n".join(vis_lines)
    # density series
    lams = sorted({o["spec"]["lam_eff"] for o in res.values() if o.get("status") == "ok" and o["spec"]["series"] == "density"})
    den_lines = ["| lam_eff | lam_eff x width | R | hyp m | hyp label kd | hyp longest edge | hyp point sees | flat m | flat label kd |"
                 " flat longest edge | flat / hyp |", "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for le in lams:
        Hd, Fd = ok_rows(res, "density", "hyp", le), ok_rows(res, "density", "flat", le)
        for R in Hd:
            h = Hd[R][0]; f = Fd.get(R, [None])[0]
            fl = (f"{f['m']:,} | {f['kd']['label_mean']:.0f} | {f['max_edge']:.2f} | {f['kd']['label_mean'] / h['kd']['label_mean']:.2f}"
                  if f else " | | | ")
            den_lines.append(f"| {le:g} | {le * WIDTH:.2f} | {R:g} | {h['m']:,} | {h['kd']['label_mean']:.0f} | {h['max_edge']:.2f} | "
                             f"{h['point_vis_mean']:.1f} | {fl} |")
    den_tab = "\n".join(den_lines)
    # plot
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for name, D, c in (("hyperbolic", H, "C0"), ("flat", F, "C3")):
        if not D:
            continue
        m = [agg(D[R], lambda o: o["m"])[0] for R in D]
        for mode, ls in (("kd", "-"), ("radial", "--")):
            y = [agg(D[R], lambda o: o[mode]["label_mean"])[0] for R in D]
            ax[0].plot(m, y, ls, marker="o", color=c, label=f"{name}, {mode} cuts")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("corners m"); ax[0].set_ylabel("mean label length")
    ax[0].set_title("near-jammed obstacles, matched area"); ax[0].legend(fontsize=8)
    for i, le in enumerate(lams):
        for name, geom, ls in (("hyp", "hyp", "-"), ("flat", "flat", ":")):
            D = ok_rows(res, "density", geom, le)
            if D:
                ax[1].plot([D[R][0]["m"] for R in D], [D[R][0]["kd"]["label_mean"] for R in D], ls, marker="o",
                           color=f"C{i}", label=f"{name}, lam={le:g}")
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("corners m"); ax[1].set_title("density series (kd cuts)")
    ax[1].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(OUT_PNG, dpi=130)
    text = open(OUT_MD).read() if os.path.exists(OUT_MD) else MD_TEMPLATE
    for tag, body in (("main", main_tab), ("fits", fit_tab), ("visibility", vis_tab), ("density", den_tab)):
        text = splice(text, tag, body)
    open(OUT_MD, "w").write(text)
    print(main_tab); print(fit_tab); print(vis_tab); print(den_tab)


WIDTH = 0.248   # mean width of a regular geodesic triangle of circumradius 0.15 (perimeter / pi, nearly flat at this size)


def splice(text, tag, body):
    b, e = f"<!-- BEGIN {tag} -->", f"<!-- END {tag} -->"
    i, j = text.index(b), text.index(e)
    return text[:i + len(b)] + "\n" + body + "\n" + text[j:]


MD_TEMPLATE = """# Exact two-point labels among obstacles: hyperbolic versus flat discs

Prototype for `notes/hyperbolic.md` section 7 (step 2). Script `bench/hyperbolic_exact.py` (`run`, then `report`), raw data
`results/hyperbolic_exact.json`, figure `results/hyperbolic_exact.png`. Tables between the markers are regenerated by `report`.

![labels versus m](../results/hyperbolic_exact.png)

## Conclusion

(written after the run)

## Main series: label length

<!-- BEGIN main -->
<!-- END main -->

## Fits

<!-- BEGIN fits -->
<!-- END fits -->

## Visibility and exactness

<!-- BEGIN visibility -->
<!-- END visibility -->

## Density series

<!-- BEGIN density -->
<!-- END density -->
"""


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("case"); c.add_argument("spec")
    r = sub.add_parser("run"); r.add_argument("--quick", action="store_true"); r.add_argument("--timeout", type=float, default=7200)
    sub.add_parser("report")
    a = ap.parse_args()
    if a.cmd == "case":
        print(json.dumps(run_case(json.loads(a.spec))))
    elif a.cmd == "run":
        cmd_run(a)
    else:
        cmd_report(a)


if __name__ == "__main__":
    main()
