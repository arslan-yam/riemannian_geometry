"""H2H tree-decomposition size on hyperbolic (H^2, curvature -1) versus Euclidean disks.

Claim under test: a balanced separator of a radius-R hyperbolic disk has O(R) = O(log n) points, a flat
one Theta(sqrt n), so H2H labels (ancestors per vertex) should grow like log^2 n against sqrt n.

  python3 bench/hyperbolic_h2h.py run     [--hyp-R 3,4,...] [--h 1.0]   # -> results/hyperbolic_h2h.json (resumable)
  python3 bench/hyperbolic_h2h.py report                                 # -> results/hyperbolic_h2h.png, notes/hyperbolic_h2h.md
  python3 bench/hyperbolic_h2h.py case '<json spec>'                     # one case, prints a json line (used by `run`)

Every case is a child process (clean peak memory, killed after --timeout seconds). Points are stored as
(r, theta); the hyperbolic distance uses the cancellation-free form
    sinh^2(d/2) = sinh^2((r1-r2)/2) + sinh r1 sinh r2 sin^2((t1-t2)/2),
which is exact for any R in double precision. Obstacle tests use the Klein model (geodesics = chords) after an
isometry that moves the obstacle centre to the origin, so a hyperbolic disk becomes a Euclidean disk of radius tanh(rho).
"""
import argparse, json, os, resource, subprocess, sys, time
from itertools import chain
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sp_dijkstra, connected_components

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT_JSON = os.path.join(ROOT, "results", "hyperbolic_h2h.json")
OUT_PNG = os.path.join(ROOT, "results", "hyperbolic_h2h.png")
OUT_MD = os.path.join(ROOT, "notes", "hyperbolic_h2h.md")

TWO_PI = 2.0 * np.pi
DMIN_OVER_H = 0.85      # Poisson-disk minimum separation = 0.85 h gives about one point per h^2 of area
KAPPA = 3.0             # dart-throwing candidates per (1/dmin^2) of area
MEAN_DEGREE = 8.0
RHO_MAX_OVER_H = 3.0    # candidate edge radius; the degree-8 radius must fall below it
MEM_LIMIT = 2.4e9       # abort a case whose H2H arrays are predicted to exceed this (bytes); real RSS is ~0.75x of it


# ---------------------------------------------------------------- geometry
class Flat:
    name = "flat"

    @staticmethod
    def area(R):
        return np.pi * R * R

    @staticmethod
    def sample(R, N, rng):
        return R * np.sqrt(rng.random(N)), TWO_PI * rng.random(N)

    @staticmethod
    def dist(r1, t1, r2, t2):
        return np.hypot(r1 * np.cos(t1) - r2 * np.cos(t2), r1 * np.sin(t1) - r2 * np.sin(t2))

    @staticmethod
    def kd_xy(r, t):
        return np.column_stack([r * np.cos(t), r * np.sin(t)])

    @staticmethod
    def kd_query(r, t, rho):  # (centres, radii) of Euclidean disks that contain the rho-ball of each point
        return Flat.kd_xy(r, t), np.full(len(r), rho)


class Hyp:
    name = "hyperbolic"

    @staticmethod
    def area(R):
        return TWO_PI * (np.cosh(R) - 1.0)

    @staticmethod
    def sample(R, N, rng):  # radial density proportional to sinh r
        return np.arccosh(1.0 + rng.random(N) * (np.cosh(R) - 1.0)), TWO_PI * rng.random(N)

    @staticmethod
    def dist(r1, t1, r2, t2):
        a = np.sinh(0.5 * (r1 - r2))
        b = np.sinh(r1) * np.sinh(r2) * np.sin(0.5 * (t1 - t2)) ** 2
        return 2.0 * np.arcsinh(np.sqrt(a * a + b))

    @staticmethod
    def kd_xy(r, t):  # Poincare model, used only to find candidate neighbours
        s = np.tanh(0.5 * r)
        return np.column_stack([s * np.cos(t), s * np.sin(t)])

    @staticmethod
    def kd_query(r, t, rho):
        # in the Poincare model the hyperbolic rho-ball about a point at Euclidean radius s is exactly a Euclidean
        # disk: centre u (1-T^2)/(1-T^2 s^2), radius T (1-s^2)/(1-T^2 s^2), with T = tanh(rho/2)
        s = np.tanh(0.5 * r)
        one_m_s2 = 1.0 / np.cosh(0.5 * r) ** 2
        tt = np.tanh(0.5 * rho)
        den = 1.0 - tt * tt * s * s
        k = (1.0 - tt * tt) / den * s
        return np.column_stack([k * np.cos(t), k * np.sin(t)]), tt * one_m_s2 / den


GEOM = {"flat": Flat, "hyp": Hyp}


# ---------------------------------------------------------------- sampling and graph
def pairs_within(g, r, t, rho, chunk=2000):
    """All pairs i<j with geometric distance <= rho (arrays I, J, D)."""
    tree = cKDTree(g.kd_xy(r, t))
    xy, B = g.kd_query(r, t, rho)
    B = B * 1.0001 + 1e-12
    Is, Js, Ds = [], [], []
    for a in range(0, len(r), chunk):
        b = min(len(r), a + chunk)
        lists = tree.query_ball_point(xy[a:b], B[a:b], return_sorted=False)
        lens = np.fromiter((len(l) for l in lists), np.int64, b - a)
        flat = np.fromiter(chain.from_iterable(lists), np.int64, int(lens.sum()))
        ii = np.repeat(np.arange(a, b), lens)
        m = flat > ii
        ii, jj = ii[m], flat[m]
        d = g.dist(r[ii], t[ii], r[jj], t[jj])
        m = d <= rho
        Is.append(ii[m]); Js.append(jj[m]); Ds.append(d[m])
    return np.concatenate(Is), np.concatenate(Js), np.concatenate(Ds)


def poisson_disk(g, R, dmin, rng):
    """Dart throwing: random area-uniform candidates in random order, accept if no accepted point within dmin."""
    N = int(KAPPA * g.area(R) / dmin ** 2)
    r, t = g.sample(R, N, rng)
    I, J, _ = pairs_within(g, r, t, dmin)
    src = np.concatenate([I, J]); dst = np.concatenate([J, I])
    order = np.argsort(src, kind="stable")
    dst = dst[order]
    ptr = np.zeros(N + 1, dtype=np.int64)
    np.cumsum(np.bincount(src, minlength=N), out=ptr[1:])
    removed = np.zeros(N, dtype=bool)
    acc = []
    for i in range(N):
        if not removed[i]:
            acc.append(i)
            removed[dst[ptr[i]:ptr[i + 1]]] = True
    acc = np.array(acc)
    r, t = r[acc], t[acc]
    o = np.argsort(r, kind="stable")  # by radius: a prefix is a smaller concentric disk
    return r[o], t[o]


def obstacles_hyp(r, t, I, J, R, rho_o, frac, rng, Lmax):
    """Hyperbolic disks of radius rho_o covering about `frac` of the area, centres area-uniform in the disk.
    Returns (alive mask over points, kept mask over edges, number of obstacles). Only the points within rho_o+Lmax of an
    obstacle and their incident edges are examined (an edge meeting the disk has an endpoint that close)."""
    n_obs = max(1, int(round(frac * Hyp.area(R) / Hyp.area(rho_o))))
    rc_all, tc_all = Hyp.sample(R, n_obs, rng)
    n = len(r)
    alive = np.ones(n, dtype=bool)
    cross = np.zeros(len(I), dtype=bool)
    tree = cKDTree(Hyp.kd_xy(r, t))
    cxy, cB = Hyp.kd_query(rc_all, tc_all, rho_o + Lmax)
    # CSR of incident edges
    ends = np.concatenate([I, J])
    eid = np.concatenate([np.arange(len(I)), np.arange(len(I))])
    o = np.argsort(ends, kind="stable")
    eid = eid[o]
    ptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(np.bincount(ends, minlength=n), out=ptr[1:])
    sr, Z = np.sinh(r), np.cosh(r)
    k_rad = np.tanh(rho_o)
    for k in range(n_obs):
        rc, tc = rc_all[k], tc_all[k]
        near = np.array(tree.query_ball_point(cxy[k], cB[k] * 1.0001 + 1e-12), dtype=np.int64)
        if len(near) == 0:
            continue
        alive[near[Hyp.dist(r[near], t[near], rc, tc) < rho_o]] = False
        cand = np.unique(np.concatenate([eid[ptr[v]:ptr[v + 1]] for v in near.tolist()]))
        cand = cand[~cross[cand]]
        if len(cand) == 0:
            continue
        e = np.concatenate([I[cand], J[cand]])
        x1 = sr[e] * np.cos(t[e] - tc)
        y1 = sr[e] * np.sin(t[e] - tc)
        z = -x1 * np.sinh(rc) + Z[e] * np.cosh(rc)           # boost the centre to the origin
        x = x1 * np.cosh(rc) - Z[e] * np.sinh(rc)
        kx, ky = x / z, y1 / z                                 # Klein coordinates: geodesics are chords
        m = len(cand)
        ax, ay, bx, by = kx[:m], ky[:m], kx[m:], ky[m:]
        dx, dy = bx - ax, by - ay
        L2 = np.maximum(dx * dx + dy * dy, 1e-300)
        s = np.clip(-(ax * dx + ay * dy) / L2, 0.0, 1.0)
        px, py = ax + s * dx, ay + s * dy
        cross[cand[px * px + py * py < k_rad * k_rad]] = True
    kept = alive[I] & alive[J] & ~cross
    return alive, kept, n_obs


def diameter_cut(theta, I, J, n_dir=16):
    """Mean over n_dir geodesics through the centre of the smaller one-sided vertex cut that separates the two halves."""
    sizes, bal = [], []
    for k in range(n_dir):
        side = np.sin(theta - np.pi * k / n_dir) > 0
        c = side[I] != side[J]
        left = np.unique(np.concatenate([I[c & side[I]], J[c & side[J]]]))
        right = np.unique(np.concatenate([I[c & ~side[I]], J[c & ~side[J]]]))
        sizes.append(min(len(left), len(right)))
        bal.append(min(side.sum(), (~side).sum()) / len(theta))
    return float(np.mean(sizes)), float(np.mean(bal))


# ---------------------------------------------------------------- one case
def footprint_mb():
    """(current, lifetime peak) physical footprint in MB (macOS proc_pid_rusage; includes compressed pages)."""
    try:
        import ctypes
        buf = (ctypes.c_uint64 * 64)()
        if ctypes.CDLL("/usr/lib/libproc.dylib").proc_pid_rusage(os.getpid(), 4, ctypes.byref(buf)) == 0:
            return buf[9] / 1e6, buf[30] / 1e6
    except OSError:
        pass
    return float("nan"), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def run_case(spec):
    from rgeo.oracles import h2h as H
    g = GEOM[spec["geom"]]
    h = spec["h"]
    rng = np.random.default_rng(spec["seed"])
    out = dict(spec=spec, status="ok")
    t0 = time.perf_counter()
    if spec["geom"] == "hyp":
        R = spec["R"]
        r, t = poisson_disk(g, R, DMIN_OVER_H * h, rng)
    else:  # sample a slightly larger disk, keep the n_target points closest to the centre
        n_t = spec["n"]
        R0 = h * np.sqrt(1.15 * n_t / (0.8 * np.pi))
        while True:
            r, t = poisson_disk(g, R0, DMIN_OVER_H * h, rng)
            if len(r) >= n_t:
                break
            R0 *= 1.1
        r, t = r[:n_t], t[:n_t]
        R = float(r[-1])
    out["R"] = float(R)
    out["peak_mb_sampling"] = footprint_mb()[1]
    n0 = len(r)
    I, J, D = pairs_within(g, r, t, RHO_MAX_OVER_H * h)
    m_edge = int(MEAN_DEGREE * n0 / 2)
    if len(D) < m_edge:
        raise RuntimeError("candidate radius too small for the target degree")
    thr = float(np.partition(D, m_edge - 1)[m_edge - 1])
    if thr >= RHO_MAX_OVER_H * h * 0.999:
        raise RuntimeError("degree-8 radius reached the candidate radius")
    sel = D <= thr
    I, J, D = I[sel], J[sel], D[sel]
    out.update(n_sampled=n0, conn_radius_over_h=thr / h, density_per_h2=float(n0 * h * h / g.area(R)),
               min_sep_over_h=float(D.min() / h))
    alive = np.ones(n0, dtype=bool)
    if spec.get("obstacle_rho"):
        Lmax = float(D.max())
        alive, kept, n_obs = obstacles_hyp(r, t, I, J, R, spec["obstacle_rho"], spec["obstacle_frac"], rng, Lmax)
        I, J, D = I[kept], J[kept], D[kept]
        out.update(n_obstacles=n_obs, removed_frac=float(1 - alive.mean()), dropped_edge_frac=float(1 - kept.mean()))
    # largest connected component, reindexed
    idx_old = np.nonzero(alive)[0]
    new = -np.ones(n0, dtype=np.int64)
    new[idx_old] = np.arange(len(idx_old))
    I, J = new[I], new[J]
    nn = len(idx_old)
    nc, lab = connected_components(csr_matrix((np.ones(len(I)), (I, J)), shape=(nn, nn)), directed=False)
    big = np.argmax(np.bincount(lab))
    keep = lab == big
    new2 = -np.ones(nn, dtype=np.int64)
    new2[keep] = np.arange(int(keep.sum()))
    em = keep[I] & keep[J]
    I, J, D = new2[I[em]], new2[J[em]], D[em]
    n = int(keep.sum())
    theta = t[idx_old][keep]
    out.update(n=n, components=int(nc), edges=len(I), mean_degree=2.0 * len(I) / n,
               sample_s=time.perf_counter() - t0)
    out["peak_mb_graph"] = footprint_mb()[1]
    cut, bal = diameter_cut(theta, I, J)
    out.update(diam_cut=cut, diam_balance=bal)

    adj = [dict() for _ in range(n)]
    for a, b, w in zip(I.tolist(), J.tolist(), D.tolist()):
        adj[a][b] = w
        adj[b][a] = w

    tm = {}
    orig_td = H.tree_decomposition

    def td(adj_):
        t1 = time.perf_counter()
        pi, bags, parent, root, middle = orig_td(adj_)
        tm["td"] = time.perf_counter() - t1
        dep = np.zeros(len(adj_), dtype=np.int64)
        for v in np.argsort(-np.asarray(pi)):  # parents are eliminated later: root first
            if parent[v] is not None:
                dep[v] = dep[parent[v]] + 1
        tm["fill_edges"] = len(middle)
        est = len(adj_) * (int(dep.max()) + 1) * 18.0 + 250.0 * len(middle)
        tm["est_bytes"] = est
        if est > MEM_LIMIT:
            raise MemoryError(f"predicted {est / 1e9:.2f} GB")
        return pi, bags, parent, root, middle

    H.tree_decomposition = td
    out["peak_mb_adj"] = footprint_mb()[1]
    try:
        t1 = time.perf_counter()
        idx = H.H2HIndex(adj)
        tm["build"] = time.perf_counter() - t1
    except MemoryError as e:
        H.tree_decomposition = orig_td
        out.update(status="memory", note=str(e), td_s=tm.get("td"), est_bytes=tm.get("est_bytes"),
                   fill_edges=tm.get("fill_edges"))
        return out
    finally:
        H.tree_decomposition = orig_td
    depth = idx.depth.astype(np.int64)
    bag = np.array([len(x) for x in idx.Xv])
    out.update(td_s=tm["td"], label_s=tm["build"] - tm["td"], build_s=tm["build"], fill_edges=tm["fill_edges"],
               max_bag=int(idx.width), mean_bag=float(bag.mean()),
               anc_mean=float(depth.mean() + 1), anc_max=int(depth.max() + 1), height=int(idx.h),
               label_entries=int(idx.label_entries), pos_entries=int(sum(len(p) for p in idx.pos)))
    # query time, then validation against Dijkstra
    rng2 = np.random.default_rng(spec["seed"] + 1)
    qs = rng2.integers(n, size=(2000, 2))
    t1 = time.perf_counter()
    for s_, t_ in qs.tolist():
        idx.query(s_, t_)
    out["query_us"] = (time.perf_counter() - t1) / len(qs) * 1e6
    G = csr_matrix((np.concatenate([D, D]), (np.concatenate([I, J]), np.concatenate([J, I]))), shape=(n, n))
    srcs = rng2.choice(n, size=min(40, n), replace=False)
    ref = sp_dijkstra(G, directed=False, indices=srcs)
    worst, npairs = 0.0, 0
    for k, s_ in enumerate(srcs.tolist()):
        for t_ in rng2.integers(n, size=8).tolist():
            if t_ == s_:
                continue
            d = idx.query(s_, t_)
            worst = max(worst, abs(d - ref[k, t_]) / max(ref[k, t_], 1e-12))
            npairs += 1
    out.update(val_pairs=npairs, val_max_rel_err=worst)
    out["rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # bytes on macOS
    out["peak_mb"] = footprint_mb()[1]
    return out


# ---------------------------------------------------------------- driver
def key(spec):
    return json.dumps({k: spec[k] for k in sorted(spec)}, sort_keys=True)


def load_results():
    if os.path.exists(OUT_JSON):
        return json.load(open(OUT_JSON))
    return []


def save_results(res):
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    json.dump(res, open(OUT_JSON, "w"), indent=1)


def launch(spec, timeout):
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, os.path.abspath(__file__), "case", json.dumps(spec)],
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return dict(spec=spec, status="timeout", wall_s=time.time() - t0)
    if p.returncode != 0:
        return dict(spec=spec, status="error", note=p.stderr[-800:], wall_s=time.time() - t0)
    res = json.loads(p.stdout.strip().splitlines()[-1])
    res["wall_s"] = time.time() - t0
    return res


def cmd_run(a):
    os.nice(10)
    res = load_results()
    done = {key(r["spec"]): r for r in res}

    def go(spec):
        k = key(spec)
        if k in done and not a.force:
            return done[k]
        r = launch(spec, a.timeout)
        print(f"{spec['geom']:4s} obs={spec.get('obstacle_rho', 0)} R={r.get('R', spec.get('R', '-'))} "
              f"status={r['status']} n={r.get('n')} maxbag={r.get('max_bag')} anc_mean={r.get('anc_mean')} "
              f"height={r.get('height')} build={r.get('build_s')} wall={r['wall_s']:.0f}s", flush=True)
        res[:] = [x for x in res if key(x["spec"]) != k] + [r]
        done[k] = r
        save_results(res)
        return r

    Rs = [float(x) for x in a.hyp_R.split(",")]
    stop = set()  # (kind, rep) pairs that already failed: larger sizes are not attempted

    def spec_of(kind, R, rep, **kw):
        off = {"hyp": 7, "flat": 8, "obs": 9}[kind]
        sp = dict(geom="flat" if kind == "flat" else "hyp", h=a.h, seed=int(R * 100) + off + 10000 * rep, **kw)
        if kind != "flat":
            sp["R"] = R
        if rep:
            sp["rep"] = rep
        return sp

    for R in Rs:
        for rep in a.reps:
            rh = None
            if ("hyp", rep) not in stop:
                rh = go(spec_of("hyp", R, rep))
                if rh["status"] != "ok":
                    stop.add(("hyp", rep))
                    rh = None
            if rh is not None and ("flat", rep) not in stop:
                sp = spec_of("flat", R, rep, n=rh["n_sampled"])
                if rep == 0 or R <= a.flat_rep_max_R:
                    if go(sp)["status"] != "ok":
                        stop.add(("flat", rep))
            for rho in a.obstacles:
                if R >= a.obs_min_R and (rho, rep) not in stop:
                    if go(spec_of("obs", R, rep, obstacle_rho=rho, obstacle_frac=a.obs_frac))["status"] != "ok":
                        stop.add((rho, rep))


# ---------------------------------------------------------------- report: fits, tables, plot
SERIES = [("hyp", None, "hyperbolic H²", "#2a78d6", "o"),
          ("flat", None, "flat R²", "#eb6834", "s"),
          ("hyp", 1.0, "H² + obstacles ρ=1", "#1baf7a", "^"),
          ("hyp", 2.0, "H² + obstacles ρ=2", "#eda100", "D")]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"


NUMERIC = ["R", "n", "n_sampled", "edges", "max_bag", "mean_bag", "anc_mean", "anc_max", "label_entries", "build_s", "td_s", "label_s",
           "query_us", "diam_cut", "peak_mb", "mean_degree", "conn_radius_over_h", "density_per_h2"]


def group_R(r):
    """nominal hyperbolic R of a case (flat cases are matched to a hyperbolic R), recovered from the seed"""
    sp = r["spec"]
    off = 9 if sp.get("obstacle_rho") else (7 if sp["geom"] == "hyp" else 8)
    return round(((sp["seed"] % 10000) - off) / 100.0, 2)


def raw_rows(res, geom, obs, h=1.0):
    return [r for r in res if r["status"] == "ok" and r["spec"]["geom"] == geom and r["spec"].get("obstacle_rho") == obs and r["spec"]["h"] == h]


def series_rows(res, geom, obs, h=1.0):
    """One row per nominal R: the mean over replicate samples (`<key>_sd` = sample std, `reps` = number of replicates)."""
    groups = {}
    for r in raw_rows(res, geom, obs, h):
        groups.setdefault(group_R(r), []).append(r)
    out = []
    for g, rs in groups.items():
        agg = dict(grp=g, reps=len(rs), spec=rs[0]["spec"], val_max_rel_err=max(r["val_max_rel_err"] for r in rs))
        for k in NUMERIC:
            v = np.array([r[k] for r in rs], float)
            agg[k] = float(v.mean())
            agg[k + "_sd"] = float(v.std(ddof=1)) if len(v) > 1 else 0.0
        out.append(agg)
    return sorted(out, key=lambda r: r["n"])


def ols(x, y):
    A = np.column_stack([np.ones_like(x), x])
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = A @ c
    r2 = 1.0 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return c, float(r2), float(np.sqrt(np.mean(((pred - y) / y) ** 2)))


MODELS = {"log n": np.log, "log² n": lambda n: np.log(n) ** 2, "√n": np.sqrt}
METRICS = [("anc_mean", "mean ancestors (label length)", ["log² n", "√n"]),
           ("anc_max", "max ancestors = tree height", ["log² n", "√n"]),
           ("max_bag", "max bag |X(v)|", ["log n", "√n"]),
           ("diam_cut", "diameter-geodesic vertex cut", ["log n", "√n"])]


def fit_table(res, fit_min_n):
    lines = ["| series | metric | model y = a + b·f(n) | a | b | R² | rel. RMS error |", "|---|---|---|---:|---:|---:|---:|"]
    expo = ["| series | metric | exponent α in y ∝ nᵅ | R² (log-log) | exponent β in y ∝ (log n)ᵝ | R² |", "|---|---|---:|---:|---:|---:|"]
    for geom, obs, label, _, _ in SERIES:
        rows = [r for r in series_rows(res, geom, obs) if r["n"] >= fit_min_n]
        if len(rows) < 3:
            continue
        n = np.array([r["n"] for r in rows], float)
        for key_, mname, models in METRICS:
            y = np.array([r[key_] for r in rows], float)
            for mod in models:
                (a, b), r2, rr = ols(MODELS[mod](n), y)
                lines.append(f"| {label} | {mname} | {mod} | {a:.2f} | {b:.3f} | {r2:.4f} | {rr * 100:.1f}% |")
            (_, al), r2a, _ = ols(np.log(n), np.log(y))
            (_, be), r2b, _ = ols(np.log(np.log(n)), np.log(y))
            expo.append(f"| {label} | {mname} | {al:.3f} | {r2a:.4f} | {be:.2f} | {r2b:.4f} |")
    return "\n".join(lines), "\n".join(expo)


def pm(r, k, fmt):
    v = f"{r[k]:{fmt}}"
    return v + (f" ± {r[k + '_sd']:{fmt}}" if r["reps"] > 1 else "")


def size_tables(res):
    out = []
    for geom, obs, label, _, _ in SERIES:
        rows = series_rows(res, geom, obs)
        if not rows:
            continue
        out.append(f"**{label}**\n\n| R | reps | n | max bag | mean bag | mean anc. | max anc. (height) | label entries | build s (TD + labels) | query µs | diam. cut | val. err |\n|" + "---:|" * 12)
        for r in rows:
            out.append(f"| {r['R']:.2f} | {r['reps']} | {r['n']:,.0f} | {pm(r, 'max_bag', '.0f')} | {r['mean_bag']:.1f} | {pm(r, 'anc_mean', '.1f')} | {pm(r, 'anc_max', '.0f')} | "
                       f"{r['label_entries']:,.0f} | {r['build_s']:.1f} ({r['td_s']:.1f} + {r['label_s']:.1f}) | {r['query_us']:.1f} | {pm(r, 'diam_cut', '.1f')} | {r['val_max_rel_err']:.0e} |")
        out.append("")
    skipped = [r for r in res if r["status"] != "ok" and r["spec"]["h"] == 1.0]
    if skipped:
        out.append("Cases that did not complete:\n")
        for r in skipped:
            sp = r["spec"]
            out.append(f"- {sp['geom']} {'obstacles ρ=%g ' % sp['obstacle_rho'] if sp.get('obstacle_rho') else ''}"
                       f"{'R=%g' % sp['R'] if 'R' in sp else 'n=%d' % sp['n']}: {r['status']} {r.get('note', '')[:120]}")
    return "\n".join(out)


def matched_table(res):
    hyp = {r["grp"]: r for r in series_rows(res, "hyp", None)}
    flat = series_rows(res, "flat", None)
    lines = ["| hyp R | n | mean anc. hyp | mean anc. flat | flat/hyp | max bag hyp | max bag flat | flat/hyp | diam. cut hyp | diam. cut flat | flat/hyp |", "|" + "---:|" * 11]
    first = {}
    for f in flat:
        h = hyp.get(f["grp"])
        if h is None:
            continue
        ra, rb, rc = f["anc_mean"] / h["anc_mean"], f["max_bag"] / h["max_bag"], f["diam_cut"] / h["diam_cut"]
        for th in (1.5, 2.0, 3.0, 5.0):
            if ra >= th and th not in first:
                first[th] = (h["grp"], h["n"], ra)
        lines.append(f"| {h['grp']:.1f} | {h['n']:,.0f} | {h['anc_mean']:.1f} | {f['anc_mean']:.1f} | {ra:.2f} | {h['max_bag']:.1f} | {f['max_bag']:.1f} | {rb:.2f} | {h['diam_cut']:.1f} | {f['diam_cut']:.1f} | {rc:.2f} |")
    return "\n".join(lines), first


def local_slopes(res):
    """log-log slope between consecutive sizes of each series, reported against the geometric-mean n of the step"""
    lines = ["| n (geometric mean of the step) | hyp: mean anc. | flat: mean anc. | hyp: max bag | flat: max bag | hyp: diam. cut | flat: diam. cut |", "|---:|---:|---:|---:|---:|---:|---:|"]
    hyp = {r["grp"]: r for r in series_rows(res, "hyp", None)}
    pairs = [(hyp[f["grp"]], f) for f in series_rows(res, "flat", None) if f["grp"] in hyp]
    for (h0, f0), (h1, f1) in zip(pairs[:-1], pairs[1:]):
        if h0["n"] < 300:
            continue
        def sl(a, b, k):
            return np.log(b[k] / a[k]) / np.log(b["n"] / a["n"])
        lines.append(f"| {np.sqrt(h0['n'] * h1['n']):,.0f} | {sl(h0, h1, 'anc_mean'):.2f} | {sl(f0, f1, 'anc_mean'):.2f} | {sl(h0, h1, 'max_bag'):.2f} | {sl(f0, f1, 'max_bag'):.2f} | "
                     f"{sl(h0, h1, 'diam_cut'):.2f} | {sl(f0, f1, 'diam_cut'):.2f} |")
    return "\n".join(lines)


def make_plot(res, path):
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
    panels = [("diam_cut", "Diameter-geodesic vertex cut (balanced separator)", "log n", np.log),
              ("max_bag", "Max bag size  max|X(v)|  (width)", "log n", np.log),
              ("anc_mean", "Label length  (mean ancestors per vertex)", "log² n", lambda n: np.log(n) ** 2)]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0), facecolor=SURFACE)
    for ax, (key_, title, hyp_ref, hyp_fn) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        for geom, obs, label, col, mk in SERIES:
            rows = series_rows(res, geom, obs)
            if not rows:
                continue
            n = np.array([r["n"] for r in rows], float)
            y = np.array([r[key_] for r in rows], float)
            raw = raw_rows(res, geom, obs)
            ax.plot([r["n"] for r in raw], [r[key_] for r in raw], ".", color=col, ms=3, alpha=0.35, zorder=2)
            ax.plot(n, y, "-", color=col, lw=1.8, marker=mk, ms=6.5, mfc=col, mec=SURFACE, mew=1.0, label=label, zorder=3)
        # reference guides through the plain-hyperbolic / flat point at the largest common n
        hyp = series_rows(res, "hyp", None)
        flat = series_rows(res, "flat", None)
        xs = np.logspace(np.log10(min(r["n"] for r in hyp)), np.log10(max(r["n"] for r in hyp)), 100)
        if flat:
            f = flat[-1]
            ax.plot(xs, f[key_] * np.sqrt(xs / f["n"]), "--", color="#898781", lw=1.0, zorder=1)
            ax.text(xs[-1] * 0.93, f[key_] * np.sqrt(xs[-1] / f["n"]) * 1.0, "∝ √n", color=INK2, ha="right", va="bottom", fontsize=9)
        if hyp:
            hr = [r for r in hyp if r["n"] >= 500]
            h0 = hr[len(hr) // 2] if hr else hyp[-1]
            ax.plot(xs, h0[key_] * hyp_fn(xs) / hyp_fn(h0["n"]), ":", color="#898781", lw=1.2, zorder=1)
            ax.text(xs[-1], h0[key_] * hyp_fn(xs[-1]) / hyp_fn(h0["n"]) * 1.08, f"∝ {hyp_ref}", color=INK2, ha="right", va="bottom", fontsize=9)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("n (vertices)"); ax.set_title(title, fontsize=10.5, color=INK, loc="left")
        ax.grid(True, which="major", color=GRID, lw=0.7); ax.grid(True, which="minor", color=GRID, lw=0.3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, fontsize=10.5, bbox_to_anchor=(0.5, 1.0))
    fig.text(0.5, -0.01, "h = 1, mean degree 8. Markers: mean over up to 3 random samples (faint dots: single samples); obstacle series: size of the giant component. "
             "Dashed: ∝ √n through the largest flat point. Dotted: ∝ log n (left, middle) or ∝ log² n (right) through a mid-range hyperbolic point.",
             ha="center", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    fig.savefig(path, dpi=150, facecolor=SURFACE, bbox_inches="tight")


def normalised_table(res):
    lines = ["| series | R | n | diam. cut / R | max bag / R | mean anc. / R² | height / R² |", "|---|---:|---:|---:|---:|---:|---:|"]
    for geom, obs, label, _, _ in SERIES:
        if geom != "hyp":
            continue
        for r in series_rows(res, geom, obs):
            if r["n"] >= 300:
                R = r["grp"]
                lines.append(f"| {label} | {R:.1f} | {r['n']:,.0f} | {r['diam_cut'] / R:.2f} | {r['max_bag'] / R:.2f} | {r['anc_mean'] / R ** 2:.2f} | {r['anc_max'] / R ** 2:.2f} |")
    lines += ["", "| series | n | diam. cut / √n | max bag / √n | mean anc. / √n | height / √n |", "|---|---:|---:|---:|---:|---:|"]
    for r in series_rows(res, "flat", None):
        if r["n"] >= 300:
            q = np.sqrt(r["n"])
            lines.append(f"| flat R² | {r['n']:,.0f} | {r['diam_cut'] / q:.2f} | {r['max_bag'] / q:.2f} | {r['anc_mean'] / q:.2f} | {r['anc_max'] / q:.2f} |")
    return "\n".join(lines)


def spacing_table(res):
    lines = ["| h | R | n | diam. cut | diam. cut · h / R | max bag | mean anc. | mean anc. · h² / R² | flat mean anc. at the same n |", "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for h in (0.5, 1.0, 2.0):
        flat = {r["grp"]: r for r in series_rows(res, "flat", None, h)}
        for r in series_rows(res, "hyp", None, h):
            if r["n"] >= 300:
                R = r["grp"]
                f = flat.get(R)
                lines.append(f"| {h:g} | {R:.1f} | {r['n']:,.0f} | {pm(r, 'diam_cut', '.1f')} | {r['diam_cut'] * h / R:.2f} | {pm(r, 'max_bag', '.0f')} | {pm(r, 'anc_mean', '.1f')} | "
                             f"{r['anc_mean'] * h * h / R ** 2:.2f} | {('%.1f' % f['anc_mean']) if f else 'not run'} |")
    return "\n".join(lines)


def splice(text, tag, body):
    """Replace the auto block `<!-- BEGIN tag -->...<!-- END tag -->` in text (append it when absent)."""
    b, e = f"<!-- BEGIN {tag} -->", f"<!-- END {tag} -->"
    block = f"{b}\n{body}\n{e}"
    if b in text and e in text:
        return text[:text.index(b)] + block + text[text.index(e) + len(e):]
    return text + "\n" + block + "\n"


def cmd_report(a):
    res = load_results()
    make_plot(res, OUT_PNG)
    fits, expo = fit_table(res, a.fit_min_n)
    matched, first = matched_table(res)
    blocks = {"sizes": size_tables(res), "matched": matched, "slopes": local_slopes(res), "fits": fits, "expo": expo,
              "normalised": normalised_table(res), "spacing": spacing_table(res)}
    text = open(OUT_MD).read() if os.path.exists(OUT_MD) else ""
    for tag, body in blocks.items():
        text = splice(text, tag, body)
    open(OUT_MD, "w").write(text)
    print("wrote", OUT_PNG, "and", OUT_MD)
    print("first R with flat/hyp mean-ancestors ratio >= threshold:", {k: (round(v[0], 1), v[1], round(v[2], 2)) for k, v in first.items()})


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("run")
    pr.add_argument("--hyp-R", default="3,4,5,6,7,8,9,10")
    pr.add_argument("--h", type=float, default=1.0)
    pr.add_argument("--obstacles", type=lambda s: [float(x) for x in s.split(",") if x], default=[1.0, 2.0])
    pr.add_argument("--obs-frac", type=float, default=0.2)
    pr.add_argument("--obs-min-R", type=float, default=5.0)
    pr.add_argument("--timeout", type=float, default=600.0)
    pr.add_argument("--reps", type=lambda s: [int(x) for x in s.split(",")], default=[0], help="replicate indices (new random samples)")
    pr.add_argument("--flat-rep-max-R", type=float, default=10.0, help="largest matched hyperbolic R at which the flat replicates run")
    pr.add_argument("--force", action="store_true")
    pc = sub.add_parser("case")
    pc.add_argument("spec")
    prep = sub.add_parser("report")
    prep.add_argument("--fit-min-n", type=int, default=500)
    a = ap.parse_args()
    if a.cmd == "case":
        print(json.dumps(run_case(json.loads(a.spec))))
    elif a.cmd == "run":
        cmd_run(a)
    else:
        cmd_report(a)


if __name__ == "__main__":
    main()
