#benchmark on real terrain meshes (the datasets of the SE-Oracle / EAR-Oracle papers): every
#oracle is built on the mesh graph (vertices + mesh edges, exact 3-d edge lengths) and queried
#vertex-to-vertex; the ground truth is the exact polyhedral geodesic (MMP, via pygeodesic).
#the same vertex pairs are written in the competitors' query format so their C++ harness
#(FixedS / UnfixedS / KAlgo / SE / EAR / MMP) can be run on exactly the same queries
import os, sys, json, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path: #runnable as `python3 bench/terrain_bench.py` without installing the package
    sys.path.insert(0, ROOT)
from rgeo.terrain import TerrainMesh, read_off, repair_height_field, ExactTerrain, TerrainA2A
from rgeo.graphs import dijkstra, get_path
from rgeo.oracles import H2HIndex, SeparatorIndex, SeparatorCutIndex
from rgeo.refine import refine_polyline
from rgeo.oracles.certified import CertifiedCutIndex
from rgeo.oracles import cquery
from rgeo.terrain.submesh import FaceGrid
from rgeo.terrain.steiner import SteinerGraph, SteinerLabels


def vertex_queries(V, F, n_queries, seed=0):
    rng = np.random.default_rng(seed)
    n = len(V)
    qp = []
    while len(qp) < n_queries:
        s, t = int(rng.integers(n)), int(rng.integers(n))
        if s != t:
            qp.append((s, t))
    return qp


def write_competitor_queries(path, V, F, qp, nudge=0.05):
    #their format: "x y z fid  x y z fid" per line, one face containing the vertex. their
    #harness is float32 and locates query points on the mesh with an AABB tree, which goes
    #wrong for a point exactly at a vertex (12 of 100 MMP answers were off, one below the
    #euclidean bound). so each point is moved `nudge` metres into its face: by the triangle
    #inequality that changes the true distance by at most 2*nudge (0.003% at 3 km)
    first_face = np.full(len(V), -1, dtype=int)
    for fi, tri in enumerate(F):
        for v in tri:
            if first_face[v] < 0:
                first_face[v] = fi
    def moved(v):
        fi = first_face[v]
        c = V[F[fi]].mean(0)
        u = c - V[v]
        return V[v] + u * (nudge / max(np.linalg.norm(u), 1e-12)), fi
    with open(path, "w") as f:
        for s, t in qp:
            p, fs = moved(s); q, ft = moved(t)
            f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f} {fs} {q[0]:.6f} {q[1]:.6f} {q[2]:.6f} {ft}\n")


def a2a_queries(tm, n_queries, seed=0):
    #arbitrary surface points: uniform random face, uniform barycentric coordinates (the
    #competitors' generator does the same). returns (p, face, q, face) in chart coordinates
    rng = np.random.default_rng(seed + 1)
    F = tm.F
    out = []
    while len(out) < n_queries:
        fp, fq = rng.integers(len(F), size=2)
        if fp == fq:
            continue
        p = rng.dirichlet([1, 1, 1]) @ tm.xy[F[fp]]
        q = rng.dirichlet([1, 1, 1]) @ tm.xy[F[fq]]
        out.append((p, int(fp), q, int(fq)))
    return out


def write_competitor_a2a(path, tm, qa):
    with open(path, "w") as f:
        for p, fp, q, fq in qa:
            P = tm.lift(p[None])[0] + np.array([*tm.origin, 0.0])
            Q = tm.lift(q[None])[0] + np.array([*tm.origin, 0.0])
            f.write(f"{P[0]:.6f} {P[1]:.6f} {P[2]:.6f} {fp} {Q[0]:.6f} {Q[1]:.6f} {Q[2]:.6f} {fq}\n")


def exact_a2a(tm, qa):
    #exact geodesic between two arbitrary surface points: both inserted as steiner vertices
    import pygeodesic.geodesic as pg
    ex = ExactTerrain(tm)
    out, t0 = [], time.perf_counter()
    for p, fp, q, fq in qa:
        V2, F2, idx = ex.insert_points(np.array([p, q]))
        geo = pg.PyGeodesicAlgorithmExact(V2, F2)
        out.append(float(geo.geodesicDistance(int(idx[0]), int(idx[1]))[0]))
    return np.array(out), (time.perf_counter() - t0) / len(qa)


def measure_a2a(query, qa, truth, repeat=2):
    est = np.array([query(*x) for x in qa])
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        for x in qa:
            query(*x)
        best = min(best, (time.perf_counter() - t0) / len(qa))
    errs = est / truth - 1.0
    return dict(query_us=best * 1e6, err_mean=float(errs.mean() * 100), err_max=float(errs.max() * 100),
                err_min=float(errs.min() * 100), err_p90=float(np.percentile(errs, 90) * 100),
                abs_err_mean=float(np.abs(errs).mean() * 100))


def exact_geodesics(V, F, qp):
    import pygeodesic.geodesic as pg
    geo = pg.PyGeodesicAlgorithmExact(V, F)
    out, t0 = {}, time.perf_counter()
    for s, t in qp:
        out[(s, t)] = float(geo.geodesicDistance(s, t)[0])
    return out, (time.perf_counter() - t0) / len(qp)


def index_bytes(index):
    if isinstance(index, H2HIndex):
        return index.label_entries * 12 #one float64 distance + one int32 middle vertex per entry
    b = sum(v.nbytes for d in index.label for v in d.values()) #graph portal labels
    if hasattr(index, "cut_label"): #the cut index queries both, so it stores both
        b += sum(v.nbytes for d in index.cut_label for v in d.values())
        b += sum(8 for d in index.foot_r for _ in d) #the foot drop per (node, cut) for the model interpolation
    b += sum(D.nbytes for D in index.leaf_table if D is not None)
    if hasattr(index, "cut_xy"):
        b += sum(x.nbytes for x in index.cut_xy)
    if getattr(index, "halo_ids", None) is not None:
        b += sum(len(h) * 12 for h in index.halo_ids if h is not None) #int32 id + float64 distance
    return b


def measure(query, qp, truth, repeat=3, keep_values=False):
    vals = np.array([query(s, t) for s, t in qp], dtype=float)
    errs = vals / np.array([truth[(s, t)] for s, t in qp]) - 1.0
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        for s, t in qp:
            query(s, t)
        best = min(best, (time.perf_counter() - t0) / len(qp))
    out = dict(query_us=best * 1e6, err_mean=float(errs.mean() * 100), err_max=float(errs.max() * 100),
               err_min=float(errs.min() * 100), err_p90=float(np.percentile(errs, 90) * 100))
    if keep_values: #weighted terrains have no exact truth: the per-query values are kept so the
        out["values"] = vals.tolist() #reference can be the best path any method found
    return out


def run(mesh_path, n_queries=100, k=32, workers=1, seed=0, out_dir=None, slow_queries=30,
        segments=20, rounds=2, leaf_size=16, methods=("dijkstra", "h2h", "portals", "exact", "cut", "refine"),
        k_list=None, halo_rings=2, n_a2a=0, a2a_near=4, weighted=None, cert=(), steiner_m=(5, 10),
        steiner_ref=0):
    k_list = list(k_list) if k_list else [k]
    name = os.path.splitext(os.path.basename(mesh_path))[0]
    V, F = read_off(mesh_path)
    V, F = repair_height_field(V, F, verbose=False)
    weights = None
    if weighted is not None: #face weights uniform in (1, 5), the competitors' own distribution
        weights = np.random.default_rng(int(weighted)).uniform(1.000001, 5.0, len(F))
        name = name + "_weighted"
    tm = TerrainMesh(V, F, weights=weights)
    points, adj = tm.mesh_graph()
    n = len(points)
    qp = vertex_queries(V, F, n_queries, seed)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        write_competitor_queries(os.path.join(out_dir, f"{name}.{'W2V' if weights is not None else 'V2V'}.query"), V, F, qp)
        if weights is not None:
            with open(os.path.join(out_dir, f"{name}.face_weight.query"), "w") as f:
                f.writelines(f"{w:.6f}\n" for w in weights)
    if weights is None:
        truth, t_exact = exact_geodesics(V, F, qp)
    else:
        #no exact solver: the reference is the heaviest query-time refinement we have; every
        #method's per-query values are stored, and the report takes the best path found by
        #any method (theirs included) as the reference
        P0 = np.array(points)
        t0 = time.perf_counter()
        truth = {}
        for s, t in qp:
            dist, pred = dijkstra(adj, s, t)
            truth[(s, t)] = refine_polyline(tm, P0[get_path(pred, s, t)], (), rounds=3, segments=40)[0]
        t_exact = (time.perf_counter() - t0) / len(qp)
    keep = weights is not None
    P = np.array(points)
    rows = []
    print(f"{name}: {n} vertices, {len(F)} faces, {sum(len(a) for a in adj) // 2} edges, "
          f"exact geodesic {t_exact * 1e3:.1f} ms/query", flush=True)
    qa, truth_a = [], None
    if n_a2a:
        qa = a2a_queries(tm, n_a2a, seed)
        if out_dir:
            write_competitor_a2a(os.path.join(out_dir, f"{name}.A2A.query"), tm, qa)
        truth_a, t_a = exact_a2a(tm, qa)
        print(f"  {len(qa)} arbitrary-point queries, exact geodesic {t_a * 1e3:.1f} ms/query", flush=True)

    def save():
        #written after every row, so a crash or a killed session keeps what is finished
        if not out_dir:
            return
        path = os.path.join(out_dir, f"{name}.json")
        result = dict(mesh=name, n=n, faces=len(F), queries=len(qp), exact_ms=t_exact * 1e3, k=k, rows=rows)
        if os.path.exists(path): #merge: a re-run of some methods replaces those rows and keeps the rest
            old = json.load(open(path))
            names = {r["method"] for r in rows}
            result["rows"] = [r for r in old["rows"] if r["method"] not in names] + rows
        with open(path, "w") as f:
            json.dump(result, f, indent=1)

    def add(label, build, size, m):
        rows.append(dict(method=label, build_s=build, index_bytes=size, **m))
        extra = f" (python {m['query_us_python']:.1f} us, batch {m['query_us_batch']:.2f} us)" if "query_us_python" in m else ""
        print(f"  {label:34s} build {build:8.1f}s  size {size / 1e6:7.2f} MB  query {m['query_us']:9.1f} us{extra}  "
              f"err mean {m['err_mean']:.3f}%  p90 {m['err_p90']:.3f}%  max {m['err_max']:.3f}%  min {m['err_min']:.3f}%", flush=True)
        save()

    if "dijkstra" in methods:
        add("dijkstra on mesh graph (no index)", 0.0, 0,
            measure(lambda s, t: dijkstra(adj, s, t)[0][t], qp, truth, repeat=1, keep_values=keep))
    if "h2h" in methods:
        t0 = time.perf_counter(); h2h = H2HIndex(adj); tb = time.perf_counter() - t0
        add("h2h", tb, index_bytes(h2h), measure(h2h.query, qp, truth, keep_values=keep))
    if "portals" in methods:
        t0 = time.perf_counter(); vp = SeparatorIndex(points, adj, k=k, leaf_size=leaf_size); tb = time.perf_counter() - t0
        add(f"vertex portals k={k}", tb, index_bytes(vp), measure(vp.query, qp, truth, keep_values=keep))
    configs = []
    if "exact" in methods: #uniform k per cut is the terrain default (see notes/error_bound.md, section 6)
        configs += [dict(tag="exact MMP labels", exact=True, k=kk, scale_k=False, halo_rings=halo_rings) for kk in k_list]
    if "exact_scaled" in methods: #the cheap variant: portals proportional to cut length, k_min=4
        configs += [dict(tag="exact MMP labels, k~|cut|", exact=True, k=k, scale_k=True, halo_rings=halo_rings)]
    if "cut" in methods: #the general-manifold pipeline: refined graph paths as labels
        configs += [dict(tag=f"refined labels (rounds={rounds}, seg={segments})", exact=False, k=k, scale_k=False, halo_rings=0)]
    for cfg in configs:
        t0 = time.perf_counter()
        cut = SeparatorCutIndex(points, adj, tm, k=cfg["k"], leaf_size=leaf_size, interpolate=False,
                                rounds=rounds, segments=segments, workers=workers, pad_speed=1.0,
                                scale_k=cfg["scale_k"], halo_rings=cfg["halo_rings"],
                                exact=ExactTerrain(tm) if cfg["exact"] else None)
        mv = lambda q: measure(q, qp, truth, keep_values=keep)
        tb = time.perf_counter() - t0
        print(f"    build breakdown: {', '.join(f'{a} {b:.1f}s' for a, b in cut.timing.items())}, "
              f"{sum(len(p) for p in cut.cut_pos)} portals on {len(cut.sep_points)} cuts, {len(cut.leaves)} leaves", flush=True)
        label = f"cut portals k={cfg['k']}, {cfg['tag']}"
        add(label, tb, index_bytes(cut), mv(cut.query))
        cut.interpolate = True
        add(label + " + interp", tb, index_bytes(cut), mv(cut.query))
        cut.interpolate = "model" #hyperbola-plus-parabola fit: handles an endpoint close to the cut
        add(label + " + model interp", tb, index_bytes(cut), mv(cut.query))
        if qa and cfg["exact"]:
            for mode in ("snap", "interp"):
                a2a = TerrainA2A(tm, cut, mode=mode, near=a2a_near)
                m = measure_a2a(a2a.query, qa, truth_a)
                add(label + f" + model interp, A2A {mode}", tb, index_bytes(cut), m)
    if steiner_ref and weights is not None: #a strong reference for weighted terrains: fine Steiner graph
        t0 = time.perf_counter()
        sg = SteinerGraph(tm, steiner_ref)
        from scipy.sparse.csgraph import dijkstra as sp_dijkstra
        src = np.unique([s for s, _ in qp])
        D = sp_dijkstra(sg.G, indices=src)
        pos = {int(v): i for i, v in enumerate(src)}
        vals = [float(D[pos[s], t]) for s, t in qp]
        dt = (time.perf_counter() - t0) / len(qp)
        errs = np.array(vals) / np.array([truth[q] for q in qp]) - 1
        add(f"Steiner graph m={steiner_ref} Dijkstra (no index)", 0.0, 0,
            dict(query_us=dt * 1e6, err_mean=float(errs.mean() * 100), err_max=float(errs.max() * 100),
                 err_min=float(errs.min() * 100), err_p90=float(np.percentile(errs, 90) * 100), values=vals))
        del sg, D
    for spec in cert:
        eps, gf = (float(x) for x in spec.split("/"))
        oracles = [("exact MMP", None)] if weights is None else [(f"Steiner m={m}", m) for m in steiner_m]
        for oname, mm in oracles:
            t0 = time.perf_counter()
            if mm is None:
                ix = CertifiedCutIndex(tm, eps=eps, gap_factor=gf)
            else:
                ix = CertifiedCutIndex(tm, eps=eps, gap_factor=gf, oracle=SteinerLabels(tm, FaceGrid(tm), m=mm))
            tb = time.perf_counter() - t0
            st = ix.stats()
            print(f"    build breakdown: {', '.join(f'{a} {b:.1f}s' for a, b in ix.timing.items())}; "
                  f"{st['cuts']} cuts, {st['portals']} portals, {st['near_entries'] / n:.0f} near entries per node", flush=True)
            label = f"certified eps={eps:g} g={gf:g}e, {oname}"
            size = ix.nbytes()
            cq = cquery.CompiledQuery(ix) if cquery.available() else None
            def timed(row_fn, c_fn, ix=ix, cq=cq, st=st):
                m = measure(row_fn, qp, truth, keep_values=keep)
                if cq is not None: #same answers from the compiled query (tests check it); its latency
                    m["query_us_python"] = m["query_us"]
                    m["query_us"] = measure(c_fn, qp, truth)["query_us"]
                    S = np.array([a for a, _ in qp]); T = np.array([b for _, b in qp])
                    t0 = time.perf_counter()
                    for _ in range(20):
                        cq.bounds_batch(S, T)
                    m["query_us_batch"] = (time.perf_counter() - t0) / (20 * len(qp)) * 1e6
                m["stats"] = st; m["timing"] = ix.timing
                return m
            add(label + " (UB)", tb, size, timed(ix.query_ub, cq.query_ub if cq else None))
            add(label + " (estimate)", tb, size, timed(ix.query, cq.query if cq else None))
            if weights is None: #check the certificate on every query
                B = np.array([ix.bounds(s, t) for s, t in qp]); tr = np.array([truth[q] for q in qp])
                w = (B[:, 1] - B[:, 0]) / B[:, 0]
                rows[-1]["certificate"] = dict(sound=bool(np.all(B[:, 0] <= tr) and np.all(B[:, 1] >= tr)),
                                               max_ub_err=float((B[:, 1] / tr - 1).max() * 100),
                                               min_lb_ratio=float((B[:, 0] / tr).min()),
                                               width_median=float(np.median(w) * 100), width_max=float(w.max() * 100),
                                               eps=eps, stats=st, timing=ix.timing)
                print(f"    certificate: sound {rows[-1]['certificate']['sound']}, max UB error "
                      f"{rows[-1]['certificate']['max_ub_err']:.3f}% (eps {eps * 100:g}%), width median "
                      f"{rows[-1]['certificate']['width_median']:.3f}% max {rows[-1]['certificate']['width_max']:.3f}%", flush=True)
                save()
            if qa and mm is None: #vertex-pair answers from the compiled query when there is one (same values)
                a2a = TerrainA2A(tm, cq if cq is not None else ix, mode="interp", near=a2a_near)
                add(label + ", A2A interp", tb, size, measure_a2a(a2a.query, qa, truth_a))
            del ix
    if "refine" in methods:
        sub = qp[:slow_queries]
        def q(s, t):
            dist, pred = dijkstra(adj, s, t)
            return refine_polyline(tm, P[get_path(pred, s, t)], (), rounds=rounds, segments=segments)[0]
        add(f"dijkstra + refine at query time (rounds={rounds}, seg={segments})", 0.0, 0,
            measure(q, sub, truth, repeat=1, keep_values=keep))
        if keep: #the heavy refinement that served as reference is a method too (all queries)
            add("dijkstra + refine at query time (rounds=3, seg=40)", 0.0, 0,
                dict(query_us=t_exact * 1e6, err_mean=0.0, err_max=0.0, err_min=0.0, err_p90=0.0,
                     values=[truth[(s, t)] for s, t in qp]))
    save()
    return dict(mesh=name, n=n, faces=len(F), queries=len(qp), exact_ms=t_exact * 1e3, k=k, rows=rows)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("--queries", type=int, default=100)
    ap.add_argument("--k", type=int, default=32)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--a2a", type=int, default=0, help="number of arbitrary-point queries (0 = none)")
    ap.add_argument("--a2a-near", type=int, default=4)
    ap.add_argument("--weighted", type=int, default=None, help="seed for uniform(1,5) face weights (weighted terrain)")
    ap.add_argument("--segments", type=int, default=20)
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--leaf-size", type=int, default=16)
    ap.add_argument("--out", default=os.path.join(ROOT, "results"))
    ap.add_argument("--methods", default="dijkstra,h2h,portals,exact,exact_scaled")
    ap.add_argument("--k-list", default="", help="comma separated k values for the exact index (default: --k)")
    ap.add_argument("--halo-rings", type=int, default=2)
    ap.add_argument("--cert", default="", help="certified index configs eps/gap_factor, comma separated, e.g. 0.1/0.5,0.05/0.25")
    ap.add_argument("--steiner-m", default="5,10", help="Steiner points per edge for weighted labels")
    ap.add_argument("--steiner-ref", type=int, default=0, help="weighted: add a Dijkstra reference on a Steiner graph with this many points per edge")
    a = ap.parse_args()
    run(a.mesh, a.queries, a.k, a.workers, out_dir=a.out, segments=a.segments, rounds=a.rounds,
        leaf_size=a.leaf_size, methods=tuple(a.methods.split(",")),
        k_list=[int(x) for x in a.k_list.split(",") if x], halo_rings=a.halo_rings, n_a2a=a.a2a, a2a_near=a.a2a_near,
        weighted=a.weighted, cert=[c for c in a.cert.split(",") if c],
        steiner_m=[int(x) for x in a.steiner_m.split(",") if x], steiner_ref=a.steiner_ref)
