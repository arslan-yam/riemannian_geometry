"""Decision workload: what each method's answers decide, and how many of those decisions are wrong.

  python3 scripts/server/decisions.py PAIRS MESH.off name:kind:file [name:kind:file ...] [--json out.json]

PAIRS is an export of cpp/tools/proximity (--export): "query object distance" per line, grouped by query, with the
exact distance. Every method file answers the same pairs, one line each, in the same order:
  kind "point"     one value per line (DGG: a graph-path length; SEP, NeuroGF: an estimate);
  kind "interval"  "LB UB [estimate]" per line (ours: proximity --dump-pairs; landmarks: landmark_bench --dump-prox).

Three operators, all from the same pairs, with budgets / radii of 1, 3 and 10% of the chart diagonal:
  * go / no-go: "is d(q, o) <= B?" for every pair. A point method says go iff its value <= B. An interval method
    says go iff UB <= B, no-go iff LB > B, and otherwise refines (one exact distance, so it is then right).
    Unsafe go: go with d > B. Lost go: no-go with d <= B.
  * coverage: the objects within r of the query (a range query). Missed objects and queries with a miss.
  * dispatch: the k nearest objects (k = 1, 10). A point method returns the k smallest values; an interval method
    decides from the intervals (in if fewer than k others can be closer, out if k others are surely closer) and
    refines the rest. Wrong: a returned set that is not a true k-nearest set.
Truth ties: a decision whose true distance is within 1e-6 (relative) of its budget, or a k-th / (k+1)-th
distance pair that close, is not scored (the exported truths are float32).
"""
import json, sys
import numpy as np

TIE = 1e-6


def load_pairs(path):
    a = np.loadtxt(path, dtype=np.float64)
    q = a[:, 0].astype(np.int64)
    starts = np.flatnonzero(np.r_[True, q[1:] != q[:-1]])
    return q, a[:, 2], np.r_[starts, len(q)]


def diag_of(mesh):
    with open(mesh) as f:
        assert f.readline().strip().upper().startswith("OFF")
        line = f.readline()
        while line.startswith("#") or not line.strip():
            line = f.readline()
        nv = int(line.split()[0])
        xy = np.loadtxt(f, max_rows=nv, usecols=(0, 1))
    return float(np.hypot(*(xy.max(0) - xy.min(0))))


def go_nogo(d, lo, hi, B):
    """counts for one budget; lo == hi for a point method"""
    keep = np.abs(d - B) > TIE * B
    d, lo, hi = d[keep], lo[keep], hi[keep]
    true_go = d <= B
    go, nogo = hi <= B, lo > B
    undecided = ~go & ~nogo
    unsafe = go & ~true_go
    lost = nogo & true_go
    wrong = unsafe | lost
    margin = float(np.max(np.abs(d[wrong] - B) / B)) if wrong.any() else 0.0
    return dict(pairs=int(keep.sum()), true_go=int(true_go.sum()), unsafe_go=int(unsafe.sum()), lost_go=int(lost.sum()),
                undecided=int(undecided.sum()), worst_margin=margin)


def coverage(d, lo, hi, B, bounds):
    """range queries: objects missed, queries with a miss, objects to refine (interval methods)"""
    missed = qmiss = fp = refine = qref = 0
    for a, b in zip(bounds[:-1], bounds[1:]):
        dd, l, h = d[a:b], lo[a:b], hi[a:b]
        keep = np.abs(dd - B) > TIE * B
        t = dd <= B
        sure_in, sure_out = h <= B, l > B
        und = ~sure_in & ~sure_out & keep  # refined: exact afterwards
        sel = sure_in | (und & t)
        m = int((t & ~sel & keep).sum()); missed += m; qmiss += m > 0
        fp += int((sel & ~t & keep).sum())
        refine += int(und.sum()); qref += und.any()
    return dict(missed=missed, queries_with_miss=qmiss, false_positives=fp, refined_objects=refine, refined_queries=qref)


def dispatch(d, lo, hi, k, bounds, point):
    wrong = qref = refine = scored = 0
    ap_hit = ap_tot = ap_over = 0; ap_worst = 1.0  # interval methods: the k smallest UB, no refine
    for a, b in zip(bounds[:-1], bounds[1:]):
        dd, l, h = d[a:b], lo[a:b], hi[a:b]
        if len(dd) <= k:
            continue
        srt = np.sort(dd)
        dk, dk1 = srt[k - 1], srt[k]
        if dk1 - dk <= TIE * dk1:  # the true k-nearest set is not unique at float32 precision
            continue
        scored += 1
        if point:
            sel = np.argpartition(l, k - 1)[:k]
            wrong += bool(np.any(dd[sel] > dk * (1 + TIE)))
            continue
        ap = np.argpartition(h, k - 1)[:k]
        ap_hit += int(np.sum(dd[ap] <= dk * (1 + TIE))); ap_tot += k
        r = float(np.max(dd[ap]) / dk); ap_worst = max(ap_worst, r); ap_over += r > 1.2 * (1 + TIE)
        # interval: o is surely in if fewer than k others have LB < UB(o); surely out if >= k others have UB < LB(o)
        Ls, Us = np.sort(l), np.sort(h)
        below_ub = np.searchsorted(Ls, h, side="left") - (l < h)  # others with LB < my UB
        below_lb = np.searchsorted(Us, l, side="left")           # others with UB < my LB
        sure_in, sure_out = below_ub < k, below_lb >= k
        und = ~sure_in & ~sure_out
        n = int(und.sum()); refine += n; qref += n > 0
        # after refining the undecided ones exactly, the answer is exact: wrong only if the bounds were unsound
        sel = sure_in | (und & (dd <= dk))
        wrong += bool(np.any(dd[sel] > dk * (1 + TIE)) or sel.sum() < k)
    out = dict(queries=scored, wrong=wrong, refined_queries=qref, refined_objects=refine)
    if not point:
        out.update(approx_recall=ap_hit / max(ap_tot, 1), approx_worst_ratio=ap_worst, approx_over_1_2=ap_over)
    return out


def main():
    args = [a for a in sys.argv[1:]]
    out = None
    if "--json" in args:
        i = args.index("--json"); out = args[i + 1]; del args[i:i + 2]
    pairs, mesh, methods = args[0], args[1], args[2:]
    q, d, bounds = load_pairs(pairs)
    diag = diag_of(mesh)
    res = dict(pairs=pairs, n_pairs=len(d), queries=len(bounds) - 1, diag=diag, methods=[])
    print(f"{len(d)} pairs, {len(bounds) - 1} queries, diag {diag:.1f}")
    for spec in methods:
        name, kind, path = spec.split(":", 2)
        v = np.loadtxt(path, dtype=np.float64, ndmin=2)
        if len(v) != len(d):
            print(f"{name}: {len(v)} lines for {len(d)} pairs, skipped"); continue
        point = kind == "point"
        lo, hi = (v[:, 0], v[:, 0]) if point else (v[:, 0], v[:, 1])
        if not point:
            unsound = int(np.sum((lo > d * (1 + 1e-6)) | (hi < d * (1 - 1e-6))))
        else:
            unsound = None
        rel = (v[:, 0] if point else (v[:, 2] if v.shape[1] > 2 else hi)) / d - 1
        m = dict(name=name, kind=kind, err_mean=float(np.mean(np.abs(rel))), err_min=float(rel.min()), err_max=float(rel.max()),
                 interval_violations=unsound, budgets=[])
        for f in (0.01, 0.03, 0.1):
            B = f * diag
            m["budgets"].append(dict(frac=f, go_nogo=go_nogo(d, lo, hi, B), coverage=coverage(d, lo, hi, B, bounds)))
        m["dispatch"] = [dict(k=k, **dispatch(d, lo, hi, k, bounds, point)) for k in (1, 10)]
        res["methods"].append(m)
        g = m["budgets"]
        print(f"{name:14s} err mean {m['err_mean']:.2e} [{m['err_min']:+.1e}, {m['err_max']:+.1e}]"
              + (f", interval violations {unsound}" if unsound is not None else ""))
        for x in g:
            a, c = x["go_nogo"], x["coverage"]
            print(f"   B={x['frac']:.2f}: unsafe go {a['unsafe_go']}, lost go {a['lost_go']}, undecided {a['undecided']} of {a['pairs']} "
                  f"(worst margin {a['worst_margin']:.1e}) | coverage: missed {c['missed']} in {c['queries_with_miss']} queries, "
                  f"false pos {c['false_positives']}, refined {c['refined_objects']} objects in {c['refined_queries']} queries")
        for x in m["dispatch"]:
            print(f"   k={x['k']}: wrong {x['wrong']} of {x['queries']} queries, refined {x['refined_objects']} objects in {x['refined_queries']} queries")
    if out:
        json.dump(res, open(out, "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))


if __name__ == "__main__":
    main()
