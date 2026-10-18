"""kNN and range queries answered by a distance oracle that has no proximity operator: one distance query
per object, then the k smallest (kNN) or those within r (range). Scores each method's answers to the pairs
written by cpp/tools/proximity --export (the same queries and objects as our certified proximity queries).

usage: python3 scripts/competitors/knn_eval.py --pairs PREFIX.pairs --mesh m.off --answers NAME=file ...
                                               [--k 1,10,25] [--radius 0.01,0.03,0.1] [--json out.json]
  PREFIX.pairs: "query object distance" per line, grouped by query (the exact distances)
  file:         "answer query_ms" per line, one per pair (cpp/tools/upo_static --dump)
  --radius:     fractions of the chart diagonal, as cpp/tools/proximity
A kNN or range query costs the sum of the method's query times over all objects of that query.
"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rgeo.terrain import read_off

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True); ap.add_argument("--mesh", required=True)
ap.add_argument("--answers", nargs="+", required=True)
ap.add_argument("--k", default="1,10,25"); ap.add_argument("--radius", default="0.01,0.03,0.1")
ap.add_argument("--json")
a = ap.parse_args()

P = np.loadtxt(a.pairs).reshape(-1, 3)
qs, first = np.unique(P[:, 0].astype(np.int64), return_index=True)
order = np.argsort(first)
groups = [np.nonzero(P[:, 0] == qs[i])[0] for i in order]   # the pairs of each query, in file order
V, _ = read_off(a.mesh)
xy = np.asarray(V)[:, :2]
diag = float(np.hypot(*(xy.max(0) - xy.min(0))))
ks = [int(x) for x in a.k.split(",")]
rs = [float(x) for x in a.radius.split(",")]

out = {"pairs": a.pairs, "queries": len(groups), "objects_per_query": float(np.mean([len(g) for g in groups])), "diag": diag, "methods": []}
for spec in a.answers:
    name, f = spec.split("=", 1)
    if not os.path.exists(f) or os.path.getsize(f) == 0:  # that method did not finish (e.g. killed for memory)
        print(f"{name}: no answers ({f}), skipped"); out.setdefault("missing", []).append(name); continue
    A = np.loadtxt(f).reshape(-1, 2)
    if len(A) != len(P):
        print(f"{name}: {len(A)} answers for {len(P)} pairs, skipped"); continue
    m = {"method": name, "pair_query_ms": float(A[:, 1].mean()), "knn": [], "range": []}
    tq = np.array([A[g, 1].sum() for g in groups])                # ms per kNN / range query
    for k in ks:
        rec, worst = [], []
        for g in groups:
            d, x = P[g, 2], A[g, 0]
            kk = min(k, len(g))
            true = set(np.argsort(d, kind="stable")[:kk]); got = np.argsort(x, kind="stable")[:kk]
            rec.append(len(true & set(got)) / kk)
            dk = np.sort(d)[kk - 1]
            worst.append(d[got].max() / dk if dk > 0 else 1.0)
        m["knn"].append({"k": k, "recall": float(np.mean(rec)), "worst_ratio": float(np.max(worst)),
                         "query_ms_mean": float(tq.mean()), "query_ms_p50": float(np.median(tq))})
    for rf in rs:
        r = rf * diag; prec, rec, fp, fn = [], [], 0, 0
        for g in groups:
            d, x = P[g, 2], A[g, 0]
            t, s = d <= r, x <= r
            inter = int((t & s).sum())
            prec.append(inter / s.sum() if s.sum() else 1.0); rec.append(inter / t.sum() if t.sum() else 1.0)
            fp += int((s & ~t).sum()); fn += int((t & ~s).sum())
        m["range"].append({"radius_frac": rf, "radius": r, "precision": float(np.mean(prec)), "recall": float(np.mean(rec)),
                           "false_positives": fp, "missed": fn, "query_ms_mean": float(tq.mean())})
    out["methods"].append(m)
    print(f"{name}: {m['pair_query_ms']:.4f} ms per pair, {tq.mean():.1f} ms per kNN/range query | " +
          ", ".join(f"k={x['k']} recall {x['recall']:.3f} worst {x['worst_ratio']:.3f}" for x in m["knn"]) + " | " +
          ", ".join(f"r={x['radius_frac']} P {x['precision']:.3f} R {x['recall']:.3f}" for x in m["range"]))
if a.json:
    json.dump(out, open(a.json, "w"), indent=1)
