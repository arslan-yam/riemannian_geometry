"""Tables of the review experiments (results/server/review, the `review` phase of run_server.sh), as markdown.

  python3 scripts/server/review_summary.py [results/server] > notes/review_results.md
"""
import glob, json, os, statistics, sys

R = sys.argv[1] if len(sys.argv) > 1 else "results/server"
RV = os.path.join(R, "review")
TERR = ["TJ", "SC", "VS", "AU", "GI", "LH"]


def J(path):
    return json.load(open(path)) if os.path.exists(path) else None


def row(*c):
    print("| " + " | ".join(str(x) for x in c) + " |")


def head(*c):
    row(*c); row(*["---"] * len(c))


def var(d, variant="dense/sym", s="random"):
    return next((v for v in d["variants"] if v["variant"] == variant and v["set"] == s), None)


def strict(d):
    vs = d["variants"]
    return sum(v.get("strict_violations", 0) for v in vs), max(v.get("max_rel_violation", 0) for v in vs), sum(v["pairs"] for v in vs)


print("# Review experiments\n")
print("## 1. DGG: full Dijkstra (their query) and point-to-point, all 2,000 random pairs, 0.5M faces\n")
head("terrain", "eps", "precompute", "graph", "full sweep", "point-to-point", "speed-up", "mean err", "max err")
for T in TERR:
    for E in ("0.001", "0.0001"):
        d = J(os.path.join(RV, f"dgg.{T}pre_500000.{E}.json"))
        if d:
            row(T, E, f"{d['precompute_s'] / 60:.1f} min" if d.get("precompute_s") else "", f"{d['graph_bytes'] / 1e9:.2f} GB", f"{d['query_ms']:.0f} ms",
                f"{d['query_ms_p2p']:.0f} ms", f"{d['query_ms'] / d['query_ms_p2p']:.1f}x", f"{100 * d['err_mean']:.4f}%", f"{100 * d['err_max']:.3f}%")

print("\n## 1. DGG kNN / range over 10,000 objects (early-stopping Dijkstra), 100 queries\n")
head("terrain", "eps", "k", "recall", "ms per query", "", "radius", "precision", "recall", "ms per query")
for T in TERR:
    for E in ("0.001", "0.0001"):
        d = J(os.path.join(RV, f"dggknn.{T}pre_500000.{E}.json"))
        if not d:
            continue
        m = d["methods"][0]
        for i in range(max(len(m["knn"]), len(m["range"]))):
            k = m["knn"][i] if i < len(m["knn"]) else None; r = m["range"][i] if i < len(m["range"]) else None
            row(T, E, k["k"] if k else "", f"{k['recall']:.4f}" if k else "", f"{k['query_ms_mean']:.2f}" if k else "", "",
                r["radius_frac"] if r else "", f"{r['precision']:.4f}" if r else "", f"{r['recall']:.4f}" if r else "", f"{r['query_ms_mean']:.1f}" if r else "")

print("\n## 2. DGG + contraction hierarchy\n")
head("terrain", "eps", "build", "stopped", "contracted", "core", "shortcuts", "index", "CH query", "p2p Dijkstra", "same answers", "peak RSS")
for T in ("TJ", "GI"):
    for E in ("0.001", "0.0001"):
        d = J(os.path.join(RV, f"dggch.{T}pre_500000.{E}.json"))
        if d:
            row(T, E, f"{d['build_s'] / 3600:.1f} h", d["stopped"], d["contracted"], d["core"], d["shortcuts"], f"{d['index_bytes'] / 1e9:.2f} GB",
                f"{d['ch_query_ms']:.1f} ms", f"{d['p2p_query_ms']:.1f} ms", f"{d['ch_vs_p2p_max_rel_diff']:.1e}", f"{(d.get('run_rss_kb') or 0) / 1e6:.0f} GB")

print("\n## 2. DGG + hub labels (pruned landmark labelling, degree order)\n")
head("terrain", "eps", "complete", "hubs", "build", "entries / vertex", "index", "query", "same answers", "mean err", "peak RSS")
for T in ("TJ", "GI"):
    for E in ("0.001", "0.0001"):
        d = J(os.path.join(RV, f"dggpll.{T}pre_500000.{E}.json"))
        if d and "pll" in d:
            p = d["pll"]
            c = p["complete"]  # a stopped build: the labels so far are a lower bound on size, and its answers are not distances
            ge = "" if c else ">= "
            row(T, E, c, p["hubs_done"], f"{p['build_s'] / 3600:.1f} h", f"{ge}{p['entries_per_vertex']:.0f}", f"{ge}{p['index_bytes'] / 1e9:.2f} GB",
                f"{p['query_us']:.1f} us" if c else "-", f"{p['vs_p2p_max_rel_diff']:.1e}" if c else "-", f"{100 * p['err_mean']:.4f}%" if c else "-",
                f"{(d.get('run_rss_kb') or 0) / 1e6:.0f} GB")

print("\n## 3. Exact-landmark intervals (farthest-point landmarks, float32), 2,000 random pairs and kNN over 10,000 objects\n")
head("terrain", "L", "size", "build", "lookup", "sound", "width median", "width p90", "certified at 1.2", "refine / pair", "end-to-end", "kNN k=10: ambiguous", "refined queries", "kNN ms")
for T in ("TJ", "GI"):
    d = J(os.path.join(RV, f"landmark.{T}.json"))
    if not d:
        continue
    for x in d["landmarks"]:
        k10 = next((k for k in x["knn"] if k["k"] == 10), None)
        row(T, x["L"], f"{x['index_bytes'] / 1e6:.0f} MB", f"{x['build_s'] / 3600:.2f} h", f"{x['query_us']:.2f} us", x["sound"], f"{100 * x['width_median']:.3f}%",
            f"{100 * x['width_p90']:.2f}%", f"{x['certified']}/{d['pairs']}", f"{x['refine_ms']:.0f} ms", f"{x['end_to_end_ms']:.2f} ms",
            f"{k10['ambiguous_mean']:.1f}" if k10 else "", f"{k10['refined_queries']}/{k10['queries']}" if k10 else "", f"{k10['end_to_end_ms']:.2f}" if k10 else "")
print("\nOurs for comparison (cap 64, 8-bit / uncapped): see tab:server-ours; width median from ours.*.json:")
for T in ("TJ", "GI"):
    for tag in ("k64", "kd"):
        d = J(os.path.join(R, f"ours.{T}.{tag}.json"))
        if d:
            v = var(d, "code8/sym" if tag == "k64" else "dense/sym")
            print(f"- {T} {tag}: {v['index_bytes'] / 1e6:.0f} MB, lookup {v['query_us']:.2f} us, width median {100 * v['width_median']:.3f}%, refined {v['refined']}/{v['pairs']}")

print("\n## 5. Adversarial variants of TJ (0.5M faces)\n")
head("kind", "vertices", "aspect p99", "slope p99", "variant", "build", "size", "query", "mean err", "max err", "width median", "refined (all sets)", "strict viol.", "sound", "DGG 1e-4 p2p / err")
for K in ("steep5", "steep10", "aniso", "slivers", "cluster", "ties", "comb"):
    mj = J(os.path.join(RV, "meshes", f"TJadv_{K}.json")) or {}
    g = J(os.path.join(RV, f"dgg.TJadv_{K}.0.0001.json"))
    for tag in ("kd", "k64"):
        d = J(os.path.join(RV, f"adv.{K}.{tag}.json"))
        if not d:
            continue
        v = var(d)
        sv, mv, _ = strict(d)
        row(K, mj.get("vertices", mj.get("vertices_out", "")), f"{mj['aspect_p99']:.0f}" if "aspect_p99" in mj else "", f"{mj['slope_p99_deg']:.0f}" if "slope_p99_deg" in mj else "",
            tag, f"{d['build_s'] / 3600:.1f} h", f"{v['index_bytes'] / 1e6:.0f} MB", f"{v['query_us']:.2f} us", f"{100 * v['est_mean']:.4f}%", f"{100 * v['est_max']:.3f}%",
            f"{100 * v['width_median']:.3f}%", sum(x["refined"] for x in d["variants"] if x["variant"] == "dense/sym"), f"{sv} (max {mv:.1e})",
            all(x["sound"] for x in d["variants"]), f"{g['query_ms_p2p']:.0f} ms / {100 * g['err_mean']:.4f}%" if g and tag == "kd" else "")

print("\n## 6. Build repeats (one thread, type_d nodes)\n")
head("index", "builds (h)", "median", "spread (max / min)")
for T in ("TJ", "GI"):
    for tag in ("kd", "k64"):
        b = []
        for p in [os.path.join(R, f"ours.{T}.{tag}.json")] + [os.path.join(RV, f"rep{r}.{T}.{tag}.json") for r in (1, 2)]:
            d = J(p)
            if d and d.get("complete", True):
                b.append(d["build_s"] / 3600)
        if b:
            row(f"{T} {tag}", ", ".join(f"{x:.1f}" for x in b), f"{statistics.median(b):.1f}", f"{max(b) / min(b):.2f}")

print("\n## 4. Strict soundness (no tolerance) over every review run of our index\n")
tot = viol = 0; mx = 0.0; nfiles = 0
for p in sorted(glob.glob(os.path.join(RV, "adv.*.json")) + glob.glob(os.path.join(RV, "rep*.json")) + glob.glob(os.path.join(RV, "small.*.json"))):
    d = J(p)
    if not d or "variants" not in d or "strict_violations" not in d["variants"][0]:
        continue
    sv, mv, n = strict(d); tot += n; viol += sv; mx = max(mx, mv); nfiles += 1
print(f"{nfiles} runs, {tot} answers, {viol} outside the raw interval, largest relative violation {mx:.2e}")
