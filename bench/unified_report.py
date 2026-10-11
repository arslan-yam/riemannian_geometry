"""One table per mesh with every method: our certified index in C++ (and Python), the EAR-Oracle
authors' harness (EAR, SE, K-Algo, fixed/unfixed Steiner, MMP), the UP-Oracle authors' code (UP,
WSPD, RC, EAR, WAV-Fly, ESP-Fly; POI setting), DGG, NeuroGF and on-the-fly exact solvers.
Writes notes/unified_results.md. usage: python3 bench/unified_report.py
"""
import glob, json, os, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
MESHES = ["HorseMount", "BigMount", "HeadlightMount", "RobinsonMountain", "Gunnison10k", "GunnisonForest100k"]


def jload(p):
    try:
        return json.load(open(p))
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def fmt_t(s):
    if s is None:
        return "—"
    return f"{s * 1e3:.0f} ms" if s < 1 else (f"{s:.1f} s" if s < 600 else f"{s / 3600:.2f} h" if s >= 3600 else f"{s / 60:.1f} min")


def fmt_q(us):
    if us is None:
        return "—"
    return f"{us:.2f} µs" if us < 1000 else (f"{us / 1e3:.1f} ms" if us < 1e6 else f"{us / 1e6:.1f} s")


def log_lines(name):
    p = os.path.join(R, "logs", name)
    return open(p).read().splitlines() if os.path.exists(p) else []


def ours_cpp(mesh):
    rows = []
    for p in sorted(glob.glob(os.path.join(R, "cpp", f"{mesh}.*.json"))):
        d = jload(p)
        if not d:
            continue
        tag = f"**ours, C++**, eps={d['eps']:g}, g={d['gap_factor']:g}e, {d['solver']} sweeps"
        rows.append(dict(method=tag, build=d["build_s"], size=d["index_bytes"] / 1e6, query=d["query_us"],
                         mean=d.get("est_mean", float("nan")) * 100, max=d.get("est_max", float("nan")) * 100,
                         note=f"certified: LB≤d≤UB on all {d['pairs']} pairs ({'sound' if d.get('sound') else 'NOT sound'}), max UB err {d.get('max_ub_err', 0) * 100:.2f}%",
                         pairs=d["pairs"], key=(0, -d["eps"], d["solver"])))
    return rows


def ours_python(mesh):
    d = jload(os.path.join(R, f"{mesh}.json"))
    rows = []
    if not d:
        return rows
    for r in d["rows"]:
        m = r["method"]
        if m.startswith("certified") and m.endswith("(estimate)"):
            c = next((x.get("certificate") for x in d["rows"] if x["method"] == m), None) or r.get("certificate") or {}
            rows.append(dict(method=f"ours, Python, {m.replace('certified ', '').replace(', exact MMP (estimate)', '')}", build=r["build_s"],
                             size=r["index_bytes"] / 1e6, query=r.get("query_us"), mean=abs(r["err_mean"]) if r["err_min"] >= 0 else r.get("abs_err_mean", abs(r["err_mean"])),
                             max=max(abs(r["err_max"]), abs(r["err_min"])), note="certified" + (f", max UB err {c.get('max_ub_err', 0):.2f}%" if c else ""),
                             pairs=d["queries"], key=(1, 0, m)))
        if "A2A certified" in m:
            c = r.get("certificate", {})
            rows.append(dict(method=f"ours, A2A certified, {m.split(',')[0].replace('certified ', '')}", build=r["build_s"], size=r["index_bytes"] / 1e6,
                             query=r.get("query_us"), mean=r.get("abs_err_mean", abs(r["err_mean"])), max=max(abs(r["err_max"]), abs(r["err_min"])),
                             note=f"A2A: sound {c.get('sound')}, max UB err {c.get('max_ub_err', 0):.2f}%, width median {c.get('width_median', 0):.2f}%",
                             pairs=100, key=(2, 0, m)))
    return rows


NAMES = dict(EAR="EAR-Oracle", SE="SE-Oracle", KAlgo="K-Algo", FixedS="fixed Steiner (on the fly)",
             UnfixedS="unfixed Steiner (on the fly)", MMP="MMP exact (on the fly, their CGAL code)")


def harness(mesh, kind=""):
    d = jload(os.path.join(R, f"{mesh}_competitors{kind}.json")) or []
    rows = []
    for r in d:
        tag = NAMES.get(r["method"], r["method"])
        if r["method"] in ("EAR", "SE", "FixedS"):
            tag += f", sp={r.get('sp', 5)}"
        if r["method"] in ("EAR", "SE", "KAlgo", "UnfixedS"):
            tag += f", eps={r['eps']}"
        rows.append(dict(method=tag + " [EAR authors' C++]", build=r["build_s"] or None, size=r["index_MB"] or None, query=r["query_us"],
                         mean=abs(r["err_mean"]), max=max(abs(r["err_max"]), abs(r["err_min"])), note="no per-query bound",
                         pairs=r.get("queries", 100), key=(3, 0, tag)))
    return rows


def dgg(mesh):
    rows, pre = [], {}
    for l in log_lines("baselines.log"):
        m = re.match(r"DGG precompute (\S+) eps=(\S+): ([\d.]+) s, graph (\d+) bytes", l)
        if m and m.group(1) == mesh:
            pre[m.group(2)] = (float(m.group(3)), int(m.group(4)) / 1e6)
        m = re.match(r"DGG (\S+)_FD([\d.]+)_c5\.binary: (\d+) queries, query ([\d.]+) ms, err mean ([-\d.]+)% max ([-\d.]+)% min ([-\d.]+)%", l)
        if m and m.group(1) == mesh:
            e = f"{float(m.group(2)):g}"
            b, s = pre.get(e, pre.get(m.group(2), (None, None)))
            rows.append(dict(method=f"DGG, eps={e} [authors' C++, macOS port]", build=b, size=s, query=float(m.group(4)) * 1e3,
                             mean=abs(float(m.group(5))), max=max(abs(float(m.group(6))), abs(float(m.group(7)))),
                             note="graph + Dijkstra per query; eps bounds the mean error", pairs=int(m.group(3)), key=(4, 0, e)))
    return rows


def neurogf(mesh):
    d = jload(os.path.join(R, "neurogf", f"{mesh}.json"))
    if not d:
        return []
    return [dict(method=f"NeuroGF (distance branch, {d['epochs']}x{d['iters']} iters) [authors' network]", build=d["build_s"], size=d["index_MB"],
                 query=d["query_us"], mean=d["err_mean"], max=d["err_max"], note=f"learned, no bound; batched {d['query_us_batch']:.1f} µs/pair",
                 pairs=1500, key=(5, 0, ""))]


def onthefly(mesh):
    rows = []
    for l in log_lines("cpp_runs.log"):
        if l.startswith(f"data/terrains/{mesh}.off: n="):
            m = re.search(r"full sweep: MMP ([\d.]+) ms, MMP-early\(all targets\) [\d.]+ ms, VTP ([\d.]+) ms", l)
            if m:
                rows.append(dict(method="MMP exact, one sweep per query [Kirsanov, C++]", build=None, size=None, query=float(m.group(1)) * 1e3,
                                 mean=0.0, max=0.0, note="exact, full sweep", pairs=5, key=(6, 0, "a")))
                rows.append(dict(method="VTP exact, one sweep per query [authors' C++]", build=None, size=None, query=float(m.group(2)) * 1e3,
                                 mean=0.0, max=0.0, note="exact, full sweep", pairs=5, key=(6, 1, "b")))
    return rows


def table(rows):
    out = ["| method | build | index | query | mean err | max err | pairs | notes |", "|---|---:|---:|---:|---:|---:|---:|---|"]
    for r in sorted(rows, key=lambda r: r["key"]):
        size = "—" if r["size"] is None else f"{r['size']:.2f} MB"
        out.append(f"| {r['method']} | {fmt_t(r['build'])} | {size} | {fmt_q(r['query'])} | {r['mean']:.4f}% | {r['max']:.3f}% | {r['pairs']} | {r['note']} |")
    return "\n".join(out)


def upo_section():
    out = []
    for tag, title in (("tj_poi50", "TJ (576 vertices), their 50 POIs"), ("tj_all", "TJ, every vertex a POI (576)"),
                       ("horse_poi50", "HorseMount, 50 random POIs"), ("horse_all", "HorseMount, every vertex a POI (793)")):
        p = os.path.join(R, "upo", f"{tag}.txt")
        if not os.path.exists(p):
            continue
        rows = []
        for l in open(p):
            m = re.match(r"RESULT (\S+)\s+build\s+([\d.]+) ms\s+size\s+([\d.]+) MB\s+query\s+([\d.]+) ms\s+err mean ([-\d.]+)% max ([-\d.]+)% min ([-\d.]+)%", l)
            if m:
                rows.append(dict(method=m.group(1) + " [UP-Oracle authors' C++]", build=float(m.group(2)) / 1e3, size=float(m.group(3)) or None,
                                 query=float(m.group(4)) * 1e3, mean=abs(float(m.group(5))), max=max(abs(float(m.group(6))), abs(float(m.group(7)))),
                                 note="POI oracle" if m.group(1) not in ("WAV-Fly", "ESP-Fly") else "on the fly", pairs=200, key=(3, 0, m.group(1))))
        for q in sorted(glob.glob(os.path.join(R, "upo", f"{tag}.cci.*.json"))):
            d = jload(q)
            if d:
                rows.append(dict(method=f"**ours, C++**, eps={d['eps']:g}, g={d['gap_factor']:g}e (all vertices indexed)", build=d["build_s"],
                                 size=d["index_bytes"] / 1e6, query=d["query_us"], mean=d["est_mean"] * 100, max=d["est_max"] * 100,
                                 note=f"certified, sound={d['sound']}, max UB err {d['max_ub_err'] * 100:.2f}%", pairs=d["pairs"], key=(0, -d["eps"], "")))
        out += [f"### {title}", "", table(rows), ""]
    return "\n".join(out)


def main():
    parts = ["# All methods, one table per mesh", "",
             "Single-threaded wall time on an Apple M1 (8 GB) for every method. `build` is index construction",
             "(on-the-fly methods: none), `query` is per pair, errors are relative to exact MMP (|error| for",
             "mean and max). Our C++ index and DGG/NeuroGF use the same 1,500 random vertex pairs per mesh (seed 7);",
             "the EAR authors' harness uses 100 vertex pairs (seed 0, points 5 cm inside a face, float32 code);",
             "on-the-fly exact solvers are timed on 5 full sweeps. Tags in brackets say whose code ran.",
             "For the harness rows the logs give the signed mean, shown in absolute value (K-Algo and SE-Oracle",
             "also answer below the truth, so their mean |error| is larger than shown).", ""]
    #summary: our main configuration against the strongest method of each class, per mesh
    summ = ["## Summary: ours (C++, eps=0.1, g=0.5e, hybrid) against the best of each class", "",
            "| mesh | n | ours: build / index / query / mean / max | EAR-Oracle sp=5 | DGG eps=1e-4 | exact on the fly (VTP) |",
            "|---|---:|---|---|---|---|"]
    for mesh in MESHES:
        o = [r for r in ours_cpp(mesh) if "eps=0.1," in r["method"] and "hybrid" in r["method"]]
        e = [r for r in harness(mesh) if r["method"].startswith("EAR-Oracle, sp=5, eps=0.2")]
        g = [r for r in dgg(mesh) if "eps=0.0001" in r["method"]]
        v = [r for r in onthefly(mesh) if r["method"].startswith("VTP")]
        if not o:
            continue
        def f(r):
            size = "—" if r["size"] is None else f"{r['size']:.1f} MB"
            return f"{fmt_t(r['build'])} / {size} / {fmt_q(r['query'])} / {r['mean']:.4f}% / {r['max']:.2f}%"
        n = jload(os.path.join(R, "cpp", f"{mesh}.hybrid.0.1.json")) or {}
        summ.append(f"| {mesh} | {n.get('n', '')} | {f(o[0])} | {f(e[0]) if e else '—'} | {f(g[0]) if g else '—'} | {fmt_q(v[0]['query']) + ' / exact' if v else '—'} |")
    parts += summ + ["", "Ours is certified on every query (LB <= d <= UB, sound on all pairs); none of the others",
                     "returns a per-query bound. Columns: build / index / query / mean |error| / max |error|.", ""]
    for mesh in MESHES:
        rows = ours_cpp(mesh) + ours_python(mesh) + harness(mesh) + dgg(mesh) + neurogf(mesh) + onthefly(mesh)
        if rows:
            parts += [f"## {mesh}", "", table(rows), ""]
    parts += ["## The UP-Oracle authors' oracles (POI setting)", "",
              "Their oracles index a set of points of interest (POIs); queries are between POIs (200 random pairs).",
              "Ours indexes every vertex and answers the same pairs. Their published A2A numbers are not",
              "measured but extrapolated (POI results multiplied by vertices/POIs, `algorithms.h`); the",
              "`every vertex a POI` rows are real runs.", "", upo_section()]
    w = jload(os.path.join(R, "HorseMount_competitors_w2v.json"))
    ow = jload(os.path.join(R, "HorseMount_weighted.json"))
    if w and ow:
        rows = harness("HorseMount", "_w2v")
        for r in ow["rows"]:
            if r["method"].startswith(("certified", "Steiner graph", "dijkstra on", "h2h")):
                rows.append(dict(method="ours: " + r["method"], build=r["build_s"] or None, size=r["index_bytes"] / 1e6 or None,
                                 query=r["query_us"], mean=abs(r["err_mean"]), max=max(abs(r["err_max"]), abs(r["err_min"])),
                                 note="weighted; UB is a path cost", pairs=100, key=(0, 0, r["method"])))
        parts += ["## Weighted HorseMount (face costs in (1, 5)), reference = best valid path any method found", "", table(rows), ""]
    open(os.path.join(ROOT, "notes", "unified_results.md"), "w").write("\n".join(parts))
    print("wrote notes/unified_results.md")


if __name__ == "__main__":
    main()
