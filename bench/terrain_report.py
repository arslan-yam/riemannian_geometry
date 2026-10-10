#collect results/<mesh>.json (ours) and results/<mesh>_competitors.json (their harness, evaluated
#against the same pygeodesic truth) into one table per mesh and a pareto figure
import json, os, sys, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


def load(mesh, out_dir=RESULTS, kind="V2V"):
    #kind "V2V": vertex-pair queries; "A2A": arbitrary-point queries (our rows are the ones
    #tagged "A2A", theirs come from <mesh>_competitors_a2a.json)
    tag = mesh + ("_weighted" if kind == "W2V" else "")
    ours = json.load(open(os.path.join(out_dir, f"{tag}.json")))
    comp_path = os.path.join(out_dir, f"{mesh}_competitors{ {'V2V': '', 'A2A': '_a2a', 'W2V': '_w2v'}[kind] }.json")
    comp = json.load(open(comp_path)) if os.path.exists(comp_path) else []
    ours = dict(ours, rows=[r for r in ours["rows"] if ("A2A" in r["method"]) == (kind == "A2A")])
    rows = []
    def order(r): #graph baselines, then the cut index by k, then query-time refinement
        m = r["method"]
        if m.startswith("dijkstra on"): return (0, 0, m)
        if m == "h2h": return (1, 0, m)
        if m.startswith("vertex portals"): return (2, 0, m)
        if m.startswith("cut portals"):
            k = int(m.split("k=")[1].split(",")[0])
            return (3 if "exact" in m else 4, k + (0.5 if "k~" in m else 0) + (0.1 if "interp" in m else 0) + (0.05 if "model" in m else 0), m)
        if m.startswith("Steiner graph"):
            return (5, 0, m)
        if m.startswith("certified"): #by eps (largest first), then oracle, UB before estimate
            eps = float(m.split("eps=")[1].split()[0])
            return (6, -eps, m.split(",")[1] if "," in m else "", "estimate" in m, m)
        return (7, 0, m)
    for r in sorted(ours["rows"], key=order):
        rows.append(dict(method=r["method"], build_s=r["build_s"], index_MB=r["index_bytes"] / 1e6,
                         query_us=r["query_us"], err_mean=r["err_mean"], err_p90=r["err_p90"], err_max=r["err_max"],
                         who="cert" if r["method"].startswith("certified") else "ours", row=r))
    names = dict(EAR="EAR-Oracle", SE="SE-Oracle", KAlgo="K-Algo", FixedS="fixed Steiner (no index)",
                 UnfixedS="unfixed Steiner (no index)", MMP="MMP (exact, no index)")
    for r in comp:
        tag = f"{names.get(r['method'], r['method'])}"
        if r["method"] in ("EAR", "SE", "FixedS"):
            tag += f" sp={r.get('sp', 5)}"
        if r["method"] in ("KAlgo", "UnfixedS", "SE", "EAR"):
            tag += f" eps={r['eps']}"
        if r.get("queries", ours["queries"]) < ours["queries"]:
            tag += f" ({r['queries']} queries)"
        rows.append(dict(method=tag, build_s=r["build_s"], index_MB=r["index_MB"], query_us=r["query_us"],
                         err_mean=r["err_mean"], err_p90=r["err_p90"], err_max=r["err_max"], who="theirs"))
    return ours, rows


def table(mesh, out_dir=RESULTS, kind="V2V"):
    ours, rows = load(mesh, out_dir, kind)
    what = {"V2V": "vertex-to-vertex", "A2A": "arbitrary-point", "W2V": "vertex-to-vertex, weighted faces"}[kind]
    lines = [f"### {mesh}: {ours['n']} vertices, {ours['faces']} faces, {ours['queries']} {what} queries"
             + (f", exact geodesic (pygeodesic) {ours['exact_ms']:.0f} ms/query" if kind == "V2V" else "")
             + (" (errors relative to the best path any method found)" if kind == "W2V" else ""),
             "", "| method | build (s) | index (MB) | query (us) | mean err | p90 err | max err |", "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['method']} | {r['build_s']:.1f} | {r['index_MB']:.2f} | {r['query_us']:.1f} | "
                     f"{r['err_mean']:.3f}% | {r['err_p90']:.3f}% | {r['err_max']:.2f}% |")
    return "\n".join(lines)


def cert_table(meshes, out_dir=RESULTS):
    #the certified index on every mesh: one line per (mesh, eps/g), estimate row + its UB twin
    lines = ["| mesh | n | eps / g | build (s) | index (MB) | query C / Python (us) | mean err | p90 err | max err | max UB err | LB <= d <= UB | UB/LB width (median) |",
             "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|"]
    for mesh in meshes:
        try:
            ours = json.load(open(os.path.join(out_dir, f"{mesh}.json")))
        except FileNotFoundError:
            continue
        rows = {r["method"]: r for r in ours["rows"]}
        ests = sorted([m for m in rows if m.startswith("certified") and m.endswith("(estimate)")],
                      key=lambda m: -float(m.split("eps=")[1].split()[0]))
        for m in ests:
            r = rows[m]; ub = rows.get(m.replace("(estimate)", "(UB)"), {})
            c = r.get("certificate", {})
            cfg = m.split("certified ")[1].split(",")[0].replace("eps=", "").replace(" g=", " / ")
            qp = r.get("query_us_python")
            q = f"{r['query_us']:.1f} / {qp:.1f}" if qp is not None else f"{r['query_us']:.1f}"
            lines.append(f"| {mesh} | {ours['n']} | {cfg} | {r['build_s']:.1f} | {r['index_bytes'] / 1e6:.2f} | {q} | "
                         f"{r['err_mean']:.4f}% | {r['err_p90']:.3f}% | {r['err_max']:.3f}% | "
                         f"{ub.get('err_max', float('nan')):.3f}% | {'yes' if c.get('sound') else ('n/a' if not c else 'NO')} | "
                         f"{c.get('width_median', float('nan')):.2f}% |")
    return "\n".join(lines)


def pareto(meshes, out_dir=RESULTS, path=os.path.join(RESULTS, "terrain_pareto.png"), kind="V2V"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    markers = {"ours": "o", "theirs": "s", "cert": "D"}
    colors = {"ours": "C0", "theirs": "C3", "cert": "C2"}
    for mesh in meshes:
        ours, rows = load(mesh, out_dir, kind)
        for r in rows:
            if r["err_mean"] <= 0 or "no index" in r["method"] and "MMP" in r["method"]:
                continue
            if kind == "A2A" and r["who"] == "ours" and "snap" in r["method"]:
                continue
            if r["who"] == "cert" and "(UB)" in r["method"]:
                continue
            c = colors[r["who"]]
            lab = r["method"].split(",")[0] if r["who"] != "theirs" else r["method"].split(" sp")[0].split(" eps")[0]
            axes[0].scatter(r["query_us"], r["err_mean"], marker=markers[r["who"]], color=c, s=40)
            axes[0].annotate(lab, (r["query_us"], r["err_mean"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
            if r["build_s"] > 0:
                axes[1].scatter(r["build_s"], r["err_mean"], marker=markers[r["who"]], color=c, s=40)
                axes[1].annotate(lab, (r["build_s"], r["err_mean"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    axes[0].set_xlabel("query time (us)"); axes[1].set_xlabel("build time (s)")
    for ax in axes:
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylabel("mean relative error (%)"); ax.grid(alpha=0.3)
    fig.suptitle(f"terrain benchmark, {'vertex-to-vertex' if kind == 'V2V' else 'arbitrary-point'} queries: "
                 "certified index (green diamonds), earlier cut index (blue circles), terrain oracles (red squares)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    return path


def scaling(meshes, out_dir=RESULTS, path=os.path.join(RESULTS, "terrain_scaling.png")):
    #build time, index size and mean error against the number of vertices: the certified index per
    #eps/g, the earlier cut index (k=32, model interp) and EAR-Oracle (sp=5, eps=0.2)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    series = {}
    def put(name, n, r):
        series.setdefault(name, []).append((n, r["build_s"], r["index_MB"], r["err_mean"]))
    for mesh in meshes:
        try:
            ours, rows = load(mesh, out_dir)
        except FileNotFoundError:
            continue
        n = ours["n"]
        for r in rows:
            m = r["method"]
            if r["who"] == "cert" and m.endswith("(estimate)"):
                put("certified " + m.split("certified ")[1].split(",")[0], n, r)
            elif m == "cut portals k=32, exact MMP labels + model interp":
                put("earlier cut index k=32 (model interp)", n, r)
            elif m.startswith("EAR-Oracle sp=5 eps=0.2"):
                put("EAR-Oracle sp=5", n, r)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for name, pts in sorted(series.items()):
        pts = sorted(pts)
        style = dict(marker="D" if name.startswith("certified") else ("s" if "EAR" in name else "o"),
                     ls="-" if name.startswith("certified") else "--")
        for ax, k in zip(axes, (1, 2, 3)):
            ax.plot([p[0] for p in pts], [max(p[k], 1e-6) for p in pts], label=name, **style)
    for ax, lab in zip(axes, ("build time (s), single thread", "index size (MB)", "mean relative error (%)")):
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("vertices"); ax.set_ylabel(lab); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    return path


def write_note(meshes, out_dir=RESULTS, path=os.path.join(ROOT, "notes", "terrain_results.md")):
    #the full results note: setup, one table per mesh and query kind, the pareto figures, caveats
    fig = pareto(meshes, out_dir)
    fig_a = pareto(meshes, out_dir, path=os.path.join(out_dir, "terrain_pareto_a2a.png"), kind="A2A")
    parts = ["# Terrain benchmark: Certified Cut Index vs terrain distance oracles", "",
             "The certified index (`rgeo/oracles/certified.py`, theorem and proofs in `notes/guarantee.md`)",
             "returns with every answer an interval LB <= d <= UB with (1-eps) d <= LB and UB <= (1+eps) d;",
             "`(UB)` rows answer with the upper bound (the length of an actual path), `(estimate)` rows with",
             "the interpolated minimum clamped into the interval. Its query is timed in C (`rgeo/oracles/",
             "cquery.py`, identical answers) and in Python. Rows `cut portals k=...` are the earlier,",
             "uncertified index, kept for comparison.", "",
             "## Certified index on every mesh", "", cert_table(meshes, out_dir), "",
             f"![scaling]({os.path.relpath(scaling(meshes, out_dir), os.path.dirname(path))})", "",
             "Every method answers the same 100 random vertex-to-vertex queries and the same 100 random",
             "arbitrary-point queries per mesh (uniform random face, uniform barycentric coordinates); the truth",
             "is the exact polyhedral geodesic (MMP, `pygeodesic`, double precision, query points inserted as",
             "Steiner vertices). Errors are relative and signed: the competitors (K-Algo, SE-Oracle) can answer",
             "below the truth, ours is >= 0 except the interpolated variants. `build` is single-threaded wall",
             "time on an Apple M1 for everyone (ours: Python + pygeodesic sweeps; theirs: their C++ harness).",
             "`index` is our label bytes vs. their reported resident-memory delta (which includes their",
             "Steiner base graph). Our `k` is portals per cut, uniform across the tree; `k~|cut|` is the budget",
             "proportional to cut length (`k_min=4`). Arbitrary points are attached to the vertex index either",
             "by `snap` (through the face's vertices, an admissible upper bound, first order in the edge",
             "length) or by `interp` (bilinear interpolation of the 3x3 vertex-pair answers, second order,",
             "two-sided); pairs within 4 rings are answered exactly by a local sweep.", "",
             "## Vertex-to-vertex queries", "",
             f"![pareto]({os.path.relpath(fig, os.path.dirname(path))})", ""]
    for m in meshes:
        try:
            parts += [table(m, out_dir), ""]
        except FileNotFoundError:
            parts += [f"### {m}: no results", ""]
    parts += ["## Arbitrary-point queries", "", f"![pareto]({os.path.relpath(fig_a, os.path.dirname(path))})", ""]
    for m in meshes:
        try:
            parts += [table(m, out_dir, kind="A2A"), ""]
        except FileNotFoundError:
            parts += [f"### {m}: no results", ""]
    weighted = [m for m in meshes if os.path.exists(os.path.join(out_dir, f"{m}_weighted.json"))]
    if weighted:
        parts += ["## Weighted terrains (face costs uniform in (1, 5))", "",
                  "No exact solver exists for weighted terrains, so the reference is the best path any method",
                  "found for each query (theirs and ours; their K-Algo is excluded because it ignores the face",
                  "weights: its weighted answers match the unweighted exact distances to 2.7%). The certified",
                  "index runs on Steiner-graph labels here (`rgeo/terrain/steiner.py`, m points per edge): its UB",
                  "is the cost of a real path, UB <= (1+eta)(1+eps) d with eta the Steiner graph's own error, and",
                  "the interval lower bound is not available (LB = chord). `refined labels` rows are the earlier",
                  "index on the vertex mesh graph.", ""]
        for m in weighted:
            parts += [table(m, out_dir, kind="W2V"), ""]
    parts += ["## Caveats", "",
              "* The competitors' harness is `float32`; meshes were translated to the origin (on the raw 1e7 m",
              "  coordinates its own MMP disagreed with the double-precision truth by up to 200%), and query",
              "  points were nudged 5 cm into a face because its AABB point location fails at vertices. Both",
              "  are applied identically to every competitor; the nudge moves the truth by < 0.005%.",
              "* BigMount (1 face) and GunnisonForest (6 faces) contain folded triangles; `repair_height_field`",
              "  slides one vertex of each fold by 20-80 m (median edge 190-450 m). All methods see the repaired mesh.",
              "* `Gunnison10k` is a 46 x 35 km crop of GunnisonForest at native resolution (9,610 vertices);",
              "  `GunnisonForest100k` is the full dataset (100,511 vertices, 147 x 112 km). At 100k the EAR-Oracle",
              "  and K-Algo runs exited without writing query results on the 8 GB machine (EAR: 51,512 s wall,",
              "  6,422 s CPU for V2V before exiting), so only MMP and the fixed Steiner scheme are reported there.",
              "* SE-Oracle was built while other jobs ran (wall 6252 s, CPU 1149 s), as were the EAR rows with",
              "  eps=0.05 or sp=10/20 on HorseMount (their sp=5/eps=0.2 twin took 59 s under that load and 23.5 s",
              "  on the quiet machine, which is the number in the table). K-Algo on the larger meshes ran on the",
              "  first 20 queries only (seconds to minutes per query).",
              "* Our builds are single-threaded here; `SeparatorCutIndex(workers=N)` forks the per-cut sweeps",
              "  and gives ~3x on 4 performance cores, but that is not what the tables show.",
              "* `+ model interp` is the interpolation that removes the explicit sqrt(r^2 + (x - x_p)^2) term of",
              "  an endpoint close to the cut before fitting (notes/error_bound.md, 5b); labels are float32.",
              "* On the weighted terrain the vertex mesh graph is the wrong skeleton: weighted geodesics refract",
              "  at edges and hug cheap ones, so the earlier index (refined mesh-graph labels) inherits the",
              "  graph's ~9% error. The certified index on a Steiner base graph removes it.",
              "* EAR-Oracle's error is that of its fixed Steiner scheme (`--sp-num`), not of `--eps`: on HorseMount",
              "  sp=5/10/20 give 3.1/2.5/2.3% at 59/156/484 s build, and eps=0.05 changes nothing at sp=5.", ""]
    open(path, "w").write("\n".join(parts))
    return path


if __name__ == "__main__":
    meshes = sys.argv[1:] or sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(RESULTS, "*.json")) if "_competitors" not in p)
    for m in meshes:
        print(table(m)); print()
