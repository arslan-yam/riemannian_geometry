"""Size experiments (results/size/*.json from scripts/run_size.sh): tables and the size-versus-error
curve against EAR-Oracle. Writes notes/size_results.md and results/size_error.png.
usage: python3 bench/size_report.py
"""
import glob, json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")


def load(p):
    try:
        return json.load(open(p))
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def rows_of(d, tag):
    out = []
    for v in d["variants"]:
        out.append(dict(config=tag, variant=v["variant"] + ("+refine" if v["refine"] else ""), build=d["build_s"],
                        mb=v["index_bytes"] / 1e6, q=v["query_us"], qr=v["query_us_refine_mean"], qr99=v["query_us_refine_p99"],
                        refined=v["refined"], pairs=v["pairs"], sound=v["sound"], ub=v["max_ub_err"] * 100,
                        mean=v["est_mean"] * 100, mx=v["est_max"] * 100, ratio=v["max_ub_over_lb"], portals=d["portals"]))
    return out


def table(rows):
    out = ["| configuration | storage | build | index | query | query with refine (mean / p99) | refined | sound | max UB err | mean / max error | max UB/LB |",
           "|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|"]
    for r in rows:
        out.append(f"| {r['config']} | {r['variant']} | {r['build']:.1f} s | {r['mb']:.2f} MB | {r['q']:.2f} µs | {r['qr']:.2f} / {r['qr99']:.1f} µs | "
                   f"{r['refined']}/{r['pairs']} | {'yes' if r['sound'] else 'NO'} | {r['ub']:.3f}% | {r['mean']:.4f}% / {r['mx']:.3f}% | {r['ratio']:.4f} |")
    return "\n".join(out)


def main():
    parts = ["# Size experiments", "",
             "One build per configuration (C++, hybrid sweeps, single thread, Apple M1), every storage evaluated",
             "on it. `kd` = the current index; `kmax` = version A (portals per cut capped, certify-or-refine);",
             "`shifted` = 3 shifted quadtrees with gap max(g, eps lambda S / 6); storages: `dense` float32,",
             "`code16` / `code8` closed-loop delta coding, `coverX` = X-covers per row. Errors are |estimate/d - 1|",
             "against exact MMP; `refined` = pairs answered by the exact local sweep because UB > (1+eps) LB.",
             "The size-versus-error section has two series of (eps, g in edges): (0.1, 0.5), (0.2, 1), (0.3, 2),",
             "(0.5, 4), (1, 8), where the near radius R = g/(eps lambda) grows from 5 to 8 edges, and (0.4, 2), (0.8, 4),",
             "(1.6, 8), (3.2, 16) with R fixed at 5 edges. The plot shows every configuration and each storage's",
             "Pareto front.", ""]
    for mesh, tag, pairs in (("Gunnison10k", "g10k", "1,500 pairs, eps 0.1, g 0.5 edge"), ("GunnisonForest100k", "g100k", "200 pairs, eps 0.2, g 1 edge")):
        rows = []
        for name in ("kd", "kmax200", "kmax100", "kmax50", "shifted"):
            d = load(os.path.join(R, "size", f"{tag}.{name}.json"))
            if d:
                rows += rows_of(d, name if name != "kd" else "kd (current)")
        if rows:
            parts += [f"## {mesh} ({pairs})", "", table(rows), ""]
    curve = []
    for p in sorted(glob.glob(os.path.join(R, "size", "curve.*.json"))):
        d = load(p)
        if d:
            mesh = os.path.basename(p).split(".")[1]
            for r in rows_of(d, f"eps {d['eps']:g}, g {d['gap_factor']:g}e"):
                r["mesh"] = mesh
                curve.append(r)
    if curve:
        parts += ["## Size versus error (kd, no cap)", ""]
        for mesh in ("HorseMount", "RobinsonMountain", "Gunnison10k"):
            rs = [r for r in curve if r["mesh"] == mesh]
            if rs:
                parts += [f"### {mesh}", "", table(rs), ""]
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
            for ax, mesh in zip(axes, ("HorseMount", "RobinsonMountain", "Gunnison10k")):
                for k, (var, mk) in enumerate((("dense", "o"), ("code16", "s"), ("code8", "D"))):
                    rs = sorted([r for r in curve if r["mesh"] == mesh and r["variant"] == var], key=lambda r: (r["mb"], r["mean"]))
                    if not rs:
                        continue
                    # every configuration as a faint point, the Pareto front (smallest error for its size) as the line
                    ax.scatter([r["mb"] for r in rs], [r["mean"] for r in rs], marker=mk, color=f"C{k}", alpha=0.3, s=18)
                    front, best = [], float("inf")
                    for r in rs:
                        if r["mean"] < best:
                            front.append(r); best = r["mean"]
                    ax.plot([r["mb"] for r in front], [r["mean"] for r in front], marker=mk, color=f"C{k}", label=f"ours, {var} (Pareto front)")
                comp = load(os.path.join(R, f"{mesh}_competitors.json")) or []
                ear = sorted([c for c in comp if c["method"] == "EAR" and c["eps"] == 0.2], key=lambda c: c["index_MB"])
                if ear:
                    ax.plot([c["index_MB"] for c in ear], [abs(c["err_mean"]) for c in ear], marker="^", ls="--", color="C3", label="EAR-Oracle (sp 5/10/20)")
                ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=0.3)
                ax.set_xlabel("index size (MB)"); ax.set_ylabel("mean relative error (%)"); ax.set_title(mesh)
            axes[0].legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(os.path.join(R, "size_error.png"), dpi=200)
            parts += ["![size versus error](../results/size_error.png)", ""]
        except ImportError:
            pass
    open(os.path.join(ROOT, "notes", "size_results.md"), "w").write("\n".join(parts))
    print("wrote notes/size_results.md")


if __name__ == "__main__":
    main()
