"""Figures of the cluster evaluation (results/server). Run from the repository root:
    MPLBACKEND=Agg python3 paper/preprint/make_server_figures.py
fig_server_methods.pdf  TJ and GI at 0.5M faces: query time and index size against mean error, every method
fig_scaling.pdf         build time and index size (float32 labels) against n, 0.5M to 2.5M faces, with fitted exponents
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.path.join(ROOT, "results", "server")
OUT = os.path.join(ROOT, "paper", "preprint", "figures")
plt.rcParams.update({"font.size": 8, "font.family": "serif", "axes.linewidth": 0.6, "pdf.fonttype": 42})
# categorical order of the paper (validated: lightness band, chroma, CVD and contrast all pass); markers differ too
C = {"ours": "#1b9e77", "EAR-Oracle": "#d95f02", "DGG": "#7570b3", "NeuroGF": "#b8860b"}


def J(name):
    p = os.path.join(R, name)
    return json.load(open(p)) if os.path.exists(p) else None


def var(d, variant, s="random"):
    return next((v for v in d["variants"] if v["variant"] == variant and v["set"] == s), None)


def front(pts):
    """Pareto front (smaller x and smaller y), sorted by x."""
    pts = sorted(pts)
    out, best = [], np.inf
    for x, y in pts:
        if y < best:
            out.append((x, y)); best = y
    return out


# ---------------------------------------------------------------- methods at 0.5M faces
fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.6))
for col, T in enumerate(("TJ", "GI")):
    M = f"{T}pre_500000"
    ours = []  # (query us, size MB, mean err %, capped)
    for f in ["ours.%s.kd" % T, "ours.%s.e01" % T, "ours.%s.k64" % T] + [f"curve.{T}.{e}" for e in
                                                                         ("e0.05.g0.25", "e0.1.g1", "e0.4.g1", "e0.4.g2", "e0.8.g4", "e1.6.g8")] + \
            [f"caps.{T}.k{k}" for k in (8, 16, 32, 128, 256)]:
        d = J(f + ".json")
        if not d:
            continue
        for vv in ("dense/sym", "code8/sym"):
            v = var(d, vv)
            if v:
                ours.append((v["query_us"], v["index_bytes"] / 1e6, 100 * v["est_mean"], d["k_max"] > 0))
    comp = {"EAR-Oracle": [], "DGG": [], "NeuroGF": []}
    for r in J(f"ear.{M}.V2V.json") or []:
        if r.get("complete"):
            comp["EAR-Oracle"].append((r["query_us_mean"], r["index_MB"], r["err_mean"], r["log"].split("_sp")[1].split(".")[0]))
    for E, lab in (("0.001", r"$10^{-3}$"), ("0.0001", r"$10^{-4}$")):
        g = J(f"dgg.{M}.{E}.json")
        if g:
            comp["DGG"].append((1e3 * g["query_ms"], g["graph_bytes"] / 1e6, 100 * g["err_mean"], lab))
    for tag, lab in (("", "400 ep."), (".epochs800", "800 ep.")):
        n = J(f"neurogf.{M}{tag}.json")
        if n:
            s = {x["set"]: x for x in n["sets"]}
            comp["NeuroGF"].append((n["query_us"], n["index_bytes"] / 1e6, 100 * s["V2V"]["err_mean"], lab))
    vtp = J(f"onthefly.{T}_500000.json")
    for row, (xi, xlab) in enumerate(((0, r"query time ($\mu$s, one CPU thread)"), (1, "index size (MB)"))):
        ax = axes[row, col]
        fr = front([(p[xi], p[2]) for p in ours])
        ax.plot([p[xi] for p in ours if not p[3]], [p[2] for p in ours if not p[3]], "D", ms=3.5, color=C["ours"], alpha=0.35, mew=0)
        ax.plot([p[xi] for p in ours if p[3]], [p[2] for p in ours if p[3]], "D", ms=3.5, mfc="none", color=C["ours"], alpha=0.5, mew=0.6)
        ax.plot([x for x, _ in fr], [y for _, y in fr], "-D", ms=4, lw=1.0, color=C["ours"], label="ours (Pareto front)")
        for name, mk in (("EAR-Oracle", "^"), ("DGG", "s"), ("NeuroGF", "o")):
            pts = comp[name]
            ax.plot([p[xi] for p in pts], [p[2] for p in pts], mk, ms=6, color=C[name], label=name, ls="none")
            for p in pts:
                lab = ("sp " + p[3]) if name == "EAR-Oracle" else p[3]
                off = (5, -9) if ("800" in lab or "sp 7" in lab) else (5, 2)
                ax.annotate(lab, (p[xi], p[2]), textcoords="offset points", xytext=off, fontsize=6, color="0.25")
        if xi == 0 and vtp:
            q = 1e3 * vtp["vtp_full_ms"]
            ax.axvline(q, color="0.55", lw=0.7, ls="--")
            ax.text(q * 0.75, 1.5e-4, "exact VTP sweep", rotation=90, color="0.35", fontsize=6, va="bottom", ha="right")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylim(5e-5, 20)
        ax.set_xlim(*((0.1, 3e7) if xi == 0 else (0.1, 1e6)))
        ax.set_xlabel(xlab)
        ax.grid(alpha=0.3, lw=0.4)
        if col == 0:
            ax.set_ylabel("mean relative error (%)")
        if row == 0:
            ax.set_title(f"{T}, 0.5M faces ({J(f'ours.{T}.kd.json')['n']:,} vertices)", fontsize=8)
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.02))
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig_server_methods.pdf"), bbox_inches="tight")

# ---------------------------------------------------------------- scaling
SIZES = [500000, 1002528, 1503378, 2000000, 2504322]
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
ls = {"TJ": "-", "GI": "--"}
for T in ("TJ", "GI"):
    for tag, lab, color, mk in (("kd", "uncapped", C["ours"], "D"), ("k64", "cap 64", C["ours"], "o"), ("shifted", "shifted", "0.45", "s")):
        rows = [(J(f"scale.{T}.{S}.{tag}.json"), S) for S in SIZES]
        rows = [(d, S) for d, S in rows if d and d.get("variants")]
        if not rows:
            continue
        n = np.array([d["n"] for d, _ in rows], float)
        b = np.array([d["build_s"] / 3600 for d, _ in rows])
        sz = np.array([(var(d, "dense/sym") or d["variants"][0])["index_bytes"] / 1e9 for d, _ in rows])
        kb = np.polyfit(np.log(n), np.log(b), 1)[0] if len(n) > 2 else None
        ks = np.polyfit(np.log(n), np.log(sz), 1)[0] if len(n) > 2 else None
        mfc = color if mk == "D" else "none"
        axes[0].plot(n, b, ls[T] + mk, color=color, mfc=mfc, ms=4, lw=1, label=f"{T} {lab}" + (f" ($n^{{{kb:.2f}}}$)" if kb else ""))
        axes[1].plot(n, sz, ls[T] + mk, color=color, mfc=mfc, ms=4, lw=1, label=f"{T} {lab}" + (f" ($n^{{{ks:.2f}}}$)" if ks else ""))
g = [(J(f"dgg.TJpre_{S}.0.0001.json"), J(f"scale.TJ.{S}.kd.json")) for S in SIZES]
g = [(a, d) for a, d in g if a and d]
axes[1].plot([d["n"] for _, d in g], [a["graph_bytes"] / 1e9 for a, _ in g], ":s", color=C["DGG"], ms=4, lw=1, label=r"TJ DGG $10^{-4}$")
for ax, yl in zip(axes, ("build time (h, one thread)", "index size (GB, float32 labels)")):
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("vertices $n$"); ax.set_ylabel(yl)
    ax.grid(alpha=0.3, lw=0.4, which="both")
    ax.legend(fontsize=5.8, frameon=False, ncol=1)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_scaling.pdf"), bbox_inches="tight")
print("written", os.path.join(OUT, "fig_server_methods.pdf"), os.path.join(OUT, "fig_scaling.pdf"))
