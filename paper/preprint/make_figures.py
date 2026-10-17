"""Figures of the preprint. Run from the repository root:  MPLBACKEND=Agg python3 paper/preprint/make_figures.py
fig_methods.pdf   query time against mean error, every method we ran (notes/unified_results.md, Apple M1, one thread)
fig_obstacles.pdf the three obstacle meshes of the pilot (Gunnison10k; scripts/competitors/make_obstacles.py, seed 1)
fig_size_error    copied from results/size_error.png (bench/size_report.py)
"""
import os, shutil, subprocess, sys, tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "paper", "preprint", "figures")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, ROOT)
plt.rcParams.update({"font.size": 8, "font.family": "serif", "axes.linewidth": 0.6, "pdf.fonttype": 42})

# ---------------------------------------------------------------- query time against mean error
# (query time in microseconds, mean relative error in %), from notes/unified_results.md
M = {
    "HorseMount (793 vertices)": {
        "ours": [(0.22, 0.0073, r"$\varepsilon$=0.05"), (0.18, 0.0235, r"$\varepsilon$=0.1"), (0.20, 0.0676, r"$\varepsilon$=0.2")],
        "EAR-Oracle": [(1300, 3.10, "sp 5"), (5700, 2.53, "sp 10"), (5700, 2.32, "sp 20")],
        "SE-Oracle": [(11600, 0.814, "sp 5")],
        "DGG": [(1400, 0.0023, r"$10^{-4}$"), (724, 0.0279, r"$10^{-3}$")],
        "K-Algo": [(8.1e6, 0.0178, "")],
        "Steiner (on the fly)": [(26400, 3.11, "sp 5"), (208300, 2.31, "sp 20")],
        "NeuroGF (short training)": [(280, 0.458, "")],
        "exact": [(3300, None, "VTP"), (9000, None, "MMP")],
    },
    "Gunnison10k (9,610 vertices)": {
        "ours": [(0.61, 0.0024, r"$\varepsilon$=0.05"), (0.36, 0.0076, r"$\varepsilon$=0.1"), (0.24, 0.0225, r"$\varepsilon$=0.2")],
        "EAR-Oracle": [(10200, 3.43, "sp 5")],
        "DGG": [(22300, 0.0028, r"$10^{-4}$"), (10700, 0.0285, r"$10^{-3}$")],
        "K-Algo": [(80.7e6, 0.0134, "")],
        "Steiner (on the fly)": [(161500, 3.37, "sp 5")],
        "exact": [(55000, None, "VTP")],
    },
}
style = {"ours": ("D", "#1b9e77", 6), "EAR-Oracle": ("^", "#d95f02", 6), "SE-Oracle": ("v", "#e7298a", 6),
         "DGG": ("s", "#7570b3", 5), "K-Algo": ("P", "#66a61e", 6), "Steiner (on the fly)": ("x", "#a6761d", 6),
         "NeuroGF (short training)": ("o", "#e6ab02", 5)}
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
for ax, (mesh, rows) in zip(axes, M.items()):
    for name, pts in rows.items():
        if name == "exact":
            for q, _, lab in pts:
                ax.axvline(q, color="0.55", lw=0.7, ls="--")
                ax.text(q * 1.15, 6e-4 * 1.3, f"exact {lab}", rotation=90, color="0.4", fontsize=6.5, va="bottom")
            continue
        mk, c, s = style[name]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], mk, color=c, ms=s, mfc=c if mk not in "x" else None,
                label=name, ls="-" if name == "ours" else "none", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.08, 3e8); ax.set_ylim(6e-4, 10)
    ax.set_xlabel(r"query time ($\mu$s)"); ax.set_title(mesh, fontsize=8)
    ax.grid(alpha=0.3, lw=0.4)
axes[0].set_ylabel("mean relative error (%)")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=7, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 1.03))
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(os.path.join(OUT, "fig_methods.pdf"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- obstacle meshes of the pilot
from rgeo.terrain import read_off
tmp = tempfile.mkdtemp()
src = os.path.join(ROOT, "data", "terrains", "Gunnison10k.off")
kinds = [("blobs10", ["--kind", "blobs", "--cover", "0.1"], "disks, 10% of the area"),
         ("blobs20", ["--kind", "blobs", "--cover", "0.2"], "disks, 20% of the area"),
         ("walls", ["--kind", "walls", "--walls", "6"], "six walls with gaps")]
fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.0))
for ax, (k, args, title) in zip(axes, kinds):
    out = os.path.join(tmp, k + ".off")
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "competitors", "make_obstacles.py"), src, out, *args, "--seed", "1"],
                   check=True, capture_output=True)
    V, F = read_off(out); V = np.asarray(V); F = np.asarray(F)
    ax.tripcolor(V[:, 0], V[:, 1], F, V[:, 2], shading="gouraud", cmap="terrain", rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_obstacles.pdf"), bbox_inches="tight", dpi=250)
plt.close(fig)
shutil.rmtree(tmp)

shutil.copy(os.path.join(ROOT, "results", "size_error.png"), os.path.join(OUT, "fig_size_error.png"))
print("figures written to", OUT)
