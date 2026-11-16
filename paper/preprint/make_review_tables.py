"""LaTeX tables of the review experiments (results/server/review): adversarial terrains and landmark intervals.

  python3 paper/preprint/make_review_tables.py [results/server]   -> paper/preprint/tab_adversarial.tex, tab_landmarks.tex

Runs still missing print as "running". Reference failures (a pair the exact MMP sweep of the ground truth did not
reach, truth 1e100) are left out of our error statistics and counted separately.
"""
import json, os, sys

R = sys.argv[1] if len(sys.argv) > 1 else "results/server"
RV = os.path.join(R, "review")
OUT = os.path.dirname(os.path.abspath(__file__))


def J(p):
    return json.load(open(p)) if os.path.exists(p) else None


def var(d, variant, s="random"):
    return next((v for v in d["variants"] if v["variant"] == variant and v["set"] == s), None)


def pct(x, nd):
    return f"{100 * x:.{nd}f}\\%"


def h(s):
    return f"{s / 3600:.1f}\\,h"


def gb(b):
    return f"{b / 1e9:.2f}"


def unreached(mesh):
    """pairs of the random set whose MMP reference did not reach the target (truth 1e100)"""
    p = os.path.join(R, "comp", f"{mesh}.V2V.truth")
    return sum(1 for l in open(p) if l.strip() and float(l) > 1e90) if os.path.exists(p) else 0


def shortest(T, tag):
    """the shortest of the repeated builds (they are deterministic: same labels)"""
    ps = [os.path.join(R, f"ours.{T}.{tag}.json")] + ([os.path.join(R, f"caps.{T}.k64.json")] if tag == "k64" else []) \
        + [os.path.join(RV, f"rep{r}.{T}.{tag}.json") for r in (1, 2)]
    return min(d["build_s"] for d in map(J, ps) if d and d.get("complete", True))


# ---- adversarial terrains
KINDS = [("TJ", "original", None),
         ("steep5", "heights $\\times5$", "steep5"), ("steep10", "heights $\\times10$", "steep10"),
         ("aniso", "$x\\times10$", "aniso"), ("slivers", "needle triangles", "slivers"),
         ("ties", "heights to 5\\,m", "ties"), ("cluster", "clustered points", "cluster"),
         ("comb", "40 comb walls", "comb")]
BASE_STATS = {"slope_p50_deg": 22, "slope_p99_deg": 55, "aspect_p99": 2.6}  # TJpre_500000 (see notes)
rows = []
for name, what, K in KINDS:
    if K is None:
        mj = dict(BASE_STATS, vertices=251001)
        kd, k64 = J(os.path.join(R, "ours.TJ.kd.json")), J(os.path.join(R, "ours.TJ.k64.json"))
        g = J(os.path.join(RV, "dgg.TJpre_500000.0.0001.json"))
    else:
        mj = J(os.path.join(RV, "meshes", f"TJadv_{K}.json")) or {}
        if K == "comb":  # TJ with faces removed: the remaining faces are TJ's
            mj = dict(BASE_STATS, **mj)
        kd, k64 = J(os.path.join(RV, f"adv.{K}.kd.json")), J(os.path.join(RV, f"adv.{K}.k64.json"))
        g = J(os.path.join(RV, f"dgg.TJadv_{K}.0.0001.json"))
    nv = mj.get("vertices", mj.get("vertices_out"))
    nu = unreached(f"TJadv_{K}") if K else 0
    sl = f"{mj['slope_p50_deg']:.0f}/{mj['slope_p99_deg']:.0f}" if "slope_p50_deg" in mj else "--"
    asp = f"{mj['aspect_p99']:.1f}" if "aspect_p99" in mj else "--"
    c = [name, what, f"{nv / 1e3:.0f}k" if nv else "--", sl, asp]
    for d, cap in ((kd, False), (k64, True)):
        if not d:
            c += ["running"] + [""] * 3
            continue
        v8 = var(d, "code8/sym")
        v = v8 if cap else var(d, "dense/sym")  # errors of the stored variant: float32 uncapped, 8-bit capped
        n = v["pairs"]  # a pair the reference did not reach has error 1 here: left out of the mean
        mean = (v["est_mean"] * n - nu) / (n - nu)
        mx = "--" if nu else pct(v["est_max"], 2)
        size = gb(v8["index_bytes"]) if cap else f"{gb(v['index_bytes'])}/{gb(v8['index_bytes'])}"
        c += [h(shortest("TJ", "k64" if cap else "kd") if K is None else d["build_s"]), size, pct(mean, 4), mx]
    c += [f"{g['query_ms_p2p']:.0f}\\,ms", pct((g["err_mean"] * g["queries"] - nu) / (g["queries"] - nu), 4)] if g else ["--", "--"]
    rows.append(" & ".join(c) + "\\\\")
with open(os.path.join(OUT, "tab_adversarial.tex"), "w") as f:
    f.write("\n".join(rows) + "\n\\bottomrule\n")  # here: after an \\input, a \\bottomrule in the tabular is misplaced
print("tab_adversarial.tex:", len(rows), "rows")

# ---- landmark intervals against ours: width by distance, share certified without a refine, worst UB
COLS = ("near", "d2", "d4", "random")


def frac(c, n):
    return "100\\%" if c == n else (f"{100 * c / n:.2f}\\%" if c > 0.999 * n else f"{100 * c / n:.1f}\\%")


lines = []
for T in ("TJ", "GI"):
    d, xs = J(os.path.join(RV, f"landmark.{T}.json")), J(os.path.join(RV, f"landmark_sets.{T}.json"))
    if not d or not xs:
        continue
    for x in d["landmarks"]:
        sx = next((y for y in xs["landmarks"] if y["L"] == x["L"]), None)
        if not sx:
            continue
        st = {y["set"]: y for y in sx["sets"]}
        st["random"] = {"width_median": x["width_median"], "certified": x["certified"], "pairs": d["pairs"], "ub_err_max": x["ub_err_max"]}
        cert, tot = sum(y["certified"] for y in st.values()), sum(y["pairs"] for y in st.values())
        wub = max(y["ub_err_max"] for y in st.values())
        lines.append(" & ".join([T, f"landmarks, $L={x['L']}$", f"{x['index_bytes'] / 1e6:.0f}\\,MB", f"{x['build_s'] / 3600:.2f}\\,h", f"{x['query_us']:.2f}\\,$\\mu$s"]
                                + [pct(st[c]["width_median"], 1 if st[c]["width_median"] >= 0.1 else 3) for c in COLS]
                                + [frac(cert, tot), f"{100 * wub:,.0f}\\%"]) + "\\\\")
    for tag, label in (("k64", "ours, cap 64, 8-bit"), ("kd", "ours, uncapped")):
        o = J(os.path.join(R, f"ours.{T}.{tag}.json"))
        if not o:
            continue
        vv = "code8/sym" if tag == "k64" else "dense/sym"
        vs = {y["set"]: y for y in o["variants"] if y["variant"] == vv and y["set"] != "a2a"}
        cert, tot = sum(y["pairs"] - y["refined"] for y in vs.values()), sum(y["pairs"] for y in vs.values())
        wub = max(y["max_ub_err"] for y in vs.values())
        lines.append(" & ".join([T, label, f"{vs['random']['index_bytes'] / 1e6:.0f}\\,MB", h(shortest(T, tag)), f"{vs['random']['query_us']:.2f}\\,$\\mu$s"]
                                + [("exact" if vs[c]["width_median"] < 1e-6 else pct(vs[c]["width_median"], 1 if vs[c]["width_median"] >= 0.001 else 3)) for c in COLS]
                                + [frac(cert, tot), f"{100 * wub:.1f}\\%"]) + "\\\\")
    lines.append("\\midrule")
if lines and lines[-1] == "\\midrule":
    lines.pop()
with open(os.path.join(OUT, "tab_landmarks.tex"), "w") as f:
    f.write("\n".join(lines) + "\n\\bottomrule\n")
print("tab_landmarks.tex:", len(lines), "lines")
