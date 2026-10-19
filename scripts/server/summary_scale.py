"""Accuracy and query time of our scale runs (or any cci_exp JSON): one line per run and storage variant, for
the random, near and arbitrary-point (a2a) query sets, plus a soundness count over every set.
usage: python3 scripts/server/summary_scale.py [glob, default results/server/scale.*.json] [--variants code8/sym,dense/sym]"""
import glob, json, sys
args = [a for a in sys.argv[1:] if not a.startswith("--")]
vs = "code8/sym,dense/sym"
if "--variants" in sys.argv: vs = sys.argv[sys.argv.index("--variants") + 1]; args = [a for a in args if a != vs]
want = vs.split(",")
files = sorted(glob.glob(args[0] if args else "results/server/scale.*.json"))
print(f"{'run':26} {'variant':10} {'size MB':>8} | {'random: us':>10} {'mean%':>8} {'max%':>7} | {'near: us':>8} {'mean%':>8} {'max%':>7} | {'a2a: us':>8} {'mean%':>8} {'max%':>7} | sound")
for f in files:
    try: d = json.load(open(f))
    except Exception: print(f.split('/')[-1], "incomplete"); continue
    rows = d.get("variants", [])
    bad = sum(1 for r in rows if not r.get("sound", True)); tot = len(rows)
    name = f.split("/")[-1][:-5].replace("scale.", "")
    for v in want:
        R = {r["set"]: r for r in rows if r.get("variant") == v}
        if not R: continue
        def c(s):
            r = R.get(s)
            if not r: return f"{'-':>10} {'-':>8} {'-':>7}" if s == "random" else f"{'-':>8} {'-':>8} {'-':>7}"
            us = r.get("query_us_refine_mean", r["query_us"]) if s == "a2a" else r["query_us"]
            w = 10 if s == "random" else 8
            return f"{us:>{w}.3f} {100 * r['est_mean']:>8.4f} {100 * r['est_max']:>7.3f}"
        size = next(iter(R.values())).get("index_bytes", 0) / 1e6
        print(f"{name:26} {v:10} {size:>8.0f} | {c('random')} | {c('near')} | {c('a2a')} | {tot - bad}/{tot}{'' if d.get('complete') else ' (running)'}")
