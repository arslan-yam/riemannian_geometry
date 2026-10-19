"""One-off repair (19 Oct): dgg_one.sh's wall-time regex never matched GNU time's line, so every dgg*.json has
precompute_s = null. Refill it from the precompute.time files kept in the work dirs.
usage: python3 scripts/server/fix_dgg_walltime.py [results/server]"""
import glob, json, os, re, sys
O = sys.argv[1] if len(sys.argv) > 1 else "results/server"
for js in sorted(glob.glob(f"{O}/dgg.*.json") + glob.glob(f"{O}/obst.dgg.*.json")):
    name = os.path.basename(js)[:-5]                       # dgg.<mesh>.<eps> or obst.dgg.<mesh>.<eps>
    rest = name.split("dgg.", 1)[1]
    tf = f"{O}/dgg/{rest}/precompute.time"
    if not os.path.exists(tf): print(f"{name}: no {tf}"); continue
    m = re.search(r"Elapsed \(wall clock\) time \([^)]*\): ([0-9:.]+)", open(tf).read())
    if not m: print(f"{name}: no wall time in {tf}"); continue
    wall = 0.0
    for part in m.group(1).split(":"): wall = wall * 60 + float(part)
    d = json.load(open(js)); old = d.get("precompute_s"); d["precompute_s"] = wall
    json.dump(d, open(js, "w"))
    print(f"{name}: precompute {wall:.0f} s (was {old})")
