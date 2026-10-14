"""Progress of a server run: finished jobs (logs/jobs.jsonl), failed jobs, jobs still running (a log
without a finished line), and result JSON files that are invalid or not complete.
usage: python3 scripts/server/status.py [results/server] [--only <job prefix>]
"""
import glob, json, os, sys

args = sys.argv[1:]
only = ""
if "--only" in args:
    i = args.index("--only"); only = args[i + 1]; del args[i:i + 2]
O = args[0] if args else "results/server"
L = os.path.join(O, "logs")

done = {}
if os.path.exists(os.path.join(L, "jobs.jsonl")):
    for line in open(os.path.join(L, "jobs.jsonl")):
        try:
            j = json.loads(line)
        except ValueError:
            continue
        if j["job"].startswith(only): done[j["job"]] = j  # a rerun replaces the earlier line
running = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(L, only + "*.log"))
                 if os.path.basename(p)[:-4] not in done)

print(f"{'job':44s} {'exit':>4s} {'wall':>9s} {'peak RSS':>9s}")
for name, j in sorted(done.items(), key=lambda kv: kv[1].get("end", "")):
    w = j.get("wall_s"); r = j.get("max_rss_kb")
    print(f"{name:44s} {j['exit']:4d} {(f'{w / 3600:.2f} h' if w is not None else '?'):>9s} "
          f"{(f'{r / 2**20:.1f} GB' if r is not None else '?'):>9s}")
failed = [n for n, j in done.items() if j["exit"] != 0]
print(f"\n{len(done)} finished, {len(failed)} failed, {len(running)} running")
for n in failed: print("  FAILED ", n, "->", os.path.join(L, n + ".log"), os.path.join(L, n + ".time"))
for n in running: print("  running", n)

bad = []
for p in sorted(glob.glob(os.path.join(O, "*.json")) + glob.glob(os.path.join(O, "smoke", "*.json"))):
    rel = os.path.relpath(p, O)
    if only and not (rel.startswith(only) or (only == "smoke" and rel.startswith("smoke/"))): continue
    try:
        d = json.load(open(p))
    except ValueError:
        bad.append((rel, "invalid JSON")); continue
    if isinstance(d, dict) and d.get("complete") is False: bad.append((rel, "complete: false"))
if bad:
    print("\nresult files not complete:")
    for rel, why in bad: print(f"  {rel}: {why}")
