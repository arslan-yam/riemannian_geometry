"""Provenance manifest of the cluster results (review item 6): which result files feed which table, figure or
paragraph of the paper, with a SHA-256 of every result file and of the query sets and meshes behind it, the command
line each run recorded (argv), and the code version. Run it wherever the files are: on the cluster it also hashes
the meshes and truth caches (results/server/meshes and truth are not copied back), on the laptop the results and
query exports.

  python3 scripts/server/manifest.py [results/server] [--out results/server/MANIFEST.json] [--md notes/manifest.md]

The snapshot is the code commit plus these hashes: regenerating a table means rerunning the recorded argv with the
binaries of that commit on the hashed meshes and query sets, then the generators listed per entry.
"""
import argparse, glob, hashlib, json, os, platform, subprocess, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default="results/server")
ap.add_argument("--out", default=None)
ap.add_argument("--md", default=None)
a = ap.parse_args()
R = a.root
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# paper item -> result-file patterns (relative to R) and the script that turns them into the paper's numbers
ENTRIES = [
    ("tab:server-ours", "our index on the six terrains at 0.5M faces", ["ours.*.kd.json", "ours.*.k64.json", "neurogf.??pre_500000.json"], "paper/preprint/server_summary.py"),
    ("tab:server-comp", "methods on TJ and GI at 0.5M faces", ["ours.TJ.*.json", "ours.GI.*.json", "ear.??pre_500000.V2V.json", "dgg.??pre_500000.*.json",
                                                              "neurogf.??pre_500000*.json", "sep_warm.*.json", "sep.*.json", "onthefly_fixed/onthefly.*_500000.json"], "paper/preprint/server_summary.py"),
    ("fig:server-methods", "query time and size against error, TJ and GI at 0.5M", ["ours.TJ.*.json", "ours.GI.*.json", "curve.*.json", "caps.*.json", "ear.??pre_500000.V2V.json",
                                                                                     "dgg.??pre_500000.*.json", "neurogf.??pre_500000*.json", "sep_warm.*.json", "onthefly_fixed/onthefly.*_500000.json"], "paper/preprint/make_server_figures.py"),
    ("tab:scale / fig:scaling", "scaling 0.5M to 5M faces", ["scale.*.json", "dgg.TJ*.json", "onthefly_fixed/*.json", "big.*.json", "neurogf.TJpre_1002528.json"], "paper/preprint/make_server_figures.py, server_summary.py"),
    ("break-even", "break-even query counts", ["scale.*.json", "ours.*.json", "dgg.*.json", "onthefly_fixed/*.json", "threads.*.json", "big.*.json"], "paper/preprint/server_summary.py"),
    ("ablations", "caps, eps_c covers, solvers, interpolation", ["caps.*.json", "ablation.*.json", "solvers.*.json", "curve.*.json"], "paper/preprint/server_summary.py"),
    ("threads", "parallel build", ["threads.*.json"], "paper/preprint/server_summary.py"),
    ("tab:obstacles", "obstacles at 0.5M faces", ["obst.*.json"], "paper/preprint/server_summary.py"),
    ("weighted", "weighted terrains", ["weighted.*.json", "ear.??pre_500000.W2V.json"], "paper/preprint/server_summary.py"),
    ("proximity", "certified kNN and range", ["prox.*.json", "proxpoi.*.json"], "paper/preprint/server_summary.py"),
    ("POI setting", "P x P table, POI index, UP-Oracle", ["poi.*.json", "upo.*.json"], "paper/preprint/server_summary.py"),
    ("small terrains (cluster reruns)", "the M1 tables rerun on the cluster", ["small.*.json", "ear.GunnisonForest100k.V2V.json"], "paper/preprint/server_summary.py"),
    ("review: DGG point-to-point, CH, hub labels, kNN", "item 1 and 2", ["review/dgg.*.json", "review/dggch.*.json", "review/dggknn.*.json"], "scripts/server/review_summary.py"),
    ("review: landmark intervals", "item 3", ["review/landmark.*.json"], "scripts/server/review_summary.py"),
    ("review: adversarial terrains", "item 5", ["review/adv.*.json", "review/meshes/*.json"], "scripts/server/review_summary.py"),
    ("review: build repeats and strict checks", "items 4 and 6", ["review/rep*.json", "review/small.*.json"], "scripts/server/review_summary.py"),
]


def sha(path, cache={}):
    if path in cache:
        return cache[path]
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    cache[path] = h.hexdigest()
    return cache[path]


def run_info(path):
    """what a result file says about how it was made"""
    try:
        d = json.load(open(path))
    except Exception:
        return {}
    if isinstance(d, list):
        d = d[0] if d and isinstance(d[0], dict) else {}
    keep = ("argv", "tool", "mesh", "graph", "n", "eps", "gap_factor", "k_max", "hierarchy", "solver", "threads", "store", "build_s",
            "complete", "log", "queries", "pairs")
    info = {k: d[k] for k in keep if k in d}
    if "variants" in d:
        vs = d["variants"]
        info["answers"] = sum(v.get("pairs", 0) for v in vs)
        info["all_sound"] = all(v.get("sound", True) for v in vs)
        if any("strict_violations" in v for v in vs):
            info["strict_violations"] = sum(v.get("strict_violations", 0) for v in vs)
            info["max_rel_violation"] = max(v.get("max_rel_violation", 0) for v in vs)
    return info


def inputs_of(info):
    """mesh and query files a run used, by the conventions of run_server.sh (hashed when present here)"""
    out = []
    argv = info.get("argv", "")
    toks = argv.split()
    for i, t in enumerate(toks):
        if t.endswith(".off") or t.endswith(".pairs") or t.endswith(".truth"):
            out.append(t)
        if t in ("--truth-cache",) and i + 1 < len(toks):
            out.append(toks[i + 1])
    return out


def git(*args):
    try:
        return subprocess.run(["git", "-C", REPO] + list(args), capture_output=True, text=True).stdout.strip()
    except Exception:
        return ""


files, entries = {}, []
for label, what, pats, gen in ENTRIES:
    got = sorted({os.path.relpath(p, R) for pat in pats for p in glob.glob(os.path.join(R, pat))})
    entries.append({"item": label, "what": what, "generator": gen, "files": got})
    for f in got:
        if f not in files:
            p = os.path.join(R, f)
            info = run_info(p)
            ins = {}
            for q in inputs_of(info):
                qp = q if os.path.isabs(q) else os.path.join(REPO, q)
                cands = [qp] + ([qp + s for s in ("", ".V2V.pairs")] if "truth" in q else [])
                for c in cands:
                    if os.path.isfile(c):
                        ins[q] = sha(c); break
                else:
                    ins[q] = None  # not on this machine (meshes and truth caches stay on the cluster)
            files[f] = {"sha256": sha(p), "bytes": os.path.getsize(p), "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(os.path.getmtime(p))),
                        "run": info, "inputs": ins}
# the exported query sets every method answers (comp/), hashed as a whole
queries = {os.path.relpath(p, R): sha(p) for p in sorted(glob.glob(os.path.join(R, "comp", "*.pairs")) + glob.glob(os.path.join(R, "comp", "*.truth"))
                                                               + glob.glob(os.path.join(R, "comp", "*.query")))}
meshes = {os.path.relpath(p, R): sha(p) for p in sorted(glob.glob(os.path.join(R, "meshes", "*.off")))}
man = {"created": time.strftime("%Y-%m-%dT%H:%M:%S"), "host": platform.node(), "machine": platform.platform(),
       "commit": git("rev-parse", "HEAD"), "dirty_files": [l for l in git("status", "--porcelain").splitlines() if l.strip()],
       "root": R, "entries": entries, "files": files, "query_sets": queries, "meshes": meshes}
out = a.out or os.path.join(R, "MANIFEST.json")
json.dump(man, open(out, "w"), indent=1)
missing = [e["item"] for e in entries if not e["files"]]
nf = len(files)
print(f"manifest {out}: {nf} result files, {len(queries)} query files, {len(meshes)} meshes, commit {man['commit'][:10]}"
      + (f" (+{len(man['dirty_files'])} uncommitted changes)" if man["dirty_files"] else ""))
unsound = [f for f, v in files.items() if v["run"].get("all_sound") is False]
if unsound:
    print("result files with an unsound answer:", ", ".join(unsound))
if missing:
    print("paper items without result files here:", ", ".join(missing))
if a.md:
    with open(a.md, "w") as f:
        f.write(f"# Result manifest\n\nCommit `{man['commit']}`, created {man['created']} on {man['host']}.\n\n")
        for e in entries:
            f.write(f"## {e['item']}\n\n{e['what']}; numbers by `{e['generator']}`.\n\n")
            for x in e["files"]:
                r = files[x]["run"]
                f.write(f"- `{x}` sha256 `{files[x]['sha256'][:12]}`" + (f", {r['answers']} answers" if "answers" in r else "")
                        + (f", strict violations {r['strict_violations']}" if "strict_violations" in r else "") + "\n")
            f.write("\n")
    print("markdown:", a.md)
