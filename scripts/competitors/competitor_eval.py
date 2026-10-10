#evaluate the logs of the EAR-Oracle harness (https://github.com/ItakEjgo/weighted_distance_oracle)
#against the same exact truth terrain_bench.py uses, into results/<mesh>_competitors[_a2a].json
#
#usage: python3 scripts/competitors/competitor_eval.py <harness results dir> <mesh> [V2V|A2A]
#log names: <mesh>_<KIND>_<METHOD>_eps<eps>_sp<sp>.log (as scripts/competitors/run_competitors.sh writes them)
import numpy as np, sys, glob, os, re, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import pygeodesic.geodesic as pg
from scipy.spatial import cKDTree
from rgeo.terrain import read_off, TerrainMesh, ExactTerrain

logdir, mesh = sys.argv[1], sys.argv[2]
kind = sys.argv[3] if len(sys.argv) > 3 else "V2V"
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mesh_path = os.path.join(root, "data", "terrains", f"{mesh}.off")
if not os.path.exists(mesh_path):
    mesh_path = os.path.join(root, "data", "terrains", "large", f"{mesh}.off")
V, F = read_off(mesh_path)
geo = pg.PyGeodesicAlgorithmExact(V, F)
tag = mesh + ("_weighted" if kind == "W2V" else "") #weighted runs are stored under <mesh>_weighted
q = [l.split() for l in open(os.path.join(root, "results", f"{tag}.{kind}.query"))]
tr = cKDTree(V)
truth = []
if kind == "W2V":
    #weighted terrain: no exact solver. the reference is the best path any method found for
    #that query, ours (per-query values stored in <mesh>_weighted.json) and theirs alike
    ours = json.load(open(os.path.join(root, "results", f"{tag}.json")))
    n_q = len(q)
    best = np.full(n_q, np.inf)
    #only answers that are the cost of an actual path define the reference: interpolated answers
    #(the certified index's estimate, the earlier index's + interp / + model interp) can lie below
    #the truth, so they are scored against the reference but never lower it
    is_path = lambda m: not ("(estimate)" in m or "interp" in m)
    used = 0
    for r in ours["rows"]:
        v = r.get("values")
        if v and is_path(r["method"]):
            best[:len(v)] = np.minimum(best[:len(v)], np.array(v, dtype=float))
            used += 1
    theirs = {}
    for log in sorted(glob.glob(os.path.join(logdir, f"{mesh}_{kind}_*.log"))):
        txt = open(log).read()
        if "Query results end" not in txt or "_KAlgo_" in log:
            #their K-Algo ignores the face weights (its weighted answers match the UNWEIGHTED
            #exact distances to 2.7%), so it is neither a reference nor a competitor here
            continue
        L = txt.split("Query results begin:")[1].split("Query results end.")[0].split()
        a = np.array(L, dtype=float).reshape(-1, 2)
        theirs[log] = a
        best[:len(a)] = np.minimum(best[:len(a)], a[:, 0])
    truth = best
    #re-score our rows against the same reference and write them back
    for r in ours["rows"]:
        v = r.get("values")
        if v:
            rel = np.array(v, dtype=float) / best[:len(v)] - 1
            r.update(err_mean=float(rel.mean() * 100), err_p90=float(np.percentile(rel, 90) * 100),
                     err_max=float(rel.max() * 100), err_min=float(rel.min() * 100))
    ours["reference"] = ("best path found by any method (no exact solver for weighted terrains); "
                         "interpolated answers are scored but excluded from the reference")
    json.dump(ours, open(os.path.join(root, "results", f"{tag}.json"), "w"), indent=1)
    print("reference = pointwise best over", used, "path-valued rows of ours and", len(theirs), "of their logs")
elif kind == "V2V":
    for l in q:
        s = tr.query([float(l[0]), float(l[1]), float(l[2])])[1]; t = tr.query([float(l[4]), float(l[5]), float(l[6])])[1]
        truth.append(geo.geodesicDistance(s, t)[0])
else: #arbitrary points: both inserted as steiner vertices (same truth as terrain_bench.exact_a2a)
    tm = TerrainMesh(V, F); ex = ExactTerrain(tm)
    for l in q:
        p = np.array([float(l[0]), float(l[1])]) - tm.origin; qq = np.array([float(l[4]), float(l[5])]) - tm.origin
        V2, F2, idx = ex.insert_points(np.array([p, qq])); g = pg.PyGeodesicAlgorithmExact(V2, F2)
        truth.append(g.geodesicDistance(int(idx[0]), int(idx[1]))[0])
truth = np.array(truth)
rows = []
for log in sorted(glob.glob(os.path.join(logdir, f"{mesh}_{kind}_*.log"))):
    txt = open(log).read()
    if "Query results end" not in txt or (kind == "W2V" and "_KAlgo_" in log):
        continue
    m = re.search(r"_" + kind + r"_(\w+)_eps([\d.]+)_sp(\d+)\.log", log)
    method, eps, spn = m.group(1), m.group(2), m.group(3)
    L = txt.split("Query results begin:")[1].split("Query results end.")[0].split()
    a = np.array(L, dtype=float).reshape(-1, 2)
    if len(a) > len(truth):
        print("length mismatch", log); continue
    rel = a[:, 0] / truth[:len(a)] - 1 #a run on a prefix of the query file (slow on-the-fly methods) is evaluated on that prefix
    it = re.search(r"\[Index Time\].*?: ([\d.]+) ms", txt); mem = re.search(r"Index memory usage: (\d+) MB", txt)
    row = dict(method=method, eps=float(eps), sp=int(spn), build_s=float(it.group(1)) / 1e3 if it else 0.0,
               index_MB=float(mem.group(1)) if mem else 0.0, query_us=float(a[:, 1].mean()), queries=int(len(a)),
               err_mean=float(rel.mean() * 100), err_p90=float(np.percentile(rel, 90) * 100),
               err_max=float(rel.max() * 100), err_min=float(rel.min() * 100))
    rows.append(row)
    print(f"{method:9s} eps={eps:5s} sp={spn:3s} q={len(a):3d} build {row['build_s']:8.1f}s  index {row['index_MB']:6.0f} MB  "
          f"query {row['query_us']:10.0f} us  err mean {row['err_mean']:.3f}%  p90 {row['err_p90']:.3f}%  max {row['err_max']:.3f}%  min {row['err_min']:.3f}%")
out = os.path.join(root, "results", f"{mesh}_competitors{ {'V2V': '', 'A2A': '_a2a', 'W2V': '_w2v'}[kind] }.json")
if os.path.exists(out): #merge with rows evaluated earlier (e.g. a method run in another session)
    old = {(r["method"], r["eps"], r.get("sp", 5)): r for r in json.load(open(out))
           if not (kind == "W2V" and r["method"] == "KAlgo")}
    for r in rows:
        old[(r["method"], r["eps"], r["sp"])] = r
    rows = list(old.values())
json.dump(rows, open(out, "w"), indent=1)
