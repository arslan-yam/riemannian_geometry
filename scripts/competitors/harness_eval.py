"""Score logs of the competitors' harness (EAR authors' code, scripts/competitors/run_competitors.sh) against
the exact distances that cpp/tools/cci_exp.cpp --export wrote next to the query file.

usage: python3 scripts/competitors/harness_eval.py <truth file> <log>... [--json out.json]

V2V queries were moved 5 cm into a face (the harness cannot locate a point exactly at a vertex), so each
harness answer may differ from the vertex-to-vertex truth by up to 0.1 m (2 x 5 cm); the "below truth"
count uses that tolerance. A2A queries are the exact points of the truth.
"""
import json, os, re, sys
import numpy as np


def parse(log):
    txt = open(log).read()
    out = {"log": os.path.basename(log), "complete": "Query results end" in txt}
    m = re.search(r"\[Index Time\][^:]*:\s*([0-9.]+)\s*ms", txt)
    out["build_s"] = float(m.group(1)) / 1e3 if m else None
    m = re.search(r"Index memory usage:\s*([0-9.]+)\s*MB", txt)
    out["index_MB"] = float(m.group(1)) if m else None
    if out["complete"]:
        a = np.array(txt.split("Query results begin:")[1].split("Query results end.")[0].split(), dtype=float).reshape(-1, 2)
        out["answers"], out["query_us"] = a[:, 0], a[:, 1]
    return out


def main():
    args = sys.argv[1:]
    js = None
    if "--json" in args:
        i = args.index("--json"); js = args[i + 1]; del args[i:i + 2]
    truth = np.loadtxt(args[0])
    rows = []
    for log in args[1:]:
        r = parse(log)
        if r["complete"]:
            a = r.pop("answers")[:len(truth)]
            t = truth[:len(a)]
            rel = a / t - 1
            r.update(pairs=int(len(a)), err_mean=float(np.abs(rel).mean() * 100), err_p99=float(np.percentile(np.abs(rel), 99) * 100),
                     err_max=float(np.abs(rel).max() * 100), err_min_signed=float(rel.min() * 100),
                     below_truth=int((a < t - 0.1).sum()), query_us_mean=float(r.pop("query_us").mean()))
        rows.append(r)
        print(f"{r['log']}: " + (f"build {r['build_s']} s, index {r['index_MB']} MB, query {r['query_us_mean']:.1f} us, "
                                  f"|err| mean {r['err_mean']:.3f}% p99 {r['err_p99']:.3f}% max {r['err_max']:.3f}%, "
                                  f"{r['below_truth']} answers below the truth" if r["complete"] else "no query results (did not finish)"))
    if js:
        json.dump(rows, open(js, "w"), indent=1)


if __name__ == "__main__":
    main()
