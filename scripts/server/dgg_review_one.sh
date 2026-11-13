#!/usr/bin/env bash
# DGG for the review experiments (items 1 and 2), one job per mode: precompute the discrete geodesic graph of one
# terrain (authors' code, one thread, a few minutes), use it, delete it. Modes:
#   bench  every exported V2V pair (2,000), one full Dijkstra (their query) and one stopped at the target
#          (cpp/tools/dgg_bench: query_ms, query_ms_p2p)          -> $O/review/dgg.<mesh>.<eps>.json
#   ch     a contraction hierarchy on the graph (cpp/tools/dgg_ch; core-based when contraction stalls, time-limited)
#                                                                    -> $O/review/dggch.<mesh>.<eps>.json
#   pll    hub labels (pruned landmark labelling, time-limited)     -> $O/review/dggpll.<mesh>.<eps>.json
#   knn    kNN / range over the 10,000 objects of the proximity runs, early-stopping Dijkstra (cpp/tools/dgg_knn)
#          on $O/comp/prox.<T>.pairs (the proxexp jobs)            -> $O/review/dggknn.<mesh>.<eps>.json
# usage: scripts/server/dgg_review_one.sh <mesh, e.g. TJpre_500000> <eps, e.g. 0.0001> <mode>
# env: DGG_BIN (default cpp/bin2), CH_HOURS (40), CH_CORE_FACTOR (3: the core starts when the mean degree triples), PLL_HOURS (72)
set -eu
cd "$(dirname "$0")/../.."
M=$1; E=$2; MODE=$3; O=${RG_OUT:-results/server}; ROOT=$PWD; DB=${DGG_BIN:-$ROOT/cpp/bin2}
T=${M%%pre_*}; T=${T%%adv_*}
mkdir -p $O/review; OA=$(cd $O && pwd)
W=$O/review/work/$M.$E.$MODE; mkdir -p $W
cp $O/meshes/$M.off $W/
cd $W
${GTIME:-/usr/bin/time} -v $DB/dgg_precompute f $M.off $E 1 > precompute.log 2> precompute.time
G=$(ls ${M}_FD*_c5.binary | head -1)
case $MODE in
  bench) $DB/dgg_bench $G $OA/comp/$M.V2V.pairs $OA/comp/$M.V2V.truth --json res.json | tee res.log; OUT=$OA/review/dgg.$M.$E.json ;;
  ch)    ${GTIME:-/usr/bin/time} -v $DB/dgg_ch $G $OA/comp/$M.V2V.pairs $OA/comp/$M.V2V.truth --settle 1000 --core-factor ${CH_CORE_FACTOR:-3} --time-limit ${CH_HOURS:-40} --json res.json > res.log 2> res.time; OUT=$OA/review/dggch.$M.$E.json ;;
  pll)   ${GTIME:-/usr/bin/time} -v $DB/dgg_ch $G $OA/comp/$M.V2V.pairs $OA/comp/$M.V2V.truth --ch 0 --pll degree --pll-time-limit ${PLL_HOURS:-72} --json res.json > res.log 2> res.time; OUT=$OA/review/dggpll.$M.$E.json ;;
  knn)   $DB/dgg_knn $G $M.off $OA/comp/prox.$T.pairs --k 1,10,100 --radius 0.01,0.03,0.1 --json res.json | tee res.log; OUT=$OA/review/dggknn.$M.$E.json ;;
  *) echo "unknown mode $MODE"; exit 1 ;;
esac
# one JSON: the result + the precompute's wall time and peak memory + graph size (+ the run's own GNU time if any)
python3 - "$OUT" "$G" <<'PY'
import json, os, re, sys
d = json.load(open("res.json"))
def gt(path):
    if not os.path.exists(path): return None, None
    t = open(path).read()
    m = re.search(r"Elapsed \(wall clock\) time \([^)]*\): ([0-9:.]+)", t)
    wall = None
    if m:
        wall = 0.0
        for part in m.group(1).split(":"): wall = wall * 60 + float(part)
    r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
    return wall, int(r.group(1)) if r else None
d["precompute_s"], d["precompute_rss_kb"] = gt("precompute.time")
d["run_wall_s"], d["run_rss_kb"] = gt("res.time")
d["graph_bytes"] = os.path.getsize(sys.argv[2])
json.dump(d, open(sys.argv[1], "w"), indent=1)
PY
cd $ROOT; rm -rf $W
