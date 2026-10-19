#!/usr/bin/env bash
# DGG on the kNN / range workload of the proxpoi phase: precompute the discrete geodesic graph of one terrain
# (authors' code, one thread), then answer the exported kNN / range queries of each POI set with an early-stopping
# Dijkstra (cpp/tools/dgg_knn.cpp). Results: $O/proxpoi.dgg.<T>.<P>.<eps>.json; the graph is deleted afterwards.
# usage: scripts/server/dgg_knn_one.sh <terrain, e.g. TJ> <eps, e.g. 0.001> <POI counts, e.g. "50 500 2500">
set -eu
cd "$(dirname "$0")/../.."
T=$1; E=$2; PS=$3; O=${RG_OUT:-results/server}; ROOT=$PWD; DB=${DGG_BIN:-$ROOT/cpp/bin2}
M=${T}pre_500000; W=$O/proxdgg/$M.$E; mkdir -p $W; OA=$(cd $O && pwd)  # absolute (we cd into $W)
cp $O/meshes/$M.off $W/
cd $W
${GTIME:-/usr/bin/time} -v $DB/dgg_precompute f $M.off $E 1 > precompute.log 2> precompute.time
G=$(ls ${M}_FD*_c5.binary | head -1)
for P in $PS; do
  $DB/dgg_knn $G $M.off $OA/comp/proxpoi.$T.$P.pairs --k 1,10,25 --radius 0.01,0.03,0.1 --json knn.$P.json | tee knn.$P.log
  # one JSON per POI set: the queries + the precompute's wall time and peak memory + graph size
  python3 - "$OA/proxpoi.dgg.$T.$P.$E.json" "$G" knn.$P.json <<'PY'
import json, os, re, sys
d = json.load(open(sys.argv[3]))
t = open("precompute.time").read()
m = re.search(r"Elapsed \(wall clock\) time \([^)]*\): ([0-9:.]+)", t)
wall = None
if m:
    wall = 0.0
    for part in m.group(1).split(":"): wall = wall * 60 + float(part)
r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
d.update(precompute_s=wall, precompute_rss_kb=int(r.group(1)) if r else None, graph_bytes=os.path.getsize(sys.argv[2]))
json.dump(d, open(sys.argv[1], "w"), indent=1)
PY
done
rm -f "$G" $M.off
