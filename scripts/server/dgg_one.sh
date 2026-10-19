#!/usr/bin/env bash
# DGG on one prepared terrain: precompute the discrete geodesic graph (authors' code), then answer the
# first 200 exported V2V pairs (one Dijkstra each) against the exact truth.
# usage: scripts/server/dgg_one.sh <mesh name, e.g. TJpre_500000> <eps, e.g. 0.001>
# env: DGG_BIN (the binaries, default cpp/bin), DGG_TAG (work dir $O/<tag>/ and result $O/<tag>.<mesh>.<eps>.json,
#      default dgg), DGG_CLEAN=1 (delete the graph and the mesh copy afterwards; the graph size stays in the JSON)
set -eu
cd "$(dirname "$0")/../.."
M=$1; E=$2; O=${RG_OUT:-results/server}; ROOT=$PWD; DB=${DGG_BIN:-$ROOT/cpp/bin}; TAG=${DGG_TAG:-dgg}
W=$O/$TAG/$M.$E; mkdir -p $W
cp $O/meshes/$M.off $W/
head -200 $O/comp/$M.V2V.pairs > $W/pairs
head -200 $O/comp/$M.V2V.truth > $W/truth
cd $W
${GTIME:-/usr/bin/time} -v $DB/dgg_precompute f $M.off $E 1 > precompute.log 2> precompute.time
G=$(ls ${M}_FD*_c5.binary | head -1)
$DB/dgg_bench $G pairs truth --json bench.json > bench.log
cat bench.log
# one JSON: the queries (dgg_bench) + the precompute's wall time and peak memory (GNU time) + graph size
python3 - "$ROOT" "$O/$TAG.$M.$E.json" "$G" <<'PY'
import json, os, re, sys
d = json.load(open("bench.json"))
t = open("precompute.time").read()
m = re.search(r"Elapsed \(wall clock\) time \([^)]*\): ([0-9:.]+)", t)
wall = None
if m:
    wall = 0.0
    for part in m.group(1).split(":"): wall = wall * 60 + float(part)
r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
d.update(precompute_s=wall, precompute_rss_kb=int(r.group(1)) if r else None, graph_bytes=os.path.getsize(sys.argv[3]))
json.dump(d, open(os.path.join(sys.argv[1], sys.argv[2]), "w"))
PY
[ "${DGG_CLEAN:-}" = 1 ] && rm -f "$G" $M.off
exit 0
