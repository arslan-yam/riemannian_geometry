#!/usr/bin/env bash
# SEP baseline (our reimplementation of Gotsman and Hormann 2022, cpp/tools/sep_bench.cpp) on one prepared
# terrain: build the index single-threaded, then answer all exported V2V pairs against the exact truth.
# usage: scripts/server/sep_one.sh <mesh name, e.g. TJpre_500000> [extra sep_bench flags, e.g. --leaf-exact]
# env: SEP_TAG (result $O/<tag>.<mesh>.json, default sep), SEP_BIN (default cpp/bin2; make -C cpp BIN=bin2 bin2/sep_bench)
set -eu
cd "$(dirname "$0")/../.."
M=$1; shift; O=${RG_OUT:-results/server}; TAG=${SEP_TAG:-sep}
W=$O/$TAG/$M; mkdir -p $W
${GTIME:-/usr/bin/time} -v ${SEP_BIN:-cpp/bin2}/sep_bench $O/meshes/$M.off $O/comp/$M.V2V.pairs $O/comp/$M.V2V.truth \
  --json $W/bench.json "$@" > $W/bench.log 2> $W/bench.time
cat $W/bench.log
python3 - "$W" "$O/$TAG.$M.json" <<'PY'
import json, re, sys
w, out = sys.argv[1:]
d = json.load(open(w + "/bench.json"))
r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(w + "/bench.time").read())
d.update(rss_kb=int(r.group(1)) if r else None)
json.dump(d, open(out, "w"))
PY
