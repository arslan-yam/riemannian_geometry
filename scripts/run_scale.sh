#!/bin/zsh
# the 10^5-vertex experiment: our index (single thread) then the competitors, all on GunnisonForest100k.
# run detached (macOS has no setsid):  nohup scripts/run_scale.sh > results/logs/scale.log 2>&1 < /dev/null &
# needs EAR_DIR (clone of ItakEjgo/weighted_distance_oracle, built) for the competitor part.
set -u
cd "$(dirname "$0")/.."
MESH=GunnisonForest100k
export MPLBACKEND=Agg
echo "=== start $(date)"
python3 -u bench/terrain_bench.py data/terrains/large/$MESH.off --queries 100 --a2a 100 --k 32 --k-list 32 --workers 1 \
    --out results --methods dijkstra,exact,exact_scaled 2>&1 | grep --line-buffered -v "Warning\|cand = \|segment_lenght\|np.ceil"
echo "=== ours done $(date)"
if [ -n "${EAR_DIR:-}" ]; then
  cp data/terrains/large/$MESH.off $EAR_DIR/datasets/repaired/$MESH.off
  head -10 results/$MESH.V2V.query > results/$MESH.V2V.head10.query
  head -10 results/$MESH.A2A.query > results/$MESH.A2A.head10.query
  KIND=V2V scripts/competitors/run_competitors.sh $MESH results/$MESH.V2V.query MMP FixedS EAR
  KIND=A2A scripts/competitors/run_competitors.sh $MESH results/$MESH.A2A.query MMP FixedS EAR
  KIND=V2V scripts/competitors/run_competitors.sh $MESH results/$MESH.V2V.head10.query KAlgo
  KIND=A2A scripts/competitors/run_competitors.sh $MESH results/$MESH.A2A.head10.query KAlgo
  python3 scripts/competitors/competitor_eval.py $EAR_DIR/results $MESH V2V
  python3 scripts/competitors/competitor_eval.py $EAR_DIR/results $MESH A2A
fi
echo "=== all done $(date)"
