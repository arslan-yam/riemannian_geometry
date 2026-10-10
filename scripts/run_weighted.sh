#!/bin/zsh
# weighted-terrain experiment on one mesh: ours (refined labels, single thread) then the competitors on
# the same face weights.  usage: scripts/run_weighted.sh HorseMount [seed]
set -u
cd "$(dirname "$0")/.."
MESH=$1; SEED=${2:-1}
export MPLBACKEND=Agg
echo "=== start $(date)"
python3 -u bench/terrain_bench.py data/terrains/$MESH.off --queries 100 --k 32 --workers 1 --rounds 1 --segments 12 \
    --weighted $SEED --out results --methods dijkstra,h2h,portals,cut,refine 2>&1 | grep --line-buffered -v "Warning\|cand = \|segment_lenght\|np.ceil"
echo "=== ours done $(date)"
if [ -n "${EAR_DIR:-}" ]; then
  W=${MESH}_weighted
  export FACE_WEIGHTS=$PWD/results/$W.face_weight.query
  head -30 results/$W.W2V.query > results/$W.W2V.head30.query
  KIND=W2V scripts/competitors/run_competitors.sh $MESH results/$W.W2V.query FixedS EAR KAlgo KAlgo@0.05
  KIND=W2V scripts/competitors/run_competitors.sh $MESH results/$W.W2V.head30.query UnfixedS@0.05
fi
echo "=== all done $(date)"
