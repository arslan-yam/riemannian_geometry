#!/bin/zsh
# final runs of the certified index (single thread, nothing else running), sequentially:
#   1. the five meshes of the comparison, eps/gap configs 0.2/1, 0.1/0.5, 0.05/0.25, V2V + A2A
#   2. weighted HorseMount (face costs uniform in (1, 5), seed 1) on Steiner labels m=5, 10,
#      with a Steiner m=24 Dijkstra reference; then the competitors' logs are re-scored (W2V)
#   3. the 100k-vertex GunnisonForest, config 0.2/1
# run detached:  nohup scripts/run_certified.sh > results/logs/certified_final.log 2>&1 < /dev/null &
set -u
cd "$(dirname "$0")/.."
export MPLBACKEND=Agg
F="Warning\|cand = \|segment_lenght\|np.ceil"
echo "=== start $(date)"
for M in HorseMount BigMount HeadlightMount RobinsonMountain Gunnison10k; do
  python3 -u bench/terrain_bench.py data/terrains/$M.off --queries 100 --a2a 100 --methods none \
      --cert 0.2/1,0.1/0.5,0.05/0.25 --out results 2>&1 | grep --line-buffered -v "$F"
  echo "=== $M done $(date)"
done
python3 -u bench/terrain_bench.py data/terrains/HorseMount.off --weighted 1 --queries 100 --methods dijkstra,h2h \
    --cert 0.1/0.5,0.2/1 --steiner-m 5,10 --steiner-ref 24 --out results 2>&1 | grep --line-buffered -v "$F"
if [ -n "${EAR_DIR:-}" ]; then
  python3 scripts/competitors/competitor_eval.py $EAR_DIR/results HorseMount W2V
fi
echo "=== weighted done $(date)"
python3 -u bench/terrain_bench.py data/terrains/large/GunnisonForest100k.off --queries 100 --a2a 100 --methods none \
    --cert 0.2/1 --out results 2>&1 | grep --line-buffered -v "$F"
echo "=== all done $(date)"
