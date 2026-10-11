#!/bin/zsh
# second batch, after scripts/run_cpp.sh: DGG (macOS port of the authors' code) and NeuroGF
# (authors' network, distance branch) on the same 1,500 pairs per mesh; then the Python benchmark's
# certified A2A rows. Waits for run_cpp.sh to finish so that nothing runs concurrently.
set -u
cd "$(dirname "$0")/.."
until grep -q "=== all done" results/logs/cpp_runs.log; do sleep 30; done
export MPLBACKEND=Agg
echo "=== start $(date)"
mkdir -p results/dgg
for M in HorseMount BigMount HeadlightMount RobinsonMountain Gunnison10k; do
  cp data/terrains/$M.off results/dgg/$M.off
  for E in 0.001 0.0001; do
    ( cd results/dgg && T0=$(perl -MTime::HiRes=time -e 'print time')
      ../../cpp/bin/dgg_precompute f $M.off $E 1 > /dev/null 2>&1
      T1=$(perl -MTime::HiRes=time -e 'print time')
      G=$(ls ${M}_FD${E}*_c5.binary 2>/dev/null | head -1)
      [ -z "$G" ] && G=$(ls ${M}_FD*_c5.binary | grep -F "$(printf '%.10f' $E)" | head -1)
      echo "DGG precompute $M eps=$E: $(perl -e "printf '%.2f', $T1-$T0") s, graph $(ls -la $G | awk '{print $5}') bytes"
      ../../cpp/bin/dgg_bench $G ../pairs/$M.seed7.txt ../cpp/$M.seed7.truth 2>/dev/null | grep "^DGG" )
  done
done
echo "=== DGG done $(date)"
for M in HorseMount RobinsonMountain; do
  python3 scripts/competitors/neurogf_terrain.py $M --sources 512 --epochs 30 --iters 100 2>&1 | grep -v Warn | tail -2
done
echo "=== NeuroGF done $(date)"
for M in HorseMount BigMount HeadlightMount RobinsonMountain Gunnison10k; do
  python3 -u bench/terrain_bench.py data/terrains/$M.off --queries 100 --a2a 100 --methods none \
      --cert 0.1/0.5 --out results 2>&1 | grep --line-buffered -E "A2A|certificate"
done
echo "=== all done $(date)"
