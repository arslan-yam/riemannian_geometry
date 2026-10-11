#!/bin/zsh
# NeuroGF re-run (single forward pass per iteration, as in the authors' script), after the other batches
set -u
cd "$(dirname "$0")/.."
until grep -q "=== all done" results/logs/after.log 2>/dev/null; do sleep 30; done
export MPLBACKEND=Agg
for M in HorseMount RobinsonMountain; do
  python3 scripts/competitors/neurogf_terrain.py $M --sources 512 --epochs 30 --iters 100 2>&1 | grep -v Warn | tail -2
done
echo "=== all done $(date)"
