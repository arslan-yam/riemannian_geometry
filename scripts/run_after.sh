#!/bin/zsh
# third batch: the Gunnison10k VTP-only run (the first attempt hit a VTP abort on an inconsistently
# oriented sub-mesh, fixed since), and the hybrid threshold (VTP when sub-mesh faces <= ratio x box faces)
set -u
cd "$(dirname "$0")/.."
until grep -q "=== all done" results/logs/baselines.log 2>/dev/null; do sleep 30; done
OUT=results/cpp; B=cpp/bin
P="--pairs results/pairs/Gunnison10k.seed7.txt --truth-cache $OUT/Gunnison10k.seed7.truth"
$B/cci_bench data/terrains/Gunnison10k.off --eps 0.1 --gap-factor 0.5 --solver vtp $=P --json $OUT/Gunnison10k.vtp.0.1.json | grep -v degenerate
for RATIO in 4 8; do for M in RobinsonMountain Gunnison10k; do
  $B/cci_bench data/terrains/$M.off --eps 0.1 --gap-factor 0.5 --solver hybrid --ratio $RATIO --pairs results/pairs/$M.seed7.txt \
      --truth-cache $OUT/$M.seed7.truth --json $OUT/ratio/$M.hybrid$RATIO.0.1.json | grep -v degenerate
done; done
echo "=== all done $(date)"
