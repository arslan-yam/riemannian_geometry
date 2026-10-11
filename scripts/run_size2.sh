#!/bin/zsh
# size-versus-error curve at a fixed near radius R = g / (eps lambda) = 5 edges: coarser gaps with eps
# growing in step (the guarantee weakens to (1+eps), the certificate stays per query)
set -u
cd "$(dirname "$0")/.."
until grep -q "=== all done" results/logs/size.log; do sleep 60; done
X=cpp/bin/cci_exp; O=results/size
for M in HorseMount RobinsonMountain Gunnison10k; do
  P="--pairs results/pairs/$M.seed7.txt --truth-cache results/cpp/$M.seed7.truth"
  for C in "0.4 2" "0.8 4" "1.6 8" "3.2 16"; do E=${C% *}; GF=${C#* }
    $X data/terrains/$M.off --eps $E --gap-factor $GF $=P --variants dense,code16,code8 --json $O/curve.$M.$E.$GF.json | grep -v degenerate
  done
done
echo "=== all done $(date)"
