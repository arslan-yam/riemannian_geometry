#!/bin/zsh
# size experiments (notes/research.md §5): one build per configuration, every storage evaluated on it.
#   1. Gunnison10k: kd (dense / 16-bit / covers), k_max caps with certify-or-refine, shifted hierarchy
#   2. size-versus-error curve: 5 (eps, g) settings on HorseMount, RobinsonMountain, Gunnison10k
#   3. GunnisonForest100k: kd (dense / 16-bit / covers), k_max 50 / 100 / 200, shifted
# run detached: nohup scripts/run_size.sh > results/logs/size.log 2>&1 < /dev/null &
set -u
cd "$(dirname "$0")/.."
X=cpp/bin/cci_exp; O=results/size
V="dense,code16,code8,cover0.002,cover0.01,cover0.05"
G="--pairs results/pairs/Gunnison10k.seed7.txt --truth-cache results/cpp/Gunnison10k.seed7.truth"
echo "=== start $(date)"
$X data/terrains/Gunnison10k.off --eps 0.1 --gap-factor 0.5 $=G --variants $V --refine --json $O/g10k.kd.json | grep -v degenerate
for K in 50 100 200; do
  $X data/terrains/Gunnison10k.off --eps 0.1 --gap-factor 0.5 --kmax $K $=G --variants dense,code16,code8 --refine --json $O/g10k.kmax$K.json | grep -v degenerate
done
$X data/terrains/Gunnison10k.off --eps 0.1 --gap-factor 0.5 --hierarchy shifted $=G --variants dense,code16,code8 --refine --json $O/g10k.shifted.json | grep -v degenerate
echo "=== Gunnison10k done $(date)"
for M in HorseMount RobinsonMountain Gunnison10k; do
  P="--pairs results/pairs/$M.seed7.txt --truth-cache results/cpp/$M.seed7.truth"
  for C in "0.1 0.5" "0.2 1" "0.3 2" "0.5 4" "1 8"; do E=${C% *}; GF=${C#* }
    $X data/terrains/$M.off --eps $E --gap-factor $GF $=P --variants dense,code16,code8 --json $O/curve.$M.$E.$GF.json | grep -v degenerate
  done
done
echo "=== curve done $(date)"
R="--random 200 --seed 7 --truth-cache results/cpp/GunnisonForest100k.r200.truth"
L=data/terrains/large/GunnisonForest100k.off
$X $L --eps 0.2 --gap-factor 1 $=R --variants $V --refine --json $O/g100k.kd.json | grep -v degenerate
for K in 100 50 200; do
  $X $L --eps 0.2 --gap-factor 1 --kmax $K $=R --variants dense,code16,code8 --refine --json $O/g100k.kmax$K.json | grep -v degenerate
done
$X $L --eps 0.2 --gap-factor 1 --hierarchy shifted $=R --variants dense,code16,code8 --refine --json $O/g100k.shifted.json | grep -v degenerate
echo "=== all done $(date)"
