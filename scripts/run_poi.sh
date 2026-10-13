#!/bin/zsh
# POI mode (index only the POIs) against the plain exact POI table (one full VTP sweep per POI, P x P
# float32), Gunnison10k, eps 0.2 / g 1 edge, 1,000 random POI pairs, single thread.
# run detached: nohup scripts/run_poi.sh > results/logs/poi.log 2>&1 < /dev/null &
set -u
cd "$(dirname "$0")/.."
X=cpp/bin/cci_exp; O=results/poi; mkdir -p $O
for M in HorseMount RobinsonMountain Gunnison10k; do
  for P in 50 200 1000 3000; do
    $X data/terrains/$M.off --eps 0.2 --gap-factor 1 --poi-random $P --random 1000 --poi-table --variants dense,code8 --refine --json $O/$M.P$P.json | grep -v degenerate
  done
done
echo "=== all done $(date)"
