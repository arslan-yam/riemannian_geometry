#!/bin/zsh
# C++ runs, sequential and single-threaded on a quiet machine:
#   1. exact solvers (MMP / VTP full sweeps, MMP stopped at nearby targets) = on-the-fly exact baselines
#   2. the certified index in C++ on the five meshes: MMP, VTP and hybrid label sweeps (eps 0.1, g 0.5e),
#      hybrid also at 0.2/1 and 0.05/0.25; 1,500 random pairs (seed 7, same as the Python sweep)
#   3. the UP-Oracle authors' oracles (UP, WSPD, RC, EAR, WAV-Fly, ESP-Fly) on their TJ terrain and on
#      HorseMount, with 50 POIs and with every vertex a POI; our index on the same pairs
#   4. the competitors' weighted methods with the endpoint fix, re-scored
#   5. the 100k-vertex mesh in C++ (hybrid, 0.2/1)
# run detached: nohup scripts/run_cpp.sh > results/logs/cpp_runs.log 2>&1 < /dev/null &
set -u
cd "$(dirname "$0")/.."
OUT=results/cpp; mkdir -p $OUT results/upo
B=cpp/bin
lim() { perl -e 'alarm shift; exec @ARGV' "$@"; } # wall-clock limit (no `timeout` on macOS)
MESHES=(HorseMount BigMount HeadlightMount RobinsonMountain Gunnison10k)
UPI=external/UpdatedStructureTerrainCode/input
echo "=== start $(date)"
for M in $MESHES; do $B/solver_check data/terrains/$M.off 5; done
$B/solver_check $UPI/TJpre_1058.off 5
echo "=== solvers done $(date)"
for M in $MESHES; do
  P="--pairs results/pairs/$M.seed7.txt --truth-cache $OUT/$M.seed7.truth"
  for S in mmp vtp hybrid; do $B/cci_bench data/terrains/$M.off --eps 0.1 --gap-factor 0.5 --solver $S $=P --json $OUT/$M.$S.0.1.json; done
  $B/cci_bench data/terrains/$M.off --eps 0.2 --gap-factor 1 --solver hybrid $=P --json $OUT/$M.hybrid.0.2.json
  $B/cci_bench data/terrains/$M.off --eps 0.05 --gap-factor 0.25 --solver hybrid $=P --json $OUT/$M.hybrid.0.05.json
done
echo "=== index done $(date)"
run_upo() { # mesh tag poi_spec
  ( cd results/upo && lim 7200 ../../$B/upo_static ../../$1 $3 0.1 200 2>&1 | grep -E "RESULT|^mesh" > $2.txt; cp upo_pairs.txt $2.pairs; cat $2.txt )
  for C in "0.2 1" "0.1 0.5" "0.05 0.25"; do E=${C% *}; G=${C#* }
    $B/cci_bench $1 --eps $E --gap-factor $G --solver hybrid --pairs results/upo/$2.pairs --json results/upo/$2.cci.$E.json; done
}
run_upo $UPI/TJpre_1058.off tj_poi50 ../../$UPI/TJpre_50_poi_on_1058.txt
run_upo $UPI/TJpre_1058.off tj_all all
run_upo data/terrains/HorseMount.off horse_poi50 random50
run_upo data/terrains/HorseMount.off horse_all all
echo "=== UP-Oracle code done $(date)"
EAR_DIR=$PWD/external/weighted_distance_oracle scripts/competitors/rerun_weighted.sh
echo "=== weighted competitors done $(date)"
$B/cci_bench data/terrains/large/GunnisonForest100k.off --eps 0.2 --gap-factor 1 --solver hybrid --random 200 --seed 7 \
    --truth-cache $OUT/GunnisonForest100k.r200.truth --json $OUT/GunnisonForest100k.hybrid.0.2.json --verbose
echo "=== all done $(date)"
