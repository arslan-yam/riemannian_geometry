#!/usr/bin/env bash
# Full-scale runs for the paper on a Linux server (notes/projection.md §4). Every job is single-threaded
# except the thread-scaling ones (marked), and jobs run side by side with xargs -P.
#
# usage: JOBS=48 scripts/server/run_server.sh [phase...]     phases: build truth ours scale threads poi
#        competitors small                                    (default: all, in this order)
# needs: g++ >= 9, make, the UP-Oracle data in external/UpdatedStructureTerrainCode/input,
#        for `competitors`: EAR_DIR (see scripts/competitors/README.md) and the UP-Oracle code built
#        (make -C cpp bin/upo_static).
# memory: our index with 8-bit labels peaks at about the coded size (~2-6 GB at 2.5M faces, eps 0.2);
#         EAR-Oracle needs > 32 GB at 500k faces. Set JOBS so that JOBS x the largest job fits.
set -u
cd "$(dirname "$0")/../.."
JOBS=${JOBS:-32}
IN=external/UpdatedStructureTerrainCode/input
X=cpp/bin/cci_exp
O=results/server; L=$O/logs; mkdir -p $O $L $O/truth
TERR="AU GI LH SC TJ VS"
SIZES="500000 1002528 1503378 2000000 2504322"
PHASES=${*:-build truth ours scale threads poi competitors small}
run() { # name, command...: runs one job, logs it, records wall time and peak memory
  local name=$1; shift
  /usr/bin/time -v "$@" > $L/$name.log 2> $L/$name.time || echo "FAILED $name" >> $L/failed.txt
}
export -f run; export L O X IN

for P in $PHASES; do
echo "=== $P $(date)"
case $P in
build)
  make -C cpp -j8 all ;;
truth) # 200 exact pairs per terrain and size (MMP, early stop), cached for every later job
  for T in $TERR; do for S in $SIZES; do
    echo "run truth.${T}pre_$S $X $IN/${T}pre_$S.off --random 200 --seed 7 --truth-cache $O/truth/${T}pre_$S.r200 --truth-only"
  done; done | xargs -P $JOBS -I{} bash -c '{}' ;;
ours) # the six terrains at 500k faces: eps 0.2 / g 1, kd and capped, 8-bit labels coded during the build
  for T in $TERR; do
    F=$IN/${T}pre_500000.off; R="--random 200 --seed 7 --truth-cache $O/truth/${T}pre_500000.r200"
    echo "run ours.$T.kd $X $F --eps 0.2 --gap-factor 1 $R --store code8 --refine --json $O/ours.$T.kd.json"
    echo "run ours.$T.k50 $X $F --eps 0.2 --gap-factor 1 --kmax 50 $R --store code8 --refine --json $O/ours.$T.k50.json"
    echo "run ours.$T.e01 $X $F --eps 0.1 --gap-factor 0.5 $R --store code8 --refine --json $O/ours.$T.e01.json"
  done | xargs -P $JOBS -I{} bash -c '{}' ;;
scale) # scaling curve on TJ and GI, 0.5M -> 2.5M faces (kd and capped; shifted from 1M on)
  for T in TJ GI; do for S in $SIZES; do
    F=$IN/${T}pre_$S.off; R="--random 200 --seed 7 --truth-cache $O/truth/${T}pre_$S.r200"
    echo "run scale.$T.$S.kd $X $F --eps 0.2 --gap-factor 1 $R --store code8 --refine --json $O/scale.$T.$S.kd.json"
    echo "run scale.$T.$S.k50 $X $F --eps 0.2 --gap-factor 1 --kmax 50 $R --store code8 --refine --json $O/scale.$T.$S.k50.json"
    [ $S != 500000 ] && echo "run scale.$T.$S.shifted $X $F --eps 0.2 --gap-factor 1 --hierarchy shifted $R --variants dense,code8 --refine --json $O/scale.$T.$S.shifted.json"
  done; done | xargs -P $JOBS -I{} bash -c '{}' ;;
threads) # MULTI-THREADED: thread scaling at 2.5M faces, one job at a time (labels are identical)
  F=$IN/TJpre_2504322.off; R="--random 200 --seed 7 --truth-cache $O/truth/TJpre_2504322.r200"
  for N in 8 16 32 64; do
    run threads.$N $X $F --eps 0.2 --gap-factor 1 $R --store code8 --threads $N --json $O/threads.$N.json
  done ;;
poi) # UP-Oracle's own setting: their POI files, 500k faces; the exact P x P table vs our POI index
  for T in $TERR; do for K in 50 100 500 1000 2500; do
    F=$IN/${T}pre_500000.off; PF=$IN/${T}pre_${K}_poi_on_500000.txt
    echo "run poi.$T.$K $X $F --poi-file $PF --random 1000 --seed 7 --poi-table --eps 0.2 --gap-factor 1 --store code8 --json $O/poi.$T.$K.json"
  done; done | xargs -P $JOBS -I{} bash -c '{}' ;;
competitors) # their code: UP-Oracle family on their POIs; DGG; EAR at 500k faces (needs > 32 GB)
  for T in $TERR; do for K in 50 500 2500; do
    echo "run upo.$T.$K cpp/bin/upo_static $IN/${T}pre_500000.off $IN/${T}pre_${K}_poi_on_500000.txt 0.1 200 UP,WSPD,RC,EAR"
  done; done | xargs -P $JOBS -I{} bash -c '{}'
  echo "EAR-Oracle and DGG at 500k+: see scripts/competitors/README.md (repair the meshes first)" ;;
small) # clean single-thread timings of the M1 tables (5 meshes, 3 configs, 100k x 2)
  for M in HorseMount BigMount HeadlightMount RobinsonMountain Gunnison10k; do
    for C in "0.2 1" "0.1 0.5" "0.05 0.25"; do E=${C% *}; G=${C#* }
      echo "run small.$M.$E $X data/terrains/$M.off --eps $E --gap-factor $G --pairs results/pairs/$M.seed7.txt --truth-cache results/cpp/$M.seed7.truth --variants dense,code16,code8 --refine --json $O/small.$M.$E.json"
    done
  done | xargs -P $JOBS -I{} bash -c '{}'
  for M in GunnisonForest100k Laramie100k; do
    run small.$M $X data/terrains/large/$M.off --eps 0.2 --gap-factor 1 --random 200 --seed 7 --truth-cache $O/truth/$M.r200 --variants dense,code16,code8 --refine --json $O/small.$M.json
  done ;;
esac
done
echo "=== all done $(date)"
