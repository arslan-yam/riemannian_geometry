#!/usr/bin/env bash
# Full-scale runs for the paper on a Linux server (notes/projection.md §4).
#
# usage: JOBS=48 JOBS_BIG=4 scripts/server/run_server.sh [phase...]
#   phases (default: all, in this order): build prepare truth ours caps scale proximity threads poi
#                                         competitors small
# needs: g++ >= 9, make, GNU xargs, python3 with numpy/scipy (rgeo), the UP-Oracle data in
#        external/UpdatedStructureTerrainCode/input; for `competitors`: EAR_DIR = a clone of
#        github.com/ItakEjgo/weighted_distance_oracle with scripts/competitors/harness.patch applied and
#        built in $EAR_DIR/bld (scripts/competitors/README.md), DGG built (cpp/bin/dgg_precompute,
#        cpp/bin/dgg_bench, see scripts/run_baselines.sh), and make -C cpp bin/upo_static.
#
# Every run is single-threaded (except the threads phase) and evaluates the same versions at every size,
# all on one build: labels float32 / 16-bit / 8-bit x near tables sym (default, lossless) / sym16 / sym8
# (ratio-coded). Query sets per run (exact truth: one exact sweep per source, cached by the truth phase):
#   random  200 sources x 10 targets (2,000 vertex pairs)
#   d0..dk  distance bands [0, R], (R, 2R], (2R, 4R], ... up to the diameter, R the near radius:
#           3 targets per band from each of the first 100 sources
#   near    2,000 pairs closer than R (answered by the near tables)
#   a2a     40 surface points x 5 targets (200 arbitrary-point pairs, certified with refine)
# memory: ~2-4 GB per job at 500k faces, ~40 GB at 2.5M faces (all label storages kept); EAR-Oracle
#         > 32 GB at 500k faces. JOBS_BIG limits the concurrent 2M / 2.5M-face and competitor jobs.
set -u
cd "$(dirname "$0")/../.."
JOBS=${JOBS:-32}
JOBS_BIG=${JOBS_BIG:-4}
SRC=external/UpdatedStructureTerrainCode/input
X=cpp/bin/cci_exp
O=results/server; L=$O/logs; IN=$O/meshes; mkdir -p $O $L $O/truth $IN $O/comp $O/dgg
TERR="AU GI LH SC TJ VS"
SIZES="500000 1002528 1503378 2000000 2504322"
PHASES=${*:-build prepare truth ours caps scale proximity threads poi competitors small}
QS="--sources 200 --per-source 10 --buckets 3 --bucket-sources 100 --near-pairs 2000 --a2a-sources 40 --a2a-per-source 5"
V="--variants dense,code16,code8 --near-variants sym,sym16,sym8"
truth() { echo "--truth-cache $O/truth/$1.s200x10"; }  # one cache name per terrain (every set derives from it)
run() { # name, command...: runs one job, logs it, records wall time and peak memory
  local name=$1; shift
  /usr/bin/time -v "$@" > $L/$name.log 2> $L/$name.time || echo "FAILED $name" >> $L/failed.txt
}
par() { xargs -d '\n' -P "$1" -I{} bash -c '{}'; } # one job per input line, $1 at a time
export -f run; export L O X IN

for P in $PHASES; do
echo "=== $P $(date)"
case $P in
build) # our tools, and the competitors built from source: DGG (our port) and the UP-Oracle static driver
  make -C cpp -j8 all dgg bin/upo_static ;;
prepare) # one prepared copy of every terrain for every method: folded faces repaired, chart at the origin
  for T in $TERR; do for S in $SIZES; do
    [ -f $IN/${T}pre_$S.off ] || echo "run prepare.${T}_$S python3 scripts/competitors/prepare_mesh.py $SRC/${T}pre_$S.off $IN/${T}pre_$S.off"
  done; done | par $JOBS
  cp $SRC/*_poi_on_*.txt $IN/ ;;
truth) # exact distances of every query set (eps 0.2 / g 1 and eps 0.1 / g 0.5 have the same near radius
       # R = 5 edges, so the same near pairs and bands), plus the competitors' query files of every terrain
  for T in $TERR; do for S in $SIZES; do
    echo "run truth.${T}_$S $X $IN/${T}pre_$S.off --eps 0.2 --gap-factor 1 $QS $(truth ${T}pre_$S) --export $O/comp/${T}pre_$S --truth-only"
  done; done | par $JOBS ;;
ours) # the six terrains at 500k faces: eps 0.2 / g 1 uncapped and cap 64, eps 0.1 / g 0.5
  for T in $TERR; do
    F=$IN/${T}pre_500000.off; R="$QS $(truth ${T}pre_500000)"
    echo "run ours.$T.kd $X $F --eps 0.2 --gap-factor 1 $R $V --refine --json $O/ours.$T.kd.json"
    echo "run ours.$T.k64 $X $F --eps 0.2 --gap-factor 1 --kmax 64 $R $V --refine --json $O/ours.$T.k64.json"
    echo "run ours.$T.e01 $X $F --eps 0.1 --gap-factor 0.5 $R $V --refine --json $O/ours.$T.e01.json"
  done | par $JOBS ;;
caps) # cap sweep (at most k portals per cut + certify-or-refine), TJ and GI at 500k faces
  for T in TJ GI; do for K in 8 16 32 64 128 256; do
    echo "run caps.$T.k$K $X $IN/${T}pre_500000.off --eps 0.2 --gap-factor 1 --kmax $K $QS $(truth ${T}pre_500000) $V --refine --json $O/caps.$T.k$K.json"
  done; done | par $JOBS ;;
scale) # scaling curve on TJ and GI, 0.5M -> 2.5M faces: uncapped and cap 64 at every size, shifted at 1M and
       # 1.5M (at 2M+ it would run ~6 days single-threaded; it already loses at 100k). 2M / 2.5M with JOBS_BIG
  job() { local T=$1 S=$2; local F=$IN/${T}pre_$S.off; local R="$QS $(truth ${T}pre_$S)"
    echo "run scale.$T.$S.kd $X $F --eps 0.2 --gap-factor 1 $R $V --refine --json $O/scale.$T.$S.kd.json"
    echo "run scale.$T.$S.k64 $X $F --eps 0.2 --gap-factor 1 --kmax 64 $R $V --refine --json $O/scale.$T.$S.k64.json"
    case $S in 1002528|1503378) echo "run scale.$T.$S.shifted $X $F --eps 0.2 --gap-factor 1 --hierarchy shifted $R $V --refine --json $O/scale.$T.$S.shifted.json";; esac; }
  for T in TJ GI; do for S in 500000 1002528 1503378; do job $T $S; done; done | par $JOBS
  for T in TJ GI; do for S in 2000000 2504322; do job $T $S; done; done | par $JOBS_BIG ;;
proximity) # kNN (k = 1, 10, 100) and range queries over 10,000 random objects from 100 query vertices:
           # certified (exact answer sets, refine only at the boundary) vs one exact sweep per query;
           # the six terrains at 500k faces and TJ at every size
  { for T in $TERR; do echo "$T 500000"; done; for S in 1002528 1503378 2000000 2504322; do echo "TJ $S"; done; } |
  while read T S; do
    echo "run prox.$T.$S cpp/bin/proximity $IN/${T}pre_$S.off --eps 0.2 --gap-factor 1 --store code8 --objects 10000 --queries 100 --k 1,10,100 --radius 0.01,0.03,0.1 --truth-cache $O/truth/${T}pre_$S.prox --json $O/prox.$T.$S.json"
  done | par $JOBS_BIG ;;
threads) # MULTI-THREADED: thread scaling at 2.5M faces, one job at a time (labels are identical)
  for N in 8 16 32 64; do
    run threads.$N $X $IN/TJpre_2504322.off --eps 0.2 --gap-factor 1 $QS $(truth TJpre_2504322) $V --threads $N --json $O/threads.$N.json
  done ;;
poi) # UP-Oracle's own setting: their POI files, 500k faces; the exact P x P table vs our POI index
  for T in $TERR; do for K in 50 100 500 1000 2500; do
    echo "run poi.$T.$K $X $IN/${T}pre_500000.off --poi-file $IN/${T}pre_${K}_poi_on_500000.txt --random 1000 --seed 7 --poi-table --eps 0.2 --gap-factor 1 $V --json $O/poi.$T.$K.json"
  done; done | par $JOBS ;;
competitors) # their code on the same prepared meshes and the same queries (exported by the truth phase)
  : ${EAR_DIR:?set EAR_DIR (see the header)}
  # 1. EAR-Oracle (authors' harness), V2V and A2A, sp 5 and 10: TJ and GI at 500k faces (> 32 GB, ~10+ h each)
  for T in TJ GI; do
    M=${T}pre_500000; cp $IN/$M.off $EAR_DIR/datasets/repaired/$M.off
    echo "run ear.$M.V2V env KIND=V2V timeout 72h scripts/competitors/run_competitors.sh $M $O/comp/$M.V2V.query EAR@0.2@5 EAR@0.2@10"
    echo "run ear.$M.A2A env KIND=A2A timeout 72h scripts/competitors/run_competitors.sh $M $O/comp/$M.A2A.query EAR@0.2@5 EAR@0.2@10"
  done | par $JOBS_BIG
  for T in TJ GI; do M=${T}pre_500000
    for K in V2V A2A; do python3 scripts/competitors/harness_eval.py $O/comp/$M.$K.truth $EAR_DIR/results/${M}_${K}_EAR_*.log --json $O/ear.$M.$K.json; done
  done
  # 2. DGG (authors' code, our port), eps 1e-3 and 1e-4, first 200 V2V pairs: TJ at every size, GI at 500k
  { for S in $SIZES; do echo "TJ $S"; done; echo "GI 500000"; } | while read T S; do
    for E in 0.001 0.0001; do
      echo "run dgg.${T}_$S.$E scripts/server/dgg_one.sh ${T}pre_$S $E"
    done
  done | par $JOBS_BIG
  # 3. exact on the fly: MMP and VTP full sweeps (5 sources), every size of TJ and GI
  for T in TJ GI; do for S in $SIZES; do
    echo "run onthefly.${T}_$S cpp/bin/solver_check $IN/${T}pre_$S.off 5"
  done; done | par $JOBS_BIG
  # 4. UP-Oracle family (their code) in their POI setting, 500k faces
  for T in $TERR; do for K in 50 500 2500; do
    echo "run upo.$T.$K cpp/bin/upo_static $IN/${T}pre_500000.off $IN/${T}pre_${K}_poi_on_500000.txt 0.1 200 UP,WSPD,RC,EAR"
  done; done | par $JOBS ;;
small) # clean single-thread reruns of the M1 tables (5 meshes x 3 configs, two 100k-vertex DEMs)
  for M in HorseMount BigMount HeadlightMount RobinsonMountain Gunnison10k; do
    for C in "0.2 1" "0.1 0.5" "0.05 0.25"; do E=${C% *}; G=${C#* }
      echo "run small.$M.$E $X data/terrains/$M.off --eps $E --gap-factor $G $QS $(truth $M.$E) $V --refine --json $O/small.$M.$E.json"
    done
  done | par $JOBS
  for M in GunnisonForest100k Laramie100k; do
    run small.$M $X data/terrains/large/$M.off --eps 0.2 --gap-factor 1 $QS $(truth $M) $V --refine --json $O/small.$M.json
  done ;;
esac
done
echo "=== all done $(date)"
