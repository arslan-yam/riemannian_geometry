#!/usr/bin/env bash
# Full-scale runs for the paper on a Linux server (notes/projection.md §4).
#
# usage: JOBS=48 JOBS_BIG=4 JOBS_PROX=2 scripts/server/run_server.sh [phase...]
#   phases (default: all, in this order): build prepare truth ours caps curve ablation solvers weighted
#                                         scale proximity threads poi competitors small
#   smoke (not in the default list): every tool once on TJpre_1058, a few minutes; run it first
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
JOBS_PROX=${JOBS_PROX:-2}  # proximity jobs at once (16 threads each): keep JOBS_PROX x 16 <= cores
SRC=external/UpdatedStructureTerrainCode/input
X=cpp/bin/cci_exp
O=results/server; L=$O/logs; IN=$O/meshes; mkdir -p $O $L $O/truth $IN $O/comp $O/dgg
TERR="AU GI LH SC TJ VS"
SIZES="500000 1002528 1503378 2000000 2504322"
PHASES=${*:-build prepare truth ours caps curve ablation solvers weighted scale proximity threads poi competitors small}
SOLVER_MEM_GB=${SOLVER_MEM_GB:-150}  # threads phase: budget for the label solvers alive at once (VTP ~0.7 KB,
                                     # MMP ~1.6 KB per face of the sub-mesh: ~1.6 / 3.5 GB per thread at 2.5M faces)
QS="--sources 200 --per-source 10 --buckets 3 --bucket-sources 100 --near-pairs 2000 --a2a-sources 40 --a2a-per-source 5"
V="--variants dense,code16,code8 --near-variants sym,sym16,sym8"
truth() { echo "--truth-cache $O/truth/$1.s200x10"; }  # one cache name per terrain (every set derives from it)
run() { # name, command...: runs one job; its output in $L/name.log, GNU time in $L/name.time, and one JSON
        # line per job in $L/jobs.jsonl (exit code, wall seconds, peak RSS)
  local name=$1; shift
  /usr/bin/time -v "$@" > $L/$name.log 2> $L/$name.time; local rc=$?
  [ $rc -ne 0 ] && echo "FAILED $name" >> $L/failed.txt
  local wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {n = split($2, a, ":"); s = 0; for (i = 1; i <= n; i++) s = s * 60 + a[i]; print s}' $L/$name.time)
  local rss=$(awk -F': ' '/Maximum resident set size/ {print $2}' $L/$name.time)
  printf '{"job":"%s","exit":%d,"wall_s":%s,"max_rss_kb":%s,"end":"%s"}\n' "$name" $rc "${wall:-null}" "${rss:-null}" "$(date -Is)" >> $L/jobs.jsonl
  return $rc
}
par() { xargs -d '\n' -P "$1" -I{} bash -c '{}'; } # one job per input line, $1 at a time
export -f run; export L O X IN

for P in $PHASES; do
echo "=== $P $(date)"
case $P in
build) # our tools, and the competitors built from source: DGG (our port) and the UP-Oracle static driver
  make -C cpp -B -j8 all dgg bin/upo_static ;;  # -B: cpp/bin holds macOS objects tracked in git
prepare) # one prepared copy of every terrain for every method: folded faces repaired, chart at the origin
  for T in $TERR; do for S in $SIZES; do
    [ -f $IN/${T}pre_$S.off ] || echo "run prepare.${T}_$S python3 scripts/competitors/prepare_mesh.py $SRC/${T}pre_$S.off $IN/${T}pre_$S.off"
  done; done | par $JOBS
  cp $SRC/*_poi_on_*.txt $IN/
  for T in TJ GI; do python3 scripts/competitors/make_weights.py $IN/${T}pre_500000.off 1 $IN/${T}pre_500000.w; done ;;
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
curve) # size versus error (the matched-accuracy comparison with EAR-Oracle at scale), TJ and GI at 500k:
       # eps from 0.05 to 1.6 at a fixed near radius R = g / eps = 5 edges (g = 5 eps edges), and eps 0.1 /
       # 0.4 at a fixed g = 1 edge (R = 10 / 2.5 edges). Near pairs and bands depend on R: cached per R
  for T in TJ GI; do
    F=$IN/${T}pre_500000.off
    for C in "0.05 0.25" "0.4 2" "0.8 4" "1.6 8" "0.1 1" "0.4 1"; do E=${C% *}; G=${C#* }
      echo "run curve.$T.e$E.g$G $X $F --eps $E --gap-factor $G $QS $(truth ${T}pre_500000) $V --refine --json $O/curve.$T.e$E.g$G.json"
    done
  done | par $JOBS ;;
ablation) # query-time ablations on one build each (TJ and GI at 500k, eps 0.2 / g 1): estimate interpolation
          # none / parabola / model, and eps_c-covers (dropped at 10k and 100k: dominated by 8-bit labels)
  for T in TJ GI; do
    echo "run ablation.$T $X $IN/${T}pre_500000.off --eps 0.2 --gap-factor 1 $QS $(truth ${T}pre_500000) --variants dense,code8,cover0.01,cover0.05 --near-variants sym --interp-variants 0,1,2 --refine --json $O/ablation.$T.json"
  done | par $JOBS ;;
weighted) # weighted terrains (face weights uniform in (1, 5)), TJ and GI at 500k faces: labels from the Steiner
          # graph (m points per edge); no exact solver exists, so the reference is a finer Steiner graph
          # (m = 12; graph distances are upper bounds that converge to the weighted distance) and the LB is
          # the chord bound. m = 5 single-threaded (~7x an unweighted build); m = 10 with 16 threads (marked)
  for T in TJ GI; do
    F=$IN/${T}pre_500000.off; WT="--weights $IN/${T}pre_500000.w --solver steiner --truth-steiner 12"
    TR="--truth-cache $O/truth/${T}pre_500000.w.s200x10"
    run truthw.$T $X $F $WT --eps 0.2 --gap-factor 1 $QS $TR --export $O/comp/${T}pre_500000.w --truth-only
  done
  for T in TJ GI; do
    F=$IN/${T}pre_500000.off; WT="--weights $IN/${T}pre_500000.w --solver steiner --truth-steiner 12"
    TR="--truth-cache $O/truth/${T}pre_500000.w.s200x10"
    echo "run weighted.$T.m5 $X $F $WT --steiner-m 5 --eps 0.2 --gap-factor 1 $QS $TR $V --json $O/weighted.$T.m5.json"
    echo "run weighted.$T.m10.t16 $X $F $WT --steiner-m 10 --threads 16 --eps 0.2 --gap-factor 1 $QS $TR $V --json $O/weighted.$T.m10.t16.json"
  done | par $JOBS_BIG ;;
solvers) # ablation of the label solver: hybrid (ratio 2 = default, 4), VTP only, MMP only; TJ and GI at 500k
  for T in TJ GI; do
    F=$IN/${T}pre_500000.off; R="$QS $(truth ${T}pre_500000)"
    for C in "hybrid 2" "hybrid 4" "vtp 2" "mmp 2"; do SV=${C% *}; RT=${C#* }
      echo "run solvers.$T.$SV$RT $X $F --eps 0.2 --gap-factor 1 --solver $SV --ratio $RT $R --variants code8 --near-variants sym --json $O/solvers.$T.$SV$RT.json"
    done
  done | par $JOBS ;;
scale) # scaling curve on TJ and GI, 0.5M -> 2.5M faces: uncapped and cap 64 at every size, shifted at 1M and
       # 1.5M (at 2M+ it would run ~6 days single-threaded; it already loses at 100k). 2M / 2.5M with JOBS_BIG
  job() { local T=$1 S=$2; local F=$IN/${T}pre_$S.off; local R="$QS $(truth ${T}pre_$S)"
    echo "run scale.$T.$S.kd $X $F --eps 0.2 --gap-factor 1 $R $V --refine --json $O/scale.$T.$S.kd.json"
    echo "run scale.$T.$S.k64 $X $F --eps 0.2 --gap-factor 1 --kmax 64 $R $V --refine --json $O/scale.$T.$S.k64.json"
    case $S in 1002528|1503378) echo "run scale.$T.$S.shifted $X $F --eps 0.2 --gap-factor 1 --hierarchy shifted $R $V --refine --json $O/scale.$T.$S.shifted.json";; esac; }
  for T in TJ GI; do for S in 500000 1002528 1503378; do job $T $S; done; done | par $JOBS
  for T in TJ GI; do for S in 2000000 2504322; do job $T $S; done; done | par $JOBS_BIG ;;
proximity) # MULTI-THREADED build (16 threads: only the query times are reported; the index is identical)
           # kNN (k = 1, 10, 100) and range queries over 10,000 random objects from 100 query vertices:
           # certified (exact answer sets, refine only at the boundary) vs one exact sweep per query;
           # the six terrains at 500k faces and TJ at every size
  { for T in $TERR; do echo "$T 500000"; done; for S in 1002528 1503378 2000000 2504322; do echo "TJ $S"; done; } |
  while read T S; do
    echo "run prox.$T.$S cpp/bin/proximity $IN/${T}pre_$S.off --eps 0.2 --gap-factor 1 --store code8 --threads 16 --solver-mem-gb 60 --objects 10000 --queries 100 --k 1,10,100 --radius 0.01,0.03,0.1 --truth-cache $O/truth/${T}pre_$S.prox --json $O/prox.$T.$S.json"
  done | par $JOBS_PROX ;;
threads) # MULTI-THREADED: thread scaling at 2.5M faces, one job at a time (labels are identical)
  for N in 8 16 32 64; do
    run threads.$N $X $IN/TJpre_2504322.off --eps 0.2 --gap-factor 1 $QS $(truth TJpre_2504322) $V --threads $N --solver-mem-gb $SOLVER_MEM_GB --json $O/threads.$N.json
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
  # 1b. EAR-Oracle on the weighted terrains (W2V, the same face weights; scored against the Steiner reference)
  for T in TJ GI; do
    M=${T}pre_500000
    echo "run ear.$M.W2V env KIND=W2V FACE_WEIGHTS=$IN/$M.w timeout 72h scripts/competitors/run_competitors.sh $M $O/comp/$M.w.V2V.query EAR@0.2@5 EAR@0.2@10"
  done | par $JOBS_BIG
  for T in TJ GI; do M=${T}pre_500000
    for K in V2V A2A; do python3 scripts/competitors/harness_eval.py $O/comp/$M.$K.truth $EAR_DIR/results/${M}_${K}_EAR_*.log --json $O/ear.$M.$K.json; done
    python3 scripts/competitors/harness_eval.py $O/comp/$M.w.V2V.truth $EAR_DIR/results/${M}_W2V_EAR_*.log --json $O/ear.$M.W2V.json
  done
  # 2. DGG (authors' code, our port), eps 1e-3 and 1e-4, first 200 V2V pairs: TJ at every size, GI at 500k
  { for S in $SIZES; do echo "TJ $S"; done; echo "GI 500000"; } | while read T S; do
    for E in 0.001 0.0001; do
      echo "run dgg.${T}_$S.$E scripts/server/dgg_one.sh ${T}pre_$S $E"
    done
  done | par $JOBS_BIG
  # 3. exact on the fly: MMP and VTP full sweeps (5 sources), every size of TJ and GI
  for T in TJ GI; do for S in $SIZES; do
    echo "run onthefly.${T}_$S cpp/bin/solver_check $IN/${T}pre_$S.off 5 --json $O/onthefly.${T}_$S.json"
  done; done | par $JOBS_BIG
  # 4. UP-Oracle family (their code) in their POI setting, 500k faces
  for T in $TERR; do for K in 50 500 2500; do
    echo "run upo.$T.$K cpp/bin/upo_static $IN/${T}pre_500000.off $IN/${T}pre_${K}_poi_on_500000.txt 0.1 200 UP,WSPD,RC,EAR --json $O/upo.$T.$K.json"
  done; done | par $JOBS ;;
smoke) # every tool and the whole pipeline once on a 1,058-face terrain: results in $O/smoke
  S=$O/smoke; M=TJpre_1058; F=$IN/$M.off; mkdir -p $S
  python3 scripts/competitors/prepare_mesh.py $SRC/$M.off $F | tail -1
  python3 scripts/competitors/make_weights.py $F 1 $IN/$M.w
  P=$SRC/TJpre_50_poi_on_1058.txt
  SQ="--sources 20 --per-source 10 --buckets 3 --bucket-sources 10 --near-pairs 200 --a2a-sources 10 --a2a-per-source 5"
  WT="--weights $IN/$M.w --solver steiner --truth-steiner 8"
  run smoke.truth $X $F --eps 0.2 --gap-factor 1 $SQ --truth-cache $S/$M.s20x10 --export $O/comp/$M --truth-only
  run smoke.cci $X $F --eps 0.2 --gap-factor 1 --kmax 64 $SQ --truth-cache $S/$M.s20x10 $V --refine --json $S/cci.json
  run smoke.threads $X $F --eps 0.2 --gap-factor 1 $SQ --truth-cache $S/$M.s20x10 --variants code8 --threads 4 --solver-mem-gb 4 --json $S/threads.json
  run smoke.weighted $X $F $WT --steiner-m 3 --eps 0.2 --gap-factor 1 $SQ --truth-cache $S/$M.w.s20x10 --export $O/comp/$M.w $V --json $S/weighted.json
  run smoke.poi $X $F --poi-file $P --random 200 --seed 7 --poi-table --eps 0.2 --gap-factor 1 $V --json $S/poi.json
  run smoke.prox cpp/bin/proximity $F --eps 0.2 --gap-factor 1 --store code8 --threads 4 --solver-mem-gb 4 --objects 200 --queries 10 --k 1,10 --radius 0.05 --truth-cache $S/$M.prox --json $S/prox.json
  run smoke.onthefly cpp/bin/solver_check $F 2 --json $S/onthefly.json
  run smoke.upo cpp/bin/upo_static $F $P 0.1 50 UP,WSPD,RC,EAR --json $S/upo.json
  run smoke.dgg scripts/server/dgg_one.sh $M 0.001 && mv $O/dgg.$M.0.001.json $S/dgg.json
  if [ -n "${EAR_DIR:-}" ]; then
    cp $F $EAR_DIR/datasets/repaired/$M.off
    run smoke.ear.V2V env KIND=V2V scripts/competitors/run_competitors.sh $M $O/comp/$M.V2V.query EAR@0.2@5
    run smoke.ear.A2A env KIND=A2A scripts/competitors/run_competitors.sh $M $O/comp/$M.A2A.query EAR@0.2@5
    run smoke.ear.W2V env KIND=W2V FACE_WEIGHTS=$IN/$M.w scripts/competitors/run_competitors.sh $M $O/comp/$M.w.V2V.query EAR@0.2@5
    for K in V2V A2A; do python3 scripts/competitors/harness_eval.py $O/comp/$M.$K.truth $EAR_DIR/results/${M}_${K}_EAR_eps0.2_sp5.log --json $S/ear.$K.json > /dev/null; done
    python3 scripts/competitors/harness_eval.py $O/comp/$M.w.V2V.truth $EAR_DIR/results/${M}_W2V_EAR_eps0.2_sp5.log --json $S/ear.W2V.json > /dev/null
  else echo "EAR_DIR not set: EAR-Oracle skipped"; fi
  python3 scripts/server/status.py $O --only smoke ;;
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
