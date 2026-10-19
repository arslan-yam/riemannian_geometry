#!/usr/bin/env bash
# usage: [KIND=V2V|A2A] run_competitors.sh <MeshName> <query file> <methods...>   (method[@eps[@sp]])
set -e
EAR=${EAR_DIR:?set EAR_DIR to the clone of ItakEjgo/weighted_distance_oracle (built in $EAR_DIR/bld, repaired meshes in $EAR_DIR/datasets/repaired)}
MESH=$1; Q=$2; shift 2
Q=$(cd "$(dirname "$Q")" && pwd)/$(basename "$Q") #absolute: we cd into the work dir below
if [ -n "${FACE_WEIGHTS:-}" ]; then FACE_WEIGHTS=$(cd "$(dirname "$FACE_WEIGHTS")" && pwd)/$(basename "$FACE_WEIGHTS"); fi
KIND=${KIND:-V2V}
W=$EAR/bld/work_${MESH}_$KIND; mkdir -p $W; cd $W
cp $Q A2A.query
NF=$(sed -n 2p $EAR/datasets/repaired/$MESH.off | awk '{print $2}')
if [ -n "${FACE_WEIGHTS:-}" ]; then cp $FACE_WEIGHTS face_weight.query; else yes 1.000000 | head -n $NF > face_weight.query; fi
NQ=$(wc -l < A2A.query)
mkdir -p $EAR/results
FAILED=0
for M in "$@"; do
  EPS=0.2; SP=5; MM=$M
  case $M in *@*@*) MM=${M%%@*}; R=${M#*@}; EPS=${R%@*}; SP=${R#*@};; *@*) MM=${M%@*}; EPS=${M#*@};; esac
  OUT=$EAR/results/${MESH}_${KIND}_${MM}_eps${EPS}_sp${SP}.log
  echo "=== $MESH $KIND $MM eps=$EPS sp-num=$SP -> $OUT"
  rm -f $OUT
  # the full console output goes to $OUT.out (the progress lines are long); a run that dies (a signal, a non-zero exit,
  # or fewer answers in its log than queries) is reported here and fails the job, after the remaining settings have run
  rc=0; ${GTIME:-/usr/bin/time} -p $EAR/bld/main --generate=0 --input=$EAR/datasets/repaired/$MESH.off --output=$OUT --grid-num=16 --query-num=$NQ --method=$MM --eps=$EPS --sp-num=$SP > $OUT.out 2>&1 || rc=$?
  grep -E "^(real|user) |Index|memory|signal|exit|status|rror|bad_alloc|Abort|Segmentation|Killed" $OUT.out | grep -v "^Load" || true
  NA=$(awk '/Query results end/{f=0} f&&NF{n++} /Query results begin/{f=1} END{print n+0}' $OUT 2>/dev/null || echo 0)
  if [ $rc -ne 0 ] || [ "$NA" -lt "$NQ" ]; then
    echo "FAILED: $MM eps=$EPS sp-num=$SP (exit $rc, $NA of $NQ answers); last lines of $OUT.out:"; tail -3 $OUT.out; FAILED=1
  fi
done
exit $FAILED
