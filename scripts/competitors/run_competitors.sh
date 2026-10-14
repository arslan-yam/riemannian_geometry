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
for M in "$@"; do
  EPS=0.2; SP=5; MM=$M
  case $M in *@*@*) MM=${M%%@*}; R=${M#*@}; EPS=${R%@*}; SP=${R#*@};; *@*) MM=${M%@*}; EPS=${M#*@};; esac
  OUT=$EAR/results/${MESH}_${KIND}_${MM}_eps${EPS}_sp${SP}.log
  echo "=== $MESH $KIND $MM eps=$EPS sp-num=$SP -> $OUT"
  /usr/bin/time -p $EAR/bld/main --generate=0 --input=$EAR/datasets/repaired/$MESH.off --output=$OUT --grid-num=16 --query-num=$NQ --method=$MM --eps=$EPS --sp-num=$SP 2>&1 | grep -E "real|user|Index|memory" || true
done
