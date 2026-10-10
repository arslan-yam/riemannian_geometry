#!/bin/zsh
# re-run the competitors' weighted methods with the endpoint fix (harness.patch, item 6 of README.md)
# on weighted HorseMount, then re-score everything against the best valid path any method found
set -u
cd "$(dirname "$0")/../.."
EAR=${EAR_DIR:-$PWD/external/weighted_distance_oracle}
export EAR_DIR=$EAR
W=results/HorseMount_weighted
echo "=== start $(date)"
FACE_WEIGHTS=$W.face_weight.query KIND=W2V scripts/competitors/run_competitors.sh HorseMount $W.W2V.query \
    FixedS EAR EAR@0.2@10 EAR@0.2@20 FixedS@0.2@20
FACE_WEIGHTS=$W.face_weight.query KIND=W2V scripts/competitors/run_competitors.sh HorseMount $W.W2V.head30.query UnfixedS@0.05
python3 scripts/competitors/competitor_eval.py $EAR/results HorseMount W2V
echo "=== done $(date)"
