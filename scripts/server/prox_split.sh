#!/usr/bin/env bash
# Split a proximity export ("query object distance" per line) into the "pairs" and "truth" files that
# landmark_bench and sep_bench read. usage: scripts/server/prox_split.sh PREFIX.pairs OUT  -> OUT.pp, OUT.tr
set -eu
awk '{print $1, $2}' "$1" > "$2.pp"
awk '{print $3}' "$1" > "$2.tr"
