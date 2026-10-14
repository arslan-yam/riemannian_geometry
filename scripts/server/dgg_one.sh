#!/usr/bin/env bash
# DGG on one prepared terrain: precompute the discrete geodesic graph (authors' code), then answer the
# first 200 exported V2V pairs (one Dijkstra each) against the exact truth.
# usage: scripts/server/dgg_one.sh <mesh name, e.g. TJpre_500000> <eps, e.g. 0.001>
set -eu
cd "$(dirname "$0")/../.."
M=$1; E=$2; O=results/server; ROOT=$PWD
W=$O/dgg/$M.$E; mkdir -p $W
cp $O/meshes/$M.off $W/
head -200 $O/comp/$M.V2V.pairs > $W/pairs
head -200 $O/comp/$M.V2V.truth > $W/truth
cd $W
/usr/bin/time -v $ROOT/cpp/bin/dgg_precompute f $M.off $E 1 > precompute.log 2> precompute.time
G=$(ls ${M}_FD*_c5.binary | head -1)
$ROOT/cpp/bin/dgg_bench $G pairs truth > bench.log
cat bench.log
