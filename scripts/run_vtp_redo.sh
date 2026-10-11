#!/bin/zsh
# fourth batch: VTP-only and hybrid-threshold runs on Gunnison10k with the build that falls back to
# MMP when VTP rejects a sub-mesh (VTP's checks used to exit the process)
set -u
cd "$(dirname "$0")/.."
until grep -q "=== all done" results/logs/neurogf.log 2>/dev/null; do sleep 30; done
mv cpp/bin/cci_bench_new cpp/bin/cci_bench; cp cpp/bin/obj/certified_new.o cpp/bin/obj/certified.o; cp cpp/bin/obj/solver_vtp_new.o cpp/bin/obj/solver_vtp.o
OUT=results/cpp; B=cpp/bin
P="--pairs results/pairs/Gunnison10k.seed7.txt --truth-cache $OUT/Gunnison10k.seed7.truth"
$B/cci_bench data/terrains/Gunnison10k.off --eps 0.1 --gap-factor 0.5 --solver vtp $=P --json $OUT/Gunnison10k.vtp.0.1.json | grep -v degenerate
for RATIO in 4 8; do
  $B/cci_bench data/terrains/Gunnison10k.off --eps 0.1 --gap-factor 0.5 --solver hybrid --ratio $RATIO $=P --json $OUT/ratio/Gunnison10k.hybrid$RATIO.0.1.json | grep -v degenerate
done
echo "=== all done $(date)"
