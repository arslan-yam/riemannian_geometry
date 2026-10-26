#!/bin/bash
# EAR-Oracle (sp-num 5) on obstacle meshes, one job per node. Memory is not scheduled on the HSE cluster, and two
# runs (225-335 GB each) on one 768 GB node kill each other, so every job is pinned to a node picked by hand
# (type_d, >= 2 idle cores, enough free memory: see the sinfo line below) and asks for 2 cores (50% load).
#
#   scripts/server/ear_pinned.sh TJ.blobs05.V2V:cn-040 GI.walls.V2V:cn-031 ...
#   sinfo -N -p rocky -o "%n %f %C %e" | grep type_d | sort -k4 -n -r | head     # free memory in MB, cores A/I/O/T
#
# Each job writes the usual logs (results/server/logs/obst.ear5.<T>.<K>.<Q>.*) and its EAR log under
# $EAR_DIR/results. Score afterwards with scripts/server/ear_pinned.sh --eval <the same sets>.
set -e
cd "$(dirname "$0")/../.."
source scripts/server/hpc_env.sh
O=results/server; L=$O/logs; IN=$O/meshes
ACC=${SLURM_ACCOUNT:?set SLURM_ACCOUNT (e.g. proj_1925)}
EVAL=""; [ "$1" = --eval ] && { EVAL=1; shift; }
for spec in "$@"; do
  set_=${spec%%:*}; node=${spec#*:}; [ "$node" = "$spec" ] && node=""
  T=${set_%%.*}; rest=${set_#*.}; K=${rest%.*}; Q=${rest##*.}; M=${T}pre_500000.$K
  [ -f $IN/$M.off ] || { echo "no mesh $IN/$M.off" >&2; exit 1; }
  if [ -n "$EVAL" ]; then
    name=obst.eareval5.$T.$K.$Q
    det=""; [ $Q = V2V ] && det="--detour $O/comp/$M.V2V.detour"
    line="run $name bash -c 'python3 scripts/competitors/harness_eval.py $O/comp/$M.$Q.truth $EAR_DIR/results/${M}_${Q}_EAR_*.log $det --json $O/obst.ear.$M.$Q.json'"
    opts=("-c 1" "-t 1:00:00")
  else
    [ -n "$node" ] || { echo "$spec: give a node (set:node)" >&2; exit 1; }
    name=obst.ear5.$T.$K.$Q
    line="run $name env KIND=$Q timeout 72h scripts/competitors/run_competitors.sh $M $O/comp/$M.$Q.query EAR@0.2@5"
    opts=("-c 2" "-t 3-02:00:00" "--nodelist=$node")
  fi
  { echo "#!/bin/bash"
    for o in "-J $name" "-o $L/slurm/$name.%j.out" "${opts[@]}" "-A $ACC" "--constraint=type_d"; do echo "#SBATCH $o"; done
    echo "cd $PWD"; echo "source scripts/server/hpc_env.sh"
    echo "SLURM= scripts/server/run_server.sh __job <<'JOB'"; echo "$line"; echo "JOB"; } > $L/slurm/$name.sbatch
  echo "$name${node:+ on $node}: $(sbatch --parsable $L/slurm/$name.sbatch)"
done
