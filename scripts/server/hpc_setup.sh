#!/usr/bin/env bash
# One-time setup on the login node of the HSE cluster (it has internet; the compute nodes do not):
# the conda environment `rg`, UP-Oracle's code and terrains (clone at 0ec6aec), EAR-Oracle (clone at 2343aa9,
# our patch, build) and our tools.
# usage: scripts/server/hpc_setup.sh      (from the repository root, after copying it to the cluster)
set -eu
cd "$(dirname "$0")/../.."
if type module > /dev/null 2>&1; then module purge; module load python; fi
eval "$(conda shell.bash hook)"
if ! conda env list | grep -q '^rg '; then
  conda create -y -n rg --override-channels -c conda-forge python=3.11 numpy scipy matplotlib cmake make cgal-cpp gmp mpfr time
fi
set +u; source scripts/server/hpc_env.sh; set -u  # conda's activation scripts use unset variables
echo "python: $(which python3), GNU time: $GTIME, g++: $(g++ --version | head -1)"

# UP-Oracle: the authors' code and their terrains (the inputs of every experiment; 2.8 GB)
U=external/UpdatedStructureTerrainCode
[ -d $U/.git ] || git clone https://github.com/yanyinzhao/UpdatedStructureTerrainCode $U
(cd $U && git checkout -q 0ec6aec)
ls $U/input/TJpre_2504322.off

# EAR-Oracle: the authors' code plus our patch (weighted-endpoint fix, query harness)
E=external/weighted_distance_oracle
[ -d $E/.git ] || git clone https://github.com/ItakEjgo/weighted_distance_oracle $E
(cd $E && git checkout -q 2343aa9 &&
 if git apply --check ../../scripts/competitors/harness.patch 2> /dev/null; then git apply ../../scripts/competitors/harness.patch;
 else echo "harness.patch: already applied"; fi &&
 mkdir -p bld datasets/repaired && cd bld && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX && make -j8)
ls -l $E/bld/main

# our tools, DGG (our port) and the UP-Oracle driver
make -C cpp -B -j8 all dgg bin/upo_static
echo "setup done"
