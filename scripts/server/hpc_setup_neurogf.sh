#!/usr/bin/env bash
# One-time setup for the obstacles and neurogf phases, on the login node of the HSE cluster (after
# hpc_setup.sh): the conda environment `ngf` (PyTorch with CUDA, for NeuroGF) and our tools in cpp/bin2.
# cpp/bin, which the main run's queued jobs use, is not touched.
# usage: scripts/server/hpc_setup_neurogf.sh      (from the repository root)
set -eu
cd "$(dirname "$0")/../.."
if type module > /dev/null 2>&1; then module purge; module load python; fi
eval "$(conda shell.bash hook)"
if ! conda env list | grep -q '^ngf '; then
  conda create -y -n ngf --override-channels -c conda-forge python=3.11 numpy scipy matplotlib
fi
NGF=$(conda env list | awk '$1 == "ngf" {print $NF}')
# PyTorch wheels bundle their CUDA runtime (12.1); the GPU nodes' driver only has to be at least as new
$NGF/bin/pip install --quiet torch --index-url https://download.pytorch.org/whl/cu121
$NGF/bin/python -c "import torch, scipy, numpy; print('torch', torch.__version__, 'CUDA build', torch.version.cuda)"
set +u; source scripts/server/hpc_env.sh; set -u
make -C cpp -B -j8 BIN=bin2 all dgg
ls cpp/bin2
echo "neurogf setup done (NGF_PY=$NGF/bin/python)"
