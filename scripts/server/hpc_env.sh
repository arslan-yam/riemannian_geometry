# Environment for the HSE cHARISMa cluster (hpc.hse.ru). Every Slurm job of run_server.sh sources it first,
# and so do you on the login node before building or submitting:  source scripts/server/hpc_env.sh
# Python, GNU time and EAR-Oracle's build dependencies (CMake, CGAL, GMP, MPFR) come from the conda
# environment `rg` (made by scripts/server/hpc_setup.sh). The C++ compiler is the system GCC of Rocky Linux 9.
if type module > /dev/null 2>&1; then module purge > /dev/null 2>&1; module load python > /dev/null 2>&1; fi
eval "$(conda shell.bash hook)"
conda activate ${RG_ENV:-rg}
export EAR_DIR=${EAR_DIR:-$PWD/external/weighted_distance_oracle}
if [ -x "$CONDA_PREFIX/bin/time" ]; then export GTIME=$CONDA_PREFIX/bin/time; fi
# NeuroGF runs in its own environment (PyTorch with CUDA; scripts/server/hpc_setup_neurogf.sh)
export NGF_PY=${NGF_PY:-$(dirname "$CONDA_PREFIX")/ngf/bin/python}
