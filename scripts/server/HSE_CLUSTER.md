# Running on the HSE cluster (cHARISMa, hpc.hse.ru)

The cluster runs Slurm on Rocky Linux 9. `run_server.sh` with `SLURM=1` submits every job of the plan as a
Slurm job of its own, with its own cores, memory and time limit (`res()` in the script). Phases are linked
by dependencies, so one command submits the whole plan. These are the cluster's rules this setup follows:

| cluster rule (hpc.hse.ru/support/instructions) | what we do |
|---|---|
| no computing on the login node | the login node only builds and submits |
| compute nodes have no internet | `hpc_setup.sh` (conda env, EAR clone, builds) runs on the login node |
| 24 h default wall time | every job gets its own `-t`, at most 6 days (the shifted hierarchy at 1.5M faces) |
| a job whose CPUs stay < 10% loaded for ~3 h can be cancelled | one Slurm job per run, so no idle cores waiting for the slowest run of a phase |
| node types differ | `-C type_d` for every job: 2× Xeon Gold 6248R (48 cores), 768 GB. One CPU type for every timing |

## 0. Once: account, key, project

- Get the key `auth_key.pem` at https://lk.hpc.hse.ru/manage_access/ (it has a passphrase).
- On the Mac, add this to `~/.ssh/config` (replace `<login>` and the key path):

  ```
  Host hse-hpc
    HostName cluster.hpc.hse.ru
    Port 2222
    User <login>
    IdentityFile ~/keys/auth_key.pem
  ```

  Then make the key private: `chmod 600 ~/keys/auth_key.pem`, and check access: `ssh hse-hpc`.
- On the cluster, `mp` lists your project IDs (the `-A` account). `checkquota` shows your quota: the run
  needs ~200 GB. If home is too small, use a workspace for the whole checkout: `ws_allocate rg 30`, which
  makes `/scratch/ws/<login>-rg`, lives 30 days, and can be extended with `ws_extend rg 30`.
- Check the longest wall time each partition allows with `sinfo -o "%P %l"`. The plan asks for up to 6
  days per job (`res()` in `run_server.sh`). If the limit is lower, shorten those `t=` values.

## 1. Code and data

On the cluster (login node):

```bash
cd ~ && git clone https://github.com/arslan-yam/riemannian_geometry.git rg
```

Everything else comes from GitHub too, cloned by `hpc_setup.sh` in step 2:
- the UP-Oracle authors' code and their terrains (`yanyinzhao/UpdatedStructureTerrainCode` at 0ec6aec), which
  are the inputs of every experiment;
- EAR-Oracle (`ItakEjgo/weighted_distance_oracle` at 2343aa9).

Our DGG port, its shim headers and the two 100k-vertex DEMs are in this repository.
For later code changes: commit and push on the Mac, then `cd ~/rg && git pull` on the cluster.

## 2. Once: setup on the login node (~20–30 min)

This makes the conda environment `rg` (Python, numpy/scipy/matplotlib, CMake, CGAL, GMP, MPFR, GNU time),
clones UP-Oracle's code and terrains, clones, patches and builds EAR-Oracle, then builds our tools, DGG and
the UP-Oracle driver with the system GCC.

```bash
ssh hse-hpc
cd ~/rg && scripts/server/hpc_setup.sh 2>&1 | tee setup.log
```

It should end with `setup done`. Before every later session on the login node, load the environment:

```bash
cd ~/rg && source scripts/server/hpc_env.sh
```

## 3. Smoke test (one 4-core job, a few minutes)

```bash
SLURM=1 SLURM_ACCOUNT=<project> scripts/server/run_server.sh smoke
```

Wait for it (`squeue -u $USER`), then check its output: it should say `12 finished, 0 failed` and no
"not complete".

```bash
cat results/server/logs/slurm/smoke.*.out
```

## 4. The whole plan (260 jobs; roughly 4,000–5,000 core-hours)

```bash
SLURM=1 SLURM_ACCOUNT=<project> THREAD_COUNTS="8 16 32 48" SOLVER_MEM_GB=150 scripts/server/run_server.sh prepare truth ours caps curve ablation solvers weighted scale proximity threads poi competitors small 2>&1 | tee results/server/submit.log
```

This returns at once. Every job of `ours` … `small` waits for all prepare and truth jobs (`afterany`).
`weighted` also waits for its Steiner reference, and the EAR scoring waits for the EAR runs.
`THREAD_COUNTS` stops at 48, the core count of a type D node; those jobs take a whole node (`--exclusive`).

How long it takes depends on the queue. The longest jobs are the 2.5M-face builds (2–3 days) and EAR-Oracle
(up to its 72 h timeout).

## 5. Monitoring

```bash
squeue -u $USER                                   # pending and running jobs (or: mj)
python3 scripts/server/status.py                  # finished / failed jobs, wall time, peak RAM, incomplete JSON
sacct -u $USER -S today --format=JobName%40,State,Elapsed,MaxRSS,ReqMem
tail -3 results/server/logs/<job>.log             # what one job is doing
```

## 6. If a job fails

- Read `results/server/logs/<job>.log` and `.time`, plus the Slurm file `logs/slurm/<job>.<id>.out`. An
  `OUT_OF_MEMORY` or `TIMEOUT` in `sacct` means its `--mem` or `-t` was too small.
- Resubmit the one job: `sbatch results/server/logs/slurm/<job>.sbatch`. You can edit its `#SBATCH` lines
  first, for example to raise `--mem`. The truth caches are reused.
- Cancel everything: `scancel -u $USER`.
- The multithreaded jobs (`prox.*`, `weighted.*.t16`, `threads.*`) answer their queries on one thread after
  the parallel build. If the cluster emails a low-load warning about one of them, that is why. A
  cancelled one can be resubmitted as above.
- Maintenance is the last working day of each month (the next one is Fri 30 Oct): running jobs are killed
  and requeued from the start. Submit so that the 3–6 day jobs end before it.

## 7. Copy the results back (on the Mac)

```bash
rsync -av hse-hpc:rg/results/server/ results/server/ --exclude meshes/ --exclude 'dgg/*/*.binary' --exclude 'dgg/*/*.off'
```

The partial JSON files are valid, so you can copy while jobs are still running. Copy the results back before
the files expire: workspaces after 30 days (plus extensions), home 3 months after the project ends.
