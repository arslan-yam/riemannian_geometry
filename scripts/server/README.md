# Server run: step-by-step guide

Everything runs from `scripts/server/run_server.sh`, one phase after another. Inside a phase, independent
single-threaded jobs run side by side. The phases and their expected costs are in `notes/projection.md` §4.
Only the `threads`, `proximity` and `weighted` (m = 10) phases use threads.

Throughout, `SERVER` is your ssh host and `~/rg` is the checkout on the server; change both to taste.

## 0. Hardware

| | needed |
|---|---|
| OS | Linux x86-64; tested commands assume Ubuntu 22.04 / 24.04 |
| CPU | 32–64 cores, good single-core speed (the comparison runs are single-threaded) |
| RAM | 256 GB. Peak per job: ~2–4 GB at 500k faces, ~35–40 GB at 2.5M faces (`JOBS_BIG` of those at once), EAR-Oracle > 32 GB at 500k faces, the weighted reference (Steiner m = 12 at 500k faces) ~10 GB |
| disk | ~200 GB free (meshes 3 GB, truth caches, logs, EAR / DGG work dirs) |
| GPU | not needed (NeuroGF comes later) |

## 1. Copy the code and the data (from the Mac)

`external/` and `data/terrains/large/` are gitignored, and `cpp/bin/` holds macOS binaries, so copy the
working tree with rsync rather than `git clone`. That is ~1.6 GB: the code, the UP-Oracle code and its
`*pre_*` meshes and POI files (1.3 GB), the DGG port, and the two 100k-vertex DEMs.

```bash
cd ~/hse/riemannian_geometry
rsync -av --progress \
  --exclude .git --exclude /cpp/bin/ --exclude /cpp/bin_dev/ --exclude /results/ \
  --exclude /external/weighted_distance_oracle/ --exclude /external/NeuroGF/ --exclude /external/DGG_and_SVG/ \
  --exclude '/external/UpdatedStructureTerrainCode/input/*post*' --exclude /external/UpdatedStructureTerrainCode/output/ \
  ./ SERVER:rg/
```

To push later code changes, run the same command again; rsync sends only what changed.

## 2. Install packages (on the server, once)

```bash
sudo apt update
sudo apt install -y build-essential cmake time tmux rsync git libcgal-dev libboost-all-dev libgmp-dev libmpfr-dev \
  python3 python3-numpy python3-scipy python3-matplotlib
```

`time` is GNU time (`/usr/bin/time -v`), which the scripts use for wall time and peak memory.

## 3. Build EAR-Oracle (on the server, once)

This is the authors' code at commit 2343aa9 plus our patch (the weighted-endpoint fix and the query
harness). The patch builds with CGAL 5 or 6.

```bash
cd ~/rg
git clone https://github.com/ItakEjgo/weighted_distance_oracle external/weighted_distance_oracle
cd external/weighted_distance_oracle && git checkout 2343aa9 && git apply ../../scripts/competitors/harness.patch
mkdir -p bld datasets/repaired && cd bld && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j8
cd ~/rg && ls -l external/weighted_distance_oracle/bld/main  # must exist
```

Our tools, DGG (our port of the authors' code) and the UP-Oracle driver are built by the `build` phase.
The build is checked with GCC 15 and Apple clang, which give identical results on the smoke test.

## 4. Smoke test (~5 minutes): every tool once on a 1,058-face terrain

```bash
cd ~/rg && mkdir -p results/server/logs
export EAR_DIR=$PWD/external/weighted_distance_oracle
scripts/server/run_server.sh build smoke 2>&1 | tee results/server/smoke.log
```

The run ends with a status table. Expect `12 finished, 0 failed`, and no "result files not complete"
section. The smoke results land in `results/server/smoke/`. If a job failed, read its
`results/server/logs/smoke.<job>.log` and `.time`.

## 5. Phase 1: prepared meshes and exact truths (~3–4 h)

Use tmux, so the run survives your ssh session.

```bash
tmux new -s rg
cd ~/rg && export EAR_DIR=$PWD/external/weighted_distance_oracle
JOBS=48 scripts/server/run_server.sh prepare truth 2>&1 | tee results/server/run1.log
# detach with Ctrl-b d; reattach later with: tmux attach -t rg
```

When it ends, check it:

```bash
python3 scripts/server/status.py --only truth   # 30 finished, 0 failed
ls results/server/comp | head                   # the competitors' query files (V2V / A2A .query/.truth/.pairs)
```

## 6. Phase 2: everything else (~3–4 days)

```bash
tmux attach -t rg   # or: tmux new -s rg
cd ~/rg && export EAR_DIR=$PWD/external/weighted_distance_oracle
JOBS=48 JOBS_BIG=4 JOBS_PROX=2 SOLVER_MEM_GB=150 scripts/server/run_server.sh \
  ours caps curve ablation solvers weighted scale proximity threads poi competitors small \
  2>&1 | tee results/server/run2.log
```

How to set the knobs:

| knob | what it limits | how to set it |
|---|---|---|
| `JOBS` | single-threaded jobs at once (500k-face jobs: ≤ 4 GB each) | physical cores − 4. Not hyperthreads: the timings are single-threaded |
| `JOBS_BIG` | 2M / 2.5M-face jobs (~40 GB each), EAR-Oracle, DGG, weighted | 4 with 256 GB of RAM; 2 with 128 GB |
| `JOBS_PROX` | proximity jobs at once, 16 threads each | keep `JOBS_PROX` × 16 ≤ cores |
| `SOLVER_MEM_GB` | memory for the label solvers alive at once in the `threads` phase | ~60% of RAM |

The `threads` phase measures 8, 16, 32 and 64 threads; the 64-thread point means something only with
≥ 64 cores. Phases run one after another, so the multithreaded phases never share the machine with
the single-threaded timings.

You can also run phases one at a time, in any order after phase 1. For example:

```bash
JOBS=48 scripts/server/run_server.sh ours caps 2>&1 | tee -a results/server/run2.log
```

## 7. Monitoring

```bash
python3 scripts/server/status.py                # finished / failed / running jobs, wall time, peak RAM
tail -f results/server/run2.log                 # phase starts ("=== <phase> <date>")
cat results/server/logs/failed.txt              # failed jobs, if any
tail -3 results/server/logs/<job>.log           # what one job is doing
free -g; uptime                                 # memory and load
```

## 8. If something fails

- Read `results/server/logs/<job>.log` (its output) and `<job>.time` (its stderr and GNU time, which
  include errors and "Killed" from the OOM killer).
- To rerun, run the phase again. It reruns all of the phase's jobs, but the truth caches
  (`results/server/truth/`) are reused, so no exact sweep is repeated.
- To rerun one job, take its exact command line from the `"argv"` field of its JSON:

  ```bash
  python3 -c "import json; print(json.load(open('results/server/ours.TJ.kd.json'))['argv'])"
  ```

- Running out of memory: lower `JOBS_BIG` (the big jobs) or `JOBS`, then rerun that phase.

## 9. Results (all JSON)

| file | content |
|---|---|
| `results/server/<phase>.<...>.json` | one per job: `cci_exp`, `proximity`, EAR via `harness_eval.py`, `upo_static`, `solver_check` (the `onthefly.*` files), DGG (`dgg.*`, merged with its precompute time and memory). Each file records the command line that made it (`argv`). |
| `results/server/logs/jobs.jsonl` | one line per finished job: exit code, wall seconds, peak RSS (from GNU time) |

`cci_exp` and `proximity` rewrite their JSON after every result, and once right after the index build.
A job that dies or times out keeps what it measured, with `"complete": false`.

Each job also leaves `results/server/logs/<job>.log` (its printed lines) and `<job>.time` (GNU time).

## 10. Copy the results back (from the Mac)

```bash
cd ~/hse/riemannian_geometry
rsync -av SERVER:rg/results/server/ results/server/ \
  --exclude meshes/ --exclude 'dgg/*/*.binary' --exclude 'dgg/*/*.off'
rsync -av SERVER:rg/external/weighted_distance_oracle/results/ results/server/ear_logs/   # EAR's raw logs
```

You can do this while the run is going too: the partial JSON files are valid.

## Expected wall time (single-threaded jobs)

- `truth`: the longest is ~3 h, at 2.5M faces.
- `ours`, `caps`, `curve`, `ablation`, `solvers`, `poi`: ~5–15 h each at 500k faces. ε = 0.05 in `curve`
  is the longest, ~15 h.
- `weighted`: ~35 h for m = 5 single-threaded.
- `scale`: ~2–3 days, the 2.5M-face builds.
- `competitors`: EAR up to 72 h (timeout).

With 48 cores the whole plan takes about 3–4 days.
