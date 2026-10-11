# Projection: the full-scale runs on a server

Everything below is extrapolated from runs on an Apple M1 (8 GB, one core per method). The server
assumption is a current x86 machine (e.g. 2 x EPYC / Xeon, 256-512 GB RAM) whose single-core speed is
within +-30% of an M1 performance core; single-thread figures carry that uncertainty, and parallel
figures assume the label sweeps are spread over cores (they are independent, but our build is
single-threaded today — a parallel loop over portals is the one change needed).

## 1. What was measured (M1, single thread)

**Our certified index, C++, hybrid sweeps** (1,500 random pairs per mesh, all sound):

| mesh | n | eps / g | build | index | query | mean / max error |
|---|---:|---|---:|---:|---:|---|
| HorseMount | 793 | 0.2 / 1e | 1.8 s | 1.0 MB | 0.20 us | 0.068% / 1.50% |
| BigMount | 1,440 | 0.2 / 1e | 6.1 s | 2.3 MB | 0.22 us | 0.067% / 1.71% |
| HeadlightMount | 2,470 | 0.2 / 1e | 15.1 s | 4.4 MB | 0.22 us | 0.047% / 1.30% |
| RobinsonMountain | 3,696 | 0.2 / 1e | 25.2 s | 7.3 MB | 0.23 us | 0.063% / 1.79% |
| Gunnison10k | 9,610 | 0.2 / 1e | 132.7 s | 24.0 MB | 0.24 us | 0.023% / 0.68% |
| GunnisonForest100k | 100,511 | 0.2 / 1e | 3,961 s | 533.8 MB | 0.51 us | 0.0069% / 0.12% (200 pairs) |
| Gunnison10k | 9,610 | 0.1 / 0.5e | 201.8 s | 37.4 MB | 0.36 us | 0.0076% / 0.37% |
| Gunnison10k | 9,610 | 0.05 / 0.25e | 418.9 s | 64.2 MB | 0.4 us | 0.0024% / 0.15% |

(exact per-mesh errors for every configuration: `notes/unified_results.md`.)

Fitted growth: build $\propto n^{1.45}$ between 10k and 100k ($n^{1.6}$-$n^{1.7}$ over 0.8k-10k, where
near tables and setup weigh more), index $\propto n^{1.32}$ (the $n^{1.5}$ asymptote is approached slowly:
near tables and tails are linear). The hybrid builds 1.2-1.7x faster than MMP-only sweeps; VTP's
advantage per full sweep grows with the mesh (2.6x at 0.6k vertices, 4.1x at 10k).

**Competitors** (their code, same machine):

| method | measured | growth | blocker at scale |
|---|---|---|---|
| EAR-Oracle (EAR authors' harness, sp=5) | 23.5 s / 189 MB at 793 v; 670 s / 702 MB at 9.6k; at 100k ran 14.3 h and exited without results (8 GB) | build $\sim n^{1.34}$; UP-Oracle paper reports 9.9 h at 100k faces (50k v) | memory; hours per 100k |
| SE-Oracle (same harness) | 1.74 h / 942 MB at 793 v | — | infeasible beyond a few thousand vertices |
| UP-Oracle (its authors' code), every vertex a POI | 15.0 s / 9.9 MB at 576 v; 26.7 s / 27.8 MB at 793 v | complete graph: $n^2$ pairs (with paths) | $n^2$ memory: ~10^10 pairs at 100k |
| WSPD / RC-Oracle (UP-Oracle code), every vertex a POI | 5-13 s, 113-244 MB at 0.6-0.8k v | $n^2$ | same |
| UP / WSPD / RC / EAR (UP-Oracle code), 50 POIs | 0.25-2 s at 0.6-0.8k v | linear in POIs x mesh | answers only POI pairs |
| DGG (authors' code, macOS port), eps = 1e-4 | 0.7 s / 4.4 MB at 793 v; 13.2 s / 73 MB at 9.6k; query 1.4 -> 22 ms | precompute ~$n^{1.2}$, query ~$n^{1.1}$ (Dijkstra) | query time: ~0.3 s at 100k, ~3-5 s at 1.25M |
| exact MMP / VTP on the fly | 9 / 3.3 ms at 793 v; 226 / 55 ms at 9.6k; MMP 3.3 s at 100k | ~$n^{1.3}$ per query | seconds per query beyond 100k |
| NeuroGF (authors' network, distance branch, 3,000 iterations on the M1 GPU) | 23 min / 0.25 MB at 793 v: 0.46% mean, 4.7% max; 23 min at 3.7k v: 0.67% mean, 10.3% max; query 280 µs (17 µs batched) | training budget fixed; data = S exact sweeps | fixed capacity: error grows with the mesh; no bound |

## 2. Projected costs of our index (single thread, eps = 0.2, g = 1 edge)

With build $\propto n^{1.45\ldots1.6}$ and size $\propto n^{1.32}$ anchored at the 100k run:

| terrain size (vertices) | where it occurs | build, 1 core | build, 64 cores (eff. 0.6) | index |
|---:|---|---:|---:|---:|
| 100k | GunnisonForest, Laramie (EAR/SE data) | 1.1 h (measured) | ~2 min | 0.53 GB (measured) |
| 250k | UP-Oracle data, 500k faces | 4.2-5.0 h | 7-8 min | 1.8 GB |
| 500k | UP-Oracle, 1M faces | 11-15 h | 18-25 min | 4.4 GB |
| 1.25M | UP-Oracle, 2.5M faces | 43-65 h | 1.1-1.7 h | 15 GB |

eps = 0.1 / g = 0.5e costs about 1.5x the build and 1.5x the index; eps = 0.05 about 3x / 2.7x.
Queries stay at 0.2-0.6 us (cache misses grow with the index). The scale-dependent gap of `notes/research.md` §5.2 would make the 1.25M index only about 1.5x
smaller with its current constants (estimate, not measured); it pays off beyond ~1M vertices.

## 3. Projected costs of the competitors at those sizes

* **EAR-Oracle**: ~4-10 h at 100k on a server (UP-Oracle paper: 9.9 h at 100k faces), tens of hours
  at 250k-1.25M; index several GB; queries 1-5 ms; mean error 2-3% (unweighted, our measurements).
* **SE-Oracle**: not feasible at these sizes.
* **UP-Oracle, WSPD, RC-Oracle** in the all-vertex setting: not feasible beyond ~5k vertices ($n^2$
  complete graph). In their own POI setting (50-2500 POIs) they are cheap to build — one exact sweep
  per POI (at 250k vertices ~10 s per sweep: 50 POIs ~8 min, 2500 POIs ~7 h) — but answer POI pairs only;
  their published "A2A" numbers are POI numbers multiplied by vertices/POIs, not runs.
* **DGG** (eps = 1e-4): precompute ~5-20 min at 250k and ~1-2 h at 1.25M, graph ~2-10 GB, queries
  0.5-5 s (one Dijkstra); accuracy as measured (0.003% mean, 0.02% max), no per-query bound.
* **On-the-fly exact (MMP / VTP)**: seconds per query at 250k+ (VTP ~4x faster than MMP).
* **NeuroGF**: training on a GPU server with their CUDA code (their full 400 x 200 schedule) is
  about an hour per mesh plus the exact training sweeps (512 sources: ~1-6 h at 250k-1.25M with VTP).

## 4. Proposed server run list (single-threaded numbers, 1 core each, independent jobs)

| job | estimate |
|---|---:|
| our index, 5 small meshes x 3 configs + 2 x 100k x {0.2, 0.1} | ~6 h |
| our index, 6 UP-Oracle terrains at 500k faces, eps 0.2 (eps 0.1 on two of them) | ~30-40 h (or ~1 h on 64 cores) |
| our index, one terrain at 1M / 1.5M / 2M / 2.5M faces (scaling curve) | ~100-150 h (or 2-4 h on 64 cores) |
| EAR-Oracle on 100k (2 meshes) and 500k faces (2 terrains), sp=5 | ~40-80 h, >32 GB RAM |
| UP-Oracle code, their POI settings on the 500k-face terrains (as in their paper) | ~10-20 h |
| DGG eps 1e-3 / 1e-4 on all terrains up to 2.5M faces | ~10 h |
| NeuroGF on 2-3 terrains (GPU) | ~3-6 h |
| ground truth: 1,000 exact pairs per terrain (MMP early stop / VTP) | ~5-20 h |

Run jobs in parallel across cores (each job single-threaded, as here), so the whole plan fits in
about 1-2 days on a 64-core server; memory per job stays below 32 GB except EAR at 500k+.

## 5. Risks

* **MMP worst cases.** Window counts grow superlinearly on very steep or noisy terrains (Gunnison10k:
  max slope 14.9); the exponent may exceed 1.6 on the 2.5M terrains. VTP-heavy hybrid settings and
  parallel sweeps mitigate it.
* **Memory.** The index is held in RAM during the build (15 GB projected at 1.25M, eps 0.2).
* **Mesh quality.** VTP needs consistently oriented, non-degenerate faces (sub-meshes are oriented
  now; folded height fields must be repaired first, as `repair_height_field` does).
* **Their harness.** EAR's float32 code needs translated meshes, nudged query points and the weighted
  endpoint fix (`scripts/competitors/README.md`); at 100k it needed more than 8 GB.
