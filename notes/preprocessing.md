# Preprocessing cost of the Separator Cut Index: what was done and what it bought

Reference point: the numbers of the final presentation (slides since removed), Poincaré disk,
16 x 16 grid (188 nodes), `k=32`, `leaf_size=16`, `rounds=1`, `segments=12`, one process:
**114 s** (k=32) and **222 s** (k=64). The build was one Dijkstra plus one `refine_polyline`
per (portal, member) pair, everything scalar Python.

## Exact changes (results bit-identical to the reference labels)

| change | where | effect |
|---|---|---|
| batched free-space test: one `in_domain_batch` + `blocked_batch` over all samples of a segment/polyline instead of a Python loop per sample | `rgeo/refine.py` (`free_segment`, `free_polyline`, `_free_points`); `rgeo/obstacles.py` (`blocked_batch`) | most of the sequential gain |
| analytic metric gradient for the energy (`metric_grad_batch`) and one fused `energy_and_grad_points` evaluation per L-BFGS step instead of finite differences / separate calls | `rgeo/manifolds.py` | 2-4x fewer metric evaluations per step |
| exact `clearance_batch` (signed distance + gradient) for the obstacle barrier | `rgeo/obstacles.py`, `rgeo/refine.py` (`barrier_and_grad`) | removes the per-point Python barrier loop |
| Dijkstra trees cached per anchor node (`self.cache`), dropped after the build | `rgeo/oracles/separator.py` (`SeparatorIndex.dijkstra`) | no repeated Dijkstra when several portals share an anchor |
| batched segment lengths (`segment_lengths`) everywhere a polyline is measured or subdivided | `rgeo/refine.py` | |
| `shortcut`: all candidate chords of the window measured in one batched call, the free-space test batched over all *shorter* chords, same greedy decisions | `rgeo/refine.py` (`shortcut`) | 8x on GaussianHills+obstacle, 2x on Poincaré (per polyline) |

Sequential build after these: **60 s** at k=32 (was 114 s), same labels.

## Parallelism

`SeparatorCutIndex(workers=N)`: (cut, portal) tasks over a forked pool (the index is inherited
by fork, only small result dicts come back), largest regions first; leaf tables fanned out per
leaf the same way. Tasks are per portal, not per cut, because the root cut alone is 75-81% of
the work. Measured 3.0x with 3 workers on the 4P+4E M1 (labels 25.9 s -> 12.5 s), i.e. near
linear on the performance cores.

## Batched refinement (`batched=True`, default)

`refine_batch` resamples every member's polyline to a portal to the same segment count (the
sequential count rounded up to a multiple of 4, so never coarser), stacks them into one
`(B, m+1, 2)` array, and runs ONE L-BFGS over the stack: the energy is separable so the gradient
just concatenates, and each numpy call now does B problems. Sequential: **27 s** at k=32 (was
60 s), accuracy unchanged (0.1312% -> 0.1308%, interpolated 0.107%).

## Summary on the Poincaré benchmark (quiet machine, k=32)

| pipeline | build |
|---|---|
| deck (before) | 114 s |
| exact optimisations, sequential | 60 s |
| + batched L-BFGS, sequential | 27 s |
| + 3 workers | **13 s** |

## Terrains: exact labels instead of refined ones

On a polyhedral terrain every energy evaluation is an exact face walk (`TerrainMesh`), so a
refined label costs ~50-140 ms instead of ~5 ms and the same index on HorseMount (793 vertices,
2110 portals, single thread) takes **1527 s** for 0.085% mean error. The terrain instantiation
replaces refinement by one exact MMP sweep per portal (`ExactTerrain`, portal inserted as a
Steiner vertex, sweep restricted to the sub-mesh its geodesics can visit): **24.6 s** single-threaded
for 0.029% mean error, and 45 s at k=64 for 0.001%. The restriction to a sub-mesh per cut is what
keeps the build near O(n log n) instead of O(n^2): with one augmented mesh for all portals, k=64 on
HorseMount took 16x longer. Single-threaded builds on the other meshes (k=32 uniform): BigMount
94 s, HeadlightMount 256 s, RobinsonMountain 304 s, Gunnison10k (9,610 vertices, 34,674 portals)
1991 s; the cheap variant (`k~|cut|`, k_min=4) is 3-4x less. `workers=N` forks the per-cut sweeps
(~3x on 4 performance cores) but all reported numbers are single-threaded.
