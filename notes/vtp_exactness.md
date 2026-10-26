# VTP returned too-long distances: two defects in the upstream code, and the fix

The cluster evaluation (`notes/server_results.md`, rows "MMP vs VTP") found VTP longer than MMP at a few
vertices on many meshes. A Steiner-graph path is a real surface path, so it bounds the true distance from
above. On plain TJ at 1M faces, source 271250, vertex 273379: MMP 25.2980, Steiner (m = 8) 25.4464, VTP
26.1581. From that source 36,990 vertices are off by more than 1e-9 relative (9,994 by more than 1e-4). On
TJ with six walls (0.5M), source 107648, vertex 115079: MMP 440.47, Steiner 442.53, VTP 446.26. VTP was
never shorter than MMP in any comparison below.

Both causes are in the upstream VTP code. Both make VTP drop information it needs, so it returns real path
lengths that are not the shortest. The fix is three local changes to
`cpp/third_party/vtp/geodesic_algorithm_exact.h`. After it, MMP and VTP agree to 1e-14 relative on every
full sweep checked (105 sources on nine meshes), and a VTP full sweep takes 1.1–1.4× as long as before (still 5–7× faster than MMP).

## Not the wrapper

`cpp/third_party/vtp` is byte-identical to the authors' release (`external/VTP_source_code`, i.e.
github.com/YipengQin/VTP_source_code, whose code has not changed since 2016), except `exit(1)` → `throw`
in the mesh checks. `solver_vtp.cpp` was ruled out piece by piece:
- **Per-sweep reset.** Every sweep uses a fresh algorithm object (window lists, queues) and resets each
  vertex's distance, state and incident face. A fresh `VTPSolver` gives the same values.
- **`split_parts`.** A connected mesh is one piece; the error is the same with and without splitting.
- **`insert_points`, `close_fans`.** Used only for surface-point queries and `exact = vtp` sub-meshes. The
  vertex sweeps above use neither.
- **Source setup, float/double.** The source is a mesh vertex and VTP is double throughout. Labels are
  rounded up to float32 later (`f32_up`), separately from this.

## Defect 1: the vertex queue loses vertices at equal distances

VTP keeps its vertex priority queue as `std::set<vertex_pointer, Vertex>`. `Vertex::operator()` compares
geodesic distances only, so two vertices at the same distance are one key:
- `insert` of the second is silently dropped;
- `erase(v)` in `create_pseudo_source_windows` and `find(v)` in `erase_from_queue` hit whichever queued
  vertex has v's distance, removing the wrong one.

A vertex lost this way keeps its distance and its FRONT state but is never popped. It never becomes
INSIDE, step (3) of `propagate` never makes it a pseudo-source, and step (4) never pushes the window lists
around it. On a terrain, every non-convex vertex and every flat interior vertex is flagged
`saddle_or_boundary` (total angle > 2π − 1e-5), and so is every hole or sub-mesh boundary vertex. When such
a vertex is lost, the paths that must bend around it are missing, and everything in its shadow gets a
too-long distance. That is the 3% error at vertex 273379.

Ties are frequent on our data: TJ, SC, VS, AU, GI and LH are regular grids (10 m at 0.5M faces, 7 m at 1M)
with integer heights, so edge lengths, and sums of them, repeat exactly. Instrumented, the sweep from
source 271250 (TJ 1M) dropped 537 inserts, erased the wrong vertex 3 times, and ended with 11 vertices that
had a distance but were never popped, 10 of them saddle or boundary vertices. GunnisonForest100k and
Gunnison10k have irregular float coordinates and almost no ties, which is why earlier checks there agreed
to 1e-15. Kirsanov's MMP orders its interval queue by (distance, start, edge id) and is immune.

A second, rarer queue defect sits in the same place. `create_pseudo_source_windows` re-queues an improved
neighbour even when it is INSIDE. `find_separating_point` then changed that queued vertex's distance (the
set's key) without taking it out first, because it erased and re-inserted FRONT vertices only. This
happened once in three instrumented sweeps.

**Fix.** The queue is ordered by (distance, vertex id), and `find_separating_point` erases the top vertex
before changing its distance, whatever its state, re-inserting it if it was queued or FRONT. With that,
no vertex with a distance is left unpopped (checked on three sweeps).

## Defect 2: Rule 2 prunes windows at crossings that need not exist on the surface

With defect 1 fixed, errors of 1e-8 to 8e-6 relative remained on about half the sources (e.g. TJ 0.5M,
source 9668: 26 vertices, worst 7.7e-6; source 48248: 247 vertices). Turning VTP's rules off one at a time
traced them all to Rule 2, `pairwise_windows_checking` / `check_between_two_windows`.

Rule 2 compares two windows on the same edge where their rays cross below the edge, in the list's
unfolded frame. Then it drops one window, or cuts the end off one (cases 1–4: overlapping windows; cases
5–6: disjoint ones). The argument is the usual one: where two paths cross, the shorter one wins beyond.

But the crossing is on the surface only if both rays pass through it in the same triangle. Each window
was unfolded through its own chain of triangles. The only triangle all windows of a list are known to
share is the face they came from, just below the edge. Rule 2's checks are weaker than that:
- cases 1–4 check only that the crossing is between the pseudo-sources' depth and the edge;
- cases 5–6 check that the innermost crossing lies in that face, but a dropped window needs every crossing
  along its extent to be real, and they go deeper.

(Case 6 also has a typo, `d1 < d2 - (1 - 1e-12)` for `d1 < d2 * (1 - 1e-12)`. It only makes the case
weaker.)

Experiments on nine sources with residual errors (TJ 0.5M plain and walls):
- skipping cases 5–6, or 1–4, removes some of the residuals but not all, so every case contributes;
- requiring the crossing that each decision rests on (the outermost one along the dropped or cut part) to
  lie inside the face below removes all of them, but the sweep becomes 1.8–2.3× slower;
- not applying Rule 2 at all also removes all of them, and is faster than the checked version.

Most windows come from far away, so most crossings lie many triangles deep, and the checked rule prunes
almost nothing while still paying for the checks.

**Fix.** `propagate_one_windows_list` no longer calls `pairwise_windows_checking`. The functions stay in
the file, unused. Rule 2's other half, `check_with_vertices`, is sound and kept: it drops a window that is
longer everywhere than the path through one of the edge's vertices, and that path is real.

## Rule 1: the same gap, left as it is

Rule 1 ("one angle two sides", `propagate_windows_to_two_edges`) sends a window to one side only. For a
window entirely on one side of the list's separating point, that is sound: the competing crossing is with
the best path to the top vertex inside the face being entered.

For a window that straddles the separating point, VTP goes one-sided when the window's outer ray crosses
the list's best ray *below* the edge. That crossing has the same gap as Rule 2. Two sound replacements were
tried:
- checking the crossing against the face below, else both sides;
- splitting the window at the separating point.

Both are exact but 25–40× slower (TJ walls 0.5M, source 173982: 60 s and 40 s against 1.7 s): the
shortcut is what keeps VTP fast. No error was traced to Rule 1 in
any sweep here, so it is unchanged. The exactness claim for VTP therefore rests on the checks below, not on
a proof.

## Evidence

All full sweeps, MMP (Kirsanov, full sweep) against VTP before and after the fix, interleaved in one
process (`probe2`, sources drawn as in `solver_check`, `std::mt19937(1)`). Relative difference |VTP − MMP| / MMP.

| mesh | sources | before: sources wrong / labels > 1e-9 (> 1e-4) / worst rel | after: labels > 1e-9 / worst rel | MMP | VTP before | VTP after | after / before |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TJ 0.5M | 20 | 14 / 32,585 (30,378) / 3.2e-02 | 0 / 9.1e-15 | 17.73 s | 2.21 s | 3.08 s | 1.39 |
| TJ 1M | 20 | 17 / 79,603 (11,436) / 3.4e-02 | 0 / 4.0e-15 | 24.51 s | 3.31 s | 4.66 s | 1.41 |
| TJ walls 0.5M | 20 | 19 / 13,194 (8,396) / 1.3e-02 | 0 / 2.8e-15 | 14.04 s | 1.73 s | 2.14 s | 1.24 |
| SC 0.5M | 5 | 4 / 310 (0) / 2.7e-06 | 0 / 2.5e-15 | 14.79 s | 1.70 s | 2.26 s | 1.33 |
| VS 0.5M | 5 | 4 / 41 (0) / 1.7e-06 | 0 / 3.1e-15 | 15.39 s | 1.66 s | 1.83 s | 1.10 |
| AU 0.5M | 5 | 5 / 25,454 (13,563) / 4.5e-02 | 0 / 2.7e-15 | 15.01 s | 1.77 s | 2.46 s | 1.39 |
| GI 0.5M | 5 | 2 / 2 (0) / 6.3e-07 | 0 / 2.4e-15 | 14.71 s | 1.74 s | 2.07 s | 1.19 |
| LH 0.5M | 5 | 4 / 2,289 (1) / 1.4e-04 | 0 / 2.6e-15 | 14.51 s | 1.69 s | 2.27 s | 1.34 |
| GunnisonForest100k (float coords) | 5 | 0 / 0 (0) / 2.9e-15 | 0 / 2.9e-15 | 2.11 s | 0.48 s | 0.54 s | 1.12 |

Times are the mean per full sweep on an M1, one thread. Peak memory of a VTP sweep is unchanged (about
380 MB on TJ 0.5M before and after).

**Labels, end to end.** A 20k-face crop of TJ 0.5M (x 2500–3500, y 3500–4500), index built three times,
eps 0.1, one thread: every cut by VTP before the fix, by VTP after it, and by MMP.
- VTP before vs MMP: 5,149 of 9,221,583 float32 labels longer than MMP's, worst by 18%, none shorter.
- VTP after vs MMP: all labels equal, except 8 that differ by one float32 step (rounding up an equal
  double).
- With the default `solver hybrid`, VTP serves 7 of the crop's 1,023 cuts. Before the fix, 1,487 labels
  were longer, worst by 11%.

**`solver_check`** (Makefile build, TJ walls, 5 sources): max |MMP − VTP| 5.5e-12 (relative 2.8e-15),
VTP 6.3× faster than MMP per full sweep.


## Consequences

- All six cluster terrains are integer grids. The labels of every VTP cut (`solver hybrid`, and every cut
  with `exact vtp`) and the `exact vtp` ground truth of the obstacle runs used the defective VTP. The
  labels were real path lengths, so upper bounds built from them stayed valid. Lower bounds were not
  guaranteed, although all 5.46M checked answers were correct. The main tables need the rerun planned
  with the fixed VTP.
- VTP's speedup over MMP per full sweep drops by 1.1–1.4× (`notes/paper_draft.md` quotes 2.6–4.1×).
- Two existing notes say "VTP stays exact" or "MMP and VTP agree on every obstacle mesh" (pitfall 1 in
  `notes/guarantee.md`, `notes/paper_draft.md`). Those checks were on Gunnison10k, with float coordinates.
  On the integer grids they held for MMP only.

## Reproduce

```bash
python3 scripts/competitors/prepare_mesh.py external/UpdatedStructureTerrainCode/input/TJpre_500000.off /tmp/TJpre_500000.off
make -C cpp BIN=bin_dev bin_dev/solver_check
cpp/bin_dev/solver_check /tmp/TJpre_500000.off 20
```

`solver_check` prints the maximum |MMP − VTP| (relative) over its sources. Before the fix it was 3.2e-2
here; after it, it is at round-off.
