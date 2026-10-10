# Running the terrain-oracle competitors

The competitors (SE-Oracle, EAR-Oracle, K-Algo, fixed/unfixed Steiner scheme, CGAL MMP) come from
the EAR-Oracle authors' harness, https://github.com/ItakEjgo/weighted_distance_oracle (datasets included).

## Building on macOS with CGAL 6 (Homebrew `cgal boost gmp mpfr`)

The harness targets Linux/CGAL 5. Patches needed (all mechanical):

1. `#include <bits/stdc++.h>` -- put `stdc++.h.shim` from this folder at `shim/bits/stdc++.h` and add
   `-I<path>/shim` to `CMAKE_CXX_FLAGS`.
2. `chrono::_V2::system_clock` -> `chrono::system_clock` (every header).
3. `physical_memory_used_by_process()` in `base.h` reads `/proc/self/status`; on macOS use
   `task_info(mach_task_self(), MACH_TASK_BASIC_INFO, ...)` and return `resident_size / 1024`
   (`#include <mach/mach.h>`).
4. `getopt.hpp` reads `/proc/<pid>/cmdline`; store `argv` from `main` instead (`store_cmdline(argc, argv)`
   before `Methods::run_new`).
5. CGAL 6: `boost::optional` -> `std::optional`, `boost::get<Point>(&x)` -> `std::get_if<Point>(&x)` in `base.h`.

Then `cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-I$PWD/../shim" && make` in `bld/`.
`harness.patch` in this folder is the complete diff (`git apply` it in a fresh clone at 2343aa9); the clone
used for the results lives in the gitignored `external/weighted_distance_oracle`.

6. **Weighted endpoint fix** (part of `harness.patch`, marked `[patch]`). Every weighted query routine
   (`kSkip::queryGraphA2A` for the fixed scheme, `kSkip::unfixedOnTheFly`, and EAR-Oracle's query in
   `quad.h`) joins the query point to the Steiner points of its face at plain Euclidean length, i.e. at
   weight 1 instead of the face's weight $w_f \in (1, 5)$. Their weighted answers were therefore below
   every weighted path on some queries (on weighted HorseMount, query 72: 212.9 m = the unweighted
   geodesic, while the straight segment through faces of weight ~4 costs 664 m), and as part of a
   best-path reference they dragged it down for every other method. The patch charges these segments
   $w_f \times$ length (`Base::qw(fid)`), which only makes their answers valid path costs. Unweighted runs
   (all weights 1) are unaffected.

## Two things that make or break the comparison

* The harness is **float32** (`CGAL::Simple_cartesian<float>`). On the raw datasets (UTM coordinates
  ~1e7 m) its own MMP disagrees with double-precision MMP by up to 200%. Translate every mesh to
  the origin first; `data/terrains/*.off` are already repaired (`terrain.repair_height_field`) and
  translated -- copy them to `$EAR_DIR/datasets/repaired/`.
* Its AABB point location fails for a query point exactly at a vertex (12 of 100 wrong). Vertex
  queries are therefore written 5 cm inside a face by `bench/terrain_bench.py:write_competitor_queries`.

## Running

```bash
export EAR_DIR=/path/to/weighted_distance_oracle
KIND=V2V scripts/competitors/run_competitors.sh HorseMount results/HorseMount.V2V.query MMP FixedS EAR KAlgo SE
KIND=A2A scripts/competitors/run_competitors.sh HorseMount results/HorseMount.A2A.query MMP FixedS EAR KAlgo
python3 scripts/competitors/competitor_eval.py $EAR_DIR/results HorseMount V2V
python3 scripts/competitors/competitor_eval.py $EAR_DIR/results HorseMount A2A
# weighted terrain (face weights written by bench/terrain_bench.py --weighted SEED):
FACE_WEIGHTS=results/HorseMount_weighted.face_weight.query KIND=W2V scripts/competitors/run_competitors.sh \
    HorseMount results/HorseMount_weighted.W2V.query FixedS EAR UnfixedS@0.05
python3 scripts/competitors/competitor_eval.py $EAR_DIR/results HorseMount W2V
```

Methods take `NAME[@eps[@sp-num]]`, e.g. `EAR@0.2@20`. K-Algo and the unfixed scheme are seconds to
minutes per query; SE-Oracle's build is ~20 CPU-minutes on the smallest mesh and impractical beyond.
Their K-Algo ignores face weights (its weighted answers equal the unweighted exact distances), so it
is excluded from weighted comparisons; on the 100k-vertex GunnisonForest EAR-Oracle and K-Algo ran
(51,512 s and 339 s wall) but exited without writing query results on the 8 GB machine.
