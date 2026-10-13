# riemannian_geometry

Shortest-path distance oracles on 2-D Riemannian manifolds and polyhedral terrains.

The central object is the **Certified Cut Index** (`rgeo.oracles.certified.CertifiedCutIndex`): a
separator tree of axis-parallel cut lines, portals placed on the cut lines by arc length, exact
geodesic labels (MMP) from every portal to its region, and exact near tables for close pairs. Every
query returns an interval `LB <= d(s,t) <= UB` with `(1 - eps) d <= LB` and `UB <= (1 + eps) d`
(theorem and proofs: `notes/guarantee.md`), plus an interpolated estimate inside it. On weighted
terrains the labels come from a Steiner base graph instead (`rgeo.terrain.steiner`). It is
benchmarked against the terrain distance oracles of the literature (SE-Oracle, EAR-Oracle, K-Algo,
fixed/unfixed Steiner schemes) on their own datasets: `notes/terrain_results.md`.

The earlier **Separator Cut Index** (`rgeo.oracles.SeparatorCutIndex`) is the general-manifold
version: graph-based, continuum labels by refinement, no certificate. It works on any `Manifold`
(Poincare disk, analytic terrains, obstacles) and is kept with its gaps fixed.

## Layout

```
rgeo/                       the library (pip install -e .  or just run from the repo root)
  manifolds.py              Manifold base (metric, lengths, energy + analytic gradient), Euclidean,
                            PoincareDisk, Terrain, GaussianHills
  obstacles.py              Circle, Rectangle, `blocked`; batched tests and signed clearance
  refine.py                 shortcut / subdivide / energy relaxation of graph paths; batched L-BFGS
  visualize.py              plotting helpers
  graphs/
    build.py                stencil graphs (grid / triangular), Riemannian edge weights, obstacles
    sampling.py             difficulty field and adaptive (quadtree) sampling
    dijkstra.py             Dijkstra with predecessor tree
    _cdijkstra.c, cdijkstra.py  multi-seed Dijkstra that stops at its last target (C, scipy fallback)
  oracles/
    tree_decomposition.py   elimination ordering for H2H
    h2h.py                  H2HIndex (hierarchical 2-hop labels) -- the graph baseline
    separator.py            SeparatorIndex (vertex portals), SeparatorRefinedIndex
    cut_index.py            SeparatorCutIndex -- the general-manifold cut index
    certified.py            CertifiedCutIndex -- the certified (1+eps) terrain oracle
    _cquery.c, cquery.py    its compiled query (C via ctypes, built on first use; same answers)
    variants.py             heat / eikonal / continuous / BFS variants (recorded negative results)
  solvers/
    eikonal.py              fast marching (anisotropic eikonal)
    heat_method.py          heat-method distances
  terrain/
    io.py                   OFF read/write, repair of folded height fields
    mesh.py                 TerrainMesh: terrain as a manifold, exact face-walk lengths + gradient,
                            optional face weights
    exact.py                ExactTerrain: exact labels by MMP sweeps (pygeodesic), Steiner insertion
    a2a.py                  TerrainA2A: arbitrary-point queries on top of the vertex index
    submesh.py              face grid, sub-mesh extraction, batch Steiner insertion, early-stopping MMP
    steiner.py              SteinerGraph (m points per edge) and SteinerLabels for weighted terrains
bench/
  terrain_bench.py          the terrain benchmark (our methods; writes results/<mesh>.json)
  terrain_report.py         tables, Pareto figures, notes/terrain_results.md
scripts/
  make_terrains.py          regenerate data/terrains/ from the EAR-Oracle repository's datasets
  run_certified.sh          the final runs of the certified index (5 meshes, weighted, 100k)
  run_scale.sh              the 100k-vertex experiment of the earlier index (detach with nohup)
  run_weighted.sh           the weighted-terrain experiment
  competitors/              build + run the competitors' C++ harness on macOS, score its logs
notebooks/
  exps.ipynb                every experiment of the project, in order (manifolds -> terrains)
  experiments_july.ipynb    the first notebook (H2H, sampling, first separator index)
data/terrains/              the benchmark meshes, repaired and translated (large/ is regenerated)
results/                    benchmark outputs: <mesh>.json ours, <mesh>_competitors[_a2a|_w2v].json
                            theirs, <mesh>.{V2V,A2A,W2V}.query the shared query sets
notes/
  code_walkthrough.md       formal walk-through of every module
  guarantee.md              the (1+eps) theorem of the certified index, with proofs and complexity
  research.md               literature survey, what was adopted / dropped, future work
  error_bound.md            the error bound of the earlier cut index and what each term implies
  preprocessing.md          the account of the build-time optimisations
  terrain_results.md        all benchmark tables (generated)
  unified_results.md        every method in one table per mesh (generated)
  size_results.md           size variants: coding, caps, shifted hierarchy, covers, size vs error (generated)
  poi_results.md            POI setting: our POI index against the exact POI table (generated)
  related_work.md           literature map and positioning (incl. Djidjev-Sommer, Kaul et al. bounds)
  projection.md             full-scale server projection; run list in scripts/server/run_server.sh
  paper_draft.md            outline of the paper with the measured numbers
  hyperbolic.md, hyperbolic_h2h.md   side project: hyperbolic domains with obstacles (sketch + H2H check)
tests/test_core.py          invariants: H2H = Dijkstra, upper bounds, exact walk, gradients, MMP labels,
                            exact obstacle test, face grid, Steiner insertion, certified bounds,
                            compiled query = Python query, Steiner labels admissible
cpp/                        C++ implementation of the certified index (`make` -> cpp/bin/):
  src/certified.{h,cpp}     the index (same construction and query as certified.py)
  src/solver_{mmp,vtp}.cpp  exact solvers: Kirsanov's MMP (early stop) and VTP (Qin et al. 2016)
  src/terrain.h, graph.h    mesh, face / vertex grids, face walk, mesh graph + Dijkstra
  tools/cci_bench.cpp       build + evaluate on vertex pairs against exact MMP
  tools/cci_exp.cpp         size experiments: storages, caps, shifted hierarchy, POI mode / POI table,
                            --threads N (bit-identical labels), --store code8 (coded during the build)
  tools/solver_check.cpp    MMP vs VTP: exactness and speed
  tools/upo_static.cpp      the UP-Oracle authors' oracles behind a static many-query driver
  tools/dgg_bench.cpp       DGG (discrete geodesic graph) queries on a precomputed graph
  third_party/              Kirsanov MMP (MIT, via pygeodesic), VTP (MIT)
external/                   (gitignored) competitors' code: EAR harness, UP-Oracle, VTP, DGG (+ macOS
                            port in dgg_port, dgg_shim), NeuroGF, pygeodesic sources
```

## Setup

```bash
pip install -r requirements.txt     # numpy, scipy, matplotlib, pygeodesic, pytest, jupyter
python3 -m pytest tests -q          # ~20 s
jupyter lab notebooks/exps.ipynb
```

The notebook and the scripts put the repository root on `sys.path` themselves, so installing the
package (`pip install -e .`) is optional.

## Ground truth

On `PoincareDisk` use `PoincareDisk.exact_dist`. On terrains use exact polyhedral geodesics
(MMP, `pygeodesic`); arbitrary surface points are inserted as Steiner vertices first.

Everywhere else, **do not use Dijkstra on a dense graph as ground truth**: graph distance converges
to the stencil's Finsler metric, not the Riemannian one, so its error plateaus at
`sec(theta_max/2) - 1` at any resolution (1.48-1.59% on `GaussianHills` for 9k, 21k and 37k nodes).
Build the reference as a dense graph for the *route* plus `rgeo.refine.refine_polyline` for the
*length*. Weighted terrains have no exact solver; there the reference is the best path any method
found.

## Certified index in five lines

```python
from rgeo.terrain import read_off, TerrainMesh
from rgeo.oracles.certified import CertifiedCutIndex
from rgeo.oracles.cquery import CompiledQuery
tm = TerrainMesh(*read_off("data/terrains/HorseMount.off"))
ix = CertifiedCutIndex(tm, eps=0.1, gap_factor=0.5)     # 4.7 s, 1.3 MB
lb, ub, est = ix.bounds(10, 700)                         # LB <= d <= UB, (1-eps) d <= LB, UB <= (1+eps) d
cq = CompiledQuery(ix); cq.bounds(10, 700)               # same answer in ~1 us
```

## C++ build

```bash
make -C cpp                                   # cpp/bin/cci_bench, solver_check
cpp/bin/cci_bench data/terrains/HorseMount.off --eps 0.1 --gap-factor 0.5 --solver hybrid --random 1000
make -C cpp bin/upo_static                    # needs external/UpdatedStructureTerrainCode
scripts/run_cpp.sh                            # all C++ runs (detach with nohup); scripts/run_baselines.sh: DGG, NeuroGF
scripts/run_size.sh; scripts/run_size2.sh     # size variants and size-vs-error; python3 bench/size_report.py
scripts/run_poi.sh                            # POI setting; python3 bench/poi_report.py
JOBS=48 scripts/server/run_server.sh          # the full-scale runs on a Linux server (notes/projection.md)
python3 bench/unified_report.py               # notes/unified_results.md: every method in one table per mesh
```

## Terrain benchmark

```bash
python3 bench/terrain_bench.py data/terrains/HorseMount.off --methods none --cert 0.2/1,0.1/0.5,0.05/0.25 --a2a 100 --out results
python3 bench/terrain_bench.py data/terrains/HorseMount.off --weighted 1 --cert 0.1/0.5 --steiner-m 5,10 --steiner-ref 24
python3 bench/terrain_bench.py data/terrains/HorseMount.off --k-list 16,32,64 --a2a 100 --out results   # the earlier index
python3 bench/terrain_report.py HorseMount
python3 -c "from bench import terrain_report as r; r.write_note(['HorseMount','BigMount','HeadlightMount','RobinsonMountain','Gunnison10k','GunnisonForest100k'])"
```

`--cert eps/gap_factor,...` builds the certified index (gap in median edge lengths) and checks the
certificate on every query. Builds are single-threaded (`--workers N` forks the per-cut sweeps of the
earlier index only). `--a2a N` adds N
arbitrary-point queries; `--weighted SEED` makes face costs uniform in (1, 5). The competitors are
run from their own harness: see `scripts/competitors/README.md` (it is float32, so meshes are
translated to the origin and vertex queries nudged 5 cm into a face).
