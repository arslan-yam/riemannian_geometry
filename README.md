# riemannian_geometry

Shortest-path distance oracles on 2-D Riemannian manifolds and polyhedral terrains.

A manifold is discretised into a graph, a separator-tree index is built over it, and the answer is
taken off the graph into the continuum, so that it beats the graph's own accuracy floor. The
central object is the **Separator Cut Index** (`rgeo.oracles.SeparatorCutIndex`): portals placed on
the cut lines themselves, labels that are continuum (refined or exact) distances, and sub-sample
interpolation of the crossing point. On real terrains it is benchmarked against the terrain
distance oracles of the literature (SE-Oracle, EAR-Oracle, K-Algo) on their own datasets.

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
  oracles/
    tree_decomposition.py   elimination ordering for H2H
    h2h.py                  H2HIndex (hierarchical 2-hop labels) -- the graph baseline
    separator.py            SeparatorIndex (vertex portals), SeparatorRefinedIndex
    cut_index.py            SeparatorCutIndex -- the contribution
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
bench/
  terrain_bench.py          the terrain benchmark (our methods; writes results/<mesh>.json)
  terrain_report.py         tables, Pareto figures, notes/terrain_results.md
scripts/
  make_terrains.py          regenerate data/terrains/ from the EAR-Oracle repository's datasets
  run_scale.sh              the 100k-vertex experiment (detach with nohup)
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
  error_bound.md            the global error bound of the cut index and what each term implies
  preprocessing.md          the account of the build-time optimisations
  terrain_results.md        all benchmark tables (generated)
tests/test_core.py          invariants: H2H = Dijkstra, upper bounds, exact walk, gradients, MMP labels
```

## Setup

```bash
pip install -r requirements.txt     # numpy, scipy, matplotlib, pygeodesic, pytest, jupyter
python3 -m pytest tests -q          # ~10 s
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

## Terrain benchmark

```bash
python3 bench/terrain_bench.py data/terrains/HorseMount.off --k-list 16,32,64 --a2a 100 --out results
python3 bench/terrain_report.py HorseMount
python3 -c "from bench import terrain_report as r; r.write_note(['HorseMount','BigMount','HeadlightMount','RobinsonMountain','Gunnison10k','GunnisonForest100k'])"
```

Builds are single-threaded by default (`--workers N` forks the per-cut sweeps). `--a2a N` adds N
arbitrary-point queries; `--weighted SEED` makes face costs uniform in (1, 5). The competitors are
run from their own harness: see `scripts/competitors/README.md` (it is float32, so meshes are
translated to the origin and vertex queries nudged 5 cm into a face).
