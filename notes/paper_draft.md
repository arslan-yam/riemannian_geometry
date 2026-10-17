# Paper draft: Certified Cut Index

Working draft for a SIGMOD/VLDB submission. Numbers are from the Apple M1 runs, single thread
(`notes/unified_results.md`, `notes/size_results.md`). Items marked **[server]** need the 1M+ runs of
`notes/projection.md`, and items marked **[check]** need verification. The proofs are in
`notes/guarantee.md` and the literature map in `notes/related_work.md`.

## Title options

1. *Certified Distance Queries on Terrains in Microseconds*
2. *H2H for Surfaces: Certified Geodesic Distance Oracles for Terrains*
3. *Every Answer with a Proof: A Certified (1+ε) Distance Oracle for Terrain Surfaces*

## Abstract (draft)

Shortest surface distances on terrains underlie kNN, range and route queries. Exact geodesic
algorithms take milliseconds to seconds per pair, and the oracles proposed in the database literature
trade this for indexes whose error is fixed by a Steiner discretisation: 2–3.5% mean error, hundreds of
MB already at 10k vertices, and no indication of how wrong a particular answer is.

We present the **Certified Cut Index**, a separator-based distance oracle that does for continuous
surfaces what H2H does for road networks. Each vertex stores exact geodesic distances to sampled portals
on the cut lines of a kd-hierarchy, and each query returns an estimate together with a certified
interval:
$$\mathrm{LB} \le d \le \mathrm{UB}, \qquad (1-\varepsilon)d \le \mathrm{LB},\quad \mathrm{UB} \le (1+\varepsilon)d.$$

Four ingredients make this possible:
- a window lemma that lets portals be sampled with a $(1+\varepsilon)$ guarantee;
- exact labels computed on provably sufficient sub-meshes, with a hybrid of MMP and VTP;
- interval bounds from the 1-Lipschitz structure of the label rows;
- closed-loop 8/16-bit label coding that keeps the bounds sound.

On terrains from the SE/EAR benchmarks and a 100k-vertex DEM, the index answers in 0.2–0.5 µs:
- against EAR-Oracle, 4–5 orders of magnitude faster, with 100–400× lower mean error and an index
  20–50× smaller;
- it builds 100k vertices in under an hour on one core, where EAR-Oracle did not finish;
- every answer is certified.

**[server]** Scaling to 1M+ vertices.

## 1. Introduction

**Problem.** Given a terrain (TIN) with $n$ vertices, answer distance queries between vertices and
between arbitrary surface points.

**Why it is hard.** Exact single-source geodesics (MMP, ICH, VTP) are super-linear per query. A
distance oracle must store something like $n^2$ information or approximate it.

**State of the art.**
- *SE / EAR / RC / UP oracles.* Steiner graphs plus well-separated pair or spanner structures. Their
  guarantee is relative to the Steiner graph, and the measured error is 0.6–3.5% mean and 6–25% max.
- *DGG.* Accurate, but queries take milliseconds (Dijkstra).
- *Learned fields (NeuroGF).* No bound.
- *Theory oracles on polyhedral surfaces* (Aleksandrov et al. 2010; Djidjev–Sommer 2011). Near-linear
  space in theory, no implementation.

**Our approach.** The index is H2H's architecture moved to the continuous surface:
- cuts instead of vertex separators;
- sampled portals instead of whole bags;
- exact geodesic labels computed on small sub-meshes;
- an interval certificate on every answer.

**Contributions.**
1. The Certified Cut Index and its $(1+\varepsilon)$ theorem with per-query certificates. Proofs:
   window lemma, ellipse lemma, 1-Lipschitz interval bound.
2. Certify-or-refine, giving a capped variant with $O(n k \log n)$ labels whose answers are
   $(1+\varepsilon)$ or exact.
3. Engineering that makes exact labels cheap:
   - ellipse sub-meshes;
   - hybrid MMP/VTP sweeps, with VTP's state reset per sweep and its orientation requirements met;
   - parallel sweeps;
   - closed-loop coded labels.
4. Certified A2A queries.
5. An evaluation against every implemented competitor, run with their own code on one machine, with two
   fixes to their harness that change its published weighted errors.

## 2. Preliminaries

Terrain, chart, weights, $\lambda = w_{\min}$, the lower bounds (1) of `guarantee.md` §1. Queries: V2V,
A2A. Error measures: mean and max relative error, max UB error, interval width.

## 3. The index

### 3.1 Separator hierarchy and windows

Median kd cuts with split values that avoid every vertex coordinate. Hub $h$, radius $\mu_C$, and the
window $W_C$ (all crossings of $\ell_C$ by a shortest $s$–$t$ path lie in it). Runs are the exact pieces
of the line on the surface.

### 3.2 Portals

- Arc gap $g$ inside the box.
- Geometric tails $\max(g, 2\varepsilon\lambda c)$ outside: the window lemma says a crossing at offset $c$
  costs $d \ge 2\lambda c$.
- Run-end stretches.

### 3.3 Labels on ellipse sub-meshes

- Ellipse lemma: a path of length $\le D$ stays within $D/(2\lambda)$ of the chord.
- Sub-mesh = faces meeting the dilated box; portals are inserted as vertices.
- One exact sweep per portal: MMP with early stop when the sub-mesh is large relative to the members'
  box, VTP otherwise.

### 3.4 Near tables

Exact distances for all pairs within $R = g/(\varepsilon\lambda)$.

### 3.5 Query and certificate

$f_j = \ell(s,p_j) + \ell(t,p_j)$.
- **UB** $= \min f$ (a path cost).
- **LB** $= \min_j \frac12(f_j + f_{j+1}) - A_j$ (with end stretches), together with the chord bound.
- **Estimate:** a model interpolation, clamped into $[\mathrm{LB}, \mathrm{UB}]$.

### 3.6 Theorem

**Theorem 1.** For every vertex pair, $(1-\varepsilon)d \le \mathrm{LB} \le d \le \mathrm{UB} \le (1+\varepsilon)d$.

Proof sketch:
- the window lemma gives that some crossing lies in $W_C$;
- every gap is $\le \varepsilon d$, because non-near pairs have $d > g/\varepsilon$ and tails grow with the offset;
- the bracket of the crossing point then bounds UB;
- LB follows from 1-Lipschitz rows.

### 3.7 Capped portals and certify-or-refine

**Theorem A.** Every answer is $(1+\varepsilon)$-certified or exact. Each cut has at most $k$ box portals plus
$O(\varepsilon^{-1}\log k)$ tail portals, so labels are $O((k+\varepsilon^{-1}\log k)\,n\log n)$. A random
pair is refined with probability about $2\pi/(k\varepsilon)^3$; measured: 0 of 1,500 pairs at 10k and
0 of 200 at 100k.

### 3.8 Arbitrary points (A2A)

Corner bounds from the vertex certificates, plus exact refine below $R$ or when the interval is wider than
$1+2\varepsilon$ (**Theorem A2A**, `guarantee.md` §7).

### 3.9 Complexity

| | kd (Theorem 1) | capped (Theorem A) |
|---|---|---|
| labels | $\sim 3\bar\sigma D/g$ per vertex, $\Theta(n^{1.5})$ total | $O((k+\varepsilon^{-1}\log k)\,n\log n)$ |
| near tables | $\pi\rho g^2/(\varepsilon\lambda)^2$ per vertex | same |
| query | $O(\text{depth} + K_C)$ | same, + rare exact refine |
| build | one sub-mesh sweep per portal: $\Theta(n^{1.5})$ sweep area, $\Theta(n^{1.5}\log n)$ time | the cap also bounds the sweeps (at most $k + O(\varepsilon^{-1}\log k)$ per cut instead of $\Theta(\sqrt n)$ at the top): $O((k+\varepsilon^{-1}\log k)\,n\log^2 n)$ |

Build times assume a sweep over $m$ faces costs $O(m\log m)$; a cut at level $\ell$ sweeps sub-meshes of
$\Theta(n/2^\ell)$ faces, and there are $O(\log n)$ levels. The near-table sweeps add $O(n\,\varepsilon^{-2}\log n)$
in both versions.

Compare with Djidjev–Sommer's space $O(n\varepsilon^{-3/2}\log^2)$ (theory). We state plainly that
near-linear space is known in theory; our contribution is the certified, implemented index.

## 4. Engineering

1. **Hybrid sweeps.** VTP is 2.6–4.1× faster per full sweep; MMP stops early at the targets. The
   ratio rule picks between them per cut. VTP needs a per-sweep state reset and consistently oriented
   faces; when it rejects a sub-mesh, the cut falls back to MMP (138 of 26,327 sweeps on Gunnison10k).
2. **Coded labels.** Closed-loop deltas with float32 anchors per run. Decoded values are $\ge$ the stored
   ones, LB subtracts the decoding error, and the estimate is de-biased.
   - 100k vertices: 16-bit −38%, 8-bit −58% size;
   - mean error 0.0069% → 0.0075%;
   - query 0.38 → 1.3 µs.
3. **Parallel sweeps, coded during the build.** The top cuts spread their portals over the threads,
   the rest go one cut per thread. Labels are bit-identical to one thread; 2 threads give 2.0×. Peak
   memory is about the coded size. The evaluation uses one thread.
4. **Float32 soundness.** Labels are rounded up and LB keeps a relative margin.

## 5. Experiments

### 5.1 Setup

- **Machine and protocol:** Apple M1, single thread, random vertex pairs (seed 7) scored against exact
  MMP in double precision.
- **Meshes:** HorseMount 793, BigMount 1,440, HeadlightMount 2,470, RobinsonMountain 3,696 (the SE/EAR
  benchmark meshes), Gunnison10k 9,610, GunnisonForest100k 100,511.
- **[server]** UP-Oracle terrains at 0.5–2.5M faces.

Competitors, each run with its authors' code:

| competitor | code | notes |
|---|---|---|
| SE-Oracle, EAR-Oracle, K-Algo, Steiner schemes, MMP | EAR authors' harness | with our weighted-endpoint fix |
| UP-Oracle, WSPD-Oracle, RC-Oracle | UP-Oracle code | |
| DGG | macOS port | |
| NeuroGF | their network and training (`scripts/neurogf/neurogf.py`), trained on exact MMP distances and paths of our terrains | GPU; test-pair vertices excluded from training |
| exact MMP/VTP on the fly | Kirsanov MMP, VTP | |
| exact POI table | ours | POI setting only |

### 5.2 Main comparison (Table: `unified_results.md` summary)

| mesh | n | ours (ε 0.1): build / index / query / mean / max | EAR (sp 5): build / index / query / mean | DGG: query / mean | VTP per query |
|---|---:|---|---|---|---:|
| HorseMount | 793 | 3.1 s / 1.3 MB / 0.18 µs / 0.024% / 1.16% | 23.5 s / 189 MB / 1.3 ms / 3.10% | 1.4 ms / 0.0023% | 3.3 ms |
| RobinsonMountain | 3,696 | 44 s / 10.7 MB / 0.28 µs / 0.020% / 1.79% | 188 s / 672 MB / 3.8 ms / 3.50% | 7.9 ms / 0.0023% | 18 ms |
| Gunnison10k | 9,610 | 202 s / 37.4 MB / 0.36 µs / 0.0076% / 0.37% | 670 s / 702 MB / 10.2 ms / 3.43% | 22 ms / 0.0028% | 55 ms |
| GunnisonForest100k | 100,511 | 41–66 min / 534 MB (224 MB 8-bit) / 0.4–0.5 µs / 0.0069% / 0.12% (ε 0.2) | ran 14.3 h, no result (8 GB) | **[server]** | MMP 3.3 s |

All of our answers are sound, i.e. LB ≤ d ≤ UB on every pair. None of the competitors gives a
per-query bound.

### 5.3 Size versus error

`results/size_error.png`, `notes/size_results.md`. Our curve lies entirely below and to the left of
EAR-Oracle's: the coarsest setting (ε=1, 21 MB, 0.31% mean) is still 10× more accurate than EAR's best
and 30× smaller.

**Storage variants at 100k (ε 0.2):**

| configuration | storage | size | mean error | query |
|---|---|---:|---:|---:|
| kd | dense | 534 MB | 0.0069% | 0.38 µs |
| kd | 16-bit | 330 MB | 0.0069% | — |
| kd | 8-bit | 224 MB | 0.0075% | 1.3 µs |
| cap 100 | 8-bit | 198 MB | 0.017% | — |
| cap 50 | dense | 346 MB | 0.039% | — |
| cap 50 | 8-bit | 177 MB | 0.041% | 0.64 µs |
| kd | covers 0.05 | 230 MB | 0.79% | — (dominated by 8-bit) |

No pair needed refining. The shifted hierarchy loses at these sizes, as predicted:

| mesh | size vs kd | mean error |
|---|---|---|
| Gunnison10k | 2.1× larger | — |
| GunnisonForest100k | 1.4× larger (742 MB) | 0.12% |

**Matched accuracy.** Pushing ε to 3.2 with a fixed near radius reaches EAR's error level at a
fraction of its size:

| mesh | ours | EAR-Oracle |
|---|---|---|
| HorseMount | 3.7% mean at 0.66 MB | 3.1% at 189 MB |
| RobinsonMountain | 2.4% at 3.9 MB | 3.5% at 672 MB |

Every answer stays certified.

### 5.4 Scalability **[server]**

- Build and size against $n$, measured at 0.8k–100k: build $\propto n^{1.45}$, index $\propto n^{1.32}$.
  Projections are in `projection.md`.
- 1M+ runs with 8-bit labels and caps.
- Thread scaling (labels are identical across thread counts).

### 5.5 Weighted terrains

Steiner base graph with $m=5$, against EAR on the fixed harness:

| method | mean error | build | size | query |
|---|---:|---:|---:|---:|
| ours | 0.35% | 10.9 s | 2.9 MB | 1.3 µs |
| EAR, sp 5 | 2.42% | 25.9 s | 166 MB | 1.14 ms |

### 5.5b Obstacles **[server]**

Holes in TJ and GI at 500k faces: disks covering 5/10/20% of the area, and six walls with gaps that force
detours. Ours (uncapped and cap 64), EAR-Oracle, DGG and NeuroGF answer the same queries. Errors are
also reported by detour ratio $d/|st|$, in the bands [1, 1.01), [1.01, 1.05), [1.05, 1.2) and ≥ 1.2.
- **Soundness is unaffected:** `guarantee.md` §8, "Obstacles".
- **Gunnison10k pilot:** every answer sound; mean error 0.03–0.04%, as without obstacles.
- **Exact solvers:** MMP and VTP agree on every obstacle mesh (checked per mesh).

### 5.6 A2A

Certified on 5×100 random surface pairs. Refine bounds the worst case: BigMount's maximum UB error goes
from 20.2% to 10.2%.

### 5.7 The POI setting

In UP-Oracle's setting (50–2,500 POIs), the exact $P\times P$ table built by one VTP sweep per POI is
exact, sub-microsecond and the cheapest to build. Neither UP-Oracle nor our index beats it below a
crossover of about 1,000 POIs on 3.7k–10k-vertex meshes (`results/poi/`, **[pending]** Gunnison10k).
Our index is for every-vertex and A2A queries.

Measured crossovers (`notes/poi_results.md`):

| mesh | our POI index smaller from | our POI index faster to build from |
|---|---|---|
| RobinsonMountain | about 1,000 POIs | about 1,000 POIs |
| Gunnison10k | about 1,000 POIs | between 1,000 and 3,000 POIs |

At 3,000 POIs on Gunnison10k our index is 3.1 MB against the table's 36 MB.

**[check]** whether UP-Oracle's paper compares against this table; its selling point is updates.

### 5.8 Ablations

| ablation | measured effect |
|---|---|
| hybrid vs MMP vs VTP sweeps | 1.2–1.7× build |
| ellipse sub-meshes vs whole-mesh sweeps | — |
| model interpolation on/off | — |
| ε and g | Gunnison10k ε 0.05/0.1/0.2: 64/37/24 MB, mean 0.0024/0.0076/0.023% |
| 16/8-bit coding | see 5.3 |
| caps | see 5.3 |
| uniform vs scale-k defaults | `error_bound.md` §6 |

## 6. Related work

Summary of `notes/related_work.md`:
- exact and approximate geodesics;
- DB terrain oracles and the Kaul et al. bounds;
- theory oracles on polyhedral surfaces;
- planar portal oracles and labelings (Thorup, Klein, H2H, PLL, hub labels).

## 7. Limitations

- **Asymptotic size.** The guaranteed index stores $n^{1.5}$ labels. The capped variant is near-linear
  but relies on refine. The shifted hierarchy fixes the exponent only beyond about 1M vertices.
- **Build.** The build is $n^{1.5}$ sweep area, dominated by the top cuts. It is parallel, but the
  exponent stays.
- **Slopes.** $\lambda$ enters the near radius and the ellipse; steep or heavily weighted terrains cost more.
- **Static.** No updates; UP-Oracle's setting.
- **Ground truth at 1M+** relies on VTP/MMP sweeps (feasible on a server).

## 8. Conclusion

(after the server results)

## TODO list for the submission

1. Server runs (`projection.md` §4) with 8-bit labels and caps, single-thread timings plus a thread
   scaling curve.
2. A clean timing rerun of the M1 size experiments on an idle machine: the 100k build times varied from
   41 to 66 min with concurrent jobs.
3. Verify the Djidjev–Sommer technique and full text, and whether UP-Oracle compares against the exact
   POI table.
4. Figures: the size-error curve, scaling, an illustration of the cut, window, portals and ellipse, and
   a certificate width histogram.
5. Paper-quality proofs from `guarantee.md` §3, with the shifting-lemma variant in the appendix.
