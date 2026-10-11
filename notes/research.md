# Research notes: what moves accuracy and speed, and what the literature offers

This note records the survey behind the Certified Cut Index (`rgeo/oracles/certified.py`), which
ideas were adopted and what each one bought (measured, single thread, Apple M1), which were tried
and dropped, and what is left as future work, ranked by expected impact.

## 1. The competition and where it stands

| oracle | guarantee | build / size / query (HorseMount, 793 vertices) | mean error |
|---|---|---|---|
| SE-Oracle — Wei, Wong, Long, Mount, SIGMOD 2017; TODS 2022 | $(1+\varepsilon)$ on Steiner points, POI-to-POI | 6,252 s / 942 MB / 11.6 ms | 0.81% (min −7.1%) |
| EAR-Oracle — Huang, Wei, Wong, Tang, PACMMOD 2023 | $(1+\varepsilon)$ w.r.t. its Steiner graph, A2A | 23.5 s / 189 MB / 1.28 ms | 3.10% (max 24.5%) |
| K-Algo — Kaul, Wong, Jensen, PVLDB 2015 | on-the-fly, no index | — / — / 8.1 s | −0.02% |
| fixed / unfixed Steiner schemes (same harness) | on-the-fly | — / — / 26 ms … 23 s | 2.3–3.1% |
| MMP exact (same harness) | exact, no index | — / — / 34 ms | 0 |
| UP-Oracle — Yan, Wong, Jensen, TKDE 2025 (its authors' code, every vertex a POI) | $(1+\varepsilon)$ on a POI spanner, updatable | 26.7 s / 27.8 MB / 0.75 ms | 0.80% (max 5.8%) |
| WSPD-Oracle (SE line) / RC-Oracle (UP-Oracle code, every vertex a POI) | $(1+\varepsilon)$ | 12-13 s / 220-244 MB / 7 µs | 0.59% (max 10.2%) |
| DGG — Adikusuma, Fang, He, TOG 2020 (authors' code, macOS port), eps $=10^{-4}$ | mean error target, no bound | 0.7 s / 4.4 MB / 1.4 ms | 0.0023% (max 0.018%) |
| NeuroGF — Zhang et al., NeurIPS 2023 (authors' network, distance branch, 3,000 iterations) | none (learned) | 23 min / 0.25 MB / 280 µs (17 µs batched) | 0.46% (max 4.7%) |
| GeodesicEmbedding — Xia et al., arXiv:2108.13821 | none | no public code; not run | — |
| exact on the fly: MMP (Kirsanov) / VTP (Qin et al. 2016) | exact | — / — / 9.0 / 3.3 ms | 0 |

Numbers are the competitors' own C++ harness ([ItakEjgo/weighted_distance_oracle]) on our machine,
same queries, scored against double-precision MMP (`scripts/competitors/`). The newest paper in the
line (UP-Oracle, TKDE 2025) reports the cost of the static baselines on larger terrains: EAR-Oracle
needs 35,500 s ≈ 9.9 h on a terrain with 100k faces, SE-Oracle-style P2P oracles 35,100 s ≈ 9.8 h on
0.5 M faces with 250 POIs, and the best on-the-fly algorithms 2.4–3.2 h for 100 queries on 0.5 M faces.
On our 8 GB machine EAR-Oracle ran 51,512 s on the 100k-vertex GunnisonForest and exited without
query results.

What these oracles share: the guarantee is relative to a Steiner graph (SE, EAR) and the error that
is actually delivered is set by the number of Steiner points per edge, not by $\varepsilon$ — EAR's mean
error is 3.1% at sp=5 and still 2.3% at sp=20 while its index grows to 683 MB. None reports a
per-query bound.

## 2. Ideas adopted, and what each bought

| idea | source / reasoning | effect (measured) |
|---|---|---|
| **Early termination of MMP sweeps** | pygeodesic (Kirsanov's MMP) ignores its stop-at-targets test unless `max_distance` is finite; `max_distance=0` stops when every target is final | identical distances (bit-equal), 10–200× faster sweeps for nearby targets |
| **Provable, tight sub-meshes** | ellipse lemma: a geodesic of length $\le D$ stays within $D/(2\lambda)$ of the chord; $D$ from hub upper bounds | on Gunnison10k (max slope 14.9) the old dilation by $\lambda_{\max}$ swept the whole mesh at depths 2 and 5; now each cut touches only its ellipse region |
| **Portals by arc length + geometric tails** | the window lemma ($d \ge 2\lambda\,c$ for a crossing at offset $c$) lets the gap grow like $2\varepsilon\lambda c$ outside the box | tails cost $O(\frac{\sigma}{\varepsilon\lambda}\log)$ portals instead of a fixed padding budget that was not provable |
| **Near tables of radius $g/(\varepsilon\lambda)$** | the endpoint-near-a-cut gap: pairs closer than the radius are stored exactly, all others have $d > g/\varepsilon$ | turns "no bound near the cut" into the $(1+\varepsilon)$ theorem; ~90 entries per vertex |
| **Certified intervals** | $f$ is a sum of two functions that are 1-Lipschitz in arc length, so $\tfrac12(f_j+f_{j+1}) - A_j \le d$ | every answer comes with $\mathrm{LB} \le d \le \mathrm{UB}$; sound on every one of the 1,500 random pairs per mesh; UB/LB median width 1.6% |
| **Exact float32 handling** | labels rounded *up*, LB margin $2^{-21}$ | the first "unsound" pairs were 1e-8 relative rounding; fixed |
| **Steiner base graph for weighted terrains** | Aleksandrov–Maheshwari–Sack (JACM 2005), Lanthier–Maheshwari–Sack; the standard skeleton of weighted-region oracles | plain $m=5$ Dijkstra is 1.3% shorter on average (up to 7.5%) than the best path found before by any method; labels from it give $\mathrm{UB} \le (1+\eta)(1+\varepsilon)d$ |
| **Compiled query** | the Python query is ~15 small numpy calls | 30 µs → about 1 µs per pair (C, ctypes), results identical |

Net effect on the four meshes of the EAR/SE papers plus a 10k crop of GunnisonForest (ε = 0.1,
g = 0.5 edge lengths, 1,500 random pairs, all sound):

| mesh | n | build | size | max UB error | mean error | max error |
|---|---:|---:|---:|---:|---:|---:|
| HorseMount | 793 | 4.7 s | 1.3 MB | 3.41% | 0.024% | 1.16% |
| BigMount | 1,440 | 15.5 s | 3.3 MB | 2.12% | 0.022% | 0.86% |
| HeadlightMount | 2,470 | 38.1 s | 6.5 MB | 2.82% | 0.016% | 1.30% |
| RobinsonMountain | 3,696 | 61.0 s | 10.7 MB | 1.79% | 0.020% | 1.79% |
| Gunnison10k | 9,610 | 328 s | 37.4 MB | 0.94% | 0.008% | 0.37% |

against EAR-Oracle (sp=5): 23.5 / 52.4 / 107 / 188 / 670 s, 189–702 MB, 3.1–3.5% mean error.

## 3. Tried and dropped

* **Interpolating the label vectors before the portal minimum (fast A2A).** One portal scan instead
  of nine, but the minimum of an average is not the average of minima (Jensen gap, first order in
  the spread of the vertices' argmins): 0.8% vs 0.18% A2A error. Removed.
* **Thorup-style $\varepsilon$-covers** (Thorup, JACM 2004). Space $O(\frac{n}{\varepsilon}\log n)$ is provable (see
  `notes/guarantee.md` §6: $\ln 2/\varepsilon_c$ cover portals per side), but the answer is only as good as
  the cover tolerance $\varepsilon_c$, and matching the measured accuracy of full rows (mean $\sim 10^{-4}$)
  needs more cover portals than the cut has. Kept as the route to near-linear space when a coarser
  answer suffices (≥ 10^6 vertices).
* **Sweeping from the members instead of the portals.** At deep levels $|M| \approx |P_C|$, and a sweep
  from a member must reach the whole window (tails included), so it is the more expensive side.
* **16-bit labels.** Residual coding against the 3-d chord would halve the label bytes at ~$10^{-5}$
  relative error. Not needed at the sizes measured (the index is 20–150× smaller than EAR-Oracle's).

## 4. Future work, ranked (updated 2026-10-11)

Done since the first version of this list: **VTP** is a label solver (`cpp/src/solver_vtp.cpp`, exact to
$10^{-15}$, 2.6-4.1x per full sweep; the hybrid builds 1.2-1.7x faster), **DGG**, **UP-Oracle** and
**NeuroGF** run as competitors (`notes/unified_results.md`), and **A2A answers are certified**
(`rgeo.terrain.CertifiedA2A`). DGG was *not* adopted as a label source: its parameter bounds the
mean error, so the labels would be $(1+\eta)$ with $\eta$ unproven and the interval LB (which needs exact
labels) would be lost; with the hybrid exact sweeps the build is already within a small factor.

1. **Scale-dependent gap on shifted hierarchies** (§5.2): the change that alters the $n^{1.5}$ exponent.
2. **Parallel build**: label sweeps are independent; one parallel loop over portals.
3. **POI mode** for the UP-Oracle / SE-Oracle setting (§5.2, item 5).
4. **Tighter A2A upper bounds**: the corner route costs $|pa|(1-\cos\theta)$; storing the distance
   gradient at vertices (the direction of the last geodesic segment) would make the corner terms second
   order.
5. **GeodesicEmbedding**: reimplement from the paper (no public code) for completeness.

## 5. Breaking the $n^{1.5}$ growth (ideas, not implemented)

**Where it comes from.** A vertex stores one label per portal of every ancestor cut. The portals
of a cut are spaced at a fixed arc gap $g$, so a cut of side $S$ has $\approx \bar\sigma S/g$ of them,
and along a root-to-leaf path the cut lengths are $D, D/2, D/2, D/4, \dots$: about $3\bar\sigma D/g$
labels per vertex, $D \propto \sqrt n$ edges, hence $n^{1.5}$ labels. The build is one sweep per
portal over its cut's region: level $k$ has $2^k$ cuts $\times\ \sqrt{n}/2^{k/2}$ portals $\times\ n/2^k$
area, a geometric series dominated by the top levels, again $n^{1.5}$. Measured at 100k vertices
($\varepsilon=0.2$, $g=1$ edge): 534 MB, of which labels are about 90%; 1.9 h, of which label sweeps
are 97%.

### 5.1 Keeping the accuracy (same exact labels, same theorem, same certificates)

| idea | size | build | query | comment |
|---|---|---|---|---|
| faster exact sweeps: VTP / hybrid (done), sweeps in parallel, GPU exact solvers (PVTP, Du et al., CAD 2021) | — | constant ($1.4\times$ measured hybrid, $2.6$–$4.1\times$ per full sweep; $\times$cores) | — | sweeps are independent; the exponent stays |
| label coding: residual against the hyperbola $\sqrt{r_u^2+\sigma^2(x-x_u)^2}$ of each row, 16-bit quantum rounded up | $2$–$3\times$ | — | $+$ decode | UB stays a path cost; LB subtracts 2 quanta; error $<10^{-5}$ |
| shorter separators: choose axis and position minimising the cut's arc length (flat valleys) instead of the median | $\sim 10$–$30\%$ | same | — | constant; balance must stay within a factor |
| boundary-only labels + dense distance graphs (Fakcharoenphol–Rao style $r$-division): a vertex stores labels to its leaf region's boundary portals only, every region stores boundary-to-boundary tables, the query composes up the hierarchy with Monge min-plus | $O(n\log n)$ | $O(n \log n)$ sweep area (sweeps stay inside regions) | $O(\sqrt n \log^2 n)$, ms | the answer crosses $O(\log n)$ boundaries, each at a portal: to keep $(1+\varepsilon)$ the gap shrinks by $\log n$ |
| Voronoi-diagram oracles (Cabello 2017; Gawrychowski et al. 2018; Charalampopoulos, Gawrychowski, Mozes, Weimann 2019) transferred to cut portals | $n^{1+o(1)}$ | $n^{1+o(1)}$ | polylog / $n^{o(1)}$ | the strongest exact route on planar graphs; research-grade, large constants |

What none of these can give is $\mu$s queries **and** near-linear space **and** exact labels:
exact distance labels on planar graphs need $\Omega(n^{1/3})$ bits per vertex (Gavoille, Peleg,
Pérennes, Raz) and the best known use $\tilde O(\sqrt n)$, which is our regime. Below $n^{1.5}$ with a
flat label scan, accuracy has to become a parameter — which our theorem already allows.

### 5.2 Trading accuracy for size and time, keeping the $(1+\varepsilon)$ theorem and the certificates

1. **Scale-dependent gap on shifted hierarchies (the main idea).** The fixed gap $g$ is needed only
   for pairs that are close but separated by a high cut. Use regular (quadtree-like) splits and
   2–3 shifted copies of the hierarchy (Chan's shifting lemma: every pair at distance $r$ shares a
   cell of side $\le C r$ in one copy), and answer each pair in the copy where its separating cell is
   smallest. Then every pair answered at a cell of side $S$ has $r \ge S/C$, so the gap can be
   $\varepsilon\lambda S/C$: **$O(C/\varepsilon)$ portals per cut at every level.** Labels per vertex
   $O(\frac{C}{\varepsilon}\log n)$ per copy, build $O(\frac{C}{\varepsilon}\, n\log n)$ sweep area, query as now plus
   one LCA per copy. Window lemma, interval LB and Lemma 4 carry over unchanged (gap $\le \varepsilon d$ for
   every pair answered at that cut).
   *Hybrid:* keep the fixed gap below the scale $S^\ast \approx 40$ edges (where both rules give the same
   count per line) and switch to $\varepsilon\lambda S/C$ above it. *Constants, counted carefully* ($\varepsilon=0.2$, $g=1$
   edge, $C=6$ from the 3-copy shifting lemma, $\bar\sigma\approx1.2$): a line gets $\approx 36$ box portals + $\approx 26$
   tail portals, a quadtree cell has 2 lines, so $\approx124$ labels per vertex per level per copy; below $S^\ast$
   each copy keeps a fixed-gap kd-tree, $\approx250$ labels per vertex. Per vertex, against today's $\approx1{,}100$
   label entries at 100k growing like $\sqrt n$:

   | vertices | today | shifted hybrid (3 copies) | ratio |
   |---:|---:|---:|---:|
   | 100k | ~1,100 | ~1,900 | 0.6x (worse) |
   | 1.25M | ~3,900 | ~2,500 | 1.5x |
   | 10M | ~11,000 | ~3,100 | 3.5x |
   | 100M | ~35,000 | ~3,600 | ~10x |

   So the exponent changes ($n^{1.5}\to n\log n$, provably, build included) but the three copies, two lines
   per cell, $C=6$ and the tails cost about $\times 1.7$ at our scales; it pays off only beyond ~1M vertices
   unless the constants come down (one copy plus a band structure for the pairs that copy cuts badly;
   a smaller $C$; a larger top-level $\varepsilon$). *(An earlier version of this note estimated 2x at 100k and
   8x at 2.5M faces; that omitted the copies, the second line per cell and the tails.)* *Accuracy:*
   top-level gaps grow from $g$ to $\varepsilon S/C$; with the model interpolation the error is second order in
   gap$/d \le \varepsilon/C$, so we expect mean errors of $10^{-4}$–$10^{-3}$ instead of $10^{-5}$–$10^{-4}$ — to be measured.
2. **Thorup $\varepsilon$-covers** (space only, `notes/guarantee.md` §6): $\ln 2/\varepsilon_c$ portals per side per cut,
   $O(\frac{n}{\varepsilon_c}\log n)$ labels; combines with 1; LB from lower reconstructions
   $\ell(s,q) - A(q,p)$, so the certificate survives but widens.
3. **Sparse rows for far vertices** (space only): a vertex at distance $r_u \gg g$ from a cut has a
   smooth row; keep every $(r_u/cg)$-th label, reconstruct by the hyperbola model for the estimate and
   by $\min_q \ell(u,q)+A(q,p)$ for the UB. Labels per vertex become $O(\log n)$ per cut level for most
   vertices; the guarantee weakens to $(1+\varepsilon)(1+\varepsilon_c)$.
4. **Approximate label sources for the build** (build only): DGG (measured on HorseMount: 0.2 s
   precompute, 0.028% mean / 0.16% max error at its $\varepsilon=10^{-3}$, Dijkstra $\sim$1 ms per
   source), AVTP, Steiner graphs. Sweeps become graph Dijkstra; UB stays a path cost but
   $\le (1+\eta)(1+\varepsilon)d$ with $\eta$ unproven (DGG's $\varepsilon$ bounds the *mean* error), and the
   interval LB needs exact labels, so the certificate falls back to the chord.
5. **POI mode.** When queries are among $k \ll n$ points of interest (the UP-Oracle / SE-Oracle
   setting), build the tree over the POIs and tie $g$ to the POI spacing: labels $O(k\cdot\text{portals})$,
   sweeps stop at the POIs.

**Recommendation.** 1 (hybrid scale-dependent gap) is the change that alters the exponent while
keeping everything the paper claims (theorem, certificates, $\mu$s queries); 2 and the label
coding stack on it; VTP/hybrid sweeps and parallelism take the constant factors. A first
experiment: implement 1 with one extra shifted copy on Gunnison10k and GunnisonForest100k and
measure size, build time and error against the current index.
