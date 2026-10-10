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
| UP-Oracle — Yan, Wong, Jensen, TKDE 2025 | $(1+\varepsilon)$, updatable, POI-based | not run here | — |

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

## 4. Future work, ranked

1. **A faster exact solver.** VTP (Qin, Han, Yu, Yu, Zhang, SIGGRAPH 2016;
   github.com/YipengQin/VTP_source_code) reports 4–15× over MMP and ICH (Xin & Wang, TOG 2009). The
   build is >95% label sweeps, so this carries over almost one to one; it needs a small C API with
   the same stop-at-targets rule.
2. **Discrete geodesic graphs as label sources.** SVG (Ying, Wang, He, TOG 2013) and DGG (Wang, Fang,
   Wu, Xin, He, CAGD 2017; github.com/GeodesicGraph/DGG_and_SVG) give $O(n/\sqrt\varepsilon)$-edge graphs with
   error below $0.6\varepsilon$; Dijkstra on them is C-speed and would replace MMP where a $(1+\eta)$ label is
   acceptable (then UB $\le (1+\eta)(1+\varepsilon)d$ as with Steiner labels). Also a natural extra baseline
   (H2H on DGG).
3. **Certified A2A.** For a point $p$ in face $f$, $|d(p,\cdot) - d(v,\cdot)| \le |p - v|$ for the face's
   vertices, so the vertex bounds give admissible A2A bounds; with the near tables for close pairs this
   is a few lines and would make the A2A answers certified too.
4. **UP-Oracle as a competitor** (github.com/yanyinzhao/UpdatedStructureTerrainCode): the newest
   $(1+\varepsilon)$ oracle of the line; it is POI-based (P2P), so the comparison needs a POI set.
5. **Learned distance fields** (NeuroGF, NeurIPS 2023; GeodesicEmbedding, arXiv:2108.13821) answer in
   µs with no guarantee; useful as a contrast for the certified answers, not as a component.
6. **Parallel build.** Cuts at the same depth are independent; reported separately from the
   single-thread numbers.

[ItakEjgo/weighted_distance_oracle]: https://github.com/ItakEjgo/weighted_distance_oracle
