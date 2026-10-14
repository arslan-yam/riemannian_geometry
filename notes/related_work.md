# Related work and positioning

This note gathers the literature around the Certified Cut Index. For each line of work it says what the
methods guarantee and how we differ, and it marks what is not yet verified. Measured numbers for the
methods we ran are in `notes/unified_results.md`; the competitor table is in `notes/research.md` §1.

## 1. Our index in one paragraph: a continuous H2H

H2H (Ouyang et al., SIGMOD 2018) takes a tree decomposition of a road network. Every vertex stores
its distance to each vertex in its ancestors' separators (bags), and a query takes a minimum over the
separator of the lowest common ancestor. The Certified Cut Index has the same structure on a
continuous surface:

| H2H on a graph | Certified Cut Index on a terrain |
|---|---|
| vertex separator (bag) of a nested dissection | axis-parallel cut line of a kd tree |
| every bag vertex | portals sampled on the cut at arc gap $g$ (tails grow like $2\varepsilon\lambda c$) |
| exact graph distances to ancestor bags | exact MMP/VTP geodesic distances to ancestor portals, on ellipse sub-meshes |
| query: min over the LCA bag | query: min over the LCA cut, plus an LB from the 1-Lipschitz rows |
| exact on the graph | $(1-\varepsilon)d \le \mathrm{LB} \le d \le \mathrm{UB} \le (1+\varepsilon)d$ on the surface, per query |

H2H must keep whole separators to stay exact. The window lemma lets us sample the cut instead and still
certify. This also answers "why not H2H on a Steiner graph?":
- With 5 Steiner points per edge the graph has about $16n$ nodes.
- Its top separators have $\Theta(\sqrt N)$ nodes, so H2H stores about $3$–$4\sqrt N$ distances per
  node: roughly 1 GB at 10k terrain vertices and 25–30 GB at 100k (estimate, not measured), against
  our 14–37 MB and 0.53 GB.
- Accuracy is capped by the Steiner graph's own error (0.5–2% at 5 points per edge).

## 2. Exact and approximate geodesics on meshes (single-source, no index)

* **Exact.**
  - MMP (Mitchell, Mount, Papadimitriou, SIAM J. Comput. 1987).
  - Chen–Han (SoCG 1990), $O(n^2)$.
  - Surazhsky et al. (SIGGRAPH 2005), the practical MMP behind Kirsanov's code that we use.
  - ICH (Xin & Wang, TOG 2009).
  - VTP (Qin et al., TOG 2016), our second label solver.
  - GPU variants (e.g. Du et al., CAD 2021).
  - These are our label sources and the "exact on the fly" baselines: 3.3 ms (VTP) to 3.3 s (MMP) per
    query at 0.8k–100k vertices.
* **Approximate.**
  - Fast marching (Kimmel & Sethian, PNAS 1998) and the heat method (Crane et al., TOG 2013).
  - Saddle vertex graph (Ying, Wang, He, TOG 2013) and its successor DGG (Adikusuma, Fang, He,
    TOG 2020), which we ran: 0.002–0.003% mean error, 1–22 ms queries, no per-query bound.
* **Learned.** NeuroGF (NeurIPS 2023), which we ran, and GeodesicEmbedding (arXiv:2108.13821), which has
  no code. No bound, and a fixed capacity, so the error grows with the mesh.

## 3. Terrain distance oracles and bounds in databases (mostly Wong's group, HKUST)

| work | setting | guarantee | per-query bound | notes |
|---|---|---|---|---|
| SE-Oracle — Wei, Wong, Long, Mount, SIGMOD 2017; TODS 2022 | POI-to-POI | $(1+\varepsilon)$ w.r.t. a Steiner graph | no | ran it: hours at 0.8k vertices |
| EAR-Oracle — Huang, Wei, Wong, Tang, SIGMOD 2023 (PACMMOD) | arbitrary points | $(1+\varepsilon)$ w.r.t. its Steiner graph, plus additive $2\delta$ | no | ran it: 2–3.5% mean error, 189–702 MB; failed at 100k on 8 GB |
| RC-Oracle — Yan, Wong, SIGMOD 2024 ("Proximity Queries on Point Clouds using Rapid Construction Path Oracle") | POIs on point clouds | $(1+\varepsilon)$ | no | ran it from UP-Oracle's code |
| UP-Oracle — Yan, Wong, Jensen, TKDE 37 (2025) | POIs, updatable terrain | $(1+\varepsilon)$ on a POI spanner | no | ran it: 0.8% mean, 5.8% max (every vertex a POI) |
| Kaul, Wong, Jensen, PVLDB 9(3) 2016, "New Lower and Upper Bounds for Shortest Distance Queries on Terrains" | pairs, on the fly | UB = a face-crossing path through cut vertices on the visited faces' edges; LB = $\lambda'\cdot$UB with $\lambda' \ge 0.9$ | **yes, but a fixed ratio $\le 1/0.9$, computed per pair by a graph search** | prunes kNN/range candidates; earlier bounds used the network distance and $\lambda = $ a function of the minimum face angle; experiments on Eagle Peak / Bearhead sub-regions (up to ~20k vertices) |
| Kaul, Wong, Yang, Jensen, PVLDB 2013/14, "Finding Shortest Paths on Terrains by Killing Two Birds with One Stone" (K-Algo) | pairs, on the fly | exact | — | ran it (8.1 s per query at 0.8k vertices) |
| Yan, Wong, MDM 2024 (best paper), "Efficient Shortest Path Queries on 3D Weighted Terrain Surfaces for Moving Objects" | weighted terrains, on the fly | approximate (Snell's law) | — | to check whether this is one of the on-the-fly algorithms in UP-Oracle's code |
| Yan, Wong, SIGMOD 2026, "Efficient Proximity Queries on Simplified Height Maps" (DOI 10.1145/3749162) | height maps; kNN / range | error up to 10% after simplification | no | on the fly after simplification, not a stored oracle. To check: is its path the continuous surface geodesic or a grid-graph path? |
| Wei, Xie, Wang, SIGMOD 2026, "Constrained Shortest Path Finding on Terrain Surfaces" (DOI 10.1145/3802031) | s–t path that visits given POIs | $2(1+\varepsilon)$, NP-hard problem | — | an application of an oracle, not a competitor (from an auto-generated summary; to check) |

**How we differ from the Kaul et al. bounds.** Our certificate claim has to be stated carefully,
because lower and upper bounds on terrain distance already exist.
- Theirs is a constant-ratio interval: LB $\ge 0.9\cdot$UB, i.e. up to about 11% wide. It is computed
  per pair by a search on a face-crossing graph, and its purpose is pruning.
- Ours is read from an index in about 1 µs. Its width is controlled by $\varepsilon$: median 0.7–2%,
  and $\mathrm{UB} \le (1+\varepsilon)d$ is guaranteed.
- Ours also comes with certify-or-refine, which makes every answer $(1+\varepsilon)$ even when an
  interval is wider.

## 4. Approximate oracles on polyhedral surfaces with theory (the closest prior work)

* **Aleksandrov, Djidjev, Guo, Maheshwari, Nussbaum, Sack** (MFCS 2006; DCG 44 (2010) 762–801),
  "Algorithms for Approximate Shortest Path Queries on Weighted Polyhedral Surfaces".
  - Answers $(1+\varepsilon)$ queries in $O(q)$ time for a query-time parameter $q$.
  - Queries from one fixed point take $O(\log\frac1\varepsilon)$.
  - According to Djidjev–Sommer, its space is roughly quadratic in $n$ for a query time comparable to theirs.
  - Uses Steiner points (Aleksandrov–Maheshwari–Sack, JACM 2005) and separators.
* **Djidjev, Sommer** (ESA 2011), "Approximate Distance Queries for Weighted Polyhedral Surfaces".
  Bounds from the abstract (planar polyhedral surface, $n$ faces, positive face weights):
  - preprocessing $O(n\varepsilon^{-2}\log^3\frac n\varepsilon\log^2\frac1\varepsilon)$
  - space $O(n\varepsilon^{-3/2}\log^2\frac n\varepsilon\log\frac1\varepsilon)$
  - query $O(\varepsilon^{-1}\log\frac1\varepsilon + \log\log n)$
  - the constants depend on the weight ratio and on the smallest face angle

  Checked against the Springer page (ESA 2011, LNCS 6942, doi 10.1007/978-3-642-23719-5_49; the full
  text is paywalled). The bounds above are confirmed. Their reference list points to the technique:
  - Steiner points (Aleksandrov–Maheshwari–Sack, JACM 2005);
  - planar distance oracles with path separators and portals (Thorup, JACM 2004; Klein, SODA 2002;
    Kawarabayashi–Klein–Sommer, ICALP 2011; Abraham–Gavoille, PODC 2006);
  - r-divisions (Fakcharoenphol–Rao).

  The page lists no experiments, and we know of no implementation.

**Where they beat us, and what we offer.** Asymptotically, Djidjev–Sommer's near-linear space beats our
guaranteed index ($\sim n^{1.5}$ labels). Near-linear space by itself is therefore not a claim we can make
as new. What neither paper offers, as far as we know:
- an implementation or measurements;
- exact continuous labels (theirs are Steiner-graph distances, so their constants include the Steiner density);
- a per-query interval;
- microsecond queries measured against the database oracles.

Our capped variant (at most $k_{\max}$ portals per cut, certify-or-refine) stores $O(n\,k_{\max}\log n)$
labels, which is near-linear like theirs. At 100k vertices it reached 177–198 MB with 8-bit labels
(`notes/size_results.md`), but its $(1+\varepsilon)$ guarantee relies on the refine step.

The comparison to make in the paper: their theory bound against our measured constants, and their
query-time term $\varepsilon^{-1}\log\frac1\varepsilon$ against our $O(K)$ scan over one cut.

## 5. Planar-graph oracles and labelings (where the portals come from)

* **Portal oracles.**
  - Klein (SODA 2002), Thorup (JACM 2004; $\varepsilon$-covers),
    Kawarabayashi–Klein–Sommer (ICALP 2011; linear space).
  - Our cuts are their separators; our certificates use that the label rows are 1-Lipschitz along
    the cut.
  - Thorup's covers are the route to near-linear space. Measured on our labels they lose to 8-bit
    coding at every size (`notes/size_results.md`).
* **Exact planar oracles.** Cabello (2017), Gawrychowski et al. (2018), Charalampopoulos, Gawrychowski,
  Mozes, Weimann (2019): Voronoi diagrams on separators. They are research-grade, with large constants.
* **Distance labeling and road-network indexes.**
  - Gavoille, Peleg, Pérennes, Raz (J. Algorithms 2004): planar graphs need labels of order $\sqrt n$
    in the worst case for exact distances; this is the graph analogue of our $n^{1.5}$.
  - Hub labels (Abraham, Delling, Goldberg, Werneck, SEA 2011), pruned landmark labeling (Akiba,
    Iwata, Yoshida, SIGMOD 2013), H2H (SIGMOD 2018).
* **Obstacles in the plane.** Arikati, Chen, Chew, Das, Smid, Zaroliagis (ESA 1996): planar spanners plus
  separators for approximate queries among polygonal obstacles. This is the Euclidean analogue of a
  Steiner graph with separators.

## 6. What to cite and claim

1. Present the index as "H2H for continuous surface geodesics, with certificates".
2. Contribution list:
   - exact continuous labels from sub-mesh sweeps (ellipse lemma, hybrid MMP/VTP);
   - the window lemma, which lets portals be sampled with a $(1+\varepsilon)$ guarantee;
   - per-query certified intervals and certify-or-refine;
   - closed-loop coded labels that keep the bounds sound;
   - certified A2A;
   - experiments against every implemented competitor.
3. In theory, stand next to Djidjev–Sommer honestly: near-linear space is theirs. Ours is the first
   implemented certified index, and its capped variant is near-linear too, with refine.
4. Cite Kaul et al. 2016 for terrain distance bounds and explain the difference: fixed ratio on the
   fly, versus $\varepsilon$-controlled bounds read from an index.

## Sources

- Djidjev & Sommer, ESA 2011: https://public.lanl.gov/djidjev/papers/esa11-quries.pdf
- Aleksandrov et al., MFCS 2006: https://doi.org/10.1007/11821069_9 (journal version DCG 44 (2010) 762–801)
- Kaul, Wong, Jensen, PVLDB 9(3) 2016: https://home.cse.hkust.edu.hk/~raywong/paper/vldb16-NewLowerUpperBound-Terrain.pdf
- UP-Oracle (TKDE 2025): https://vbn.aau.dk/da/publications/an-efficiently-updatable-path-oracle-for-terrain-surfaces/
- EAR-Oracle (SIGMOD 2023): https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf
- SE-Oracle (SIGMOD 2017): https://www.cs.umd.edu/~mount/Papers/sigmod-2017-dist-oracle.pdf
- Height maps (SIGMOD 2026): https://api.crossref.org/works/10.1145%2F3749162
- Constrained shortest paths (SIGMOD 2026): https://rmarcus.info/dbscholar/papers/hd0d93927d28322e2
- MDM 2024 best paper: https://cse.hkust.edu.hk/News/MDM2024/
