# Database-community prior art on terrain/surface distance oracles and surface spatial queries (novelty check for the Certified Cut Index)

Method note: full texts of Kaul et al. 2013 and 2015/16, SE-Oracle (SIGMOD'17), EAR-Oracle (SIGMOD'23), UP-Oracle (TKDE'24/'25), the TODS'26 point-cloud oracles paper, Shahabi et al. (VLDB'08), Deng et al. (ICDE'06, VLDBJ'08) and Yan's HKUST thesis (2025) were downloaded and text-searched. Items marked **[verified]** come from the paper text. Items marked **[abstract/secondary]** come only from abstracts, search snippets or other papers' related-work sections. **[inference]** marks our own reasoning.

## Q1. Terrain distance oracles in DB venues (SE, EAR, UP, RC, SP/WSPD): what they index, what they guarantee, and whether they give a lower bound or certificate

### Takeaway
Every DB terrain oracle stores a single **point estimate** and proves an **a-priori, one-sided or (1±ε) worst-case** bound. In practice the bound is usually measured against a **Steiner-graph distance** or an approximate sub-solver, not the exact continuous geodesic. None of them returns a per-query lower bound or interval certificate. None uses kd-style straight cut lines with portals, and none uses certify-or-refine. EAR-Oracle is the closest structurally. It cuts the 2D projection into a uniform grid of ζ² boxes, uses vertices of faces crossing box boundaries as "highway nodes", and gives each highway node a distance map to the Steiner points in its box. That is a one-level portal-like scheme, but its distances are approximate (Steiner/A-SSAD), it has no hierarchy of cuts, and it has no lower bound.

### Cited Findings
**SE-Oracle**: Wei, Wong, Long, Mount, "Distance Oracle on Terrain Surface", SIGMOD 2017, pp. 1211–1226, DOI 10.1145/3035918.3064038 — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod17-ProximityTerrain.pdf), [DBLP](https://dblp.org/rec/conf/sigmod/WeiWLM17.html)
- **Indexes [verified]:** a set P of n POIs, using a "compressed partition tree" of surface disks over the POIs plus a well-separated node-pair set (a WSPD-style decomposition with separation factor (2/ε+2)). The distance stored for each node pair is the geodesic distance between the disk centres, computed with an SSAD algorithm (Chen–Han cited) — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod17-ProximityTerrain.pdf)
- **Guarantee [verified]:** Theorem 1 says that for any p, q in P exactly one node pair contains ⟨p,q⟩, and its stored distance is an "ε-approximate distance of d_g(p,q)". The proof uses the triangle inequality d_g(c1,c2) − r terms ≤ d_g(p1,p2) ≤ d_g(c1,c2) + r terms, which gives a symmetric (1±ε) a-priori bound against the geodesic d_g. The node-pair set has size O(nh/ε^{2β}), with β∈[1.5,2] in practice, and the query takes O(h) — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod17-ProximityTerrain.pdf)
- **Arbitrary points [verified, appendix C]:** SE builds the oracle over **Steiner points**. A query for s,t returns min over neighbourhood Steiner points p∈N(s), q∈N(t) of d(s,p)+ẽd(p,q)+d(q,t), with |N(s)|·|N(t)| = 1/(sinθ·ε). Size is O(Nh/(sinθ√ε·ε^{2β})·log 1/ε) — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod17-ProximityTerrain.pdf). The TODS extension is described in search snippets as giving (1+2ε) for arbitrary points — [TODS'22 PDF](https://cse.hkust.edu.hk/~raywong/paper/tods22-proximityQueries-terrain.pdf) **[abstract/secondary]**
- **Experiments [verified]:** 48 GB memory budget, with oracle sizes reported in MB — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod17-ProximityTerrain.pdf). EAR-Oracle later reports that SE adapted to arbitrary points "exceeds 256GB of memory and takes more than 3 hours to build in a dataset with only 10,243 vertices" — [EAR PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf)
- **No lower bound or per-query certificate is output [verified by reading the query algorithm].**

**SE TODS extension**: Wei, Wong, Long, Mount, Samet, "Proximity Queries on Terrain Surface", ACM TODS 47(4), 2022 — [PDF](https://cse.hkust.edu.hk/~raywong/paper/tods22-proximityQueries-terrain.pdf)
- Adds a path oracle and uses SE for k-nearest/farthest-neighbour and top-k closest/farthest-pair queries. Claims construction time, size and query time up to 2, 3 and 5 orders of magnitude better than the best-known prior method — [PDF](https://cse.hkust.edu.hk/~raywong/paper/tods22-proximityQueries-terrain.pdf) **[abstract/secondary]**

**EAR-Oracle**: Huang, Wei, Wong, Tang, "EAR-Oracle: On Efficient Indexing for Distance Queries between Arbitrary Points on Terrain Surface", PACMMOD 1(1) (SIGMOD 2023), Art. 14, DOI 10.1145/3588694 — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf), [slides](https://home.cse.ust.hk/~raywong/ppt/sigmod24-EAR-Oracle-terrain-slides.pdf)
- **Indexes [verified]:** arbitrary query points. The method:
  1. Builds a base graph G_B with m Steiner points on each face's angle bisectors (O(mN) Steiner points).
  2. Imposes a uniform ζ×ζ grid on the x–y plane.
  3. Takes as highway nodes all vertices of faces that overlap box boundaries.
  4. Gives each highway node a distance map of *approximate* distances d̃_g to every Steiner point in its box, computed with the "A-SSAD-Box" approximate SSAD.
  5. Connects highway nodes by highway edges from a well-separated disk-pair decomposition (WSDPD) adopted from SE-Oracle.

  An inner-box query runs an online local search. An inter-box query builds a query graph from the highway network and the distance maps and runs Dijkstra on it — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf)
- **Guarantee [verified]:** Theorem 3 states d(s,t) ≤ (1+ε)(d_g(s,t)+2δ), "where δ is a small additive error". This is an **upper-side, multiplicative-plus-additive** bound, where the additive part comes from bisector discretisation (k faces × partition length, with an L_max term). Lemma 4 bounds highway-network distance by (1+ε) times the **approximate** distance d̃_g, not the exact one. No lower bound is stated or returned — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf)
- **Complexity [verified, Thm 5]:** build O(λγζ·mN log(mN) + N log N/ε^{2β} + Nh log N + Nh/ε^{2β}), space O(mλγN/ζ + Nh/ε^{2β}), query O(λγζ log(λγζ)) — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf)
- **Experiments [verified]:** 256 GB machine. Datasets go up to 1M faces ("EAR-Oracle could scale up to a large dataset with 1 million faces", while other indexes handled only 0.2M). Ground truth is the CGAL surface-mesh shortest path (exact, unweighted). On **weighted** terrains the reference is "Fixed Scheme" (a Steiner method) "since there is no algorithm could find exact shortest distances on weighted terrain surfaces". Reported practical relative errors are "less than 0.03" against a default theoretical ε of 0.2. Only mean relative error is plotted; no max error is reported — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf)

**UP-Oracle**: Yan, Wong, Jensen, "An Efficiently Updatable Path Oracle for Terrain Surfaces", IEEE TKDE 37(2):557–571, online Oct 2024, issue Feb 2025 — [PDF](https://yanyinzhao.github.io/files/TKDE2024-UpdatedStructureTerrain-paper.pdf)
- **Indexes [verified]:** POIs, using pairwise POI-to-POI **exact** shortest paths on the terrain before the update, a complete graph, and a (1+ε)-sub-graph (spanner) regenerated after updates. It can be adapted to arbitrary points (A2A) when no POIs are given — [PDF](https://yanyinzhao.github.io/files/TKDE2024-UpdatedStructureTerrain-paper.pdf)
- **Guarantee [verified]:** |Π′(s,t|T_aft)| ≤ (1+ε)|Π(s,t|T_aft)| for s,t in P, an upper-side bound only. The thesis version uses a lower bound internally for pruning, of the form "min_{v∈ΔV}|Π(u,v|T_bef)| − l_max", but only to decide which paths need recomputation, not as a query output — [Yan thesis](https://yanyinzhao.github.io/files/Thesis.pdf)
- **Figures [verified]:** on 0.5M faces with 250 POIs, update time is 400 s and output 22 MB, versus 35,100 s and 250 MB for the best-known oracle. In the A2A setting the figures are 480 s, 3 MB and 0.05 s query, versus 7,100 s, 150 MB and 5 s — [PDF](https://yanyinzhao.github.io/files/TKDE2024-UpdatedStructureTerrain-paper.pdf)

**RC-Oracle and TI-Oracle on point clouds**: Yan, Wong, (Sun), "Proximity Queries on Point Clouds using Rapid Construction Path Oracle", PACMMOD 2(1) (SIGMOD 2024) — [PDF](https://cse.hkust.edu.hk/~raywong/paper/sigmod24-proximityQueriesPointCloud.pdf), [DBLP](https://dblp.org/rec/conf/sigmod/YanWS24). Extended as Yan & Wong, "Efficient Path Oracles for Proximity Queries on Point Clouds", ACM TODS 51(1), Art. 5, Jan 2026, DOI 10.1145/3770577 — [PDF](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf)
- **Indexes and guarantee [verified]:** POIs (P2P), or arbitrary point to POI (TI-Oracle, A2P), on a **point cloud graph** in which each point connects to 8 neighbours with Euclidean edge weights. The guarantee is |Π_A(s,t|C)| ≤ (1+ε)|Π*(s,t|C)|, where Π* is the exact shortest path **on the point-cloud graph**, not the continuous surface geodesic. FastFly is plain Dijkstra on that graph — [TODS PDF](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf)
- **Figures [verified]:** best-known adapted TIN oracle (SE): 78,000 s, 1.5 GB, 150 s query. EAR adapted to A2P on 250k points with 500 POIs: 300 GB and 600 s query. RC-Oracle-A2A: 42,000 s, 100 GB, 12.5 s — [TODS PDF](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf)

**Height maps**: Yan & Wong, "Efficient Proximity Queries on Simplified Height Maps", PACMMOD 3(4) (SIGMOD 2026), Art. 244, DOI 10.1145/3749162 — [PDF](https://yanyinzhao.github.io/files/SIGMOD2026-HeightMap-paper.pdf)
- Simplifies a height map and answers shortest-path, kNN and range queries. Reported as up to 5× faster than a simplified point cloud and 1,340× faster than a simplified TIN, "with an error at most 10%" — [PDF](https://yanyinzhao.github.io/files/SIGMOD2026-HeightMap-paper.pdf) **[abstract/secondary]**. The distance is on the height-map grid graph, not a continuous TIN geodesic **[inference from the abstract framing; not verified in the body]**.

**Exact on-the-fly search**: Wei, Wong, Long, Mount, Samet, "On Efficient Shortest Path Computation on Terrain Surface: A Direction-Oriented Approach", IEEE TKDE 2024 — [PDF](https://cse.hkust.edu.hk/~raywong/paper/tkde24-shortestPathTerrainDirectionOriented.pdf). This is exact on-the-fly search, not an oracle **[abstract/secondary]**.

**Learned, non-DB**: "De-coupled NeuroGF for Shortest Path Distance Approximations on Large Terrain Graphs", ICML 2025 — [poster page](https://icml.cc/virtual/2025/poster/45376). A learned approximator on terrain DEM graphs with no certified bound reported **[abstract/secondary]**.

### Inferences
- **[inference]** Novelty axis 1 (certificate). DB terrain oracles (SE, EAR, UP, RC, TI) return only an estimate, together with a global a-priori guarantee. SE's is two-sided (1±ε). EAR's, UP's and RC's are upper-side and relative to an approximate (Steiner) or graph metric. A per-query interval [LB, UB] with UB/LB ≤ (1+ε)/(1−ε) against the exact continuous geodesic does not appear in any of them.
- **[inference]** Novelty axis 2 (reference metric). EAR's own distance maps hold Steiner-based approximate distances, and its weighted-terrain evaluation uses a Steiner method (Fixed Scheme) as ground truth. So EAR's 2–3% error claims are measured against approximate references on weighted inputs. The authors' observation of 2–3% mean and up to 25% max error against the exact geodesic is consistent with the paper reporting only mean relative errors (< 0.03); the max error is never reported.
- **[inference]** Closest structural precedent. EAR's "grid boxes + highway nodes on box boundaries + per-node distance maps to in-box points" is a flat, single-level analogue of "cuts + portals + per-vertex distances to portals". The differences:
  - The Certified Cut Index uses a recursive kd hierarchy of straight cuts with ancestor-cut portals; EAR has one level.
  - Portals are sampled within a bounded window with geometrically growing spacing; EAR's highway nodes are all mesh vertices of boundary-crossing faces.
  - Stored distances are exact MMP/VTP; EAR's come from approximate A-SSAD.
  - The cut index combines portal distances by min-plus over a separator, not Dijkstra over a highway graph.
  - It outputs a lower bound.

  The paper should cite EAR and say this explicitly.
- **[inference]** SE is the only DB oracle whose stored distances are exact geodesics (SSAD from disk centres). Its WSPD structure (disk pairs), however, differs from separator or portal labelling.

### Gaps
- No paper text for an "RC-Oracle in PVLDB 2024" was found. RC-Oracle appeared in SIGMOD/PACMMOD 2024 and TODS 2026, not PVLDB. No "WSPD-Oracle" as a separately named DB paper was found; WSPD/WSDPD is internal to SE and EAR.
- EAR-Oracle's absolute oracle sizes in GB at 0.5M or 1M faces were not extracted (they are in figures only). The authors' measurement of 333 GB at 0.5M faces exceeds the paper's 256 GB budget statement, which may reflect different parameters (m, ζ, ε). This should be reconciled before citing.
- The SIGMOD'26 height-map paper's error definition (max or mean, against which reference) was not verified in the body.

## Q2. Kaul et al. lower/upper bounds (PVLDB 2013 and PVLDB 2015/16): what bounds, how computed per query, ratio

### Takeaway
These are the only DB works found that compute a provable per-query **[LB, UB] interval on the exact surface distance**, and they are therefore the most important prior art for the "certified interval" claim. They are **index-free and on-the-fly**: a Dijkstra-style shortest-path computation on a (face-crossing) augmented graph per query. The LB/UB ratio is a constant fixed by a parameter. The 2016 paper's λ is about 0.9 at its default and tighter at smaller ε, but it is not presented as (1−ε)-tight with ε→0 at microsecond cost. The bounds are used only to prune candidates before exact Chen–Han computations in kNN, RkNN and range queries.

### Cited Findings
**Kaul, Wong, Yang, Jensen, "Finding Shortest Paths on Terrains by Killing Two Birds with One Stone", PVLDB 7(1):73–84, 2013 (VLDB 2014)** — [PDF](https://www.vldb.org/pvldb/vol7/p73-kaul.pdf)
- **[verified]** Both bounds come from one network (Delaunay-graph) shortest path Π_G: λ·|Π_G(s,t)| ≤ |Π(s,t)| ≤ |Π_G(s,t)|, with λ = min{sin θ_min/2, sin θ_min·cos θ_min}, where θ_min is the minimum interior angle of the triangulation. The lower bound is "upto 2.8 times larger than the previous lower bound" (Euclidean) on the largest dataset.
- **[verified]** A smaller graph G′ trades tightness for speed with guaranteed looseness ω: (1/ω)L ≤ L̃ ≤ L and U ≤ Ũ ≤ ωU.
- Reported speedups are up to 43× for spatial queries — [secondary summary](https://vbn.aau.dk/en/publications/new-lower-and-upper-bounds-for-shortest-distance-queries-on-terra)

**Kaul, Wong, Jensen, "New Lower and Upper Bounds for Shortest Distance Queries on Terrains", PVLDB 9(3):168–179, 2015 (VLDB 2016)** — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/vldb16-NewLowerUpperBound-Terrain.pdf), [DBLP](https://dblp1.uni-trier.de/rec/journals/pvldb/KaulWJ15.html)
- **[verified]** The paper criticises the 2013 bound because it depends on θ_min. Users were forced to use a constrained Delaunay triangulation with θ_min > 45°, and on skinny triangles the bound becomes "worse than even the Euclidean distance".
- **[verified]** New bound (Theorem 1): λ·|Π_FC(s,t)| ≤ |Π(s,t)| ≤ |Π_FC(s,t)|, with λ = 1 − 1/K. Π_FC is the shortest "face-crossing" path through cut-vertices (Steiner points) placed on the edges of faces **visited** by an A*-like expansion. With the error parameter, λ = 1 − ε/(K+c), K ≥ 10, and λ = 0.9 at ε = 1. The bound is computed per query by Algorithm 1, which places cut vertices on the fly. There is no precomputed index.
- **[verified]** Empirical tightness on the EP dataset at ε = 0.25: mean exact distance D_s = 9,075 m, OurLB = 8,862 m, previous LB = 5,075 m. That gives LB/D_s ≈ 0.977 on average (our arithmetic). The upper bound is about 0.92× the network-distance upper bound at default settings. Pruning removes up to 450% more objects, with "nearly 82 times faster" surface kNN. Datasets reach 1 million vertices in the scalability study.
- **[verified]** The related-work section cites Varadarajan et al., who partition the polyhedron with the Lipton–Tarjan planar separator theorem to get a 7-approximation. This is cited only as a computational-geometry approximation, not used in a DB index.

### Inferences
- **[inference]** The certified-interval idea per se, LB ≤ d ≤ UB against the exact geodesic, is **not new in DB**. Kaul et al. 2013 and 2016 do exactly this, for pruning. The Certified Cut Index differs in four ways:
  1. **Speed.** Kaul's bounds need an on-the-fly Dijkstra per pair, which is far from microsecond queries. Our intervals come from a precomputed index in about 1 µs.
  2. **Tightness.** Kaul's ratio is a fixed constant (λ ≈ 0.9–0.98 in practice), with no ε→0 index-based guarantee of (1−ε)d ≤ LB and UB ≤ (1+ε)d.
  3. **Domain.** Kaul's bounds are vertex-to-vertex, on unweighted TINs.
  4. **Use of the interval.** Kaul's bounds are never exposed as a query answer and are not combined with certify-or-refine. They are a filter before Chen–Han refinement, which is close to our capped certify-or-refine in spirit.

  The paper should position certify-or-refine as "the Kaul filter-and-refine paradigm, but with index-backed (1±ε) intervals".
- **[inference]** Kaul's LB proof works by bounding each face-crossing segment (λ|o_i o_{i+1}| ≤ |p_i p_{i+1}|). Our LB comes from separator and portal reasoning instead (a path must cross a cut, and portal spacing bounds the detour). These are different proof techniques.

### Gaps
- Query-time numbers for the bound computation itself (ms per pair) in Kaul 2016 were not extracted. Figure 10(b) shows that PrevLB is faster to compute than OurLB.
- Kaul 2016's handling of weighted terrains and arbitrary (non-vertex) points was not verified. The paper's examples are vertex-based.

## Q3. Surface kNN and range queries (Deng ICDE'06 and VLDBJ'08 MR3, Shahabi VLDB'08 SIR-tree, Xing GIS'10, Liu & Wong SIGMOD'11, skyline): how they use bounds

### Takeaway
Early surface-kNN work all uses **filter-and-refine with lower and upper bounds**: Euclidean distance as the lower bound, network distance or multi-resolution distance as the upper bound, and multi-resolution "support networks" for tighter lower bounds. The bounds serve only for pruning and ranking candidates. They come from multi-resolution meshes or Voronoi-like cells, not from an ε-tight distance index, and no work returns certified intervals as answers.

### Cited Findings
- **Deng, Zhou, Shen, Xu, Lin, "Surface k-NN Query Processing", ICDE 2006** — [PDF](https://infolab.usc.edu/csci599/Fall2009/papers/Surface%20k-NN%20Query%20Processing.pdf) **[verified]**. It ranks objects by a distance range [lb, ub] computed on a multi-resolution terrain (Direct Mesh). It observes that Euclidean distance "can always be used as a lower bound" and network distance as an upper bound, but that these "do not satisfy" monotonicity across resolutions. Chen–Han is cited as the O(n²) exact method.
- **Deng, Zhou, Shen, Liu, Xu, Lin, "A multi-resolution surface distance model for k-NN query processing", VLDB Journal 17(5):1101–1119, 2008** — [Springer](https://link.springer.com/article/10.1007/s00778-007-0053-2), [PDF](https://cgi.cse.unsw.edu.au/~lxue/paper/VLDBJ2007.pdf) **[verified]**. Algorithm MR3 uses a Distance Multi-resolution Mesh (DMTM) for progressive upper bounds and a Multi-resolution Support Distance Network (MSDN) for lower bounds. It estimates d_S at resolution r by lb_r(q,a) ≤ d_S(q,a) ≤ ub_r(q,a) and ranks by distance range.
- **Shahabi, Tang, Xing, "Indexing Land Surface for Efficient kNN Query", PVLDB 1(1):1020–1031, 2008** — [PDF](https://infolab.usc.edu/csci599/Fall2009/papers/Indexing%20land%20surface%20for%20efficient%20kNN%20query.pdf) **[verified]**.
  - Index: a Tight Surface Index (TSI) and a Loose Surface Index (LSI) stored in an SIR-tree. They form Voronoi-like tight and loose cells around **sites**. A query inside a tight cell is answered immediately; sites whose loose cell excludes the query are pruned. Exact Chen–Han computations handle the rest.
  - Bounds: 3D Euclidean distance as the lower bound and network distance as the upper bound.
  - The paper criticises Deng et al. because MR3 "only provides approximate answers… neither guarantees the correctness of the computed distance nor provides the actual shortest surface path".
- **Xing, Shahabi, Pan, "Continuous monitoring of nearest neighbors on land surface", PVLDB 2(1):1114–1125, 2009** — cited in [TODS'26 refs](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf) **[secondary; content not verified]**.
- **Liu & Wong, "Finding shortest path on land surface", SIGMOD 2011, pp. 433–444** — cited in [TODS'26 refs](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf) **[secondary; content not verified]**.
- **Deng & Zhou, "Expansion-based algorithms for finding single pair shortest path on surface", W2GIS 2004** — cited in [TODS'26 refs](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf) **[secondary]**.
- **[verified]** Yan's thesis and TODS'26 call Shahabi et al.'s index "SU-Oracle (SUrface Oracle)" with O(N log² N) kNN. They call Deng et al.'s multi-resolution work O(N²). Both are adapted as baselines: SU-Oracle-Adapt kNN took 1,520 s versus 25 s for RC-Oracle on their workload — [Yan thesis](https://yanyinzhao.github.io/files/Thesis.pdf)
- **Kaul 2013 and 2016 [verified]** plug their bounds into surface kNN, reverse NN (the MSRNN variants) and range queries as filters before exact refinement. See Q2 for the links.
- **SE TODS'22 [abstract]** answers kNN and top-k pair queries from ε-approximate oracle distances, not certified intervals — [PDF](https://cse.hkust.edu.hk/~raywong/paper/tods22-proximityQueries-terrain.pdf)
- **RC-Oracle and TI-Oracle [verified]** answer "(1+ε)-approximate kNN and range queries" relative to point-cloud graph distance — [TODS PDF](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf)

### Inferences
- **[inference]** Certified kNN and range queries decided from (1±ε) intervals, with exact refinement only for ambiguous candidates, are a direct descendant of the Deng/Shahabi/Kaul filter-and-refine line. Our contribution there is (a) intervals that are ε-tight and index-backed, rather than Euclidean, network or multi-resolution bounds, and (b) the ability to decide results exactly, with an exactness certificate, rather than return "(1+ε)-approximate kNN" as RC-Oracle and SE do. Claims should be phrased against both.

### Gaps
- "Scalable shortest paths browsing on land surface" (Xing & Shahabi, ACM GIS 2010) and a "skyline on terrain" paper by Yan, Liu or Xiong were not located or verified in this pass. They are not listed in TODS'26's reference list excerpts that were checked.
- The exact bound formulas of MSDN in MR3 were not extracted beyond lb_r ≤ d_S ≤ ub_r.

## Q4. Recent work (2023–2026), separator/portal labelling on terrains, and certified or interval-valued distance answers in DB

### Takeaway
2023–2026 DB terrain and surface work comes almost entirely from the HKUST Wong group: EAR (SIGMOD'23), UP (TKDE'24/25), RC and TI (SIGMOD'24, TODS'26), height maps (SIGMOD'26), and direction-oriented exact search (TKDE'24). None uses hierarchical separators with portals, and none returns certified intervals. Outside terrains, DB road-network work (ALT, PLL, landmark estimation) gives either exact answers or uncertified estimates. No DB paper was found that returns a guaranteed (1±ε) per-query interval as the answer with certify-or-refine.

### Cited Findings
- **[verified]** The terrain-oracle lineage and figures are as listed in Q1 — [EAR](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf), [UP](https://yanyinzhao.github.io/files/TKDE2024-UpdatedStructureTerrain-paper.pdf), [TODS'26](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf), [HeightMap SIGMOD'26](https://yanyinzhao.github.io/files/SIGMOD2026-HeightMap-paper.pdf)
- **[verified]** Yan's HKUST PhD thesis, "Shortest Path Queries on Terrain Surfaces and Point Clouds" (2025), surveys this line; text search found no separator- or portal-based terrain oracle — [Thesis](https://yanyinzhao.github.io/files/Thesis.pdf), [defence page](https://cse.hkust.edu.hk/pg/defenses/S25/yyanas-26-02-2025.html)
- **[verified]** EAR-Oracle cites Djidjev & Sommer 2011, "Approximate distance queries for weighted polyhedral surfaces", which is separator-based work from computational geometry, not DB. EAR does not adopt a separator approach — [EAR PDF](https://home.cse.hkust.edu.hk/~raywong/paper/sigmod23-EAR-Oracle-terrain.pdf)
- **[verified]** Kaul 2016 cites Varadarajan et al.'s planar-separator-based 7-approximation, as computational-geometry context only — [PDF](https://home.cse.hkust.edu.hk/~raywong/paper/vldb16-NewLowerUpperBound-Terrain.pdf)
- **[abstract/secondary]** On spatial networks, Sankaranarayanan & Samet, "Distance oracles for spatial networks" (ICDE 2009) and "Path oracles for spatial networks" (PVLDB 2(1), 2009) give WSPD-based ε-approximate oracles. These are the conceptual ancestors of SE — cited in [TODS'26 refs](https://yanyinzhao.github.io/files/TODS2025-PointCloudOracles-paper.pdf)
- **[abstract/secondary]** Qiao et al., "Approximate Shortest Distance Computing: A Query-Dependent Local Landmark Scheme" (ICDE 2012), is a landmark-based estimate. It notes that global landmark embedding can have worst-case relative error 68 on US roads and gives no certified interval — [PDF](https://www1.se.cuhk.edu.hk/~hcheng/paper/ICDE12_conf_full_262.pdf)
- **[abstract/secondary]** Goldberg & Harrelson's ALT landmarks give triangle-inequality lower bounds for A* search (exact search, not interval answers) — [PDF](https://www.cs.princeton.edu/courses/archive/spring06/cos423/Handouts/GH05.pdf)
- **[abstract/secondary]** The learned terrain distance approximator De-coupled NeuroGF (ICML 2025, not a DB venue) gives no certificates — [ICML](https://icml.cc/virtual/2025/poster/45376)

### Inferences
- **[inference]** No DB work found combines (i) a recursive straight-cut kd hierarchy over the 2D chart, (ii) portal sampling with geometric spacing in a bounded window, (iii) exact continuous-geodesic (MMP/VTP) labels from vertices to ancestor-cut portals, and (iv) per-query certified (1±ε) intervals. The combination appears novel for the DB community. Individually:
  - Separator and portal labelling is known in computational geometry and in road-network theory (Thorup-, Klein- or Djidjev–Sommer-style).
  - Grid partitions with boundary hub nodes exist in EAR.
  - Per-query LB/UB for pruning exist in Kaul 2013 and 2016.
  - WSPD oracles exist in SE.
- **[inference]** Strongest likely reviewer objections: "Kaul et al. already give LB/UB on terrain distance" (answer: on-the-fly, constant-factor, not index-backed or ε-tight), and "EAR already partitions the plane and stores boundary-node distance maps" (answer: single-level, approximate Steiner distances, no lower bound, Steiner-graph-relative guarantee with an additive δ).

### Gaps
- No DBLP-wide sweep of PVLDB, SIGMOD or ICDE 2025–2026 proceedings was performed for other groups' terrain papers. Searches returned only the HKUST line, so a missing paper from another group cannot be excluded.
- No DB paper on "exact-or-certified" (certify-or-refine) distance answers in road networks was found. This is absence of evidence from a limited search (2 queries), not proof of absence.
