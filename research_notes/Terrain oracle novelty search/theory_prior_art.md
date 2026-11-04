# Theory / computational-geometry prior art for the Certified Cut Index (CCI)

Scope: approximate geodesic distance oracles on polyhedral surfaces and terrains; separator- and portal-based oracles for planar graphs and polygonal domains; per-query certified [LB, UB].
How to read the evidence labels:
- **[V]** means verified from an abstract, a publisher or LIPIcs/arXiv page, or the paper text fetched in this session.
- **[R]** means recalled from background knowledge of the paper and not re-checked in this session. Treat it as inferred and check it before citing.

Access note: public.lanl.gov (Djidjev's PDFs) refused connections, Springer redirected to a login page, and the sommer.jp/shortestpaths.org PDFs had certificate errors. So the internal construction details of Djidjev–Sommer 2011 and Aleksandrov et al. 2010 are **[R]**. Their headline bounds are **[V]**.

## Q1. Planar-graph portal oracles: Thorup 2004, Klein 2002, Kawarabayashi–Klein–Sommer, Gu–Xu, Chan–Skrepetos. What is the separator, and what is the portal argument?

### Takeaway
These oracles all cut with **shortest-path separators**: O(1) shortest paths of a shortest-path tree, applied recursively. Each vertex stores O(1/ε) "connections" (portals) on the separator paths of every ancestor level. Correctness uses the fact that the separator is itself a shortest path: distance along the separator is a true graph distance, so the triangle inequality gives an ε-cover. The guarantee is relative to the graph metric. The CCI cuts the continuous surface with straight axis-parallel lines in the xy chart. Those lines are not geodesics, so the shortest-path ε-cover lemma does not apply. The CCI instead bounds where a crossing can occur (the window argument) and places portals more sparsely with geometric spacing. It stores exact continuous geodesic distances and derives a per-query lower bound from the 1-Lipschitz property of each label row.

### Cited Findings
- **Thorup 2004.** "Compact oracles for reachability and approximate distances in planar digraphs", J. ACM 51(6):993–1024 (FOCS 2001). [V] It answers (1+ε)-approximate distances in weighted planar digraphs in O(log log Δ + 1/ε) time with near-linear space. [V] It uses a separator made of a constant number of dipaths (a "dipath decomposition") instead of an O(√n) vertex separator. [V] Each vertex gets an O(log n)-space label. — [KU pure record](https://test01.ku.dk/:obvius/pureproxy/3974268/en/publications/compact-oracles-for-reachability-and-approximate-distances-in-pla-2); [unpaywall/DOI 10.1145/1039488.1039493](https://unpaywall.org/10.1145%2F1039488.1039493)
- **Thorup 2004, details [R].** In the undirected case each vertex stores distances to O(1/ε) connection vertices on each separator path of each of O(log n) recursion levels. That gives space O(n ε⁻¹ log n) and query time O(1/ε). The ε-cover lemma needs the separator path to be a shortest path. — [Thorup JACM 2004 (DOI)](https://unpaywall.org/10.1145%2F1039488.1039493)
- **Klein 2002 [R].** "Preprocessing an undirected planar network to enable fast approximate distance queries" (SODA 2002) is an independent undirected analogue with the same shortest-path-separator plus portal structure. I did not fetch it this session.
- **Kawarabayashi–Klein–Sommer [R].** "Linear-space approximate distance oracles for planar, bounded-genus and minor-free graphs", ICALP 2011. It uses O(n) space with roughly O(ε⁻² log² n) query time. Kawarabayashi–Sommer–Thorup (SODA 2013) later gave "more compact" oracles. Both still use shortest-path separators and portals.
- **Gu & Xu.** "Constant query time (1+ε)-approximate distance oracle for planar graphs", ISAAC 2015. [V] Its preprocessing and space depend exponentially on 1/ε, per Chan–Skrepetos. — [Chan–Skrepetos ESA 2017 LIPIcs](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ESA.2017.25)
- **Chan & Skrepetos.** "Faster approximate diameter and distance oracles in planar graphs", ESA 2017 and Algorithmica 81(8), 2019. [V] It makes the 1/ε dependence of Gu–Xu's preprocessing and space polynomial. [V] The journal version corrects the query time to O(log(1/ε)) instead of the O(1) claimed earlier. — [LIPIcs ESA 2017](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ESA.2017.25); [Illinois Experts](https://experts.illinois.edu/en/publications/faster-approximate-diameter-and-distance-oracles-in-planar-graphs/)
- **Muller & Zachariasen.** "Fast and compact oracles for approximate distances in planar graphs", ESA 2007, LNCS 4698. [V, bibliographic only] It is an engineered variant of the Thorup/Klein approach. — [listed alongside Thorup ESA 2007](https://dblp.dagstuhl.de/rec/conf/esa/Thorup07.html)
- **Lineage of separators [V].** Early planar oracles relied on the Lipton–Tarjan separator and its refinements by Miller and Frederickson. — [Almost optimal exact distance oracles for planar graphs (J. ACM)](https://dl.acm.org/doi/10.1145/3580474)

### Inferences
- **The "single crossing" argument in this literature.** A shortest u–v path that crosses a separator path Q meets Q at some vertex w. Because Q is a shortest path, d(u,w) and d(w,v) can be approximated by portals spaced ε·d(u,Q) apart along Q. This works on any separator that is a geodesic.
- **What the CCI uses instead.**
  - A straight xy-line on a terrain is in general not a geodesic.
  - The CCI therefore replaces "distance along the separator is a shortest path" with a metric-distortion (λ) bound. Crossing the cut at offset c outside the bounding box costs at least 2λc, so crossings are confined to a window of radius μ/λ, and spacing can grow as max(g, 2ελc).
  - It replaces graph distances with exact continuous geodesics computed by MMP/VTP.
  - These are genuine differences in the separator type and in the proof device.
- Space comparison: Thorup/Klein labels are O(ε⁻¹ log n) per vertex, while the capped CCI uses O(k + ε⁻¹ log k) per level times O(log n) levels. The two are of similar shape. A reviewer will see the capped CCI as "Thorup-style labels on a geometric separator".

### Gaps
- I did not fetch the Klein 2002, KKS 2011 or Gu–Xu full texts this session, so their exact bounds are [R].

## Q2. Djidjev & Sommer 2011, "Approximate distance queries for weighted polyhedral surfaces" (ESA)

### Takeaway
This is the closest theoretical relative of the CCI on surfaces. It gives a (1+ε) oracle on a **weighted planar polyhedral surface**: space O(n ε^{-3/2} log²(n/ε) log(1/ε)) and query O(ε⁻¹ log(1/ε) + log log n). It is built on the Aleksandrov-style Steiner-point discretization, plus planar-style separator and portal machinery on that discretization **[R]**. The CCI differs in four ways:
- it uses straight axis-parallel xy cuts, not graph or Steiner-graph separators;
- it stores exact continuous MMP/VTP distances, not Steiner-graph distances;
- it is unweighted, with terrain-specific λ;
- it returns an adaptive certified interval.

### Cited Findings
- **Setting [V].** P is a planar polyhedral surface with n triangular faces, each with a positive weight. Path cost is the weighted Euclidean length within each face. — [Springer chapter (LNCS 6942, pp. 579–590)](https://rd.springer.com/chapter/10.1007/978-3-642-23719-5_49); [author PDF (LANL)](https://public.lanl.gov/djidjev/papers/esa11-quries.pdf)
- **Bounds [V].** For every ε ∈ (0,1) there is a distance oracle of size O(n ε^{-3/2} log²(n/ε) log(1/ε)) that answers (1+ε)-approximate queries in O(ε⁻¹ log(1/ε) + log log n) time. The authors say this improves significantly in n over the earlier tradeoff of Aleksandrov et al. (DCG 2010). — [search abstract of ESA'11 paper / LANL PDF](https://public.lanl.gov/djidjev/papers/esa11-quries.pdf)
- **Hidden constants [V].** The big-O hides constants that are logarithmic in the ratio of the largest to the smallest face weight and reciprocal in the sine of the smallest angle of P. — [Sommer survey / citing text](https://dl.acm.org/doi/10.1145/2530531)
- **Construction [R].**
  - Steiner points are placed on the edges in geometric progression away from the vertices, as in Aleksandrov–Maheshwari–Sack (J. ACM 2005), giving a graph G_ε.
  - For a planar surface, G_ε inherits a planar-like structure. Recursive separators of the surface triangulation, with Thorup-style portals (connections) on them, give the label structure.
  - The approximation is relative to the **exact weighted surface distance**: G_ε paths are surface paths, and the discretization is (1+ε)-faithful.
  - I could not open the PDF to confirm the exact separator type (vertex or cycle separator versus shortest-path separator).

### Inferences
- **Guarantee [R/inferred].** The guarantee is against the exact continuous metric, not only against the Steiner graph. Any estimate of this kind is the length of a realizable path, so it is an upper bound UB. A lower bound est/(1+ε) follows automatically. **A reviewer can say that every such oracle already "certifies" [est/(1+ε), est].** The CCI's novelty on the interval must therefore rest on two things:
  - the LB is per-query and adaptive, often much tighter than est/(1+ε), and computed from data (Lipschitz rows), not from the worst-case ε;
  - the interval drives an exact-refine fallback.
- Geometric spacing of Steiner points away from vertices (Aleksandrov–Maheshwari–Sack) is a close conceptual cousin of the CCI's portal spacing max(g, 2ελc) away from the box. Cite it and contrast: theirs is along mesh edges and around vertices, for discretization; the CCI's is along a cut line and around a region's box, for portal sparsification.
- Space comparison: Djidjev–Sommer is near-linear in n with ε^{-3/2}. The uncapped CCI is Θ(n^{1.5}), which is worse in theory. The capped CCI is O((k + ε⁻¹ log k) n log n), comparable in n but with an exact-refine fallback. Do not claim asymptotic superiority over Djidjev–Sommer.

### Gaps
- I could not retrieve the separator type and the per-vertex stored data for Djidjev–Sommer, because the PDF hosts were unreachable. This needs a manual check.

## Q3. Aleksandrov, Djidjev, Guo, Maheshwari, Nussbaum, Sack 2010, "Algorithms for approximate shortest path queries on weighted polyhedral surfaces" (DCG)

### Takeaway
This paper gives ε-approximate distance queries on **weighted polyhedral surfaces of arbitrary genus**. Space trades against query time through a parameter q, with O(q) query time, and single-source queries take O(log 1/ε). It is built on the Aleksandrov–Maheshwari–Sack Steiner discretization with separators of the discretization graph **[R for internals]**. The guarantee is relative to the exact weighted surface distance.

### Cited Findings
- **Publication [V].** DCG 44(4):762–801, 2010. An earlier version appeared at MFCS 2006. — [MaRDI portal](https://portal.mardi4nfdi.de/wiki/Publication:603866); [DOI 10.1007/11821069_9 (MFCS)](https://doi.org/10.1007/11821069_9); [LANL PDF](https://public.lanl.gov/djidjev/papers/SPquries.pdf)
- **Results [V].**
  - A data structure answers ε-approximate distance queries in O(q) time, for a tradeoff parameter q.
  - With one query point fixed, queries take O(log 1/ε).
  - The surface may be weighted and of arbitrary genus.
  - Source: same pages as above.
- **Relative weakness [V].** Djidjev–Sommer 2011 describe their own space/query tradeoff as a significant improvement in n over this one. — [ESA'11](https://rd.springer.com/chapter/10.1007/978-3-642-23719-5_49)
- **Construction [R].**
  - Steiner points are placed in geometric progression along the edges, about O(n ε^{-1/2} log(1/ε)) of them.
  - A separator theorem is used for the resulting graph or surface. For genus g, the separator is a small set of Steiner nodes and edges found by a planar/genus separator on the face graph.
  - Each node stores distances to separator nodes in a recursive decomposition.

### Inferences
- This is the same paradigm as Djidjev–Sommer: a combinatorial or graph separator on a mesh-derived discretization, with labels to separator nodes. The CCI's straight chart cuts with a window and geometric portal spacing, together with exact continuous labels, are different. The "store distances from every vertex to the separator portals of ancestor regions" idea is shared and must be credited.

### Gaps
- The precise space formula in terms of q and ε could not be verified.

## Q4. Agarwal–Har-Peled–Sharir–Varadarajan, Varadarajan–Agarwal, Har-Peled, WSPD for geodesic metrics, and recent SoCG/ESA/SODA work

### Takeaway
The classic AHSV, VA and Har-Peled works are **approximate single-pair path algorithms** on convex or nonconvex polyhedra. They are not oracles. The exception is Har-Peled 1999, which is about approximate shortest-path maps and the geodesic diameter on convex polytopes. Recent SoCG/ESA work on query structures (2016–2026) targets **exact** or approximate two-point queries in **planar polygonal domains**, not terrains. I found no SoCG/SODA/ESA geodesic distance oracle specific to terrains. The terrain oracle line is in the DB community (SIGMOD/VLDB/TODS).

### Cited Findings
- **Agarwal, Har-Peled, Sharir, Varadarajan.** "Approximating shortest paths on a convex polytope in three dimensions", J. ACM 1997. [V] It gives a simple, efficient approximate shortest-path algorithm on convex polytopes. — [JACM bibliography](https://projects.csail.mit.edu/jacm/References/agarwalhsv1997:567.html); [BibBase](https://bibbase.org/network/publication/agarwal-harpeled-sharir-varadarajan-approximatingshortestpathsonaconvexpolytopeinthreedimensions-1997)
- **Varadarajan & Agarwal.** "Approximating shortest paths on a nonconvex polyhedron", SIAM J. Comput. 30(4), 2001 (FOCS 1997). [V] Its paths are within 7(1+ε) of optimal, or 15(1+ε) in O(n^{8/5} log^{8/5} n) time. It is a single-pair algorithm. — [Iowa repository](https://iro.uiowa.edu/esploro/outputs/journalArticle/Approximating-shortest-paths-on-a-nonconvex/9984259490902771); [MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:2706111)
- **Har-Peled.** "Approximate shortest paths and geodesic diameter on convex polytopes in three dimensions", DCG 21:216–231, 1999. [V, bibliographic] [R] It builds approximate shortest-path maps from a source, giving (1+ε)-approximate queries from a fixed source in O(log n) time, and approximates the geodesic diameter. — [JACM author page](https://projects.csail.mit.edu/jacm/Authors/harpeledsariel.html)
- **WSPD for geodesic surface metrics.**
  - [V] I found no SoCG/SODA paper on a WSPD for terrain geodesic metrics.
  - [V] Har-Peled 2025, "Well-separated pairs decomposition revisited", gives WSPDs for general finite metrics and notes the application to distance oracles. — [arXiv 2509.05997](https://arxiv.org/html/2509.05997v1)
  - [V] Gao & Zhang (SICOMP 2005) gave a WSPD for unit-disk graph metrics. — (via search summary; same arXiv page cites it)
- **DB terrain oracles (for completeness; probably covered by another researcher).**
  - [V] Wei, Wong, Long, Mount, "Distance oracle on terrain surface", SIGMOD 2017. The SE oracle is ε-approximate and relative to exact terrain geodesics between POIs. — [Mount's PDF](https://www.cs.umd.edu/%7Emount/Papers/sigmod-2017-dist-oracle.pdf); [DBLP](https://dblp.dagstuhl.de/rec/conf/sigmod/WeiWLM17.html)
  - [R] SE uses a partition tree with a well-separated node-pair set, which is WSPD-like. Verify this.
  - [V] The TODS 2022 extension is "Proximity queries on terrain surface". — [PDF](https://cse.hkust.edu.hk/~raywong/paper/tods22-proximityQueries-terrain.pdf)
  - [V] EAR-Oracle (SIGMOD 2023) handles arbitrary points on the terrain. — [poster](https://itakejgo.github.io/static/poster/EAR-Oracle-poster.pdf)
- **Haitao Wang, ESA 2026.** "Shortest path map equivalence decompositions and applications", LIPIcs ESA.2026.160 (full version arXiv 2607.03416). [V] It gives new complexity bounds for SPM-equivalence decompositions (for example |Ψ_B(B)| = O(n^{4+ε}) and |Ψ_B(P)| = O(n^5)), with applications to exact two-point queries, geodesic diameter and center in polygonal domains. These are exact and use large space. — [LIPIcs PDF](https://drops.dagstuhl.de/storage/00lipics/lipics-vol388-esa2026/LIPIcs.ESA.2026.160/LIPIcs.ESA.2026.160.pdf)

### Inferences
- In terrain-specific CG theory, the main oracle papers remain Aleksandrov et al. 2010 and Djidjev–Sommer 2011 (weighted surfaces). Beyond these, I found no theory paper that targets terrains (xy-monotone surfaces) and exploits their 2D chart.

### Gaps
- The search may have missed niche EuroCG or CCCG abstracts on terrain oracles. I did not search Oh Eunjin's geodesic-query papers in depth. They concern simple polygons and polygonal domains (Voronoi diagrams, NN), not terrains.

## Q5. Is there prior work combining a hierarchy of straight cut lines over a terrain or polygonal domain with sampled points on the cut and a window/offset argument for sparse sampling?

### Takeaway
I found **no** work using straight axis-parallel (kd) cut lines on a terrain chart with portals for a geodesic oracle. The closest works use different separators:
- Thorup's ESA 2007 polygonal-domain oracle is the closest structural analogue. It uses recursive **shortest-path separators** of the free space (at most 6 shortest paths, 3 suffice), O(1/ε) anchor points per vertex per separator, space O((n/ε) log n), and is built on a Clarkson cone graph.
- van der Laan–Staals–Theunissen (arXiv 2026) follows the same "recursive separator plus anchor points" pattern.
- Gaur–Inkulu (arXiv 2026) handles convex obstacles.

I found nothing with a geometric-progression window on the cut whose bound comes from the cost of offset crossings.

### Cited Findings
- **Thorup.** "Compact oracles for approximate distances around obstacles in the plane", ESA 2007, LNCS 4698, pp. 383–394. — [DBLP](https://dblp.dagstuhl.de/rec/conf/esa/Thorup07.html); [researchr](https://researchr.org/publication/Thorup07%3A0)
  - Details as described by van der Laan et al. 2026 [V]:
    - recursive balanced separators made of shortest paths (≤6; 3 suffice);
    - O(k) "anchor points" or "connections" per polygon vertex on each separator, with k = ⌈1/ε⌉;
    - distances come from a cone graph with ε/3 error;
    - space O((n/ε) log n) and query O(1/ε³ + log n/(ε log log n));
    - a non-Real-RAM bit-manipulation model.
  - Source: [arXiv 2603.11775](https://arxiv.org/html/2603.11775v1)
- **van der Laan, Staals, Theunissen**, "Approximate dynamic nearest neighbor searching in a polygonal domain", arXiv 2603.11775, March 2026. [V]
  - It gives a (1+ε) two-point structure in a polygon with space O((n/ε) log n) and query O(ε⁻² log n).
  - Its separators are built from shortest paths using the Lipton–Tarjan fundamental-cycle separator on a triangulated shortest-path-tree graph: at most 3 shortest paths plus one segment.
  - The authors describe their approach as the same as Thorup's: recursively build a separator and estimate distances via discrete anchor points.
  - Source: [arXiv HTML](https://arxiv.org/html/2603.11775v1). A separate search result places it at SoCG 2026; the fetched page did not state a venue, so this is unconfirmed.
- **Gaur & Inkulu**, "Two-point approximate shortest path queries among convex polygonal obstacles in the plane", arXiv 2608.09677, August 2026. [V]
  - It gives (1+ε) multiplicative plus 13ℓ additive stretch.
  - Space O(n + (h/√ε)(h + 1/ε)); query O(ε^{-1/2} lg(h/√ε) + (h/ε^{2.5}) lg lg(h/√ε)).
  - The abstract page does not state the separator or portal technique.
  - Source: [arXiv abs](https://arxiv.org/abs/2608.09677)
- **Arikati, Chen, Chew, Das, Smid, Zaroliagis**, "Planar spanners and approximate shortest path queries among obstacles in the plane", ESA 1996, LNCS 1136:514–528. [V]
  - It builds planar spanners among obstacles in O(n log n) time.
  - The query structures run planar-graph separator machinery on the spanner.
  - It cites Chen's earlier (6+ε)-approximate L2 structure with O(n log n) space.
  - Source: [preprint PDF](https://ranger.uta.edu/~gdas/websitepages/preprints-papers/arikati.pdf); [researchr](https://researchr.org/publication/ArikatiCCDSZ96)
- **Chiang & Mitchell**, "Two-point Euclidean shortest path queries in the plane", SODA 1999, pp. 215–224. [V]
  - Exact queries with about O(n^11) space and O(log n) query time, with tradeoffs such as O(n^{5+10δ+ε}) space and O(n^{1-δ} log n) query.
  - The paper also notes that approximate two-point results predated it.
  - Sources: [DBLP](https://dblp.org/rec/conf/soda/ChiangM99.html); [de Berg–Miltzow–Staals SoCG 2024](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2024.17)
- **Guo, Maheshwari, Sack**, "Shortest path queries in polygonal domains", AAIM 2008, pp. 200–211. [V] Exact queries with O(n²) space and O(h log n) query time. — [DBLP](https://dblp.uni-trier.de/rec/conf/aaim/GuoMS08.html); cited in [de Berg–Miltzow–Staals 2024](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2024.17)
- **de Berg, Miltzow, Staals**, "Towards space efficient two-point shortest path queries in a polygonal domain", SoCG 2024. [V] Exact two-point queries with lower space. — [LIPIcs](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2024.17)

### Inferences
- **No prior separator–portal oracle uses straight cuts.** Every one I found (Thorup 2004/2007, Klein, KKS, Djidjev–Sommer, Aleksandrov et al., van der Laan et al. 2026) uses **shortest-path or graph separators**. The reason is the portal lemma: it needs the separator to be a shortest path, or at least needs distances along it to be controlled.
- **The CCI's straight-cut device.** Straight xy cuts plus a terrain distortion constant λ is a different device. It trades the generality of shortest-path separators (which work on any planar metric) for terrain-specific geometry. That makes it plausibly novel as a technique, but only for terrains or surfaces with bounded slope or distortion. Present λ honestly as an assumption-dependent constant, since the window μ/λ and the spacing 2ελc depend on it.
- **Geometric separators in the CG literature [R].** Straight-line or quadtree separators appear in Euclidean spanners, WSPD and Arora/Mitchell PTAS "portals" (Arora 1998 TSP uses quadtree dissection lines with portals). Arora-style portals on straight dissection lines are the best-known precedent for "portals on straight cuts". Cite and contrast them: Arora's portals approximate tours in Euclidean space under patching with m portals per side, with no oracle and no geodesic surface.
- **The geometric-progression window.** It is closest in spirit to the geometric-progression Steiner points of Aleksandrov–Maheshwari–Sack, and to Thorup's anchor spacing proportional to ε·distance. Neither confines crossings to a window derived from a lower bound on the cost of offset crossings.

### Gaps
- I did not verify the Arora 1998 portal details, or Mitchell's guillotine subdivisions, this session [R].
- I could not confirm whether any EuroCG or CCCG abstract uses kd cuts for surface geodesic oracles.

## Q6. Does any theory result give a per-query certified interval [LB, UB] around the true geodesic?

### Takeaway
I found no theory oracle paper that outputs an explicit, adaptive per-query [LB, UB] against the exact surface geodesic. However, all the path-realizing (1+ε) oracles above implicitly certify [est/(1+ε), est]. The CCI's claim must be framed as an adaptive, data-dependent LB (Lipschitz rows) that is usually much tighter than est/(1+ε), drives an exact-refine fallback, and is valid for the exact continuous geodesic. It should not be framed as "the first certified interval".

### Cited Findings
- The (1+ε) guarantees in Djidjev–Sommer and Aleksandrov et al. are stated relative to the weighted surface distance. — [ESA'11](https://rd.springer.com/chapter/10.1007/978-3-642-23719-5_49); [DCG'10](https://portal.mardi4nfdi.de/wiki/Publication:603866)
- In Thorup-style oracles the estimate is d(u,w) + d(w,v) for a portal w, which is a realizable walk length and so ≥ d [R, standard].
- Gaur–Inkulu 2026 even mixes multiplicative and additive error, so its implicit interval is not purely relative. — [arXiv 2608.09677](https://arxiv.org/abs/2608.09677)

### Inferences
- **Lipschitz lower bounds.** These appear in DB/AI practice: landmark or ALT lower bounds |d(L,u) − d(L,v)| by the triangle inequality, as in Goldberg–Harrelson ALT (SODA 2005) [R]. The CCI's LB is a cousin. It is a min over cut positions of the lower envelope of 1-Lipschitz label rows, which is a triangle-inequality bound over a continuum of virtual landmarks on the cut. Reviewers may connect it to ALT/landmark bounds, so cite and differentiate.
- **Recommended novelty claim.** "To our knowledge, the first terrain geodesic oracle with
  1. straight chart-aligned cut separators, using a window/offset lemma instead of shortest-path separators;
  2. exact continuous geodesic labels rather than Steiner-graph distances;
  3. a per-query adaptive certified interval used to trigger exact refinement."

  The prior-art set that must be cited:
  - Thorup 2004 and 2007, and Klein 2002;
  - Aleksandrov et al. 2010 and Djidjev–Sommer 2011;
  - Aleksandrov–Maheshwari–Sack 2005 (geometric Steiner spacing);
  - Wei et al. 2017 (SE) and EAR-Oracle 2023;
  - Arora-style portals and ALT landmarks.

### Gaps
- I did not search the certified-interval angle in numerical geometry (interval arithmetic for geodesics) or in exact-refinement "filter and refine" oracles in depth.
