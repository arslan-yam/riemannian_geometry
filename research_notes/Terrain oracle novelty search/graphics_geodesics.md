# Graphics / geometry-processing prior art on precomputed, fast geodesic distance queries (novelty check for the Certified Cut Index)

Verification legend: **[V]** means I read it in the paper's full text during this search. **[V-abs]** means it comes from the abstract or a search snippet only. **[2nd]** means a later paper reports it. **[M]** means it comes from my background knowledge and was not re-checked in this session; treat it as unverified.

## Q0. Closest prior art overall: Gotsman & Hormann, separator-based geodesic query database ("SEP"), SIGGRAPH Asia 2022

### Takeaway
Gotsman & Hormann 2022 is the closest prior art I found. It is a graphics paper that already does three things: it recursively bisects the mesh with separator curves, it stores per-vertex distances to the separators of every ancestor node, and it answers a query as min over shared separators of min_t (d_x(t) + d_y(t)), with exact leaf pairs as an option. Its conclusion even suggests replacing graph-cut separators with axis-parallel kd-style cuts. What it does **not** have is any per-query certified interval, any lower bound, exact continuous label values, or a window or geometric-spacing portal scheme. Its relative-error lemma only holds under a sup-norm fit assumption that the implementation does not enforce. Novelty claims for the Certified Cut Index (CCI) have to be repositioned around **certification** and **exact, sparse, provably sufficient portal labels**, not around the separator hierarchy.

### Cited Findings
- Citation: C. Gotsman, K. Hormann, "Compressing Geodesic Information for Fast Point-to-Point Geodesic Distance Queries", SIGGRAPH Asia 2022 Conference Papers, Article 4, DOI 10.1145/3550469.3555412 [V] — [author PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf); [dblp](https://dblp.org/rec/conf/siggrapha/GotsmanH22)
- Complexity claims: database size O(m n log n), built in O((f(n) + m³ n)√n) preprocessing time, where f(n) is the cost of one single-source field and m is "a constant that depends on the geometric complexity". Query time is O(m log n) worst case and O(m) on average [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- Theorem 2.1, the core identity: for vertices x, y in different leaves, d(x,y) = min over c in C of min_t [d^c_x(t) + d^c_y(t)]. Here C is the set of separator curves on the common ancestor path, and d^c_x(t) is the distance from x to c(t) **restricted to the patch P(N) of the node that owns c**. All ancestor curves must be scanned, not just the lowest one, because the true shortest path can leave the lowest common patch (paper Fig. 3) [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- What is stored for each vertex [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf):
  - the index of its leaf patch;
  - the indices of all ancestor separator curves;
  - for each such curve, a polynomial approximation of d^c_x(t): degree up to m_max = 12, fitted to either d or d², picking whichever is better;
  - optionally, the exact distances between vertex pairs inside the same leaf, stored as a sparse matrix D.
- In the experiments, leaves of V_max = 30 vertices were "sufficiently close to planar", so same-leaf pairs used **Euclidean distance** and D was not stored [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- The separator distances are not computed exactly at the separator samples [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf):
  - Single-source fields come from VTP (exact), fast DGG, fast marching or the heat method, run from the vertices of the triangles that contain separator samples.
  - The distance from a vertex to a separator sample inside a triangle is then interpolated with the barycentric formula of Rustamov et al. 2009, which is exact only in the plane.
- How the polynomial is fitted [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf):
  - weighted least squares that targets **average** relative error;
  - the degree is raised until the error drops below 0.5% or the degree reaches 12;
  - the authors note they "tacitly ignored" the fact that d^c_x is non-smooth where the separator meets the cut locus.
- The query minimises over t by brute force: the curve parameter is sampled at a spacing of 0.005 (mesh scaled to the unit cube) [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- Lemma A.3: if the fitted d̂^c_x has relative error at most ε for **every** t (a sup-norm condition), then the query has relative error at most ε. This is a conditional two-sided bound, not a per-query certificate. The proof assumes exact minimisation over continuous t and ignores the Euclidean fallback in the leaves [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- Experiments [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf):
  - only small meshes, with n = 900 to 5,051 for exact VTP compression;
  - mean relative error 0.1–0.4% and RMSE 0.2–1.3%;
  - 17–86 doubles stored per vertex;
  - "several thousand queries per second" in unoptimised MATLAB;
  - preprocessing took about 30 minutes for n = 20k with DGG or FM, and VTP was limited to about 5k vertices in the same time.
  - No maximum error is reported.
- Same paper, comparison with landmarks: about 10× more accurate than the landmark method of Burghard & Klein (LAND) at the same storage, which had ME 3.7–14.3% and RMSE 18.9–40.8% [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- From the conclusion, verbatim: "it is possible to partition the mesh using axis-parallel planes to cut through the bounding box of the mesh (similar to a 3D-tree construction)". This anticipates CCI's kd cuts as a design option [V] — [PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- NeuroGF describes SEP as handling point-to-point queries efficiently, "limited to query points located on mesh vertices", unable to produce paths, and possibly weak on anisotropic meshes [V] — [NeuroGF arXiv](https://arxiv.org/abs/2306.00658)
- I found no Gotsman–Hormann geodesic follow-up in 2023–2026. Hormann's later talks (2024, 2025) present the same work [V-abs] — [ZČU talk page](https://info.zcu.cz/clanek.jsp?id=6199); [Hormann publication list](https://www.inf.usi.ch/hormann/pub_chrono.html)

### Inferences
- Several CCI ingredients cannot be claimed as novel against SEP. The recursive separator hierarchy, labels from each vertex to its ancestor separators, the "min over shared separators of min over the cut parameter" query, exact near pairs, and even axis-parallel cuts are all anticipated, at least as an explicit suggestion.
- CCI differs from SEP on the points below; these are the defensible novelty axes.
  1. **Certified per-query interval.** CCI returns LB ≤ d ≤ UB with (1−ε)d ≤ LB and UB ≤ (1+ε)d, and the LB comes from 1-Lipschitz label rows along the cut. SEP returns only a point estimate and has no lower bound at all. Its Lemma A.3 rests on an unenforced sup-norm hypothesis: the fit criterion is average error with a capped degree, plus Rustamov interpolation, t-grid minimisation and a Euclidean fallback in the leaves.
  2. **Exact continuous labels.** CCI stores exact MMP/VTP distances computed on small ellipse sub-meshes. SEP stores lossy polynomial fits of exact or approximate fields, and its separator-point distances are interpolated.
  3. **Portal design.** CCI uses sparse portals inside a bounded window with geometrically growing spacing, which is what keeps the labels small while still certifying (1±ε). SEP samples the whole separator uniformly and compresses with polynomials.
  4. **Locality of preprocessing.** CCI runs exact sweeps on small ellipse sub-meshes. SEP computes fields over the whole patch, which costs O(f(n)√n) and, with VTP, limited it to about 5k vertices.
  5. **Scale and speed.** CCI answers in about 1 µs on large terrains, against SEP's thousands of queries per second on meshes of at most about 5k vertices (MATLAB). This is inferred from the two papers' numbers, not from a head-to-head run.
  6. **Domain.** CCI works on a 2.5D terrain chart with straight cuts. SEP targets general manifold meshes with curved, smoothed separators.
- **Recommendation (inference):** cite SEP prominently as the closest graphics prior art, and consider adding it as a baseline. Its algorithm is simple to reimplement on top of VTP.

### Gaps
- I did not find whether SEP has public code.
- I could not confirm whether any 2023–2026 paper extends SEP with bounds. Neither search nor the LiteGE and NeuroGF related-work sections turned one up.

## Q1. Xin, Ying, He, "Constant-time all-pairs geodesic distance query on triangle meshes" (I3D 2012): landmark-based, any bound?

### Takeaway
The method samples landmarks, precomputes all landmark-to-landmark distances, builds a coarse intrinsic Delaunay triangulation, and interpolates a query by unfolding quadrilaterals into the plane. Queries are O(1) with O(m² + n) storage. There is **no per-query bound**: Burghard & Klein state that it "has no bounds on the approximation error", and its output may not be a metric and can be discontinuous.

### Cited Findings
- Citation: S.-Q. Xin, X. Ying, Y. He, I3D 2012, pp. 31–38 [V-abs] — [dblp](https://dblp.uni-trier.de/rec/conf/si3d/XinYH12.html); [NTU record](https://libmast.utm.my/Record/ntu-10356%2F96690)
- Method, as described by Burghard & Klein: from evenly spaced landmarks, build a coarse intrinsic Delaunay triangulation, precompute pairwise landmark distances, and interpolate vertex distances by projecting quadrilaterals into the plane, since their side lengths and one diagonal are known [2nd, V] — [Burghard & Klein 2015](https://diglib.eg.org/handle/10.2312/vmv20151253)
- Same source: "They have no bounds on the approximation error and their approximated distance might not be a distance metric." The fields can also be non-continuous, because the coarse triangulation changes discontinuously [2nd, V] — [Burghard & Klein 2015](https://diglib.eg.org/handle/10.2312/vmv20151253)
- Gotsman & Hormann summarise it as O(m² + n) storage, O(1) query time and "O(1/√m) absolute accuracy", where m is the number of samples [2nd, V] — [Gotsman & Hormann 2022](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)

### Inferences
- XYH12 is a coarse-interpolation scheme with an asymptotic, absolute, non-certified accuracy. It differs from CCI in every respect that matters: no cuts or portals, no exact labels, and no interval.

### Gaps
- I did not read the I3D paper's full text. The ACM DL is paywalled, so I could not check whether it states any theorem of its own.

## Q2. Landmark / triangle-inequality bounds: Burghard & Klein, VMV 2015, the main (b)-type prior art

### Takeaway
Burghard & Klein give per-query **lower and upper bounds** from landmark distance fields:
- d_min(p,q) = max over r of |d(r,p) − d(r,q)|
- d_max(p,q) = min over r of [d(r,p) + d(r,q)]

The bounds are valid when the landmark fields are exact, and a query costs O(#landmarks). They prove an **absolute** gap bound, e_min + e_max ≤ min over r and z ∈ {p,q} of d(z,r), but **no worst-case (1±ε) relative guarantee**. Relative error is argued only locally or empirically, and in their experiments the fields came from the heat method, so the "bounds" were not certified against the exact polyhedral geodesic. This is the ALT idea transplanted to surfaces, and it is the graphics precedent for "certified LB/UB per query".

### Cited Findings
- Citation: O. Burghard, R. Klein, "Simple, Robust, Constant-Time Bounds on Surface Geodesic Distances using Point Landmarks", VMV 2015, pp. 17–24, DOI 10.2312/vmv.20151253 [V] — [EG DL](https://diglib.eg.org/handle/10.2312/vmv20151253); [PDF](https://diglib7.eg.org/bitstream/handle/10.2312/vmv20151253/017-024.pdf)
- The bounds d_min and d_max are defined as above (Def. 1). d_min is a pseudo-metric, and a metric if no point is equidistant from all landmarks (Thm 1) [V] — [PDF](https://diglib7.eg.org/bitstream/handle/10.2312/vmv20151253/017-024.pdf)
- Error results (Eq. 13 and Eqs. 16–17) [V] — [PDF](https://diglib7.eg.org/bitstream/handle/10.2312/vmv20151253/017-024.pdf):
  - absolute: e_max + e_min ≤ min over r and z of d(z,r);
  - relative, for d_min in the Euclidean plane: e_min/d ≤ 1 − cos α_r. On a smooth manifold this holds only "for some small neighborhood". Globally, the paper only says largest relative errors appeared locally in experiments, and there is no proven global ε. The upper bound's relative error is labelled "unbounded" (Fig. 10).
- Practical setup [V] — [PDF](https://diglib7.eg.org/bitstream/handle/10.2312/vmv20151253/017-024.pdf):
  - landmarks are chosen by farthest-point sampling;
  - landmark fields come from the heat method (Crane et al. 2013), and MMP gave "similar results";
  - storage is O(kn) and the query is O(k), where k = 30–1000 landmarks;
  - the figures show relative-error colour scales up to 35%.
- Gotsman & Hormann's experiment on SEP's meshes: LAND gives mean relative error 3.7–14.3% at a storage equal to SEP's [V] — [Gotsman & Hormann 2022](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- Gotsman & Hormann characterise LAND as "O(mn) storage with O(m) time queries and O(1/m) relative accuracy". This is their summary; B&K do not prove a worst-case relative rate [2nd, V] — [Gotsman & Hormann 2022](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)

### Inferences
- CCI's LB/UB idea has precedent in graphics through B&K's triangle-inequality landmark bounds, and in road networks through ALT. CCI's distinguishing claims are:
  1. the LB comes from separator portals, using 1-Lipschitz continuity of the label row along the cut, rather than from global landmarks;
  2. the interval is **guaranteed** to be within (1±ε)d for every query, whereas B&K certify only an absolute width tied to landmark spacing and give no relative guarantee for d_max;
  3. the labels are exact polyhedral geodesics, whereas B&K used the heat method in practice.

### Gaps
- None material.

## Q3. Graph-based exact or approximate methods (SVG, DGG, fast DGG, DGG-VTP) and exact single-source solvers (MMP, ICH, VTP, FWP)

### Takeaway
None of these is a pair-query oracle. Each needs a graph Dijkstra or a single-source sweep per query, so query time is well above O(1) per pair. DGG-type methods have an accuracy parameter ε. DGG (CAGD 2017) proves an O(ε) relative bound only under assumptions such as uniform vertex distribution. The DGG-VTP paper (CAD 2022) states that its ε is **not** a guaranteed bound. None returns a certified per-query interval.

### Cited Findings
- **SVG** (Ying, Wang, He, TOG 32(6):170, 2013): a sparse graph whose user parameter K trades size for accuracy; distances are a metric [V-abs] — [hgpu](https://hgpu.org/?p=11214); [researchr](https://researchr.org/publication/YingW013/bibtex)
- **DGG** (Wang, Fang, Wu, Xin, He, CAGD 52–53:262–284, 2017):
  - it criticises SVG because SVG cannot directly specify accuracy;
  - it proves that O(n/√ε) edges give relative error O(ε) under uniform vertex distribution;
  - observed error was below 0.6ε [V-abs]
  - Sources: [NTU PDF](https://personal.ntu.edu.sg/yhe/papers/2017-cagd_dgg.pdf); [Gotsman & Hormann ref 25](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- **Fast DGG** (Adikusuma, Fang, He, TOG 39(2):14, 2020): computes an approximate field in O(n) time with a tailored Dijkstra, and the user specifies the accuracy [2nd, V] — [Gotsman & Hormann 2022](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf); [NeuroGF](https://arxiv.org/abs/2306.00658)
- **DGG-VTP** (Adikusuma, Du, Fang, He, CAD 150:103333, 2022): "accuracy controllable and memory efficient". A search summary says its ε is not a guaranteed error bound; I did not read the paper itself [V-abs] — [NTU record](https://libmast.utm.my/Record/ntu-10356%2F161756/Details)
- NeuroGF on graph methods in general: results are "guaranteed to be a true distance metric and have controllable accuracy", but high accuracy needs large graphs [V] — [NeuroGF](https://arxiv.org/abs/2306.00658)
- **Geodesic Tracks** (Meng, Xin, Tu, Chen, He, Wang, TVCG 28(12), 2022): source-independent Steiner points plus a track partition. It is inefficient for point-to-point queries because it needs a Dijkstra-like sweep [V, via NeuroGF] — [NeuroGF](https://arxiv.org/abs/2306.00658)
- **Exact single-source solvers** [V, via Gotsman & Hormann and Burghard & Klein refs]:
  - MMP (Mitchell, Mount, Papadimitriou 1987); Surazhsky et al. (TOG 24(3), 2005); Xin & Wang ICH (2009); VTP (Qin, Han, Yu, Yu, Zhang, TOG 35(4):125, 2016).
  - These compute one exact field per source in about O(n²) worst case and are what CCI uses on its ellipse sub-meshes.
  - Sources: [G&H refs](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf); [Surazhsky et al. PDF](https://hhoppe.com/geodesics.pdf)
  - FWP (Xu et al. 2015) appears as [XWL*15] in B&K's references; its details were not verified here [V-ref only] — [B&K PDF](https://diglib7.eg.org/bitstream/handle/10.2312/vmv20151253/017-024.pdf)

### Inferences
- For a SIGMOD/VLDB audience, these are "online" methods. CCI uses them as subroutines and does not compete with them on pair queries, except with DGG as a pair-query baseline, which the authors already include.

### Gaps
- I did not verify the exact theorem statements in the SVG and DGG-VTP papers.

## Q4. Embedding / spectral / PDE methods (Panozzo et al., Rustamov, biharmonic and diffusion, heat method, edge flips, Kimmel spectral and Nyström)

### Takeaway
These produce smooth surrogate metrics or approximate fields: a Euclidean embedding queried in O(D), an O(n) field per source, or geodesic **paths** through edge flips. None stores distances to cut or separator samples, and none gives a worst-case or per-query bound against the exact polyhedral geodesic. They are not competitors for certified queries.

### Cited Findings
- **Panozzo, Baran, Diamanti, Sorkine-Hornung, "Weighted averages on surfaces"** (TOG 32(4):60, 2013) [V] — [PDF](https://cims.nyu.edu/gcl/papers/weighted-averages-siggraph-2013-panozzo-et-al.pdf); [dblp](https://dblp1.uni-trier.de/rec/journals/tog/PanozzoBDS13.html)
  - It decimates the mesh into a sample set S, computes **fast-marching** distances between samples, embeds them with metric MDS into R⁸ (stress minimisation), and extends to all vertices with least-squares meshes.
  - Euclidean distance in R⁸ is used as a geodesic-mimicking metric, and no error bound is stated.
  - It cites Rustamov et al.'s full all-pairs MMDS as not scaling ("just under 20 GB" for 100k vertices).
- **Rustamov, Lipman, Funkhouser, "Interior distance using barycentric coordinates"** (CGF 28(5), 2009): the barycentric interpolation of distances is "exact in the plane". Gotsman & Hormann reuse it [V, via G&H] — [G&H PDF](https://www.inf.usi.ch/hormann/papers/Gotsman.2022.CGI.pdf)
- **Diffusion, commute-time and biharmonic distances** (Lipman et al. 2010) are Euclidean-embedding metrics, different from geodesic distance by design [V, via Panozzo] — [Panozzo PDF](https://cims.nyu.edu/gcl/papers/weighted-averages-siggraph-2013-panozzo-et-al.pdf)
- **Heat method** (Crane, Weischedel, Wardetzky, TOG 32(5):152, 2013): O(n) per field after a one-time factorisation, a first-order approximation. NeuroGF reports degradation on anisotropic meshes [V, via B&K and NeuroGF] — [NeuroGF](https://arxiv.org/abs/2306.00658)
- **Sharp & Crane, "You can find geodesic paths in triangle meshes by just flipping edges"** (TOG 39(6):249, 2020): an iterative path method that needs an initial path. Its cost depends on the swept area, and it is not a distance oracle [V, via NeuroGF] — [NeuroGF](https://arxiv.org/abs/2306.00658)
- **Geodesic embedding with cascaded optimisation** (arXiv 2108.13821, TVCG version):
  - Euclidean embedding plus extra coordinates to reduce residual error;
  - bunny query time about 0.02 ms, mean relative error about 0.74%;
  - the embedding step takes minutes to hours, and there is no bound [V-abs] — [arXiv](https://arxiv.org/abs/2108.13821v1); [Cardiff repository](https://orca.cardiff.ac.uk/143769/3/NEWGeodesicEmbedding-TVCG.pdf)
- **Earth mover's distances on discrete surfaces** (Solomon et al.), **spectral MDS** (Aflalo & Kimmel), **Nyström / sparse approximations of geodesic matrices** (Shamai & Kimmel), and **geodesic distance descriptors**: from background knowledge, these approximate distance matrices for shape analysis without worst-case guarantees against the exact polyhedral geodesic [M]

### Inferences
- Embedding methods give O(D) queries with no bound, so they occupy the same "fast but uncertified" niche as NeuroGF.

### Gaps
- I did not fetch the Solomon, Aflalo–Kimmel and Shamai–Kimmel papers in this session. Their exact claims are unverified, so treat that line as [M].
- I could not confirm the authors of the cascaded geodesic embedding paper.

## Q5. Learned geodesic distance queries (NeuroGF, GeGnn, LiteGE, 2023–2026): any guarantees?

### Takeaway
None provides an error guarantee. They report mean relative error around 0.5–2.6%. NeuroGF explicitly admits that maximum relative error reaches about 5% and that there is "no strict guarantee".

### Cited Findings
- **NeuroGF** (Q. Zhang, J. Hou, Y. Y. Adikusuma, W. Wang, Y. He, NeurIPS 2023) [V] — [arXiv](https://arxiv.org/abs/2306.00658); [NeurIPS PDF](https://proceedings.neurips.cc/paper_files/paper/2023/file/3e22abb329d44080460b0eb11bf21da1-Paper-Conference.pdf)
  - An overfitted per-shape MLP. Ground truth comes from DGG-VTP and fDGG with ε = 1e-7.
  - Mean relative error below 0.5%, and a single-source all-destination (SSAD) query under 1 ms on GPU.
  - Max relative error "can rise up to approximately 5%"; the paper says "no strict guarantee".
- **GeGnn** (B. Pang, Z. Zheng, G. Wang, P.-S. Wang, "Learning the Geodesic Embedding with Graph Neural Networks", SIGGRAPH Asia 2023 / TOG, DOI 10.1145/3618317) [V-abs] — [arXiv](https://arxiv.org/abs/2309.05613v2); [GitHub](https://github.com/IntelligentGeometry/GeGnn)
  - A GNN embeds the mesh once, and a lightweight decoder maps pairs of embeddings to distances, giving constant-time queries.
  - Reported average relative error is below 1.4% on ShapeNet and 2.6% on unseen categories, empirical only.
- **LiteGE** (Y. Y. Adikusuma, Q. Huang, Y. He, arXiv 2512.17781, Dec 2025) [V] — [arXiv](https://arxiv.org/html/2512.17781)
  - Uses a PCA of unsigned distance fields plus an MLP.
  - About 1.9 ms per query against NeuroGF's 88 ms, and about 600K parameters.
  - Reports mean L1 error on normalised distances, with no max error and no guarantee stated.
  - Its baselines are NeuroGF and GeGnn; SEP and landmark methods are not mentioned.

### Inferences
- Learned methods are fast but uncertified, so CCI's certification is a clean differentiator. NeuroGF's own admission of about 5% max error is a quotable contrast.

### Gaps
- I did not find any 2024–2026 learned geodesic oracle with a provable bound.

## Q6. Separator hierarchies and bounded-error oracles outside graphics (context only)

### Takeaway
The (1+ε) separator-based oracles with **worst-case** guarantees come from computational geometry and graph algorithms, not graphics. Two examples are Thorup-style planar oracles and a 2026 polygonal-domain oracle. In graphics, Gotsman & Hormann is the separator work and Burghard & Klein the landmark-bounds work.

### Cited Findings
- J. van der Laan, F. Staals, L. Theunissen, "Approximate Dynamic Nearest Neighbor Searching in a Polygonal Domain" (arXiv 2603.11775, March 2026) [V] — [arXiv](https://arxiv.org/html/2603.11775v1)
  - It includes a (1+ε) two-point distance oracle for polygonal domains built on hierarchical shortest-path separators with anchor points: O((n/ε) log n) space and O((1/ε²) log n) query time.
  - It compares against Thorup's oracle.
  - It covers polygons, not terrains.
- Gotsman & Hormann also wrote "On landmark distances in polygons" (CGF 2021) [V-abs] — [NJIT record](https://digitalcommons.njit.edu/fac_pubs/3898)

### Inferences
- Terrain and polygon oracle prior art in computational geometry is presumably covered by another researcher. The graphics-side finding is that SEP plus B&K together anticipate "separator labels" and "per-query LB/UB" separately. **No graphics work combines them into a certified (1±ε) interval with exact labels.**

### Gaps
- The content of "On landmark distances in polygons" was not read.
