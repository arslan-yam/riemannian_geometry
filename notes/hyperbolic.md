# Side project: certified distance oracles on hyperbolic domains with obstacles

This is a separate idea from the terrain paper. It is not implemented, and the empirical check of its
core claim is running separately (`notes/hyperbolic_h2h.md`, written by that session).

## 1. Why hyperbolic geometry changes the exponent

The terrain index stores about $n^{1.5}$ labels because a balanced separator of a flat region with $n$
samples has length $\Theta(\sqrt n)$ in sample spacings. Terrains cannot escape this: planar distance
is a lower bound on terrain distance, so a terrain disk of radius $r$ lies over a planar disk of
radius $r$ and has area at most $\lambda_{\max}\pi r^2$. Local saddles do not help.

The hyperbolic plane $\mathbb H^2$ (curvature $-1$) is different. A disk of radius $R$ has area
$2\pi(\cosh R - 1) \approx \pi e^R$, but a geodesic through its centre has length $2R$. With $n$ samples
at fixed spacing, $R \approx \ln n$. So a balanced geodesic separator holds $O(\log n)$ samples, and
recursing gives $O(\log n)$ levels of $O(\log n)$ portals each.

**Known results to build on** (to read in full before writing):
- Kisfaludi-Bak (SODA 2020, arXiv:1812.03960): a new separator theorem for ball graphs in $\mathbb H^d$.
  They behave like Euclidean ones one dimension lower, i.e. $n^{1-1/(d-1)}$ for $d\ge3$. In the
  hyperbolic plane, noisy uniform disk graphs get small treewidth and quasi-polynomial algorithms.
- Bläsius, von der Heydt, Kisfaludi-Bak, Wilhelm, van Wordragen (SoCG 2025, arXiv:2407.09362): disk
  graphs of radius $r$ in $\mathbb H^2$ have balanced separators of $((1+1/r)\log n)$ cliques. For $n$
  points at pairwise distance $\ge 2r$, the Delaunay complex has outerplanarity $1+O(\log n/r)$, hence
  $O(\log n)$-type separators and treewidth.
- Kisfaludi-Bak, van Wordragen (SoCG 2024, arXiv:2305.01356): a hyperbolic quadtree, an L-order, and a
  Steiner $(1+\varepsilon)$-spanner with $O(n)$ edges. Without Steiner points, near-linear $(1+\varepsilon)$
  spanners are impossible in $\mathbb H^d$. This is for point sets *without obstacles*, where distances
  have a closed form.
- A quick search found no distance oracle for hyperbolic domains with obstacles. The nearest hits are
  Euclidean polygonal-domain oracles (Gudmundsson et al.; Arikati et al., ESA 1996) and work on
  quasihyperbolic metrics.

**The effect is two-dimensional.** In $\mathbb H^3$ separators are again $\sqrt n$.

**Topology matters.** Negative curvature alone is not enough. Triply periodic minimal surfaces have
$K<0$ everywhere but huge genus, fill space, and their area grows like $r^3$, so their separators are
about $n^{2/3}$, worse than terrains. The setting needs exponential area growth: a domain of $\mathbb H^2$,
or a simply connected surface with $K \le -\kappa < 0$.

## 2. Setting

$\Omega = \mathbb H^2 \cap$ (a disk of radius $R_0$) minus obstacles. Obstacles are polygons with geodesic
edges, or are triangulated if curved. Points and queries are at fixed spacing $h$, and $n \approx$ area$/h^2$.
Distances are the shortest-path metric of $\Omega$; without obstacles they would have the closed form
$\cosh d = \cosh r_1\cosh r_2 - \sinh r_1 \sinh r_2 \cos\Delta\theta$.

## 3. The index, transferred

* **Chart: the Beltrami–Klein model.** Geodesics are straight chords, so a straight cut is a geodesic,
  and visibility among geodesic polygons is Euclidean visibility in Klein coordinates. The distance is
  $\cosh d(x,y) = (1-\langle x,y\rangle)/\sqrt{(1-|x|^2)(1-|y|^2)}$.
* **Exact labels.** Shortest paths among geodesic polygons are taut: they are geodesic pieces that bend
  only at obstacle vertices (uniquely geodesic, CAT(−1)). A label sweep is therefore Dijkstra on the
  visibility graph of obstacle vertices, with each target connected by visibility. This is exact and
  replaces MMP/VTP; for triangulated negatively curved surfaces MMP/VTP still apply unchanged.
* **Cuts.** Use centre-based recursion rather than median cuts in Klein coordinates, which compress
  near the rim:
  - split the disk by a diameter (length $2R$);
  - then by the perpendicular radii (length $R$);
  - then bisect sectors from the centre, or use the hyperbolic quadtree.

  Every cut has hyperbolic length $\le 2R = O(\log n)$. A cut through an obstacle is only shorter and
  still separates.
* **Window lemma.** It gets stronger. For a point at distance $a$ from the cut line, a crossing at offset
  $c$ from its foot costs $\cosh d = \cosh a\cosh c$, so tails thin out exponentially.
* **Near tables, certificates, certify-or-refine, coding.** Unchanged. The 1-Lipschitz argument along the
  cut uses only the triangle inequality.

## 4. Theorem candidate

> **Superseded (18 Oct 2026).** Section 7 shows that the $(1+\varepsilon)$ version below gains nothing
> over flat oracles, and replaces it with an exact two-point theorem. The text below is kept for the record.

**Claim.** For a domain of $\mathbb H^2$ with geodesic-polygon obstacles sampled at fixed spacing, the index
answers certified $(1\pm\varepsilon)$ queries using **exact full label rows**, with
$O(n\log^2 n\cdot c(\varepsilon))$ space and $O(\log n\cdot c(\varepsilon))$ query time.

Proof outline:
- separator length $O(R) = O(\log n)$ per level, over $O(\log n)$ levels;
- the portal count per cut is (length)/$g$ plus tails;
- the window lemma and the $(1+\varepsilon)$ argument carry over with $\lambda = 1$ in the hyperbolic
  metric itself; the chart lower bound becomes the hyperbolic distance.

**Why "full rows" is the point.** Near-linear space is already known in flat settings:
- Thorup-style covers for planar graphs;
- Djidjev–Sommer for polyhedral surfaces;
- a Steiner spanner (Kisfaludi-Bak–van Wordragen) followed by a planar oracle would also give it in
  $\mathbb H^2$, through covers.

All of these pay in accuracy. Our cover experiment at 100k vertices measured 0.79% mean error against
0.0075% for 8-bit full rows at the same size. In $\mathbb H^2$ the geometry gives near-linear space
without covers. The paper's comparison is therefore: (Steiner spanner + Thorup) versus (geodesic cuts +
full rows), on accuracy at equal size.

## 5. Applications (the weak spot)

A theory venue (SoCG) needs none. Candidates otherwise:
- hyperbolic maps of networks with forbidden regions;
- motion planning in hyperbolic configuration spaces;
- simply connected surfaces with $K\le -\kappa$, such as saddle-shaped shells and hyperbolic paraboloid
  roofs (these are bounded, so the regime $R\gg1$ is unclear).

## 6. Plan

1. **[done]** Empirical check of H2H separator and label growth on hyperbolic-uniform samples with
   $R$ up to about 12, flat versus hyperbolic, with obstacles (`notes/hyperbolic_h2h.md`). It holds.
2. **[done, 18 Oct 2026]** Literature and theorem check: section 7. The plan from here is in 7.6;
   steps 3–4 below apply only if the prototype in 7.6 holds up.
3. Implement in C++ by reusing the cut machinery: Klein coordinates, visibility-graph labels, and
   centre-based cuts.
4. Compare against a Steiner spanner + Dijkstra (exact on the spanner) and the closed-form distance (no
   obstacles), for sanity.

## 7. Step 1: literature and theorem check (18 Oct 2026)

### 7.1 What is known

| work | setting | result | relevance |
|---|---|---|---|
| Thorup, ESA 2007 | flat plane, polygonal obstacles, $m$ corners, any two points | $(1+\varepsilon)$: $O(m\log m/\varepsilon)$ space, $O(1/\varepsilon^3 + \log m/(\varepsilon\log\log m))$ query | the $(1+\varepsilon)$ flat problem is already near-linear |
| Gudmundsson et al. 2008 | flat, $t$-rounded domains | $(1+\varepsilon)$: $O(m\log m)$ space, $O(\log m)$ query | the same, under a bounded-detour assumption |
| Chiang–Mitchell 1999; de Berg–Miltzow–Staals, SoCG 2024 | flat, **exact** two-point | $O(m^{11})$ → $O(m^{10+\varepsilon})$ space for $O(\log m)$ query; $O(m^2)$ space with $O(h\log m)$ query (Guo–Maheshwari–Sack) | exact queries are very expensive |
| de Berg–Miltzow–Staals, Thm 5 | flat, exact | polylog query needs $\Omega(m^2)$ space, conditional on Integer Hopcroft | the hard instance uses tiny holes and long sightlines across a square |
| Mitchell, Handbook ch. 31, open problem 2 | flat, exact | asks for sublinear queries in subquadratic space | open |
| Kisfaludi-Bak 2020; Bläsius et al. 2025 | $\mathbb H^2$ ball and disk graphs | $O(\log n)$ separators and treewidth | gives the separators below |
| Kisfaludi-Bak et al., SoCG 2024 (planar δ-hyperbolic graphs) | Gromov-hyperbolic planar graphs | geodesic separators, an MIS approximation scheme, TSP | no oracles, labels or obstacles; open problems (§8) do not mention them |
| Kisfaludi-Bak–van Wordragen, SoCG 2024 | point sets in $\mathbb H^d$, no obstacles | Steiner $(1+\varepsilon)$-spanner with $n\,d^{O(d)}\log(1/\varepsilon)/\varepsilon^d$ edges; approximate nearest neighbour | no distance oracle, no obstacles |
| Kisfaludi-Bak–Poon–van Wordragen 2025 | $\{p,q\}$ tilings | treewidth $O(\log n)$ of convex hulls; exact graph shortest paths | no oracles, no obstacles |
| Kopczyński–Celińska-Kopczyńska 2021 | hyperbolic tilings | exact distances in $O(\text{dist})$ time; dynamic distance counting | a graph setting, no obstacles |
| Chepoi–Maftuleac; Lubiw–Maftuleac–Owen | CAT(0) planar complexes | two-point queries, $O(n^2)$ space, $O(d(p,q))$ query | the closest "nonpositive curvature + queries" work |

Not found: any shortest-path algorithm or distance oracle for obstacles in $\mathbb H^2$. Not checked
(paywalled or rate-limited): Gavoille–Ly 2005 and Chepoi et al. 2008, which give additive-error labels for
δ-hyperbolic graphs, and the body of Thorup 2007.

### 7.2 The $(1+\varepsilon)$ version gains nothing

- Thorup already gives near-linear $(1+\varepsilon)$ space in the flat plane, for any two points.
  Section 4's argument ("near-linear without covers") therefore has no flat baseline to beat.
- **Portal count per cut is the same in both geometries.** Take a point at distance $a$ from a cut. A
  crossing at offset $c$ needs portal spacing about $\varepsilon\,d(c)$. Flat: $d \approx \sqrt{a^2+c^2}$.
  Hyperbolic: $\cosh d = \cosh a\cosh c$, so $d \approx a + c$ for large $c$. Either way
  $\int dc/(\varepsilon d(c)) = O(\varepsilon^{-1}\log(L/a))$.
- **Additive-error labels give the same bound.** With full rows at spacing $\delta$, $\mathbb H^2$ costs
  $O(R^2/\delta)$ per point. A flat $(1+\varepsilon)$ oracle with $\varepsilon = \delta/\mathrm{diam}$ costs
  $O(R\log n/\delta)$ in a region of the same diameter.
- So the hyperbolic advantage measured in `notes/hyperbolic_h2h.md` only applies to **exact** labels.

### 7.3 New candidate: exact two-point queries

**Setting.**
- $\Omega$ is the hyperbolic disc $D_R$ minus disjoint open polygons with geodesic edges. The outer boundary
  is a geodesic polygon.
- $m$ is the number of corners.
- $G$ is the tangent visibility graph on the corners: free geodesic segments between corners that can lie
  on a shortest path. Shortest corner-to-corner paths in $\Omega$ are shortest paths in $G$.

**Assumptions** (constants $c$, $\rho$):
- **(D) density:** every hyperbolic disc of radius 1 contains at most $c$ corners.
- **(V) bounded visibility:** every free geodesic segment has length at most $\rho$.

**Consequences:**
- every edge of $G$ has length at most $\rho$;
- every free point sees $k = O(c\,e^{\rho})$ corners;
- every geodesic segment of length 1 is crossed by $\kappa = O(c^2 e^{2\rho})$ edges of $G$.

**Theorem (candidate).** Under (D) and (V) there is a structure with:
- $O(\kappa\, m R\log m)$ space in words (real RAM);
- $O(k\kappa\, m R\log^2 m)$ preprocessing, up to log factors;
- exact $d_\Omega(p,q)$ for **any** $p,q\in\Omega$ in $O(\log m + k^2\kappa R\log m)$ time.

When obstacle edges have bounded length, (V) forces $m = \Omega(e^{R-\rho})$, so $R = O(\log m)$. Space is
then $O(m\log^2 m)$ and query time $O(\log^2 m)$.

**Proof sketch.**
1. **Radial cuts.** Split $D_R$ by a diameter into two half-discs. Split each wedge (angle $\le\pi$, hence
   convex) by the radius at the median angle of its corners. Every cut has length $R$. An edge of $G$ between
   two child wedges crosses their common radius, so the endpoints of the crossing edges form a separator of
   at most $\kappa R + O(1)$ corners. The corner count halves at each level, giving depth $\log m$.
2. **Labels.** Each corner stores, for every ancestor piece $P$ and every $s\in S(P)$, the distance $d_P(u,s)$
   inside the subgraph of $G$ induced by $P$. That is $O(\kappa R\log m)$ entries per corner. Correctness:
   take the deepest piece $P$ that contains a whole shortest $u$–$w$ path. Then the path meets $S(P)$, so
   $d(u,w) = \min_{P,\,s} d_P(u,s) + d_P(s,w)$ over the common ancestor pieces. Preprocessing is one Dijkstra
   per separator vertex inside its own piece.
3. **Arbitrary query points.** A shortest $p$–$q$ path is either the segment $pq$ or bends only at corners.
   Its first bend is visible from $p$ and its last bend from $q$. So
   $d(p,q) = \min\big(|pq|\text{ if free},\ \min_{u\in V(p),\,w\in V(q)} |pu| + d_G(u,w) + |wq|\big)$.
   This is the corner-SPM idea of Chen–Daescu–Klenk, but with labels in place of $O(m^2)$ space.
   - $V(p)$: point location, then a local scan, costing $O(\log m + k)$.
   - Is $pq$ free? If $|pq| > \rho$, no. Otherwise test the $O(k)$ obstacle edges near $p$.
4. **Klein model.** Geodesics are straight chords, so visibility, point location and triangulation are the
   Euclidean algorithms in Klein coordinates.

**The flat counterpart is a result too.**
- With balanced straight cuts of length $O(\sqrt{\text{area}})$, the same proof gives exact two-point
  queries in $O(\kappa\, m^{1.5})$ space with $O(k^2\kappa\sqrt m)$ query time.
- That is subquadratic space and sublinear query: Mitchell's open problem 2, restricted to inputs
  satisfying (D) and (V).
- It does not contradict the conditional $\Omega(m^2)$ bound, because that instance violates (V).
- Chen–Daescu–Klenk give $O(m^2)$ space and $O(\log m)$ query under (V), so the flat result is a trade-off,
  not an improvement on every axis.

**Lower bounds to prove** (they would make it a separation rather than two upper bounds):
- **Flat:** any exact labeling scheme under (D) and (V) needs $\Omega(\sqrt m)$ reals per label.
  - Idea: a corridor maze that realises a weighted $\sqrt m\times\sqrt m$ grid.
  - The left-to-right distance matrix should have Jacobian rank $\Omega(m)$ in the corridor lengths.
  - The $2\sqrt m$ labels on the two sides must determine it.
  - This is for labels only. Lower bounds for centralised oracles are open even in general.
- **Hyperbolic:** $\Omega(\log m)$ by the same counting on tilings. That would leave a $\log m$ gap.

### 7.4 Main risks

1. **Random obstacles violate (V) at worst case.** Under a Poisson model the longest free segment is
   $\Theta(\log m)$, so $e^{\rho}$ is polynomial in $m$ and the worst-case constants are useless.
   - What is needed is an expected-case version: expected crossings per unit of cut, and expected $|V(p)|$.
   - In $\mathbb H^2$ the number of corners within distance $\ell$ grows like $e^{\ell}$. A segment of length
     $\ell$ is free with probability $e^{-\lambda\ell}$, so the expected number of crossings is finite only
     when $\lambda > 1$.
   - This should be the visibility phase transition of Benjamini–Jonasson–Schramm–Tykesson, "Visibility to
     infinity in the hyperbolic plane, despite obstacles" (2009, to verify). Above a critical obstacle
     intensity, sightlines are short and the theorem applies on average. Below it, long sightlines exist and
     labels should grow.
   - That is a clean experiment, and a nice structural point: flat geometry has no such transition.
2. **"Separators plus labels" is standard.** The new parts are:
   - the model;
   - the reduction for arbitrary query points;
   - the radial-cut lemma;
   - the lower bounds;
   - the random-obstacle analysis.
   
   Reviewers may still find the upper bound easy.
3. **Real RAM.** Labels store exact reals. Fine for computational geometry; it should be stated.
4. **Realistic-input literature not fully checked.** One search found nothing on exact two-point queries
   under fatness, low density or bounded visibility. Still to read: de Berg–Katz–van der Stappen–Vleugels on
   realistic input models, and van der Stappen on low obstacle density.

### 7.5 Novelty and venue estimate (revised)

| version | venue | estimate |
|---|---|---|
| upper bounds only (flat $m^{1.5}$, hyperbolic $m\log^2 m$, exact) | EuroCG or CCCG | likely |
| plus the flat $\Omega(\sqrt m)$ labeling lower bound and the random-obstacle phase transition (theory + experiment) | ESA, ISAAC or SWAT | about 30–40% |
| SoCG | — | about 10–15%; needs a stronger lower bound or (V) removed |

Database venues are out. The terrain paper is unaffected.

### 7.6 Prototype (step 2) — done 19 Oct 2026, results in `notes/hyperbolic_exact.md`: passes the stop/go

Reuse the sampler and Klein geometry in `bench/hyperbolic_h2h.py` and `H2HIndex`:

1. **Obstacles.** Small geodesic polygons (triangles and squares of hyperbolic radius $r_o$) with Poisson
   centres in $D_R$, at several intensities $\lambda$. Flat discs with matched $m$.
2. **Tangent visibility graph** in Klein coordinates: Euclidean visibility, using a grid or KD-tree of
   obstacle edges.
3. **Measure:**
   - the longest free segment, $|V(p)|$ and crossings per unit length along radii;
   - the radial-cut separator sizes and label length against $R$, at each $\lambda$;
   - the same for the flat discs.
4. **Validate** exact two-point answers against Dijkstra on $G\cup\{p,q\}$, for random free $p,q$.
5. **Stop/go:** labels about $\log^2 m$ above the critical $\lambda$, and visibly worse below it.

## Sources

- Kisfaludi-Bak, SODA 2020: https://arxiv.org/abs/1812.03960v2
- Bläsius et al., SoCG 2025: https://arxiv.org/html/2407.09362v2
- Kisfaludi-Bak & van Wordragen, SoCG 2024: https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2024.68
- Euclidean polygonal-domain oracles (Gudmundsson, Levcopoulos, Narasimhan, Smid): https://people.scs.carleton.ca/~michiel/distoraclespanner.pdf
- Quasihyperbolic shortest paths in planar domains: https://www.utupub.fi/handle/11111/61700
- Thorup, Compact oracles for approximate distances around obstacles in the plane, ESA 2007, LNCS 4698: https://link.springer.com/chapter/10.1007/978-3-540-75520-3_35
- Mitchell, Shortest paths and networks, Handbook of DCG ch. 31 (2017): https://www.csun.edu/~ctoth/Handbook/chap31.pdf
- de Berg, Miltzow, Staals, Towards space efficient two-point shortest path queries in a polygonal domain, SoCG 2024: https://arxiv.org/abs/2303.00666
- Kisfaludi-Bak, Masaříková, van Leeuwen, Walczak, Węgrzycki, Separator theorem and algorithms for planar hyperbolic graphs, SoCG 2024: https://arxiv.org/abs/2310.11283
- Kisfaludi-Bak, Poon, van Wordragen, Shortest paths, convexity, and treewidth in regular hyperbolic tilings (2025): https://arxiv.org/abs/2510.26110
- Kisfaludi-Bak & van Wordragen, quadtree and Steiner spanner in hyperbolic space: https://arxiv.org/abs/2305.01356
- Kopczyński & Celińska-Kopczyńska, Dynamic distances in hyperbolic graphs (2021): https://arxiv.org/abs/2111.01019
- Chepoi & Maftuleac, Shortest path problem in rectangular complexes of global nonpositive curvature: https://arxiv.org/abs/1010.0852
- Lubiw, Maftuleac, Owen, Shortest paths and convex hulls in 2D complexes with non-positive curvature: https://arxiv.org/abs/1603.00847
