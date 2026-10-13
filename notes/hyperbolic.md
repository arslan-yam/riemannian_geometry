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

1. **[running]** Empirical check of H2H separator and label growth on hyperbolic-uniform samples with
   $R$ up to about 12, flat versus hyperbolic, with obstacles (`notes/hyperbolic_h2h.md`).
2. If it holds: write the proof (sections 3–4) and read the three papers above for overlap.
3. Implement in C++ by reusing the cut machinery: Klein coordinates, visibility-graph labels, and
   centre-based cuts.
4. Compare against a Steiner spanner + Dijkstra (exact on the spanner) and the closed-form distance (no
   obstacles), for sanity.

## Sources

- Kisfaludi-Bak, SODA 2020: https://arxiv.org/abs/1812.03960v2
- Bläsius et al., SoCG 2025: https://arxiv.org/html/2407.09362v2
- Kisfaludi-Bak & van Wordragen, SoCG 2024: https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2024.68
- Euclidean polygonal-domain oracles (Gudmundsson, Levcopoulos, Narasimhan, Smid): https://people.scs.carleton.ca/~michiel/distoraclespanner.pdf
- Quasihyperbolic shortest paths in planar domains: https://www.utupub.fi/handle/11111/61700
