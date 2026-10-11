# A global error bound for the Separator Cut Index

> This note is about the earlier, uncertified `SeparatorCutIndex` on general manifolds. The terrain
> oracle with a proved $(1+\varepsilon)$ guarantee and per-query certificates is the Certified Cut Index;
> its theorem and proofs are in `notes/guarantee.md` (the window lemma of §2 here is its Lemma 2, and
> the tail spacing that resolves crossings outside the box is its gap rule (2)).

Everything the index returns is the length of an admissible curve, so the query is an
upper bound on the true distance. This note bounds how far above it can be, and ties the
bound to the three knobs of the construction: the portal spacing $h$, the label
(refinement) tolerance $\eta$, and the padding of the cut window.

## Setting

* $(M, g)$ a 2-D Riemannian manifold with obstacles removed; $d(\cdot,\cdot)$ its geodesic
  distance. $L(\gamma)$ is the Riemannian length of a curve.
* $\lambda_{\min}, \lambda_{\max}$: the smallest and largest metric speed over the domain,
  $\lambda_{\min}^2 = \min_p \lambda_{\text{min}}(g_p)$, so that
  $\lambda_{\min}\,|\gamma|_{\text{eucl}} \le L(\gamma) \le \lambda_{\max}\,|\gamma|_{\text{eucl}}$.
* A **cut** $C$ is a straight coordinate line $\{x = c\}$ (or $\{y = c\}$) restricted to a
  window $W \subset C$. A **region** $B$ is the bounding box of the members of a separator;
  the cut splits $B$ into two halves. $\operatorname{diam}_g(B) = \sup_{u,v\in B} d(u,v)$.
* **Portals** $x_1 < \dots < x_k$ on $W$ (ordered along the line). $h$ is the largest gap
  between consecutive portals (in Riemannian arc length along $C$) inside every obstacle-free
  run of $W$; runs always get their two endpoints as portals.
* **Labels** $\ell(u, x_j)$: the refined length of an admissible curve from node $u$ to
  portal $x_j$. Refinement never lengthens a curve, so
  $d(u, x_j) \le \ell(u, x_j) \le (1+\eta)\, d(u, x_j)$ for some $\eta \ge 0$ discussed below.
* **Query.** For nodes $s, t$ with lowest common separator $C$:
  $$\hat D(s,t) \;=\; \min\Big\{\,B_{\text{leaf}}(s,t),\; \min_j \big[\ell(s,x_j) + \ell(x_j,t)\big]\Big\}.$$

## 1. Soundness (lower bound)

$\ell(s,x_j) + \ell(x_j,t)$ is the length of the concatenation of two admissible curves
through $x_j$, and $B_{\text{leaf}}$ is a refined admissible curve, so

$$ d(s,t) \le \hat D(s,t). $$

The interpolated variant (section 5) is the one exception: a parabola vertex is not a
curve, so it is two-sided.

## 2. Where the shortest path crosses the cut (padding)

Let $\gamma$ be a shortest path from $s$ to $t$, both in $B$ and on opposite sides of the
cut. By the intermediate value theorem $\gamma$ meets the line of $C$ at some point $x^*$.
Since $s, t \in B$,
$$ d(s,t) = d(s,x^*) + d(x^*,t) \;\ge\; 2\, d(x^*, B), \qquad d(s,t) \le \operatorname{diam}_g(B), $$
so $d(x^*, B) \le \tfrac12 \operatorname{diam}_g(B)$. For $x^*$ on the line outside the
region interval $I = C \cap B$, the Euclidean distance from $x^*$ to $B$ is its coordinate
offset $c$ from $I$, and $d(x^*, B) \ge \lambda_{\min}\, c$. Hence

**Lemma (padding).** If the window $W$ extends $I$ by at least
$\operatorname{diam}_g(B) / (2\lambda_{\min})$ in coordinate on both sides, every shortest
path between two members crosses $C$ inside $W$.

On a terrain $g = I + \nabla f \nabla f^{\!\top}$ has eigenvalues $1$ and $1 + |\nabla f|^2$,
so $\lambda_{\min} = 1$ exactly and the padding is $\tfrac12 \operatorname{diam}_g(B)$ in
plain coordinates (`pad_speed=1`). The default implementation (`pad="diam"` without
`pad_speed`) inverts the arc length *along the cut* instead, i.e. it uses the cut's own
speed as a local $\lambda_{\min}$; that is exact when the metric is isotropic and no detour
off the cut is cheaper than the cut itself, and a heuristic otherwise.
$\operatorname{diam}_g(B)$ itself is estimated by the Riemannian length of the box diagonal.

## 3. Portal discretisation

Let $f(x) = d(s,x) + d(x,t)$ for $x$ on the free run of $W$ containing $x^*$; $f$ attains
its minimum $d(s,t)$ at $x^*$. Let $x_j$ be the portal nearest to $x^*$ along the run, so
its arc distance is at most $h/2$ (runs end in portals, so the nearest one is never across
an obstacle).

*Lipschitz bound.* $d(s,\cdot)$ and $d(\cdot,t)$ are 1-Lipschitz and the arc of $C$ between
$x^*$ and $x_j$ is an admissible curve, so
$$ f(x_j) - f(x^*) \le 2\, d(x^*, x_j) \le 2 \cdot \tfrac h2 = h. $$

*Second-order bound.* When $x^*$ is not on the cut locus of $s$ or $t$ and the
run is straight, $f$ is $C^2$ near $x^*$ with $f'(x^*) = 0$ (transversal or tangential
crossing alike, as $x^*$ is a minimum), so
$$ f(x_j) - f(x^*) \le \tfrac12 \Lambda \,(h/2)^2 = \tfrac{\Lambda h^2}{8},
   \qquad \Lambda = \sup_{|x - x^*| \le h/2} f''(x). $$
By the Hessian comparison theorem, if the sectional curvature is $\ge -\kappa^2$ then
$\nabla^2 d(s,\cdot) \le \kappa \coth(\kappa\, r) \le 1/r + \kappa$ with $r = d(s,x)$
(in the flat case $\nabla^2 d = (I - uu^{\!\top})/r$). Restricting to the unit tangent of the
cut (length $\lambda_C \le \lambda_{\max}$ in coordinates is irrelevant since we measure $h$
in arc length) gives
$$ \Lambda \;\le\; \frac{1}{r_s} + \frac{1}{r_t} + 2\kappa, \qquad r_s = d(s,x^*),\; r_t = d(x^*,t),\; r_s + r_t = d(s,t). $$

Together:

$$ \boxed{\;\delta_h \;:=\; \min_j f(x_j) - d(s,t) \;\le\; \min\Big\{\, h,\; \frac{h^2}{8}\Big(\frac{1}{r_s} + \frac{1}{r_t} + 2\kappa\Big) \Big\}\;} $$

The second term is the generic case ($O(h^2)$, matching the observed $\approx 2.5\times$
error drop per doubling of $k$ on the Poincaré disk); the first takes over when the crossing
is very close to one endpoint or on a cut-locus kink.

## 4. Labels and leaves

The labels are refined graph paths. Two sources of error:

* *Homotopy.* Refinement descends inside the homotopy class of the graph path. If the graph
  path is in the wrong class (obstacles, or a valley the graph skipped), the label converges to
  a locally shortest curve, not the shortest one. This is the only term that does not vanish
  with more resolution; it is bounded by the graph's own error $E_{\text{graph}}$, since the
  refined curve is never longer than the graph path.
* *Discretisation of the polyline.* A polyline with $m$ chart-straight segments through
  points of the geodesic $\gamma$ of length $L$: each segment costs at most
  $L(\text{chord}) - d(\text{ends})$, which is $O(|\partial g|^2 (L/m)^3)$ for a smooth
  metric (chords are geodesic to first order), so the total is $O(|\partial g|^2 L^3/m^2)$.
  On a polyhedral terrain a chord that crosses a face edge with bend angle $\theta$ exceeds
  the unfolded straight line by $O(\theta^2 \ell)$, so the total is $O(\theta_{\max}^2 L/m)$
  times the number of crossed edges — first order in $1/m$, which is why terrains want a
  larger `segments`.

Writing $\eta = \eta_{\text{hom}} + \eta_{\text{disc}}$ for the relative label error
(the optimiser tolerance adds a constant $\epsilon_{\text{opt}}$ to it),
$\ell(u,x) \le (1+\eta)\, d(u,x)$ for every label and every leaf-table entry (the leaf tables
are refined by the same routine, so same-leaf queries carry the error $\eta$ alone).

## 5. The bound

**Theorem.** For nodes $s,t$ with lowest common separator $C$, window padded as in §2,
portal gap $h$ and label error $\eta$,

$$ d(s,t) \;\le\; \hat D(s,t) \;\le\; (1+\eta)\,\big(d(s,t) + \delta_h\big), $$

i.e. the relative error satisfies

$$ \boxed{\; \varepsilon(s,t) \;\le\; \eta + (1+\eta)\,\frac{\delta_h}{d(s,t)}
   \;\le\; \eta + (1+\eta)\min\Big\{\frac{h}{d},\; \frac{h^2}{8d}\Big(\frac1{r_s}+\frac1{r_t}+2\kappa\Big)\Big\} \;} $$

*Proof.* Lower bound: §1. Upper bound: by §2 the crossing $x^*$ lies in $W$; take the
nearest portal $x_j$; then
$\hat D \le \ell(s,x_j) + \ell(x_j,t) \le (1+\eta)(d(s,x_j) + d(x_j,t)) = (1+\eta) f(x_j) \le (1+\eta)(d(s,t) + \delta_h)$ by §3. $\square$

With $k$ portals on a cut of Riemannian length $|C|$ we have $h \approx |C|/k$ and, for a
crossing near the middle of a path of length $d \sim |C|$, the generic term is
$\approx |C|^2/(2k^2 d^2) \sim 1/(2k^2)$: doubling $k$ quarters the portal term until the
label floor $\eta$ is reached. Budgeting $k \propto |C|$ (`scale_k`) makes $h$ uniform
across the tree, so one $h$ serves every level.

**Interpolation.** Fitting the parabola through $(x_{i-1}, f_{i-1}), (x_i, f_i), (x_{i+1}, f_{i+1})$
around the discrete argmin and reporting its vertex value is standard Lagrange interpolation:
if $f \in C^3$ near $x^*$, $|p(x) - f(x)| \le \tfrac{1}{9\sqrt3}\, \|f'''\|_\infty h^3$ on the
middle interval, so $|\min p - d(s,t)| \le 0.065\,\|f'''\|_\infty h^3 + 1.25\,\eta\, d$
(the Lebesgue constant of three equispaced nodes is $1.25$, which is how label noise enters).
This is two-sided; the guard in `cut_bound` caps the downward correction at a quarter of
the smaller neighbour gap, so the reported value is never below
$f_i - \tfrac14\min(f_{i-1}-f_i,\, f_{i+1}-f_i)$ and the worst case is no worse than the
discrete bound above.

## 5b. Interpolation when an endpoint is close to the cut (model interpolation)

The parabola fit assumes $f$ is smooth at the scale $h$. It is not when an endpoint sits within
a few $h$ of the cut: with $r_s = d(s, C)$ the drop and $x_s$ the foot,
$d(s,x) \approx \sqrt{r_s^2 + (x - x_s)^2}$ along the cut, whose curvature $1/r_s$ is of order
$1/h$ exactly there. Measured on HorseMount (k=32, 2000 pairs), *every* error above 1% had
$r_s$ or $r_t$ below a few portal gaps, at any $d$ (the worst: $d = 566$ m, $r_s = 12$ m; $d =
3063$ m, $r_s \approx 30$ m). So the sharp part is known in closed form and can be removed:

$$ g(x) = f(x) - \sqrt{r_s^2 + (x-x_s)^2} - \sqrt{r_t^2 + (x-x_t)^2} $$

is smooth (residuals of a few metres where $f$ varies by hundreds), gets the parabola through
the same three samples, and the model sum $m(x) = \sqrt{\cdot} + \sqrt{\cdot} + g(x)$ is minimised
in 1-D by bracketed Newton on $[x_{i-1}, x_{i+1}]$. $r_s, x_s$ are precomputed per (node, cut)
at build (one batched walk per cut, +8 bytes per label vector). When $r_s, r_t \gg h$ both
square roots are straight lines over the window and the model is the parabola.
HorseMount, k=32, 2000 pairs: mean 0.037% → 0.019%, p99 0.73% → 0.27%, max 4.2% → 1.2%;
k=16: max 11% → 2.2%. Query 9 → 20 µs. Error is now bounded by the smoothness of $g$
(second derivative of the *remainder*, not of $d(s,\cdot)$), which is $O(1/d)$ off the cut locus.

## 5c. Arbitrary points

The index stores labels at vertices. For a surface point $p$ in face $f$ with barycentric
weights $w$, two attachments:

* *Snap*: $\hat D(p,q) = \min_{a \in f, b \in f'} d_f(p,a) + \hat D(a,b) + d_{f'}(b,q)$ with $d_f$ the
  in-face straight distance. Admissible (an upper bound), but first order: a corner at $a$ costs
  $|p-a|(1-\cos\theta)$, tens of metres on 190 m edges, and exact vertex distances do not help
  (same 2% mean). Every terrain oracle we ran attaches points this way in effect: K-Algo's worst
  arbitrary-point case on HorseMount (160.8%) is the same pair as ours.
* *Interpolate*: $\hat D(p,q) = \sum_{a,b} w_a w_b\, \hat D(a,b)$. Off the cut locus the distance
  field is $C^2$ inside a face with Hessian $\lesssim 1/d$, so linear interpolation over a face of
  edge length $e$ errs by $O(e^2/d)$ per end (measured for a point source: 0.03% at 2–4 km,
  0.11% at 1–2 km, 0.5% below 1 km on HorseMount, $e \approx 190$ m). Two-sided. Pairs whose faces
  are within a few rings are answered exactly by one local sweep with both points inserted, since
  $e^2/d \to e$ there. Interpolating the *label vectors* first and taking one portal minimum would
  be nine times cheaper but has a Jensen gap, first order in the spread of the three vertices'
  argmins (0.8% vs 0.18%): the minimum must be taken per vertex pair.

Total arbitrary-point error $\lesssim \varepsilon_{V2V} + C\,e^2\,(1/d_p + 1/d_q)/d$ away from the cut
locus, with the near-pair radius $R$ chosen so that $e^2/R$ is below the target.

## 6. Where the bound is loose, and what it says to do

1. **Close pairs across a cut, and endpoints close to a cut.** For $s,t$ within a few node
   spacings but on opposite sides of the cut, $d \lesssim h$ and the Lipschitz term $h/d$ is
   $O(1)$; and for an endpoint within a few $h$ of the cut at *any* $d$, the parabola fails
   (§5b). On terrains the second case turned out to be the whole tail. Three remedies, all
   implemented (the third is §5b):
   * **Uniform $k$ across the tree** (`scale_k=False`). With the budget proportional to cut
     length the deep cuts get `k_min` portals, and the pairs they separate are close pairs;
     the generic term $\sim (h/d)^2$ is scale-free only if $h$ shrinks with the region, i.e.
     if every cut has the same $k$. Measured on BigMount (400 random pairs, exact labels,
     interpolated), error by depth of the lowest common separator:

     | depth | pairs | mean $d$ | $k\propto\lvert C\rvert$ (portals) | uniform $k=32$ |
     |---|---|---|---|---|
     | 0 | 183 | 4.7 km | 0.073% (32) | 0.073% |
     | 1 | 113 | 3.5 km | 0.124% (18) | 0.044% |
     | 2 | 56 | 2.5 km | 0.166% (16) | 0.016% |
     | 3 | 25 | 2.0 km | 0.277% (11) | 0.000% |
     | 4 | 12 | 1.2 km | 0.407% (9) | 0.020% |
     | 5–6 | 11 | 0.7 km | 1.1% (6) | 0.003% |

     Mean 0.150% → 0.049% for 4.5× the portals; the root level is untouched, as the bound says.

     **Revision (2026-10-11): uniform $k$ was the wrong fix, and it is not the default any more.**
     The table measures *interpolated* answers before the model interpolation and the halo existed.
     With both in place the proportional budget is nearly as accurate at a fraction of the cost, and
     uniform $k$ spends its extra portals on deep cuts whose pairs are better answered exactly:

     | GunnisonForest100k (100 pairs, exact labels, model interp.) | build | index | mean err |
     |---|---:|---:|---:|
     | uniform $k=32$ | 83,464 s | 252 MB | 0.063% |
     | $k\propto\lvert C\rvert$ ($k_{\min}=4$) | 16,545 s | 130 MB | 0.076% |
     | certified index, $\varepsilon=0.2$, $g=1$ edge (arc-length gap + near tables) | 6,854 s | 534 MB | 0.005% |

     | HorseMount (100 pairs) | build | index | mean err |
     |---|---:|---:|---:|
     | uniform $k=32$ + model interp. | 19.7 s | 1.14 MB | 0.008% |
     | $k\propto\lvert C\rvert$ + model interp. | 8.7 s | 0.85 MB | 0.016% |
     | certified, $\varepsilon=0.2$, $g=1$ edge | 2.8 s | 1.00 MB | 0.036% |
     | certified, $\varepsilon=0.1$, $g=0.5$ edge | 4.7 s | 1.34 MB | 0.011% |

     What the close pairs need is not more portals on every cut but an exact answer below a
     radius tied to the spacing: the certified index places portals at a fixed arc-length gap $g$
     (so the budget is proportional to cut length, like `scale_k`), and stores exact distances for
     every pair closer than $R = g/(\varepsilon\lambda_{\min})$; every other pair has $d > g/\varepsilon$, so the gap
     is at most $\varepsilon d$ at every depth (`notes/guarantee.md`, Lemma 4). Uniform $k$ is kept only as an
     option of the earlier index (`scale_k=False`); the benchmark's default for it is now
     `scale_k=True`.
   * **A halo table** (`halo_rings`): for every node, distances to all vertices within $r$
     hops regardless of the tree, answered before the portal minimum. With exact labels it
     is free (same sweeps, more targets); HorseMount, uniform $k=32$: max error 3.2% → 1.4%,
     mean 0.047% → 0.032%, +0.4 MB. In the refined pipeline it would cost $O(n\,\deg^r)$
     batched refines.
2. **Label floor $\eta$.** Once $h$ is small the error is all $\eta$: more portals cannot
   help (measured: k=48 on Poincaré gives 0.097% vs the 0.06% of query-time refinement).
   Lower it with more `segments` (quadratic gain on smooth metrics, linear on terrains) or a
   tighter optimiser tolerance, both linear in build time.
3. **Padding.** The provable window needs $\lambda_{\min}$ (`pad_speed`); without it the
   arc-length window is a heuristic that can miss a crossing when a detour off the cut is
   cheaper than the cut itself. On terrains $\lambda_{\min} = 1$, so the proof applies.
   The same argument bounds where a geodesic between two points of a box $B$ can go
   (within $\lambda_{\max}(B)\cdot\operatorname{diag}(B)$ of the box), which is what lets
   the exact terrain labels be computed on a sub-mesh per cut without losing exactness
   (`ExactTerrain.region_geo`).
4. **Curvature constant $\kappa$.** Enters only through $2\kappa h^2/8$; on the Poincaré
   disk $\kappa = 1$ and $h \ll 1$, so it is negligible next to $1/r_s + 1/r_t$.
5. **Homotopy term $\eta_{\text{hom}}$ is not small on weighted terrains.** With face costs in
   $(1, 5)$ the weighted geodesic refracts at every edge (Snell) and hugs cheap edges; the vertex
   mesh graph, even with edges charged at the cheaper face, routes 9% longer than a Steiner
   graph with 5 points per edge (HorseMount, 100 pairs, reference = best path found by any
   method), with one pair at 3x. Refinement descends inside the graph path's corridor and
   cannot cross to the right one, so the labels inherit the 9% (8.6% after refinement).
   Portals at the mesh-edge crossings of the cut (`edge_portals`) change 1.26% → 1.14% against
   our own heavy refinement and nothing against theirs. The bound is not wrong — it says
   $\eta_{\text{hom}} \le E_{\text{graph}}$ — the graph is. The remedy is the competitors' own
   base graph (Steiner points on edges) under the same index.
