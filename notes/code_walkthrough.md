# Code walk-through

This document explains every module of the repository: what problem it solves, the mathematics it
implements, the algorithm in the order the code executes it, its cost, and the properties the rest of
the code relies on. Section numbers follow the dependency order (each section only uses earlier ones).
Notation is fixed in §0; file paths are relative to the repository root.

---------------------------------------------------------------------------------------------------

## 0. The problem and the notation

**Manifold.** A 2-D Riemannian manifold is given in one chart: a domain $U \subset \mathbb R^2$ and a
metric field $g : U \to \mathrm{SPD}(2)$. Obstacles are closed sets $O_1, \dots, O_m \subset U$; the free
space is $U_{\text{free}} = U \setminus \bigcup O_i$.

**Length and distance.** For a piecewise-$C^1$ curve $\gamma : [0,1] \to U_{\text{free}}$,
$$ L(\gamma) = \int_0^1 \sqrt{\dot\gamma(t)^{\top} g(\gamma(t))\, \dot\gamma(t)}\; dt, \qquad
   d(p,q) = \inf\{ L(\gamma) : \gamma(0)=p,\ \gamma(1)=q \}. $$
$\lambda_{\min}(p), \lambda_{\max}(p)$ denote the square roots of the eigenvalues of $g(p)$, so
$\lambda_{\min}|v| \le |v|_g \le \lambda_{\max}|v|$.

**Energy.** For the same curve $E(\gamma) = \int_0^1 |\dot\gamma|_g^2\,dt$. By Cauchy–Schwarz
$L(\gamma)^2 \le E(\gamma)$, with equality iff $\gamma$ has constant speed; minimisers of $E$ are
constant-speed geodesics. The code always minimises a discrete energy and reports a length.

**Oracle.** A finite node set $V \subset U_{\text{free}}$, $n = |V|$, is preprocessed into an index that
answers $\hat D(s,t) \approx d(s,t)$ for $s, t \in V$ (and, on terrains, for arbitrary surface
points). An answer is *admissible* if it is the length of an actual free curve, hence $\hat D \ge d$.

**Graph.** $G = (V, E, w)$ with $w(u,v) = L(\overline{uv})$, the length of the chart-straight segment.
$d_G$ is the graph distance. `adj` is always a list of dicts, `adj[u][v] = w(u,v)`.

---------------------------------------------------------------------------------------------------

## 1. `rgeo/manifolds.py` — metrics, lengths, energies

### 1.1 The `Manifold` interface

| member | meaning |
|---|---|
| `metric(p)`, `metric_batch(P)` | $g(p)$, $(N,2,2)$ for $N$ points |
| `in_domain(p)`, `in_domain_batch(P)` | $p \in U$ |
| `segment_lenght(p, q, n=4)` | $L(\overline{pq})$ by $n$-point Gauss–Legendre quadrature |
| `metric_grad_batch(P)` | $(\partial_x g, \partial_y g)$; central differences by default, analytic in subclasses |
| `energy`, `energy_and_grad_points(P, m)` | discrete energy and its gradient (below) |
| `geodesic_path(p, q, m)` | straight-line initialisation + L-BFGS on the discrete energy (no obstacles) |

**Segment length.** With $\Delta = q - p$ and Gauss–Legendre nodes/weights $(t_i, \omega_i)$ on
$[0,1]$: $L(\overline{pq}) \approx \sum_i \omega_i \sqrt{\Delta^\top g(p + t_i \Delta)\, \Delta}$. Four
nodes integrate polynomials of degree 7 exactly; for smooth metrics and short segments the error is far
below every error measured in the project. (On terrains this is replaced by an exact face walk, §8.2.)

**Discrete energy.** For a polyline $P_0, \dots, P_m$ with $\Delta_i = P_{i+1} - P_i$ and midpoints
$c_i = \tfrac12(P_i + P_{i+1})$,
$$ E(P) = m \sum_{i=0}^{m-1} \Delta_i^\top g(c_i)\, \Delta_i \;\approx\; m \sum_i L_i^2 \;\ge\; \Big(\sum_i L_i\Big)^2 , $$
the last step by Cauchy–Schwarz, with equality iff all $L_i$ are equal. Minimising $E$ over the
interior vertices therefore both shortens the polyline and equalises its pieces. Its gradient with
respect to $P_{i+1}$ and $P_i$ (one term per segment, accumulated):
$$ \partial_{P_{i+1}} = m\big(2 g(c_i)\Delta_i + \tfrac12 \nabla_c(\Delta_i^\top g\, \Delta_i)\big), \qquad
   \partial_{P_{i}} = m\big(-2 g(c_i)\Delta_i + \tfrac12 \nabla_c(\Delta_i^\top g\, \Delta_i)\big), $$
with $\nabla_c(\Delta^\top g \Delta) = (\Delta^\top \partial_x g\, \Delta,\ \Delta^\top \partial_y g\, \Delta)$.
`energy_and_grad_points` computes both from one metric evaluation.

### 1.2 The test manifolds

* **`Euclidean`**: $g = I$.
* **`PoincareDisk(rmax)`**: $g = \lambda^2 I$, $\lambda = 2/(1-|p|^2)$, Gaussian curvature $K \equiv -1$,
  domain $|p| < r_{\max}$. Closed-form distance (used as ground truth)
  $d(u,v) = \operatorname{arccosh}\!\big(1 + \frac{2|u-v|^2}{(1-|u|^2)(1-|v|^2)}\big)$, and geodesics by a
  Möbius transform. Analytic gradient $\partial_k(\lambda^2) = 2\lambda^3 p_k$.
* **`Terrain(amp, freq)`** and **`GaussianHills(hills)`**: the surface $z = f(x,y)$ with the induced metric
  $g = I + \nabla f \nabla f^\top$; $f = a\sin(\omega x)\cos(\omega y)$, resp.
  $f = \sum_i a_i \exp(-|p - c_i|^2/2\sigma_i^2)$. `build_metric_grad` uses
  $\partial_k g = h_k \nabla f^\top + \nabla f h_k^\top$ with $h_k$ the $k$-th Hessian column, so the
  energy gradient is analytic. Here $\lambda_{\min} = 1$ and $\lambda_{\max} = \sqrt{1 + |\nabla f|^2}$.

---------------------------------------------------------------------------------------------------

## 2. `rgeo/obstacles.py` — obstacles

`Obstacle` has `blocked(p)` (point test), `blocked_batch(P)` (vectorised) and `clearance_batch(P)`,
which returns the signed distance $c(p)$ to the boundary (negative inside) and $\nabla c$. `Circle` and
`Rectangle` implement all three; the rectangle's interior clearance is minus the distance to the nearest
side with that side's outward normal as gradient. `blocked(p, obstacles)` is the union test used
everywhere. Clearance is chart-Euclidean: it feeds a penalty (§5), not a length.

---------------------------------------------------------------------------------------------------

## 3. `rgeo/graphs/` — from a manifold to a graph

### 3.1 `build.py`

`build_graph(manifold, bounds, nx, ny, obstacles, type)` samples a regular lattice (rectangular, or
triangular with odd rows shifted by $\tfrac12 dx$), drops nodes outside $U_{\text{free}}$, and connects
each node to the offsets in `ds` = $\{(1,0),(0,1),(1,\pm1),(2,\pm1),(1,\pm2)\}$ (16 directions with their
reverses) or, for the triangular lattice, to its 6 lattice neighbours. An edge is kept if
`segment_free` (4 interior samples) finds no obstacle; its weight is `segment_lenght`.

**Why the graph has an accuracy floor.** As the lattice is refined at a fixed stencil, $d_G$ converges to
the distance of the *Finsler* metric $F(v) = \min\{\sum_i a_i |e_i|_g : v = \sum_i a_i e_i,\ a_i \ge 0\}$
over the stencil directions $e_i$, not to $d$. If $v$ lies between adjacent directions $e_1, e_2$ at
angle $\theta$ (measured in $g$), writing $v = a e_1 + b e_2$ gives $a|e_1| + b|e_2| \le |v| \sec(\theta/2)$,
with equality when $v$ bisects them. Hence
$$ 1 \le d_G / d \;\lesssim\; \sec(\theta_{\max}/2) \quad (h \to 0), $$
independently of the number of nodes. For the 16-direction stencil in the Euclidean metric
$\theta_{\max} = \arctan 2 - \pi/4 \approx 26.6^\circ$, a 2.7% worst case; measured means are 1.5–6%. This is
why no amount of graph resolution reaches the accuracy targets, and why §5 exists.

`connected_comps(adj)` labels components by DFS (used by `MultiH2HIndex`).

### 3.2 `dijkstra.py`

Binary-heap Dijkstra on `adj`, optional early exit at a target, returning distance and predecessor
arrays; `get_path(pred, s, t)` walks predecessors (returns `None` if unreachable). $O(|E|\log n)$.

### 3.3 `sampling.py` — adaptive node sets

A **difficulty field** on a $\mathrm{res}\times\mathrm{res}$ grid,
$$ D(p) = w_1\, \widehat{\|\nabla g(p)\|_F} + w_2\, \widehat{|K(p)|} + w_3\, e^{-\mathrm{dist}(p, O)/\sigma}, $$
where hats are normalisations by the 95th percentile clipped to $[0,1]$, $K$ is the Gaussian curvature
from the Brioschi formula with finite differences, $\sigma = 0.04 \times$ (domain diagonal), and $D$ is
finally divided by its maximum. **Quadtree sampling** (`sample_points_adaptively`): the domain is split
into $n_0 \times n_0$ cells; a cell gets refinement level
$\ell = \min(\ell_{\max}, \lfloor \max_{5\ \text{control points}} D / \tau \rfloor)$ and contributes the
$(2^\ell + 1)^2$ lattice points of its subdivision; points are deduplicated on the finest lattice, so
neighbouring cells of different levels share nodes. `boundary` points can be added on an offset ring
around each obstacle (`samlpe_boundary_points`). `graph_from_points` builds a $k$-NN graph with the same
free-segment test and Riemannian weights; `build_adaptive_graph` / `build_uniform_graph` are the two
constructions compared in the notebook at matched node counts.

---------------------------------------------------------------------------------------------------

## 4. `rgeo/oracles/` — graph distance oracles

### 4.1 `tree_decomposition.py` and `h2h.py` — the baseline (H2H)

**Tree decomposition by minimum-degree elimination.** Repeatedly eliminate the node $v$ of minimum
current degree; its bag is $X(v)$ = its current neighbours with their current edge weights; every pair
$u, x \in X(v)$ gets the fill-in edge $w(u,x) \leftarrow \min(w(u,x), w(u,v) + w(v,x))$ (`middle`
records $v$ for path unpacking); $v$ is removed. The parent of $v$ is the member of $X(v)$ eliminated
first after $v$; the last node is the root. This yields a tree $T$ of height $h$ and width
$\omega = \max |X(v)|$, with the property that $X(v) \cup \{v\}$ separates the subtree of $v$ from the rest
of $G$ along shortest paths.

**H2H labels** (Ouyang et al., SIGMOD 2018). For every node $v$ and every ancestor $a_i$ at depth $i$,
`dis[v, i]` $= d_G(v, a_i)$, computed top-down:
$\mathrm{dis}[v, i] = \min_{x \in X(v)} \big(w(v,x) + \mathrm{dis}[x, i]\big)$ (with the symmetric lookup when
$x$ is deeper than $i$). **Query**: $c = \mathrm{LCA}_T(s,t)$ by Euler tour and sparse table ($O(1)$), and
$$ d_G(s,t) = \min_{i \in \mathrm{pos}(c)} \mathrm{dis}[s,i] + \mathrm{dis}[t,i], \qquad \mathrm{pos}(c) = \{\mathrm{depth}(x) : x \in X(c) \cup \{c\}\}. $$
Exact for $d_G$; $O(\omega)$ per query, $O(nh)$ label entries. Its error is therefore exactly the graph's
(§3.1). `query_path` unpacks fill-in edges through `middle`. `MultiH2HIndex` builds one index per
connected component.

### 4.2 `separator.py` — `SeparatorIndex` and `SeparatorRefinedIndex`

**Separator tree.** `build(R, depth, B)` on a node set $R$ with bounding box $B$:
1. if $|R| \le$ `leaf_size` (or the depth limit is hit) $R$ becomes a leaf;
2. otherwise choose the axis $a$ of largest extent of $R$, the median $v$ of the $a$-coordinates, and the
   cut line $C = \{p : p_a = v\}$; $R$ splits into $R_< = \{p_a < v\}$, $R_\ge$ (ties broken by sorting);
3. **vertex portals**: the members within the band $|p_a - v| \le \delta$ ($\delta = 1.1\times$ median
   nearest-neighbour spacing), subsampled to $k$ evenly in the other coordinate;
4. record the separator (axis, value, depth, box, members, parent) and recurse on both halves with the
   box split at $v$.

Every node $u$ stores the list `anc[u]` of separators containing it; `rep[u]` is the deepest one.
**Lowest common separator** `lcs(s,t)` is the LCA of `rep[s]`, `rep[t]` in the separator tree
(Euler tour + sparse table, $O(1)$). For $s, t$ in different leaves it is the cut that separates them.

**Labels.** For every separator $C$, member $u$ and portal $p$: $\ell(u,p) = d_G(u,p)$ (one Dijkstra per
portal, cached). **Leaf tables**: for each leaf, exact $d_G$ among its members, computed in two passes:
a Dijkstra restricted to the leaf plus a `halo` of graph neighbours gives an upper bound
$U_s = \max_{t} d^{\text{local}}(s,t)$; a second, global Dijkstra truncated at $U_s$ is exact for every
target within $U_s$, so the table does not depend on the halo being large enough.

**Query.** $\hat D(s,t) = \min\big(\min_p \ell(s,p) + \ell(t,p),\ \text{leaf table}(s,t)\big)$, an
admissible answer (concatenation of two graph paths). It is *not* exact for $d_G$: the graph's shortest
path crosses the cut line, but not necessarily at a band portal, and possibly outside $B$.

**`SeparatorRefinedIndex`** keeps the same portals and replaces every label and leaf-table entry by
$\min(d_G, L(\mathrm{refine}(\text{graph path})))$ (§5), so labels become continuum lengths.

### 4.3 `cut_index.py` — the Separator Cut Index (the contribution)

`SeparatorCutIndex(points, adj, manifold, ...)` builds the separator tree of §4.2 and then, for every
cut, a set of **virtual portals on the cut line itself** with continuum labels. Three changes against
§4.2, each with a reason:

* portals on the line, not at nearby vertices — the crossing cost $f(x) = d(s,x) + d(x,t)$ becomes a
  function of **one** variable along the cut, so its minimum can be interpolated;
* labels are continuum distances (refined, or exact on terrains) — $f$ is then smooth, which
  interpolation needs and graph labels do not provide;
* the window and the portal budget are set by Riemannian quantities — the same coordinate spacing can
  differ in length by an order of magnitude (Poincaré disk near the boundary).

**Step 1 — cut geometry** (`cut_arclength`). For cut $C$ with region interval $[r_{lo}, r_{hi}]$ (the
box $B$ projected on the cut) the code samples the line, computes the metric speed along it,
$\sigma(t) = \sqrt{g_{bb}}$ ($b$ the other axis), and its cumulative arc length. The **window**
$W = [lo, hi]$ pads the region interval by a margin $\mu = \tfrac12 L(\mathrm{diag}\, B)$:
* with `pad_speed` $= \lambda_{\min}$ (terrains: $\lambda_{\min}=1$): $W = [r_{lo} - \mu/\lambda_{\min},\, r_{hi} + \mu/\lambda_{\min}]$;
* otherwise the arc length along the cut is inverted (a heuristic: a detour off the cut can be cheaper
  than the cut itself).

The padding lemma (notes/error_bound.md §2): if $s,t \in B$ are on opposite sides and the shortest path
crosses the cut line at $x^\ast$ outside $B$ by coordinate offset $c$, then
$d(s,t) \ge 2\, d(x^\ast, B) \ge 2 \lambda_{\min} c$, while $d(s,t) \le \mathrm{diam}_g(B)$; hence
$c \le \mathrm{diam}_g(B) / (2\lambda_{\min})$. A placement density $w(t) = \sigma(t)^{-\alpha}\,\mathrm{feas}(t)$
($\alpha=0$ by default: uniform in the chart) is integrated to a placement measure; `feas` decays linearly
from 1 at the region boundary to 0 at the end of the padding.

**Step 2 — portal budget and placement** (`build_cut_labels`).
* Budget $k_C = k$ for every cut (`scale_k=False`, the terrain default), or
  $k_C = \mathrm{clip}(k\, |C| / |C_{\max}|, k_{\min}, k_{\max})$ with $|C|$ the Riemannian length of the
  region part (`scale_k=True`).
* The window is split into **free runs** (maximal intervals outside obstacles and inside $U$) and each
  run into its in-region piece and its padded tails. The $k_C$ portals are shared between region pieces
  in proportion to their placement measure and spaced equally in it; each tail gets `k_pad` (=1)
  interior portals; optionally portals at the crossings with mesh edges (`edge_portals`, terrains).
* Portals are sorted along the cut; `cut_seg` records the free run of each, because interpolation must
  never fit across an obstacle gap.

**Step 3 — labels** $\ell(u, x_j)$ for every member $u$ of the region and every portal $x_j$:
* *refined* (`_label_portal`): Dijkstra from the node $a_j$ nearest to $x_j$; for every member $s$ the
  polyline $(x_j, a_j, \dots, s)$; all of them refined in **one** batched optimiser call (`refine_batch`,
  §5.3); the label is the refined length. Graph-portal labels of §4.2 are kept as a safety net.
* *exact* (`_label_separator_exact`, when `exact=ExactTerrain(...)`): one exact solver per cut on the
  sub-mesh its geodesics can visit, with the cut's portals inserted as Steiner vertices, and one MMP
  sweep per portal to all members (§8.3). Labels are stored as float32 (`label_dtype`).
* `workers > 1` forks a pool: (cut, portal) tasks for refined labels (the root cut alone is ~80% of the
  work), one task per cut for exact labels; children inherit the index by `fork`, only small dicts
  return. All reported numbers are single-threaded.

**Step 4 — leaf tables and halo** (`_leaf_table`). Refined: the two-pass Dijkstra of §4.2, then all
member-pair polylines of the leaf refined in one batched call. Exact: one local sweep per member on the
leaf's sub-mesh, which also returns the exact distance to every vertex within `halo_rings` graph hops of
the leaf (`halo_ids`, `halo_d`): pairs a few edges apart on opposite sides of a cut are answered there.

**Step 5 — feet.** For every member $u$ and cut, the foot $x_u$ (the chart projection on the line) and
the drop $r_u = L(\overline{u\, x_u})$, one batched length call per cut.

**Query** (`query`, `cut_bound`, `portal_min`). With $C = \mathrm{lcs}(s,t)$ and
$f_j = \ell(s, x_j) + \ell(t, x_j)$, $i = \arg\min_j f_j$:
* `interpolate=False`: $f_i$ — admissible.
* `interpolate=True`: the parabola through $(x_{i-1}, f_{i-1}), (x_i, f_i), (x_{i+1}, f_{i+1})$ in Newton
  form; if convex and its vertex lies in $[x_{i-1}, x_{i+1}]$, report
  $\max\big(f_v,\ f_i - \tfrac14 \min(f_{i-1}-f_i,\, f_{i+1}-f_i)\big)$. The guard is exact for a
  symmetric parabola and rejects fits through samples from different basins. $O(h^3)$ instead of $O(h^2)$
  for smooth $f$; two-sided.
* `interpolate="model"`: when an endpoint is within `model_reach` windows of the cut, $f$ is not smooth at
  scale $h$ because $d(s, x) \approx m_s(x) = \sqrt{r_s^2 + \sigma^2 (x - x_s)^2}$ has curvature $1/r_s$. The
  remainder $g_j = f_j - m_s(x_j) - m_t(x_j)$ is fitted by the parabola (or a line), and
  $m_s + m_t + \hat g$ is minimised on $[x_{i-1}, x_{i+1}]$ by bracketed Newton (bisection whenever a step
  leaves the bracket); the result is floored at $m_s + m_t$. Measured: maximum error 4.2% → 1.2% on
  HorseMount, k=32.
* $\hat D = \min(\text{cut bound}, \text{graph-portal bound (refined only)}, \text{leaf table}, \text{halo})$.

**Cost.** Build $= \sum_C k_C \times$ (one label solve on $C$'s region): refined, a Dijkstra plus one
batched refinement of $|R_C|$ polylines; exact, one MMP sweep on the sub-mesh of $C$. Space
$= \sum_u \sum_{C \ni u} k_C$ float32 labels $+$ leaf tables $+$ halo $+$ one foot per (node, cut): with
uniform $k$ and $O(\log n)$ depth, $O(n k \log n)$. Query: $O(1)$ for the LCS, $O(k_C)$ for the minimum,
$O(1)$ for the interpolation (≈10 Newton steps in the model case): 3–30 µs measured.

### 4.4 `variants.py` — recorded negative results

Kept because the notebook records why each fails (cells 33, 43, 48, 49, 53):
`SeparatorHeatIndex` (labels from the heat method: smooth but biased by the time step);
`SeparatorEikonalLeafIndex` / `SeparatorEikonalIndex` (fast marching on a grid: $O(h)$ absolute error,
and a grid fine enough to beat the graph is far too slow); `SeparatorContinuousIndex` (allows travel
between two portals, $\min_{i,j} \ell(s,p_i) + D_{ij} + \ell(t,p_j)$ — provably identical to §4.2 because
$\ell(s,p_j) \le \ell(s,p_i) + D_{ij}$ by the triangle inequality for graph distances, so it only doubles
the query time); `SeparatorBFSRefinedIndex` (balanced BFS frontiers, true vertex separators, so no path
escapes round a cut end — but the frontier is not a line, so no interpolation).

---------------------------------------------------------------------------------------------------

## 5. `rgeo/refine.py` — from a graph path to a continuum curve

`refine_polyline(manifold, coords, obstacles, rounds, segments, ...)` takes a free polyline (a graph path,
possibly extended to a cut point) and returns $(L, \text{curve})$ with $L$ no larger than the input's
length. Step size `step` defaults to a quarter of the median positive edge; it is the sampling
resolution of every free-space test.

### 5.1 One round

1. **Shortcut** (`shortcut`). Greedy farthest jump: from vertex $i$, move to the largest
   $j \le i + \text{window}$ whose chord $\overline{P_i P_j}$ is free and shorter than the path
   $P_i \dots P_j$. All candidate chords of all windows are measured in one batched length call, and the
   sampled free-space test runs once, batched, over the chords that are shorter. Same decisions as the
   pairwise loop.
2. **Subdivide** (`subdivide`). Split segment $i$ into $\lceil L_i / (L / \text{segments}) \rceil$ equal chart
   pieces, so the polyline can bend where it needs to (lengths are Riemannian, so the large-metric parts
   get more vertices).
3. **Relax** (`relax`). Minimise over the interior vertices, with L-BFGS-B and the analytic gradient,
   $$ E(P) + w \sum_{q \in \text{vertices} \cup \text{midpoints}} \max\big(0,\ \mu - c(q)\big)^2 , $$
   the discrete energy of §1.1 plus a quadratic barrier on the signed clearance $c$ (margin $\mu = 0$ by
   default, $w = 10^4$). If the optimum is not free, backtrack towards the input by factors
   $\tfrac12, \tfrac14, \tfrac18, \tfrac1{16}$ with a projection out of the obstacles; otherwise keep the input.
4. Keep the shortest *free* polyline seen so far (the input counts if it is free).

### 5.2 Properties the index relies on

* **Admissibility**: the returned curve is free at the sampling resolution `step` and its length is the
  quadrature length of an actual polyline, so labels are upper bounds of $d$ up to (i) obstacles thinner
  than `step` and (ii) Gauss–Legendre error. On terrains (ii) vanishes (exact walk).
* **Homotopy**: shortcut accepts only free chords and relaxation is a local descent with a barrier, so
  the result stays in the homotopy class of the graph path. This is why the graph route matters (it
  decides which side of every obstacle the curve passes) and why refinement cannot repair a graph path
  in the wrong corridor (the weighted-terrain result, notes/error_bound.md §6.5).

### 5.3 Batched refinement (`refine_batch`, `relax_stack`)

`relax` spends its time in interpreter overhead (~50 L-BFGS steps on a 12-point array). `refine_batch`
runs the same rounds for $B$ polylines at once: each polyline is resampled (`resample_uniform`) to
$4 \lceil n_i / 4 \rceil$ equal-arc pieces ($n_i$ the sequential subdivision count, so never coarser),
polylines with equal counts are stacked into a $(B, n+1, 2)$ array, and one L-BFGS minimises the summed
energy plus barrier; because the objective is separable the gradient is the concatenation of the
per-polyline gradients. Per-polyline free checks and fallbacks as in §5.1. 2.2–2.6× faster than the
sequential path with the same accuracy (0.1312% vs 0.1308% on the Poincaré benchmark).

---------------------------------------------------------------------------------------------------

## 6. `rgeo/solvers/` — continuum baselines

* **`eikonal.py`** — fast marching for the anisotropic eikonal equation $|\nabla u|_{g^{-1}} = 1$ on a
  regular grid: a heap of trial values, one- and two-point upwind updates (`two_point` solves the
  quadratic for the anisotropic metric $(e, f, g)$ on a triangle), obstacles as impassable cells;
  `sample_field` bilinearly interpolates the field, `anisotropy_ratio` measures $\lambda_{\max}/\lambda_{\min}$.
  First-order accurate; its error is a fixed fraction of the grid step, which is the reason it loses.
* **`heat_method.py`** — Crane et al.'s heat method on a Delaunay triangulation of the nodes (long edges,
  faces in obstacles and faces crossing them dropped): P1 finite elements for the Laplace–Beltrami
  operator in the chart, stiffness $L_{ab} = \sum_T |T|_g\, \nabla\varphi_a^\top g_T^{-1} \nabla\varphi_b$ (the metric
  cotangent Laplacian) and lumped mass $M$; one implicit heat step $(M + tL)u = \delta_s$ with
  $t = t_{\text{factor}}\,\bar h^2$, normalised field $X = -\nabla u/|\nabla u|_g$, Poisson solve $L\phi = \nabla\cdot X$.
  Smooth but biased by the time step $t$.

---------------------------------------------------------------------------------------------------

## 7. `rgeo/visualize.py`

Plotting helpers: metric-scale background (`plot_background`), obstacles, graphs, components, paths,
fields, the separator tree with portals (`plot_separators`), and `visualize_query` which composes them.

---------------------------------------------------------------------------------------------------

## 8. `rgeo/terrain/` — triangulated terrains

A terrain is a triangulation $\mathcal T$ of points $(x_v, y_v, z_v)$ whose $(x,y)$ projection is a planar
triangulation (a height field). It is the object of the terrain-oracle literature; its geodesics are the
exact polyhedral ones computed by MMP (Mitchell–Mount–Papadimitriou 1987).

### 8.1 `io.py`

`read_off`, `write_off` (OFF format). `repair_height_field(V, F)`: if most faces are clockwise, flip
all; then for every face with non-positive projected area, try each of its vertices: slide it toward the
centroid of its 1-ring by the smallest fraction (bisection, 30 steps) that makes all incident faces
positively oriented, and take the vertex needing the smallest displacement; $z$ is unchanged. BigMount
(1 face) and GunnisonForest (6 faces) need it; moves are 20–80 m on 190–450 m edges.

### 8.2 `mesh.py` — `TerrainMesh`, the terrain as a manifold

On face $f$ the surface is the plane $z = a_f x + b_f y + c_f$, so the metric is constant per face:
$$ g_f = w_f^2 \,(I + \nabla f\,\nabla f^\top), \qquad \nabla f = (a_f, b_f), $$
with optional face weights $w_f$ (weighted-region terrains; $w \equiv 1$ otherwise). Point location uses
matplotlib's trapezoid-map finder, with a barycentric-tolerance fallback around the nearest vertex
(points on the boundary).

**Exact segment length** (`segment_lengths_batch`). For a chart segment $A \to B$, $d = B - A$, the
segment is walked face by face in lock-step for all segments of a batch. In face $k$ the lifted speed is
$s_k = w_k |M_k d|$ with $M_k = \begin{psmallmatrix} 1 & 0 \\ 0 & 1 \\ a_k & b_k\end{psmallmatrix}$, and
$$ L(A,B) = \sum_k (t_{k+1} - t_k)\, s_k , $$
$t_k$ the parameters where the segment crosses face edges, found by solving
$A + t d = v_0 + \tau e$ per edge (Cramer). Quadrature would be wrong by $O(\text{face jump})$ at every
edge, which is the size of the errors being measured. Robustness, all needed on real meshes: the walk
runs on the segment shrunk by $10^{-7}$ at both ends (endpoints are usually vertices) and adds the cut-off
pieces back; edges parallel to the segment are not crossings; a segment *on* an edge is charged at the
cheaper adjacent face; a walk stuck at a vertex is relocated by point location just ahead; a boundary
edge is confirmed by point location beyond it before the segment is declared outside ($L = \infty$).
Verified: mesh edges reproduce their 3-D lengths to $10^{-14}$; chords agree with 20,000-point brute
force to $10^{-4}$.

**Gradient.** $L$ depends on $B$ through the speeds and through the crossing parameters:
$$ \frac{\partial L}{\partial B} = \sum_k (t_{k+1} - t_k)\, \frac{w_k^2 M_k^\top M_k d}{s_k}
   \;+\; \sum_{\text{crossings } k} (s_k - s_{k+1})\, \frac{\partial t_k}{\partial B},\qquad
   t_k = \frac{(v_0 - A) \times e}{d \times e}, $$
and symmetrically for $A$. Verified against central differences to $10^{-8}$. With it the polyhedral
energy $E = m \sum_i L_i^2$ is exact and continuous in the vertices (chart quadrature jumps when a
sample crosses an edge), so `refine.py` works on terrains unchanged
(`energy_and_grad_points`, `energy_and_grad_stack`; segments that leave the surface get a large linear
penalty instead of $\infty$, which would stop L-BFGS for the whole stack).

**`mesh_graph()`**: nodes = vertices, edges = mesh edges, weight = 3-D length times the cheaper adjacent
face weight — the graph every terrain oracle starts from. **`edge_crossings(axis, v, lo, hi)`**: where a
cut line crosses mesh edges.

### 8.3 `exact.py` — `ExactTerrain`, exact labels

* `insert_points(X)`: each point is inserted as a Steiner vertex by splitting its (sub-)face in three
  (repeatedly, if several points share a face); the surface is unchanged.
* `distances_from_points(X, targets)`: exact geodesic distances (pygeodesic's MMP) from each point of $X$.
* `region_geo(box, portals)`: an exact solver on the part of the terrain that geodesics between points of
  `box` (extended to contain the portals) can visit. **Lemma.** On a terrain $\lambda_{\min} = 1$, so a
  curve of length $\ell$ stays within Euclidean distance $\ell$ of its start; for $p, q$ in the box,
  $d(p,q) \le \lambda_{\max}(\text{box}) \cdot \mathrm{diag}$; hence every shortest path between them lies in
  the box dilated by $R = \lambda_{\max}\,\mathrm{diag}$, and distances on the sub-mesh of faces touching that
  dilation are exact. Deep cuts get small sub-meshes; this is what keeps the exact build near
  $O(n \log n)$ sweeps-size instead of $O(n^2)$ (16× on HorseMount at k=64; verified identical labels).
* `submesh`, `local_table`: the ring-based variant used before `region_geo`.
* Weighted terrains are refused: MMP does not handle face weights.

### 8.4 `a2a.py` — `TerrainA2A`, arbitrary surface points

For $p$ in face $f$ (barycentric weights $\beta$) and $q$ in face $f'$:
* faces within `near` rings of each other: both points are inserted into the local sub-mesh and one
  exact sweep answers (no interpolation is valid when $d \sim$ edge length);
* `mode="snap"`: $\min_{a \in f, b \in f'} |pa| + \hat D(a,b) + |bq|$ — admissible but **first order**: the
  corner at $a$ costs $|pa|(1-\cos\theta)$, independently of how accurate $\hat D(a,b)$ is (exact vertex
  distances give the same 2%);
* `mode="interp"` (default): $\hat D(p,q) = \sum_{a,b} \beta_a \beta'_b\, \hat D(a,b)$, bilinear
  interpolation of the $3\times 3$ vertex-pair answers. Off the cut locus the distance field is $C^2$ in
  a face with Hessian $\lesssim 1/d$, so the error is $O(e^2/d)$, second order. Interpolating the label
  vectors first and taking one portal minimum would be 9× cheaper but has a Jensen gap first order in the
  spread of the vertices' argmins (measured 0.8% vs 0.18%), so the minimum is taken per vertex pair.

---------------------------------------------------------------------------------------------------

## 9. `bench/` — the benchmark

**`terrain_bench.py`** — `run(mesh_path, ...)`:
1. read and repair the mesh, build `TerrainMesh` and the mesh graph; with `--weighted SEED` draw face
   weights uniform in $(1, 5)$ (the competitors' distribution);
2. query sets, seeded: 100 vertex pairs (`vertex_queries`) and, with `--a2a N`, $N$ arbitrary points
   (uniform face, uniform barycentric — `a2a_queries`); both are also written in the competitors' format,
   vertex points nudged 5 cm into a face (`write_competitor_queries`, see §10);
3. truth: exact MMP for vertex pairs; both points inserted for arbitrary pairs; on weighted terrains a
   heavy refinement whose values are later replaced by the best path any method found;
4. methods: Dijkstra, H2H, vertex portals, the cut index in several configurations (`exact` with a list of
   $k$, `exact_scaled`, `cut` = refined labels), each with no / parabola / model interpolation and the
   two A2A modes;
5. `measure`: signed relative errors and the best of three timing loops per query set;
   `index_bytes`: label arrays + leaf tables + portal coordinates + feet + halo (H2H: 12 bytes per label
   entry). Rows are appended to `results/<mesh>.json` after every method, merged with earlier rows.

**`terrain_report.py`** — `load` merges our rows with the competitors' (`<mesh>_competitors[_a2a|_w2v].json`),
`table` renders one Markdown table per mesh and query kind, `pareto` plots error against query time and
build time, `write_note` regenerates `notes/terrain_results.md` with its caveats.

---------------------------------------------------------------------------------------------------

## 10. `scripts/` — data and long runs

* `make_terrains.py`: regenerates `data/terrains/` from the EAR-Oracle repository's datasets (repair,
  translation to the origin, the 10k and 1k crops of GunnisonForest, `--large` for the 100k meshes).
* `run_scale.sh`, `run_weighted.sh`: the 100k and weighted experiments, ours then theirs; meant to run
  detached (`nohup ... &`).
* `competitors/`: `README.md` (macOS/CGAL-6 patches of their harness, the float32 and point-location
  pitfalls), `run_competitors.sh` (runs `method[@eps[@sp]]` on a query file; `KIND`, `FACE_WEIGHTS`),
  `competitor_eval.py` (scores their logs against the same truth; on weighted terrains it builds the
  best-path reference and excludes their K-Algo, which ignores face weights).

---------------------------------------------------------------------------------------------------

## 11. `notebooks/`, `tests/`, `results/`, `notes/`

* `notebooks/exps.ipynb` — the experiments in project order: discretisation, graph vs exact distance,
  H2H, obstacles, components, adaptive sampling, edge impact, the separator tree, error analysis
  ($E = E_{quad} + E_{curve} + E_{metr}$, the stencil floor), refinement, the portal-budget study, the
  full comparison of index variants, adaptive sampling with the cut index, and the terrain benchmark
  (cells 63–67). Cell 1 changes into the repository root, so it runs from `notebooks/` or the root.
* `tests/test_core.py` — invariants: H2H equals Dijkstra; the cut index without interpolation is an upper
  bound and beats the graph; the exact walk reproduces mesh edges; walk gradients match finite
  differences; exact labels match MMP; weighted terrains refuse exact labels; refinement avoids
  obstacles.
* `results/` — see §9; `results/logs/` (gitignored) holds the run logs.
* `notes/` — this file, the error bound, the preprocessing account, the generated benchmark tables.

---------------------------------------------------------------------------------------------------

## 12. Gaps found while writing this walk-through

These do not change any reported number by much, but a paper's proofs need them fixed.

1. **Padding margin.** `cut_arclength` uses $\mu = \tfrac12 L(\mathrm{diag}\,B)$, the length of one chart
   segment, as a stand-in for $\tfrac12\mathrm{diam}_g(B)$. The lemma needs an *upper* bound of
   $d(s,t)$ over member pairs, e.g. $\tfrac12 \lambda_{\max}(B)\,|\mathrm{diag}\,B|$ as `region_geo` already
   uses. On terrains the two differ by the slope variation inside $B$.
2. **Tail resolution.** Each padded tail gets `k_pad = 1` portal, so for a crossing in a tail the portal gap
   in the bound is the tail length. The window makes the crossing *covered*, not *well resolved*;
   empirically such crossings are rare.
3. **Model-interpolation units.** `_model_min` scales cut offsets by the slope of the placement measure
   (`cut_arc`), which for $\alpha = 0$ is 1 inside the region — chart units, not the metric speed
   $\sigma$ the comment describes. Harmless on terrains (the fitted remainder absorbs it), worth fixing to
   the true arc length for strongly anisotropic metrics.
4. **Admissibility of refined labels** holds at the sampling resolution `step` of the free-space test and
   up to Gauss–Legendre error on analytic metrics (§5.2). Exact on terrains.

Minor: `foot_r` falls back to a drop of 0 when the foot segment leaves the domain (the model then sees a
cone); interpolated answers are two-sided by design; for two nodes in the same leaf `lcs` returns the
leaf's parent cut, which still gives an admissible portal answer, but the leaf table is what answers them.
