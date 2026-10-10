# A (1+ε) guarantee for the Certified Separator Cut Index

This note states and proves the guarantee implemented by `rgeo/oracles/certified.py`
(`CertifiedCutIndex`). Every definition below is exactly what the code computes; the line of code
that realises each step is named in brackets.

## 1. Setting

A terrain $T$ is a triangulated height field: a triangulation of points $(x_v, y_v, z_v)$ whose
$(x,y)$-projection is a planar triangulation of a chart domain $U \subset \mathbb R^2$. Each face
$f$ carries a cost $w_f \ge w_{\min} > 0$ (unweighted terrains: $w \equiv 1$). The cost of a curve
is $\sum_f w_f \times$ (3-D length inside $f$); $d(p,q)$ is the infimum over curves on $T$. Two
elementary lower bounds are used throughout. For any curve $\gamma$ from $p$ to $q$:

$$ L(\gamma) \;\ge\; \lambda \, L_{\text{chart}}(\gamma) \;\ge\; \lambda\,|p - q|, \qquad
   L(\gamma) \;\ge\; w_{\min}\,|P - Q|_{3D}, \qquad \lambda := w_{\min}, \tag{1}$$

because a face's metric is $w_f^2(I + \nabla f \nabla f^\top) \succeq w_f^2 I$ in the chart, and
3-D length dominates the length of the straight chord. Here $|p-q|$ is chart distance and $P, Q$
are the lifted points. Shortest paths exist (compactness).

**Parameters:** $\varepsilon \in (0,1)$ and a portal gap $g > 0$ (an arc length). Put
$R := g/(\varepsilon\lambda)$ [`self.Rn`].

## 2. The index

**Separator tree** [`_tree`]. Regions are vertex sets. For a region $M$, let $B_M$ be the
bounding box of its chart points. If $\operatorname{diag}(B_M) \le R$ (or $|M| \le 1$), $M$ is a
leaf. Otherwise choose the axis $a$ of the larger extent and a value $v$ that lies strictly between
the two middle member coordinates **and differs from the $a$-coordinate of every vertex of $T$**
[`_split_value`]. The cut is the line $\ell_C = \{x_a = v\}$; the children are
$M_< = \{u : u_a < v\}$ and $M_> = \{u : u_a > v\}$ — both non-empty, and no vertex of $T$ lies on
$\ell_C$. Let $I_C = \ell_C \cap B_M$, a segment (because $v$ lies strictly inside $B_M$'s range).

**Window** [`_cut_geometry`]. Pick a hub $h \in M$ and let $\mu_C = \max_{u \in M} \tilde U(h,u)$,
where $\tilde U(h,u) \ge d(h,u)$ is the mesh-graph distance (a graph path is a surface curve; an
edge costs its length times the cheaper adjacent face weight, which is the cost of walking along
it). The window is
$$ W_C = \{\, x \in \ell_C \;:\; \operatorname{dist}(x, I_C) \le \mu_C / \lambda \,\}. $$
$W_C \cap T$ is a finite union of closed segments, the **runs**, computed exactly from the faces the
line crosses [`_line_pieces`]. On each run the lifted cut is a polyline (straight inside every
face), and its arc length $A(\cdot,\cdot)$ is exact and additive.

**Portals** [`_place`]. A finite set $P_C \subset W_C \cap T$, sorted along the line, containing
both endpoints of every run and both endpoints of $I_C \cap$ run. For two consecutive portals of the
same run with arc distance $A_j$:
$$ A_j \le g \ \text{ if the interval lies in } I_C, \qquad
   A_j \le \max(g,\; 2\varepsilon\lambda\, c_j) \ \text{ otherwise}, \tag{2}$$
where $c_j$ is the chart distance from the interval's near end to $I_C$. (Implementation detail: a
portal within $10^{-7}$ of a face edge is moved along the line by $\sim 10^{-5}$ edge lengths before it is
inserted; all arc positions and gaps are recomputed from the inserted positions [`_finalize_cut`], and
a run end that is no longer a portal leaves an *end stretch* $[\text{run end}, p_j]$ of arc $e_j$.)

**Labels** [`_labels`, `ExactLabels.cut`]. $\ell(u,p) = d(u,p)$ for every $u \in M$, $p \in P_C$.

**Near tables** [`_near`]. For every vertex $u$: $N(u) = \{(t, d(u,t)) : |u - t| \le R\}$.

**Query** [`bounds`]. For $s \ne t$:
* if $|s-t| \le R$: return $d(s,t)$ from $N(s)$;
* otherwise let $C$ be the deepest cut whose region contains both [`_cut_of`] — $s$ and $t$ lie in
  different leaves (a leaf has diagonal $\le R$), hence on opposite sides of $\ell_C$. With
  $f_j = \ell(s,p_j) + \ell(t,p_j)$:
$$ \mathrm{UB} = \min_j f_j, \qquad
   \mathrm{LB} = \max\Big( w_{\min} |S - T|_{3D},\ \min\big\{\ \tfrac12 (f_j + f_{j+1}) - A_j,\ \ f_k - 2e_k\ \big\} \Big), $$
  the inner minimum over consecutive portals of the same run and over the end stretches. The
  returned estimate is any value in $[\mathrm{LB}, \mathrm{UB}]$ (the code clamps an interpolated minimum
  of $f$ into it).

## 3. Proofs

**Lemma 1 (crossing).** Every curve on $T$ from $s$ to $t$ meets $\ell_C$.
*Proof.* Its chart projection is continuous with $s_a < v < t_a$ (intermediate value theorem). ∎

**Lemma 2 (window).** Let $\gamma$ be a shortest path from $s$ to $t$ and $x^\ast \in \gamma \cap \ell_C$.
Then $x^\ast \in W_C$; moreover $d(s,t) \ge 2\lambda\, \operatorname{dist}(x^\ast, I_C)$.
*Proof.* By (1), $d(s,t) = L(\gamma_{s x^\ast}) + L(\gamma_{x^\ast t}) \ge
\lambda(|s - x^\ast| + |x^\ast - t|) \ge 2\lambda\, \operatorname{dist}(x^\ast, B_M)$, and for a point of
$\ell_C$ the nearest point of the box $B_M$ lies on $I_C$, so
$\operatorname{dist}(x^\ast, B_M) = \operatorname{dist}(x^\ast, I_C)$. On the other side,
$d(s,t) \le d(s,h) + d(h,t) \le 2\mu_C$. ∎

**Lemma 3 (sandwich).** $\mathrm{LB} \le d(s,t) \le \mathrm{UB}$ for every pair.
*Proof.* Near pairs are exact. Otherwise every $f_j$ is the cost of a curve through $p_j$, so
$f_j \ge d(s,t)$. Let $x^\ast$ be as in Lemma 2; it lies on a run, either in an interval $[p_j, p_{j+1}]$
of consecutive portals or in an end stretch. In the first case let $a = A(p_j, x^\ast)$. The lifted cut
between $p_j$ and $x^\ast$ is a curve, so $d(x^\ast, p_j) \le a$ and
$$ f_j \le \big(d(s,x^\ast) + a\big) + \big(a + d(x^\ast,t)\big) = d(s,t) + 2a, \qquad
   f_{j+1} \le d(s,t) + 2(A_j - a). $$
Adding, $\tfrac12(f_j + f_{j+1}) - A_j \le d(s,t)$. In the second case, the same argument with the
single portal $p_k$ of the stretch gives $f_k \le d(s,t) + 2A(x^\ast, p_k) \le d(s,t) + 2e_k$. Either way
the inner minimum is $\le d(s,t)$; the chord term is $\le d(s,t)$ by (1). (If a run has no portal at all —
possible only if a label failed — the code drops the interval term for that cut, LB = chord.) ∎

**Lemma 4 (accuracy).** If $|s - t| > R$, then $\mathrm{UB} \le (1+\varepsilon)\, d(s,t)$ and
$\mathrm{LB} \ge (1-\varepsilon)\, d(s,t)$.
*Proof.* By (1), $d := d(s,t) \ge \lambda |s-t| > \lambda R = g/\varepsilon$.
*Upper bound.* From the proof of Lemma 3, $\min(f_j, f_{j+1}) \le d + \min(2a, 2(A_j - a)) \le d + A_j$
(and $f_k \le d + 2e_k$ in an end stretch, $2e_k \le 2\cdot 8\cdot 10^{-5}$ edge lengths $\ll g < \varepsilon d$).
If the interval lies in $I_C$, $A_j \le g < \varepsilon d$. Otherwise every point of the interval,
$x^\ast$ included, has $\operatorname{dist}(\cdot, I_C) \ge c_j$, so $d \ge 2\lambda c_j$ by Lemma 2,
and (2) gives $A_j \le \max(g, 2\varepsilon\lambda c_j) \le \varepsilon d$.
*Lower bound.* An end-stretch term is $f_k - 2e_k \ge d - 2e_k > (1-\varepsilon)d$. Each interval term
satisfies $\tfrac12(f_j + f_{j+1}) - A_j \ge d - A_j$ (as $f \ge d$). If
$A_j \le \varepsilon d$ we are done. Otherwise $A_j > \varepsilon d > g$, so the interval is a tail interval
with $A_j \le 2\varepsilon\lambda c_j$, hence $c_j > d/(2\lambda)$. Both of its portals $p$ are at chart
distance $\ge c_j$ from $B_M$, so by (1) $f(p) \ge \lambda(|s - p| + |p - t|) \ge 2\lambda c_j$, and the
term is $\ge 2\lambda c_j - 2\varepsilon\lambda c_j = 2\lambda c_j (1 - \varepsilon) > (1-\varepsilon) d$. ∎

**Theorem.** For every pair of vertices $s,t$:
$$ (1-\varepsilon)\, d(s,t) \;\le\; \mathrm{LB} \;\le\; d(s,t) \;\le\; \mathrm{UB} \;\le\; (1+\varepsilon)\, d(s,t), $$
with equality throughout when $|s-t| \le R$. Consequently the returned estimate has relative error
at most $\varepsilon$, $\mathrm{UB}$ is the cost of an actual path (admissible), and every query carries the
certificate $\mathrm{UB}/d - 1 \le \mathrm{UB}/\mathrm{LB} - 1$, usually far below $\varepsilon$.

**Approximate labels.** If labels satisfy only $d \le \ell \le (1+\eta)d$ (Steiner-graph labels on
weighted terrains), the same argument gives $d \le \mathrm{UB} \le (1+\eta)(1+\varepsilon)\,d$; the
interval lower bound is not valid then and the code keeps only the chord (`oracle.exact = False`).

## 4. Exactness of the restricted computations

**Lemma 5 (sub-mesh).** Let $x$ be a portal and $u$ a member, and suppose
$d(x,u) \le D$. Every point $z$ of a shortest path from $x$ to $u$ lies within chart distance
$D/(2\lambda)$ of the segment $[x, u]$.
*Proof.* By (1), $|z - x| + |z - u| \le L_{\text{chart}}(\gamma) \le D/\lambda$: $z$ lies in the
ellipse with foci $x, u$ and major axis $D/\lambda$, whose points are within its semi-minor axis
$\tfrac12\sqrt{(D/\lambda)^2 - |x-u|^2} \le D/(2\lambda)$ of the segment. ∎

With $D = \max_x \tilde U(x,h) + \mu_C$ (through the hub), every label geodesic stays inside the box
of $M \cup P_C$ dilated by $D/(2\lambda)$; the labels are computed by MMP on the faces meeting that
dilated box [`_labels`], with the portals inserted as Steiner vertices (each new vertex lies on its
face's plane, so the surface is unchanged). A geodesic of $T$ inside the sub-mesh is a curve of the
sub-mesh, and sub-mesh curves are curves of $T$, so the sub-mesh distance equals $d$. The near tables
use the same argument with graph upper bounds of the near distances [`_near`].

**Early termination.** MMP (Mitchell–Mount–Papadimitriou; Kirsanov's implementation in pygeodesic)
processes windows in non-decreasing order of their minimum distance. Once that minimum exceeds the
current value at every target, no later window can lower any target value. With `max_distance=0`
pygeodesic stops exactly then (its default, an infinite bound, disables the stop condition), so the
distances are unchanged — verified bit-identical — while a sweep only explores the neighbourhood its
targets require (10–200× faster for nearby targets).

**Floating point.** Stored distances are float32 rounded *up* (`f32_up`), so labels and near
distances remain upper bounds bit-exactly; lower bounds subtract the matching relative margin
($2^{-21}$).

## 5. Size, query time, construction

Let $\sigma$ be the maximum metric speed along the cut lines ($\sigma = \max_f w_f\sqrt{1+|\nabla f|^2}$
in the worst case, far smaller in practice) and $\rho$ the vertex density.

* **Portals per cut.** In the box, $\lceil A(I_C)/g \rceil + 1 \le \sigma |I_C|/g + 2$. In each tail,
  arc gaps equal $g$ until the offset reaches $R/2$ — at most $\sigma/(2\varepsilon\lambda) + 1$ portals —
  and then the chart offset grows by a factor $\ge 1 + 2\varepsilon\lambda/\sigma$ per portal up to
  $\mu_C/\lambda$, so a tail has $O\big(\frac{\sigma}{\varepsilon\lambda}\log\frac{\mu_C}{\lambda R}\big)$ portals.
* **Labels per vertex** $= \sum_{C \ni u} |P_C|$. For balanced cuts of a square domain of side $D$ the
  cut lengths along a root-to-leaf path are $D, D/2, D/2, D/4, D/4, \dots$, so the box parts sum to about
  $3\bar\sigma D/g$ ($\bar\sigma$ the mean speed along the cuts), plus
  $O\big(\frac{\sigma}{\varepsilon\lambda}\log\frac{D}{R}\big)$ per level for the tails.
* **Near table** per vertex: $\approx \pi R^2 \rho = \pi \rho\, g^2/(\varepsilon\lambda)^2$ entries.
* **Query:** $O(\text{depth})$ to find $C$, $O(|P_C|)$ for $f$, UB and LB; near pairs $O(\log |N(u)|)$.
* **Construction:** one early-terminated MMP sweep per portal on its sub-mesh, one per vertex for the
  near tables, one graph Dijkstra per cut for $\mu_C$.

With $g$ a fixed multiple of the edge length, the labels grow like $n^{3/2}$ in total (the top cuts
have $\Theta(\sqrt n)$ portals and every vertex stores them). Thorup-style $\varepsilon$-covers — each vertex
keeps only the $O(1/\varepsilon)$ portals needed to reach the rest *along the cut* — reduce this to
$O(\frac{n}{\varepsilon}\log n)$ without changing the construction; see §6.

## 6. Covers (space reduction, same construction)

For a vertex $s$ and cut $C$, call $Q \subseteq P_C$ an $\varepsilon_c$-cover if for every portal $p$ of the
same run, $\min_{q \in Q} \ell(s,q) + A(q,p) \le (1+\varepsilon_c)\,\ell(s,p)$. Storing only $\ell(s,\cdot)$ on
covers and answering with
$\widetilde{\mathrm{UB}} = \min_{q \in Q_s,\, q' \in Q_t} \ell(s,q) + A(q,q') + \ell(q',t)$ (a merge along the
cut in $O(|Q_s| + |Q_t|)$) is the cost of a real path (s to $q$, along the cut, $q'$ to $t$), and at the
portal $p^\ast$ that attains UB the covers of $s$ and $t$ give $q, q'$ with
$\ell(s,q) + A(q,q') + \ell(q',t) \le \ell(s,q) + A(q,p^\ast) + A(p^\ast,q') + \ell(q',t) \le (1+\varepsilon_c)\,\mathrm{UB}$,
so $d \le \widetilde{\mathrm{UB}} \le (1+\varepsilon_c)(1+\varepsilon)\,d$. In flat geometry, with $s$ at distance $r$ from
the cut and $\theta$ the angle between the cut and the ray from $s$ to the point $x$ of the cut, a greedy
cover adds a portal whenever $\int (1-\cos\theta)\,dx$ since the last one reaches $\varepsilon_c\,\ell(s,x)$; the count is
$\int (1-\cos\theta)/(\varepsilon_c \ell)\,dx = \int_0^{\pi/2} \frac{\cos\varphi}{1+\sin\varphi}\,\frac{d\varphi}{\varepsilon_c} = \frac{\ln 2}{\varepsilon_c}$
per side, independent of $r$ and $n$ — hence $O(\frac{n}{\varepsilon_c}\log n)$ labels in total. The number of
MMP sweeps does not change (every portal still needs its sweep), and the answer is only as good as
$\varepsilon_c$: matching the measured accuracy of full rows (mean error $\sim 10^{-4}$) needs $\varepsilon_c \sim 10^{-4}$,
i.e. more cover portals than the cut has. The implementation therefore keeps full rows; covers are the
route to provably near-linear space when a coarser answer suffices.
