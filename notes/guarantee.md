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

## 7. Arbitrary surface points (`rgeo.terrain.CertifiedA2A`)

For $p$ in face $f$ and $q$ in face $f'$ with corners $a_i$, $b_j$: the straight segment inside a face
is a path, so $d(p,a_i) \le w_f|P-A_i|$ and likewise for $q$. From the vertex certificates
$[\mathrm{LB}(a_i,b_j), \mathrm{UB}(a_i,b_j)]$ and the triangle inequality:
$$ \mathrm{UB}(p,q) = \min_{i,j}\, w_f|P-A_i| + \mathrm{UB}(a_i,b_j) + w_{f'}|Q-B_j| \;\ge\; d(p,q), $$
$$ \mathrm{LB}(p,q) = \max\Big(\lambda|P-Q|,\ \max_{i,j}\, \mathrm{LB}(a_i,b_j) - w_f|P-A_i| - w_{f'}|Q-B_j|\Big) \;\le\; d(p,q). $$
Pairs with $|p-q| \le R$ — and, with the certify-or-refine rule, every pair whose interval is wider
than $\mathrm{UB} > (1+w)\,\mathrm{LB}$ — are computed exactly: one MMP sweep on the faces meeting
$\mathrm{bbox}(p,q)$ dilated by $\mathrm{UB}/(2\lambda)$ (ellipse lemma), with both points inserted (a point within
$10^{-7}$ of an edge is moved by $\sim10^{-5}$ edge lengths; the move is added to both bounds).
**Theorem (A2A).** Every answer satisfies $\mathrm{LB} \le d(p,q) \le \mathrm{UB} \le (1+w)\,\mathrm{LB} \le (1+w)\,d(p,q)$;
the default $w = 2\varepsilon$. Without refinement the bound is additive:
$\mathrm{UB} \le (1+\varepsilon)\,d + 2(1+\varepsilon)(w_f r_p + w_{f'} r_q)$, $r$ the distance to the nearest corner. Measured
(100 random surface-point pairs per mesh, $\varepsilon = 0.1$): sound everywhere; on BigMount one pair in 100
was refined, the largest UB error went from 20.2% to 10.2%.

## 8. Size variants (C++: `cpp/src/certified.cpp`, driver `cpp/tools/cci_exp.cpp`)

**A. Capped portals with certify-or-refine** (`Options::k_max`, `bounds_refined`). A cut puts at most
$k$ portals inside its box: gap $g_C = \max(g, A_C/k)$, $A_C$ the arc length of the box part (tails step
$\max(g_C, 2\varepsilon\lambda c)$). Lemma 3 (soundness) does not use the gap, so $\mathrm{LB} \le d \le \mathrm{UB}$ for every
pair. Lemma 4 needs $g_C \le \varepsilon d$; when that fails the interval may be wide, and the query then
computes the pair exactly (one MMP sweep from $s$ stopped at $t$ on the faces meeting
$\mathrm{bbox}(s,t)$ dilated by $\mathrm{UB}/(2\lambda)$ — ellipse lemma). **Theorem A.** Every answer satisfies
$\mathrm{UB} \le (1+\varepsilon)\,\mathrm{LB} \le (1+\varepsilon)\,d$, or is exact. Per cut: at most $k$ box portals, plus the tails.
With $g_C \approx \sigma D_C/k$, a tail steps $g_C$ until the offset reaches $g_C/(2\varepsilon\lambda)$, which takes
$O(1/(\varepsilon\lambda))$ portals. It then grows geometrically up to $\mu_C/\lambda \approx D_C/\lambda$, which takes
$O\big(\frac{\sigma}{\varepsilon\lambda}\log\frac{2\varepsilon k}{\sigma}\big)$ portals. So a cut has
$k + O(\varepsilon^{-1}\log k)$ portals, labels per vertex are $\le (k + O(\varepsilon^{-1}\log k))\cdot$ depth, and the
total is $O((k+\varepsilon^{-1}\log k)\,n\log n)$. The cap bounds the label sweeps the same way: a level has
$2^\ell$ cuts of $k + O(\varepsilon^{-1}\log k)$ sweeps over $\Theta(n/2^\ell)$ faces, so the build is
$O((k+\varepsilon^{-1}\log k)\,n\log^2 n)$ for $O(m\log m)$ sweeps, against $\Theta(n^{1.5}\log n)$ uncapped; a uniform random pair needs the exact route with probability about
$2\pi/(k\varepsilon)^3$ (it must straddle a capped cut closer than $g_C/\varepsilon$).

**Shifted hierarchy** (`Options::hierarchy = "shifted"`). Three regular quadtrees over the chart, copy $j$
shifted by $j/3$ of the top cell (offset by a tiny irrational amount so that no vertex lies on a grid
line); a cell splits along its two midlines (a line is built only if members lie on both sides).
*Shifting lemma.* At every level the three copies' lines lie at offsets $0, S/3, 2S/3$ modulo the cell
side $S$ (because $2^\ell \not\equiv 0 \bmod 3$); an interval shorter than $S/3$ meets at most one of them per
axis, so for $\|s-t\|_\infty = r$ and the level with $3r < S \le 6r$, one copy keeps $s$ and $t$ in one cell:
the smallest common cell over the copies has side $S \le 6r$. The query takes that copy and the cell's
line that separates $s, t$; then $d \ge \lambda r \ge \lambda S/6$ and (non-near) $d > g/\varepsilon$, so the gap
$\max(g, \varepsilon\lambda S/6)$ is $\le \varepsilon d$ and Lemma 4 holds unchanged. Portals per line are bounded
independently of $S$ once $\varepsilon\lambda S/6 > g$: $O(n\log n/\varepsilon)$ labels and sweep area per copy.

**16-bit coding** (`encode16`). Along a run a row is 1-Lipschitz in arc length, so its step from portal
$j-1$ to $j$ lies in $[-(A_j + e_{j-1}), A_j]$, $e_{j-1}$ the previous decoding error. Each step is stored as
$\lceil \text{step}/q_j\rceil$ in an int16 with $q_j = (A_j + e_{j-1})/32767$ (shared by the rows of the cut),
decoded closed-loop (the decoder's previous value is the reference), with the first portal of each
run stored as a float32 anchor. Every decoded label is $\ge$ the stored one and exceeds it by $< q_j$
($q_j \approx g/32767$, i.e. $\lesssim 2\varepsilon/32767$ relative to $d$): UB stays a path cost, and every interval
term of the LB subtracts the decoding error of its two labels.

**Covers** (`build_covers`). Per row an $\varepsilon_c$-cover $Q$ (§6), chosen greedily outwards from each run's
minimum. The query reconstructs $\ell^+(s,p) = \min_{q\in Q} \ell(s,q) + A(q,p)$ (a path cost, so UB stays
admissible and $\le (1+\varepsilon_c)\,\mathrm{UB}_{\text{dense}}$) and $\ell^-(s,p) = \max_{q\in Q} \ell(s,q) - A(q,p) \le d(s,p)$,
whose interval terms are still valid lower bounds (they only need labels $\le$ distances). The
certificate stays sound; it widens with $\varepsilon_c$.

**POI mode** (`Options::poi`). The tree, cuts, labels and near tables are built over a set $P$ of vertices
only: regions are subsets of $P$, members are POIs, and near tables hold POI pairs within $R$. Nothing in
§3 uses vertices outside the members:
- the window lemma needs only $d(s,t) \le d(s,h) + d(h,t) \le 2\mu_C$ for members $s, t$;
- Lemma 4 needs only the gaps;
- the near-table radius argument is per pair.

So Theorem 1 holds for every POI pair. Measured: sound on 1,000 random POI pairs per run
(`results/poi/`).

**Parallel build** (`Options::threads`). Each portal's sweep, each cut's geometry and each leaf's near
table is computed independently of the others and written to its own slots, so the index is
bit-identical for any thread count (checked on HorseMount: same portals, labels, near entries and
per-pair bounds for 1 and 2 threads). **Coding during the build** (`Options::store`) applies the
coding above to each cut as soon as its labels are final; it gives the same coded index as coding
afterwards.

**Compact near tables** (`Options::near_store`, `encode_near`; default `sym`). A near pair is stored
once, in the row of its smaller vertex, as (varint id delta, distance). The two directions of a pair
come from two exact sweeps, each rounded up to float32. We keep the larger, which is still $\ge d$ and
within float32 of it, so the dense LB $d(1-2^{-22})$ stays valid (`sym`, lossless).

With `sym16` / `sym8` the distance is a 16 / 8-bit code $c$ of its ratio to the chord. The decoded value
is $D(c) = \chi\,(1 + c\,q)$, where:
- $\chi = \lambda|S-T|_{3D} \le d$ is the chord;
- $q$ is a per-row relative quantum (the row's largest $d/\chi - 1$ over 60000 or 240 steps, rounded up to
  float32);
- $c$ is the least code with $D(c) \ge$ the stored value.

One out-of-line function computes $D$ for the encoder and the query, so both round identically. Then
- $\mathrm{UB} = D(c) \ge d$;
- $D(c-1) <$ the stored value (or $c = 0$ and $D = \chi$), so $\mathrm{LB} = \max(\chi, D(c-1)(1-2^{-22})) \le d$;
- the estimate is $\tfrac12(D(c) + D(c-1))$.

Every pair is within one relative quantum of its distance. Measured on 3,000 near pairs on Gunnison10k:

| mode | max UB error | mean / max estimate error |
|---|---:|---|
| `sym16` | 0.001% | — / 0.001% |
| `sym8` | 0.37% | 0.005% / 0.18% |

**Blockwise decoding of coded labels** (`eval_cut`, `block_scan`). The query sums a run's steps in
eight independent blocks and then shifts each block by the totals before it. The additions are
therefore reordered relative to the encoder's running sum, and the decoded value can differ from the
encoder's by a few ulps per step. To keep "decoded ≥ stored", `encode_cut` codes every label for
$L + \delta$ with $\delta = L_{\max}2^{-30}$, which is far larger than $K$ ulps of $L_{\max}$ for any
realistic $K$. Its error bound becomes $q_j + 2\delta$, which the LB subtracts.

**A2A in C++** (`CertifiedIndex::a2a_bounds`). The bounds are those of §7. The exact refine runs MMP from
the surface point itself (Kirsanov's `SurfacePoint` on its face), with no point inserted into the
sub-mesh, so there is no move and no slack term. The early stop stays exact for surface points:
Kirsanov's rule takes each target's closest vertex $v$ and stops only once the queue is beyond
$d(v) + |{\rm target} - v| \ge d({\rm target})$.

**Certified proximity queries** (`cpp/tools/proximity.cpp`). For a query vertex $q$ and objects $O$:
- **Candidates.** Every object with $d \le \tau$ satisfies $\lambda|q - o| \le \tau$, because
  $d \ge \lambda \cdot$ chart distance. So a grid query of radius $\tau/\lambda$ finds all candidates.
- **kNN.** $\tau$ is the $k$-th smallest UB among the objects first found. A candidate is
  - surely in when fewer than $k$ other candidates have LB below its UB;
  - surely out when at least $k$ candidates have UB below its LB.

  The remaining candidates get exact distances from one MMP sweep from $q$ on the faces meeting their
  bounding box dilated by $\max \mathrm{UB}/(2\lambda)$ (ellipse lemma). The top-$k$ is then the
  surely-in set plus the nearest of those exact distances.
- **Range $r$.** $\tau = r$. In if $\mathrm{UB} \le r$, out if $\mathrm{LB} > r$, the rest exact.

Both answer sets are exact, so a wrong answer is a bug. Measured: every answer matched the exact one
(HorseMount, TJpre_1058).

**Certified approximate kNN and range (no refine).**
- **kNN.** Return the $k$ candidates of smallest UB. Every true top-$k$ object $o$ has
  $\mathrm{UB}(o) \le (1+\varepsilon)d(o) \le (1+\varepsilon)d_k$, so the $k$-th smallest UB is
  $\tau \le (1+\varepsilon)d_k$. Every returned object therefore has $d \le \mathrm{UB} \le \tau \le (1+\varepsilon)d_k$:
  the standard $\varepsilon$-approximate kNN.
- **Range $r$.** Return the objects with $\mathrm{UB} \le r$. There are no false positives, and every
  object with $d \le r/(1+\varepsilon)$ is returned.

Both cost one candidate pass, with no exact sweep. Measured on HorseMount:
- kNN, $k = 50$: 27 µs against 0.5 ms for the exact answer and 0.85 ms for the pruned on-the-fly search;
  recall 99.4%, worst $d/d_k = 1.027$.
- Range: no false positives, and nothing missed within $r/(1+\varepsilon)$.

**Weighted terrains in C++** (`--solver steiner`, `src/steiner.h`). One Steiner graph is built for the
terrain (m points per edge, chains at the cheaper adjacent weight, straight in-face segments at
$w_f\times$ length) and shared by every search.
- **Labels.** A portal's labels come from a Dijkstra seeded at the boundary nodes of its face (at
  $w_f\times$ the straight distance), stopped once every member is settled. They are costs of real
  paths, so $\ell(u,p) \ge d(u,p)$ and $\mathrm{UB} \ge d$ holds as before.
- **LB.** The labels are not exact, so the interval LB falls back to the chord bound $\lambda|S-T|$
  (as in Python), and nothing is refined exactly.
- **A2A.** Near pairs take a Steiner search between the two points, which only lowers the UB.
- **Reference.** There is no exact solver, so the truth is a finer Steiner graph (`--truth-steiner M`),
  itself an upper bound. Soundness is checked for the LB only, and a UB below the reference is counted
  as a path the reference missed.
- **Measured.** Weighted HorseMount (the competitors' weights), ε 0.1, against Steiner m = 24: estimate
  mean error 0.38% (m = 5) and 0.23% (m = 10). The Python version gave 0.35% / 0.20% against the best
  path of any method.

**Obstacles as holes** (`scripts/competitors/make_obstacles.py`, server phase `obstacles`). An obstacle
removes the faces whose chart centroid lies inside it, so the surface gets a hole and paths go around it.
The kinds are disks covering a given fraction of the area, and walls with gaps. The certificate needs no
change, for four reasons:
- **Lower bounds.** A hole only lengthens paths, so every chart lower bound (chord, $\lambda|S-T|$)
  still holds.
- **Windows.** $\mu_C$ is a graph distance on the holed mesh, so a detour widens the window by itself,
  and the tails cover the crossing.
- **Cut lines.** A path from one side of a cut line to the other crosses it at a surface point. A line
  broken by holes gives several runs (`line_pieces`), each bracketed by portals at its ends.
- **Near tables.** A table holds every pair within chart radius $R$, with the exact distance on a
  sub-mesh dilated by a graph upper bound (ellipse lemma). A pair close in the chart but separated by a
  wall is therefore exact too, and a pair not in the table has $d \ge \lambda|S-T| > g/\varepsilon$.

Measured on Gunnison10k (disks covering 10% and 20%, six walls), ε 0.2:
- every answer is sound, uncapped and with cap 64, over random, near, band and A2A pairs, and over 208
  pairs with detour ratio $d/|st| \ge 1.2$;
- estimate mean error 0.03–0.04%, the same as without obstacles.

*Two pitfalls found on the way (both handled):*
1. **Touching holes.** Two holes touching at one vertex make it non-manifold (two fans of faces).
   Kirsanov's MMP then returns distances up to 5% too long, longer than a Steiner-graph path, while VTP
   stays exact. The generator therefore splits such vertices, keeping the largest fan, and keeps the
   largest connected piece. On the cleaned meshes MMP and VTP agree to $10^{-15}$, and the phase checks
   that on every obstacle mesh (`obst.onthefly`).
2. **Disconnected sub-meshes.** A sub-mesh cut from a holed terrain can fall into several pieces. VTP
   links its window lists from the first edge only and rejects such a mesh, so the cut falls back to
   MMP. `--exact vtp` (`Options::exact`) takes VTP for every exact computation (truth, labels, near
   tables, refine, A2A; surface points become inserted vertices). It closes the fans at non-manifold
   sub-mesh vertices, which only adds faces, and sweeps each piece on its own; MMP fallbacks are
   counted. On the cleaned meshes it gives the same results as the default solvers.
