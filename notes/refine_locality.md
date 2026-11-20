# Refine locality for the capped index

This note develops a theorem that says which pairs the capped index can refine, what a refine costs,
and how many such pairs there are. It is meant to replace "we measure how often the refine runs" in
the paper (Section 4.4, Theorem 2). Notation is that of `notes/guarantee.md` and
`paper/preprint/main.tex`.

**Status.** Proved on paper and checked against the proof of Theorem 1 and the portal placement in
`CertifiedIndex::place` (`cpp/src/certified.cpp`). Not yet run in a build. Section 6 is a sketch only.

## 1. Setting

The capped index puts at most $k$ portals inside the box of each cut $C$. Its in-box gap is
$$ g_C = \max\big(g,\; A(I_C)/k\big), $$
and consecutive portals of a run satisfy
$$ A_j \le g_C \ \text{ inside } I_C, \qquad A_j \le \max(g_C,\; 2\delta\lambda\, c_j) \ \text{ in the tails}, \tag{1}$$
where $c_j$ is the chart distance from the interval to $I_C$ and $\delta$ is the **tail parameter**.

Today the code uses $\delta = \varepsilon$ [`place`: `k2 = 2.0 * eps * lam`]. The theorem needs
$$ \delta = \frac{\varepsilon}{2+\varepsilon}, \qquad\text{so that}\qquad \frac{1+\delta}{1-\delta} = 1+\varepsilon. $$
For $\varepsilon = 0.2$ this is $\delta = 0.0909$.

A query for a non-near pair is answered at the deepest cut $C$ whose region contains $s$ and $t$,
and is **accepted** if $\mathrm{UB} \le (1+\varepsilon)\,\mathrm{LB}$; otherwise it is **refined** by one
exact sweep on the faces meeting $\mathrm{bbox}(s,t)$ dilated by $\mathrm{UB}/(2\lambda)$
[`bounds_refined`, `exact_pair`].

For a cut $C$ write $b_C = |I_C|$ (the chart length of the cut inside its box), $\sigma_C$ for the
largest cost per unit chart length along the cut, and
$$ w_C = \frac{g_C}{\delta\lambda}, \qquad \beta = \max_C \frac{\sigma_C}{\lambda\,k\,\delta}. $$

## 2. Theorem

**Theorem (refine locality).** Let the capped index be built with tails (1) at
$\delta = \varepsilon/(2+\varepsilon)$, with exact labels, in exact arithmetic. Let $(s,t)$ be a
non-near pair, $C$ its answering cut and $d = d(s,t)$.

1. **Where.** If $d \ge g_C/\delta$ the pair is accepted. Hence a refined pair has $d < g_C/\delta$,
   and $s$ and $t$ lie within chart distance $w_C$ of the cut line $\ell_C$ and of each other.
   In particular only capped cuts ($A(I_C) > kg$) can refine pairs with $d \ge g/\delta$.
2. **Cost.** A refine sweeps the faces meeting a chart box of side less than $(2+\delta)\,w_C$.
3. **Count.** Suppose every chart square of side $r \ge g$ holds at most $\rho_{\max} r^2$ vertices and
   every box $B_M$ of the hierarchy holds at least $\rho_{\min}\,\mathrm{area}(B_M)$ vertices. Then at
   most
   $$ \pi \Big(\frac{\rho_{\max}}{\rho_{\min}}\Big)^{2} \beta^{3}\, n^{2} $$
   vertex pairs can be refined: a fraction $O\big((\sigma/(\lambda k \delta))^{3}\big)$ of all pairs.

**Corollary (expected query time).** For a pair drawn uniformly from the vertex pairs, under the
assumptions of part 3 and with a sweep over $m$ faces costing $S(m)$,
$$ \mathbb E[\text{time}] \;=\; O(\text{depth} + K_C) \;+\; O(\beta^{3})\cdot S\big(O(\beta^{2} n)\big). $$
The worst case over all pairs is $O(\text{depth}+K_C) + S(O(\beta^2 n))$ instead of $S(n)$.

## 3. Proof

**Part 1.** Put $m = \max(g_C, \delta d)$. We show $\mathrm{UB} \le d + m$ and $\mathrm{LB} \ge d - m$.

*Upper bound.* Let $x^\ast$ be a crossing of a shortest path, in the interval $[p_j, p_{j+1}]$ with
arc gap $A_j$. As in the sandwich lemma, $\mathrm{UB} \le \min(f_j, f_{j+1}) \le d + A_j$.
If the interval lies in $I_C$ then $A_j \le g_C$. Otherwise $x^\ast$ has offset at least $c_j$ from
$I_C$, so $d \ge 2\lambda c_j$ by the window lemma, and (1) gives
$A_j \le \max(g_C, 2\delta\lambda c_j) \le \max(g_C, \delta d)$. So $\mathrm{UB} \le d + m$.

*Lower bound.* $\mathrm{LB}$ is at least the minimum of the interval terms
$\tfrac12(f_j + f_{j+1}) - A_j$ over **all** intervals, and every $f \ge d$.

* If $A_j \le g_C$ (every in-box interval, and the tail intervals still at the box gap), the term is
  at least $d - g_C$.
* Otherwise the interval is a tail interval with $A_j \le 2\delta\lambda c_j$. Both of its portals
  have offset at least $c_j$ from $B_M$, and $s,t \in B_M$, so
  $f(p) \ge \lambda(|s-p| + |p-t|) \ge 2\lambda c_j$. With $f \ge d$ as well, the term is at least
  $\max(d, 2\lambda c_j) - 2\delta\lambda c_j \ge (1-\delta)\max(d, 2\lambda c_j) \ge (1-\delta)\,d$.

So every term is at least $d - m$. (End-stretch terms $f_k - 2e_k$ are at least $d - 2e_k$ with
$2e_k \ll g$, as in the appendix of the paper.)

*Conclusion.* If $d \ge g_C/\delta$ then $m = \delta d$, so
$\mathrm{UB}/\mathrm{LB} \le (1+\delta)/(1-\delta) = 1+\varepsilon$ and the pair is accepted.
For a refined pair, $d < g_C/\delta$. By the chart bound $\lambda|s-t| \le d$, so
$|s-t| < w_C$; and $s,t$ lie on opposite sides of $\ell_C$, so
$\operatorname{dist}(s,\ell_C) + \operatorname{dist}(t,\ell_C) \le |s-t| < w_C$.
For an uncapped cut $g_C = g$, so a refined pair there has $d < g/\delta$. $\qquad\blacksquare$

**Part 2.** For a refined pair, $\mathrm{UB} \le d + g_C < g_C(1 + 1/\delta)$. The sweep box is
$\mathrm{bbox}(s,t)$ dilated by $\mathrm{UB}/(2\lambda)$, of side at most
$|s-t| + \mathrm{UB}/\lambda < w_C + (1+\delta)\,w_C = (2+\delta)\,w_C$. $\qquad\blacksquare$

**Part 3.** Fix a cut $C$. By part 1 a refinable pair has one endpoint, say $s$, in the strip of
width $w_C$ on one side of $\ell_C$ inside $B_M$, and $t$ in the half-disc of radius $w_C$ around $s$
on the other side. There are at most $\rho_{\max}\, w_C\, b_C$ such $s$ and at most
$\rho_{\max}\,\pi w_C^2/2$ such $t$ for each, so at most
$$ \tfrac{\pi}{2}\,\rho_{\max}^2\, b_C\, w_C^3 $$
pairs. Since $g_C = A(I_C)/k \le \sigma_C b_C/k$ on a capped cut, $w_C \le \beta\, b_C$, and the
count is at most $\tfrac{\pi}{2}\rho_{\max}^2\beta^3 b_C^4$.

Now sum over cuts. The boxes $B_M$ at one depth $\ell$ have disjoint interiors (siblings lie in the
two open half-planes of their parent's cut), so $\sum_{C \text{ at depth } \ell}\mathrm{area}(B_M) \le \mathcal A$,
the area of the chart's bounding box. The cut runs across the shorter side, so
$b_C^2 \le \mathrm{area}(B_M) \le |M|/\rho_{\min} \le n/(2^\ell \rho_{\min})$ with median cuts. Hence
$$ \sum_{C \text{ at depth } \ell} b_C^4 \;\le\; \frac{n}{2^\ell\rho_{\min}}\,\mathcal A,
   \qquad \sum_C b_C^4 \;\le\; \frac{2n\mathcal A}{\rho_{\min}} \;\le\; \frac{2n^2}{\rho_{\min}^2}, $$
using $n \ge \rho_{\min}\mathcal A$. Multiplying by $\tfrac{\pi}{2}\rho_{\max}^2\beta^3$ gives the
bound. $\qquad\blacksquare$

**Corollary.** A refinable pair at $C$ costs $S(O(\rho_{\max} w_C^2)) = S(O(\beta^2 \rho_{\max} b_C^2))$
by part 2, and both the number of pairs ($\propto b_C^4$) and the cost are largest at the top cuts,
where $\rho_{\max} b_C^2 = O(n)$.

## 4. What it gives, and how tight it is

**What it gives.**

* A statement where the paper now has a measurement: refines are confined to a band around the capped
  cuts, and a refine is a local sweep.
* A scaling law to test: the refinable fraction falls at least like $k^{-3}$ once $k\delta > 1$.
* A design rule. The expected refine cost is $O(\beta^5 n\log n)$ under $S(m) = O(m \log m)$, so it
  is $O(1)$ once $k \gtrsim (\sigma/\lambda\delta)\,(n\log n)^{1/5}$. To check before claiming: with
  such a $k$ the labels are $O(n^{6/5}\,\mathrm{polylog})$.

**How tight.** At $k = 64$, $\varepsilon = 0.2$ and $\sigma/\lambda \approx 1.1$:
$\delta = 0.091$, $k\delta = 5.8$, $\beta \approx 0.19$. On a uniform grid the exact constant is
$2\beta^3$ of all pairs (the count is $\tfrac23\rho^2 b_C w_C^3$ per cut and the sum over depths is
$1.5\,L^4$), so the bound allows about 1.4% of the pairs to refine, each on up to roughly
$(2+\delta)^2\beta^2 \approx 16\%$ of the region. The measurements show none of the 2,000 random
pairs refining at cap 64. At $k = 256$ the bound is about 0.02%.

**Why it is loose.** The interval lower bound really is about $d - g_C$ for a typical pair near a
capped cut. Such pairs are accepted anyway because $\mathrm{LB}$ also takes the chord
$\lambda|S-T|$, and on these terrains $d$ is within $1+\varepsilon$ of the chord unless the path
detours. A refined pair therefore needs both
$$ d < g_C/\delta \qquad\text{and}\qquad \mathrm{UB} > (1+\varepsilon)\,\lambda|S-T|. $$
The second condition is rigorous but not an a-priori statement about the terrain: turning it into one
needs a second-order bound on $\mathrm{UB} - d$, which does not hold cleanly on a polyhedral surface
(the distance function is not uniformly semiconcave near saddle vertices). Report it as a
measurement instead (Section 5).

## 5. What to change and what to measure

**Code.** `eps` is used in three places that the theorem separates:

| use | where | value for the theorem |
| --- | --- | --- |
| near radius $R = g/(\varepsilon\lambda)$ | `CertifiedIndex` constructor, `Rn` | unchanged |
| tail spacing | `place`, `k2 = 2.0 * eps * lam` | $\delta = \varepsilon/(2+\varepsilon)$ |
| acceptance $\mathrm{UB} \le (1+\varepsilon)\mathrm{LB}$ | `bounds_refined` | unchanged |

So add one option for the tail parameter. The tails get about $(2+\varepsilon)$ times denser near the
box; expect the capped index to grow by some tens of percent. Measure it.

With the near radius left at $g/(\varepsilon\lambda)$, pairs with $g/\varepsilon < d < g/\delta$ on
uncapped cuts are inside the theorem's refinable set (part 1 with $g_C = g$). Their sweeps cover a
box of side about $(2+\delta)g/(\delta\lambda)$, some 23 edges at $\varepsilon = 0.2$, $g = 1$.

**Experiments.**

1. Rebuild cap 64 on TJ and GI with $\delta$ tails; report size, build time and errors next to the
   current build.
2. For every refined pair, over all caps and query sets, log $d/g_C$ and $d/(\lambda|S-T|)$. Check
   $d < g_C/\delta$ on all of them (the theorem) and show the stretch distribution (the explanation).
3. Refined fraction of random pairs against $k$ for $k = 8,\dots,256$, with the $k^{-3}$ bound.
4. Report $\sigma_C/\lambda$, $\kappa_C$ and $r_C$ per terrain, including the adversarial ones, so the
   assumptions of the capped theorem are numbers rather than words.

## 6. A stronger version without sweeps (sketch, not checked in full)

The refine can be removed altogether by giving each capped cut extra rows at finer gaps, stored only
where they are needed.

**Structure.** On a capped cut, place portal levels $i = 0,1,\dots$ at gaps $\gamma_i = 2^i g$ up to
$g_C$. A vertex $u$ stores its level-$i$ labels only for the portals within chart distance
$W_i = 3\gamma_i/(\delta\lambda)$ of its foot on $\ell_C$, and only if
$\operatorname{dist}(u,\ell_C) \le 2\gamma_i/(\delta\lambda)$. That is $O(\sigma/(\delta\lambda))$
entries per vertex and level.

**Query.** Go from the capped row to finer levels. At level $i$, take $\mathrm{UB}_i$ as the minimum of
$f$ over the portals both endpoints store; it is a path length, so it is always valid. A shortest
path crosses $\ell_C$ within chart distance $d/\lambda \le \mathrm{UB}_i/\lambda$ of both $s$ and $t$,
so if the common window covers both discs of radius $\mathrm{UB}_i/\lambda$ on the line, the minimum of
the interval terms over that window is a valid $\mathrm{LB}_i$ (sandwich lemma). Accept if
$\mathrm{UB}_i \le (1+\varepsilon)\mathrm{LB}_i$, else go one level finer.

**Why it terminates.** At the level with $\gamma_i \le \delta d < 2\gamma_i$:
$\mathrm{UB}_i \le d + \gamma_i \le (1+\delta)d$, the discs have radius at most
$(1+\delta)\cdot 2\gamma_i/(\delta\lambda) \le W_i$, and the argument of Section 3 gives
$\mathrm{LB}_i \ge (1-\delta)d$.

**Cost, roughly.** Per capped cut the extra labels are about
$\rho\, b_C \int_0^{\beta b_C} \delta^{-1}\log(\beta b_C/h)\,dh = O(\beta\rho b_C^2/\delta)$, so
$O\big(\sigma\, n\log n/(\lambda k\delta^2)\big)$ in total, some tens of percent of the capped index at
$k = 64$. The extra sweeps are local (radius $W_i$), so the build grows by a similar factor.

**What it would give.** A worst-case $(1+\varepsilon)$ certificate with no exact sweep at query time,
query time $O(K_C + \delta^{-1}\log(g_C/g))$, in near-linear space under the same assumptions.

**Open points.** Runs broken by obstacles (a window must hold the ends of every run it meets); the
exact constant in $W_i$; whether the coarse-to-fine walk or a direct level choice from
$\lambda|S-T|$ is faster; and the whole thing needs implementing and measuring.
