"""Continuum post-processing of graph paths: shortcut, subdivide, energy relaxation (L-BFGS),
obstacle barrier; and the batched variant that relaxes many polylines in one optimiser call."""
import numpy as np
from scipy.optimize import minimize
from .obstacles import Circle, Rectangle
from .manifolds import gauss_legendre


def project_free(coords, obstacles, eps):
    out = np.array(coords, dtype=float)
    
    for i in range(1, len(out) - 1):
        for o in obstacles:
            p = out[i]
            
            if isinstance(o, Circle):
                d = np.hypot(p[0] - o.cx, p[1] - o.cy)
                if d < o.r + eps:
                    if d < 1e-12:
                        out[i] = np.array([o.cx + o.r + eps, o.cy])
                    else:
                        out[i] = np.array([o.cx, o.cy]) + (p - np.array([o.cx, o.cy])) * (o.r + eps) / d
                        
            elif isinstance(o, Rectangle):
                if o.sx - eps <= p[0] <= o.ex + eps and o.sy - eps <= p[1] <= o.ey + eps:
                    left, right = p[0] - (o.sx - eps), (o.ex + eps) - p[0]
                    down, up = p[1] - (o.sy - eps), (o.ey + eps) - p[1]
                    m = min(left, right, down, up)
                    if m == left:
                        out[i][0] = o.sx - eps
                    elif m == right:
                        out[i][0] = o.ex + eps
                    elif m == down:
                        out[i][1] = o.sy - eps
                    else:
                        out[i][1] = o.ey + eps
    return out


def clearance_grad(p, obstacles):
    best, grad = np.inf, np.zeros(2)
    for o in obstacles:
        
        if isinstance(o, Circle):
            v = p - np.array([o.cx, o.cy])
            d = np.hypot(v[0], v[1])
            c = d - o.r
            g = v / d if d > 1e-12 else np.array([1.0, 0.0])
            
        elif isinstance(o, Rectangle):
            dx = max(o.sx - p[0], 0.0, p[0] - o.ex)
            dy = max(o.sy - p[1], 0.0, p[1] - o.ey)
            
            if dx > 0.0 or dy > 0.0:
                c = np.hypot(dx, dy)
                sx = -1.0 if p[0] < o.sx else (1.0 if p[0] > o.ex else 0.0)
                sy = -1.0 if p[1] < o.sy else (1.0 if p[1] > o.ey else 0.0)
                g = np.array([sx * dx, sy * dy]) / c if c > 1e-12 else np.zeros(2)
            else:
                left, right = p[0] - o.sx, o.ex - p[0]
                down, up = p[1] - o.sy, o.ey - p[1]
                c = -min(left, right, down, up)
                if -c == left:
                    g = np.array([-1.0, 0.0])
                elif -c == right:
                    g = np.array([1.0, 0.0])
                elif -c == down:
                    g = np.array([0.0, -1.0])
                else:
                    g = np.array([0.0, 1.0])
                    
        else:
            continue
        if c < best:
            best, grad = c, g
    return best, grad


def clearance_batch(points, obstacles):
    #min over obstacles of the signed clearance, and the gradient of that min, vectorised
    best = np.full(len(points), np.inf)
    grad = np.zeros((len(points), 2))
    for o in obstacles:
        c, g = o.clearance_batch(points)
        take = c < best
        best = np.where(take, c, best)
        grad = np.where(take[:, None], g, grad)
    return best, grad


def barrier_and_grad(coords, obstacles, margin, weight):
    #same penalty as before, on the interior vertices and on every segment midpoint, but every
    #clearance is evaluated in one batched call per obstacle instead of a python loop per point
    n = len(coords)
    grad = np.zeros_like(coords)
    val = 0.0
    if n > 2:
        c, g = clearance_batch(coords[1:-1], obstacles)
        hit = c < margin
        if hit.any():
            d = margin - c[hit]
            val += float((d * d).sum())
            grad[1:-1][hit] += -2.0 * d[:, None] * g[hit]
    mid = 0.5 * (coords[:-1] + coords[1:])
    c, g = clearance_batch(mid, obstacles)
    hit = c < margin
    if hit.any():
        d = margin - c[hit]
        val += float((d * d).sum())
        gm = -2.0 * d[:, None] * g[hit]
        idx = np.nonzero(hit)[0]
        np.add.at(grad, idx, 0.5 * gm)
        np.add.at(grad, idx + 1, 0.5 * gm)
    return weight * val, weight * grad[1:-1].ravel()


def _free_points(manifold, pts, obstacles):
    ok = manifold.in_domain_batch(pts) if hasattr(manifold, "in_domain_batch") \
        else np.array([manifold.in_domain(q) for q in pts], dtype=bool)
    for o in obstacles:
        if not ok.any():
            break
        ok &= ~o.blocked_batch(pts)
    return bool(ok.all())


def free_segment(manifold, a, b, obstacles, step):
    m = max(4, int(np.ceil(np.linalg.norm(b - a) / step)))
    t = np.linspace(0.0, 1.0, m + 1)[:, None]
    return _free_points(manifold, (1.0 - t) * a + t * b, obstacles)


def segment_lengths(manifold, A, B, n=4):
    #lengths of chart-straight segments A[i]->B[i]. a manifold that can measure them exactly
    #(a polyhedral terrain) provides segment_lengths_batch; otherwise gauss-legendre on the metric
    A = np.asarray(A, dtype=float); B = np.asarray(B, dtype=float)
    if hasattr(manifold, "segment_lengths_batch"):
        return manifold.segment_lengths_batch(A, B)
    delta = B - A
    ts, ws = gauss_legendre(n)
    pts = (A[:, None, :] + ts[None, :, None] * delta[:, None, :]).reshape(-1, 2)
    g = manifold.metric_batch(pts).reshape(len(delta), n, 2, 2)
    q = np.einsum("ij,itjk,ik->it", delta, g, delta)
    return (np.sqrt(np.maximum(q, 0.0)) * ws[None, :]).sum(1)


def polyline_lenght(manifold, coords, n=4):
    coords = np.asarray(coords, dtype=float)
    if len(coords) < 2:
        return 0.0
    return float(segment_lengths(manifold, coords[:-1], coords[1:], n).sum())


def free_polyline(manifold, coords, obstacles, step):
    #sample every segment, then one batched free-space test over all samples
    coords = np.asarray(coords, dtype=float)
    if len(coords) < 2:
        return True
    delta = coords[1:] - coords[:-1]
    counts = np.maximum(4, np.ceil(np.linalg.norm(delta, axis=1) / step).astype(int))
    pieces = []
    for i, m in enumerate(counts):
        t = np.linspace(0.0, 1.0, m + 1)[:, None]
        pieces.append(coords[i] + t * delta[i])
    return _free_points(manifold, np.vstack(pieces), obstacles)


def shortcut(manifold, coords, obstacles, step, window=30):
    #greedy farthest shortcut: from vertex i jump to the farthest j (within the window) whose
    #chord is free and shorter than the path between them. all candidate chords are measured
    #in ONE batched call up front; the free-space test (the expensive sampled one) only runs
    #for chords that are actually shorter. same decisions as the pairwise loop, ~100x fewer
    #length calls, which matters when every length is an exact face walk on a terrain
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if n < 3:
        return coords.copy()
    cum = np.concatenate([[0.0], np.cumsum(segment_lengths(manifold, coords[:-1], coords[1:]))])
    I, J = [], []
    for i in range(n - 2):
        for j in range(i + 2, min(n - 1, i + window) + 1):
            I.append(i); J.append(j)
    I = np.array(I); J = np.array(J)
    chord = segment_lengths(manifold, coords[I], coords[J])
    cand = np.nonzero(chord < cum[J] - cum[I])[0]
    #free-space test of every shorter chord at once: the same samples free_segment would take,
    #concatenated, one batched domain/obstacle query, then reduced per chord
    usable = {}
    if len(cand):
        A, B = coords[I[cand]], coords[J[cand]]
        counts = np.maximum(4, np.ceil(np.linalg.norm(B - A, axis=1) / step).astype(int)) + 1
        owner = np.repeat(np.arange(len(cand)), counts)
        t = (np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)) / np.repeat(counts - 1, counts)
        pts = A[owner] + t[:, None] * (B - A)[owner]
        ok = manifold.in_domain_batch(pts) if hasattr(manifold, "in_domain_batch") \
            else np.array([manifold.in_domain(q) for q in pts], dtype=bool)
        for o in obstacles:
            ok &= ~o.blocked_batch(pts)
        bad = np.bincount(owner, weights=~ok, minlength=len(cand)) > 0
        for c, b in zip(cand, bad):
            if not b:
                usable.setdefault(int(I[c]), set()).add(int(J[c]))
    out = [coords[0]]
    i = 0
    while i < n - 1:
        best = i + 1
        for j in range(min(n - 1, i + window), i + 1, -1):
            if j in usable.get(i, ()):
                best = j
                break
        out.append(coords[best])
        i = best
    return np.array(out)


def subdivide(manifold, coords, max_step):
    #max_step comes from polyline_lenght, which is a RIEMANNIAN length, so the segment must be
    #measured the same way. using the coordinate norm under-resolves wherever the metric is
    #large: on the poincare disk at r=0.9 that is a factor of 10
    coords = np.asarray(coords, dtype=float)
    delta = coords[1:] - coords[:-1]
    seg = segment_lengths(manifold, coords[:-1], coords[1:])
    out = [coords[0]]
    for i in range(len(coords) - 1):
        m = max(1, int(np.ceil(seg[i] / max_step)))
        u = np.arange(1, m + 1)[:, None] / m
        out.append(coords[i] + delta[i] * u)
    return np.vstack([out[0][None, :]] + out[1:])


def relax(manifold, coords, obstacles, step, eps=1e-4, margin=0.0, weight=1e4):
    m = len(coords) - 1
    if m < 2:
        return coords
    p, q = coords[0], coords[-1]
    start = coords[1:-1].ravel()
    if not obstacles:
        res = minimize(lambda z: manifold.energy_and_grad(p, q, z, m), start, method="L-BFGS-B", jac=True)
        out = np.vstack([p, res.x.reshape(-1, 2), q])
        return out if free_polyline(manifold, out, obstacles, step) else coords

    def objective(z):
        pts = np.vstack([p, z.reshape(-1, 2), q])
        e, ge = manifold.energy_and_grad_points(pts, m) #one metric evaluation, not two
        b, gb = barrier_and_grad(pts, obstacles, margin, weight)
        return e + b, ge + gb

    res = minimize(objective, start, method="L-BFGS-B", jac=True)
    out = np.vstack([p, res.x.reshape(-1, 2), q])
    if free_polyline(manifold, out, obstacles, step):
        return out
    for a in (0.5, 0.25, 0.125, 0.0625):
        trial = project_free(coords + a * (out - coords), obstacles, eps)
        if free_polyline(manifold, trial, obstacles, step):
            return trial
    return coords


def refine_polyline(manifold, coords, obstacles=(), rounds=2, window=30, segments=20, step=None, use_relax=True, eps=None, margin=None, weight=1e4):
    #same as refine_path but starting from coordinates, so the endpoints need not be graph nodes
    return _refine(manifold, np.asarray(coords, dtype=float), obstacles, rounds, window, segments,
                   step, use_relax, eps, margin, weight)


def refine_path(manifold, points, path, obstacles=(), rounds=2, window=30, segments=20, step=None, use_relax=True, eps=None, margin=None, weight=1e4):
    if path is None or len(path) < 2:
        return 0.0, None
    coords = np.array([np.asarray(points[i], dtype=float) for i in path])
    return _refine(manifold, coords, obstacles, rounds, window, segments, step, use_relax, eps, margin, weight)


def _refine(manifold, coords, obstacles, rounds, window, segments, step, use_relax, eps, margin, weight):
    if len(coords) < 2:
        return 0.0, coords
    given = coords #the loop rebinds coords, so keep the input around for the fallback
    if step is None: #median of the POSITIVE edges: a duplicated point would give step 0
        edges = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
        positive = edges[edges > 0]
        step = 0.25 * float(np.median(positive)) if len(positive) else 1.0
    if eps is None:
        eps = 0.1 * step
    if margin is None:
        margin = 0.0

    best = coords if free_polyline(manifold, coords, obstacles, step) else None
    for _ in range(rounds):
        coords = shortcut(manifold, coords, obstacles, step, window)
        if use_relax:
            size = polyline_lenght(manifold, coords) / segments
            if size > 0 and np.isfinite(size):
                coords = subdivide(manifold, coords, size)
            coords = relax(manifold, coords, obstacles, step, eps, margin, weight)
        if free_polyline(manifold, coords, obstacles, step):
            if best is None or polyline_lenght(manifold, coords) < polyline_lenght(manifold, best):
                best = coords
                
    if best is None: #nothing was free, so hand back what we were given rather than a bad path
        best = given
    return polyline_lenght(manifold, best), best


# ----------------------------------------------------------------------------------------------
# batched refinement: many polylines, one optimiser
#
# relax() spends its time in ~50 L-BFGS evaluations per polyline, each a handful of numpy calls
# on a 12-point array, so the cost is interpreter overhead, not arithmetic. resampling every
# polyline to the same segment count lets B of them stack into one (B, m+1, 2) array; the energy
# is separable so the gradient just concatenates, and one L-BFGS over the stacked vector does
# B problems per numpy call. the result is the same fixed point up to optimiser tolerance
# ----------------------------------------------------------------------------------------------

def resample_uniform(manifold, coords, m):
    #m equal riemannian-arc-length pieces along the polyline, endpoints kept
    coords = np.asarray(coords, dtype=float)
    if len(coords) < 2:
        return np.repeat(coords[:1], m + 1, axis=0)
    seg = segment_lengths(manifold, coords[:-1], coords[1:])
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    if cum[-1] <= 0.0:
        return np.repeat(coords[:1], m + 1, axis=0)
    want = np.linspace(0.0, cum[-1], m + 1)
    return np.stack([np.interp(want, cum, coords[:, 0]), np.interp(want, cum, coords[:, 1])], axis=1)


def energy_and_grad_stack(manifold, P, m):
    #P: (B, m+1, 2). energy summed over the batch, gradient w.r.t. the interior vertices (B, m-1, 2)
    if hasattr(manifold, "energy_and_grad_stack"): #exact polyhedral energy on a terrain mesh
        return manifold.energy_and_grad_stack(P, m)
    B = P.shape[0]
    delta = P[:, 1:] - P[:, :-1]
    mid = 0.5 * (P[:, 1:] + P[:, :-1])
    g = manifold.metric_batch(mid.reshape(-1, 2)).reshape(B, m, 2, 2)
    gd = np.einsum("bijk,bik->bij", g, delta)
    e = m * float(np.einsum("bij,bij->", delta, gd))
    gx, gy = manifold.metric_grad_batch(mid.reshape(-1, 2))
    gx = gx.reshape(B, m, 2, 2); gy = gy.reshape(B, m, 2, 2)
    gmid = 0.5 * np.stack([np.einsum("bij,bijk,bik->bi", delta, gx, delta),
                           np.einsum("bij,bijk,bik->bi", delta, gy, delta)], axis=2)
    grad = np.zeros_like(P)
    grad[:, 1:] += 2.0 * gd + gmid
    grad[:, :-1] += -2.0 * gd + gmid
    return e, m * grad[:, 1:-1]


def barrier_and_grad_stack(P, obstacles, margin, weight):
    B, n, _ = P.shape
    grad = np.zeros_like(P)
    val = 0.0
    if n > 2:
        inner = P[:, 1:-1].reshape(-1, 2)
        c, g = clearance_batch(inner, obstacles)
        hit = c < margin
        if hit.any():
            d = margin - c[hit]
            val += float((d * d).sum())
            gi = np.zeros_like(inner)
            gi[hit] = -2.0 * d[:, None] * g[hit]
            grad[:, 1:-1] += gi.reshape(B, n - 2, 2)
    mid = (0.5 * (P[:, :-1] + P[:, 1:])).reshape(-1, 2)
    c, g = clearance_batch(mid, obstacles)
    hit = c < margin
    if hit.any():
        d = margin - c[hit]
        val += float((d * d).sum())
        gm = np.zeros_like(mid)
        gm[hit] = -2.0 * d[:, None] * g[hit]
        gm = gm.reshape(B, n - 1, 2)
        grad[:, :-1] += 0.5 * gm
        grad[:, 1:] += 0.5 * gm
    return weight * val, weight * grad[:, 1:-1]


def relax_stack(manifold, P, obstacles, margin=0.0, weight=1e4, maxiter=2000):
    #P: (B, m+1, 2) with endpoints fixed. returns the relaxed stack
    B, n, _ = P.shape
    m = n - 1
    if m < 2:
        return P
    ends0, ends1 = P[:, 0], P[:, -1]

    def assemble(z):
        Q = np.empty_like(P)
        Q[:, 0] = ends0; Q[:, -1] = ends1
        Q[:, 1:-1] = z.reshape(B, m - 1, 2)
        return Q

    def objective(z):
        Q = assemble(z)
        e, ge = energy_and_grad_stack(manifold, Q, m)
        if obstacles:
            b, gb = barrier_and_grad_stack(Q, obstacles, margin, weight)
            return e + b, (ge + gb).ravel()
        return e, ge.ravel()

    res = minimize(objective, P[:, 1:-1].ravel(), method="L-BFGS-B", jac=True,
                   options=dict(maxiter=maxiter, maxcor=10))
    return assemble(res.x)


def refine_batch(manifold, polylines, obstacles=(), rounds=1, window=30, segments=12, step=None,
                 use_relax=True, eps=None, margin=None, weight=1e4):
    #same rounds as _refine, but every polyline in the batch is relaxed in one optimiser call.
    #returns a list of (length, coords), one per input polyline, in order
    polylines = [np.asarray(c, dtype=float) for c in polylines]
    B = len(polylines)
    if B == 0:
        return []
    #per-polyline step, exactly as _refine picks it, so free-space sampling and shortcut
    #decisions match the sequential path polyline by polyline
    def pick_step(c):
        if step is not None:
            return step
        e = np.linalg.norm(c[1:] - c[:-1], axis=1) if len(c) > 1 else np.zeros(0)
        pos = e[e > 0]
        return 0.25 * float(np.median(pos)) if len(pos) else 1.0
    steps = [pick_step(c) for c in polylines]
    epss = [0.1 * st if eps is None else eps for st in steps]
    if margin is None:
        margin = 0.0

    best = [c if free_polyline(manifold, c, obstacles, st) else None for c, st in zip(polylines, steps)]
    cur = list(polylines)
    for _ in range(rounds):
        cur = [shortcut(manifold, c, obstacles, st, window) if len(c) > 2 else c for c, st in zip(cur, steps)]
        if use_relax:
            #the SAME subdivision as _refine, so each polyline gets the segment count it would get
            #sequentially; then stack polylines with equal counts and relax each group at once
            #the sequential subdivide count, rounded UP to a multiple of 4: never fewer segments
            #than _refine would use (so never a coarser discretisation), and 4x fewer groups
            subs = {}
            for i, c in enumerate(cur):
                if len(c) < 2 or np.allclose(c[0], c[-1]):
                    continue
                size = polyline_lenght(manifold, c) / segments
                if not np.isfinite(size): #a segment left the domain (e.g. the mesh): leave it
                    continue
                count = len(subdivide(manifold, c, size)) - 1 if size > 0 else 1
                count = max(2, -(-count // 4) * 4)
                subs.setdefault(count, []).append(i)
            for count, idx in subs.items():
                P = np.stack([resample_uniform(manifold, cur[i], count) for i in idx])
                #a resampled polyline can leave the domain (jagged mesh boundary); it would
                #make the stacked objective infinite and stop the optimiser for everyone
                fin = np.isfinite([polyline_lenght(manifold, p) for p in P])
                if not fin.all():
                    idx = [i for i, f in zip(idx, fin) if f]
                    if not idx:
                        continue
                    P = P[fin]
                R = relax_stack(manifold, P, obstacles, margin, weight)
                for k, i in enumerate(idx):
                    out = R[k]
                    if free_polyline(manifold, out, obstacles, steps[i]):
                        cur[i] = out
                    else:
                        base = P[k]
                        for a in (0.5, 0.25, 0.125, 0.0625):
                            trial = project_free(base + a * (out - base), obstacles, epss[i])
                            if free_polyline(manifold, trial, obstacles, steps[i]):
                                cur[i] = trial
                                break
                        else:
                            cur[i] = base
        for i, c in enumerate(cur):
            if free_polyline(manifold, c, obstacles, steps[i]):
                if best[i] is None or polyline_lenght(manifold, c) < polyline_lenght(manifold, best[i]):
                    best[i] = c
    result = []
    for i in range(B):
        b = best[i] if best[i] is not None else polylines[i]
        result.append((polyline_lenght(manifold, b), b))
    return result
