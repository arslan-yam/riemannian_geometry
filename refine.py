import numpy as np
from scipy.optimize import minimize
from graph_build import blocked
from obstacles import Circle, Rectangle


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


def barrier_and_grad(coords, obstacles, margin, weight):
    val = 0.0
    grad = np.zeros_like(coords)
    
    for i in range(1, len(coords) - 1):
        c, g = clearance_grad(coords[i], obstacles)
        if c < margin:
            val += (margin - c) ** 2
            grad[i] += -2.0 * (margin - c) * g
            
    for i in range(len(coords) - 1):
        c, g = clearance_grad(0.5 * (coords[i] + coords[i + 1]), obstacles)
        if c < margin:
            val += (margin - c) ** 2
            gm = -2.0 * (margin - c) * g
            grad[i] += 0.5 * gm
            grad[i + 1] += 0.5 * gm
    return weight * val, weight * grad[1:-1].ravel()


def free_segment(manifold, a, b, obstacles, step):
    m = max(4, int(np.ceil(np.linalg.norm(b - a) / step)))
    for x in np.linspace(0.0, 1.0, m + 1):
        q = (1.0 - x) * a + x * b
        if not manifold.in_domain(q) or blocked(q, obstacles):
            return False
    return True


def polyline_lenght(manifold, coords):
    return sum(manifold.segment_lenght(coords[i], coords[i + 1]) for i in range(len(coords) - 1))


def free_polyline(manifold, coords, obstacles, step):
    return all(free_segment(manifold, coords[i], coords[i + 1], obstacles, step) for i in range(len(coords) - 1))


def shortcut(manifold, coords, obstacles, step, window=30):
    cum = np.zeros(len(coords))
    for i in range(len(coords) - 1):
        cum[i + 1] = cum[i] + manifold.segment_lenght(coords[i], coords[i + 1])
    out = [coords[0]]
    i = 0
    
    while i < len(coords) - 1:
        best = i + 1
        for j in range(min(len(coords) - 1, i + window), i + 1, -1):
            if free_segment(manifold, coords[i], coords[j], obstacles, step):
                if manifold.segment_lenght(coords[i], coords[j]) < cum[j] - cum[i]:
                    best = j
                    break
        out.append(coords[best])
        i = best
        
    return np.array(out)


def subdivide(coords, max_step):
    out = [coords[0]]
    for i in range(len(coords) - 1):
        a, b = coords[i], coords[i + 1]
        m = max(1, int(np.ceil(np.linalg.norm(b - a) / max_step)))
        for u in range(1, m + 1):
            out.append(a + (b - a) * u / m)
    return np.array(out)


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
        e, ge = manifold.energy(pts, m), manifold.energy_grad(pts, m)
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


def refine_path(manifold, points, path, obstacles=(), rounds=2, window=30, segments=20, step=None, use_relax=True, eps=None, margin=None, weight=1e4):
    if path is None or len(path) < 2:
        return 0.0, None
    coords = np.array([np.asarray(points[i], dtype=float) for i in path])
    if step is None:
        edges = [np.linalg.norm(coords[i + 1] - coords[i]) for i in range(len(coords) - 1)]
        step = 0.25 * float(np.median(edges)) if edges else 1.0
    if eps is None:
        eps = 0.1 * step
    if margin is None:
        margin = 0.0

    best = coords if free_polyline(manifold, coords, obstacles, step) else None
    for _ in range(rounds):
        coords = shortcut(manifold, coords, obstacles, step, window)
        if use_relax:
            size = polyline_lenght(manifold, coords) / segments
            if size > 0:
                coords = subdivide(coords, size)
            coords = relax(manifold, coords, obstacles, step, eps, margin, weight)
        if free_polyline(manifold, coords, obstacles, step):
            if best is None or polyline_lenght(manifold, coords) < polyline_lenght(manifold, best):
                best = coords
                
    if best is None:
        best = np.array([np.asarray(points[i], dtype=float) for i in path])
    return polyline_lenght(manifold, best), best
