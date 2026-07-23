import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
from obstacles import Circle, Rectangle
from manifolds import Manifold
from graph_build import blocked, segment_free

decay_rate = 0.04 #for obstacle distance field normalization

def metric_grad_norm(manifold: Manifold, p, h=5e-3):
    dx, dy = np.array([h, 0.0]), np.array([0.0, h])
    fx = (manifold.metric(p + dx) - manifold.metric(p - dx)) / (2 * h)
    fy = (manifold.metric(p + dy) - manifold.metric(p - dy)) / (2 * h)
    return float(np.sqrt(np.sum(fx * fx) + np.sum(fy * fy)))
    
    
def obstacle_distance(p, obstacles):
    dist = np.inf
    for o in obstacles:
        if isinstance(o, Circle):
            dist = min(dist, abs(np.hypot(p[0] - o.cx, p[1] - o.cy) - o.r))
        elif isinstance(o, Rectangle):
            dx = max(o.sx - p[0], 0.0, p[0] - o.ex)
            dy = max(o.sy - p[1], 0.0, p[1] - o.ey)
            dist = min(dist, np.hypot(dx, dy))
    return dist


def gaussian_curvature(manifold, p, h=1e-2): #approximation using Brioschi formula
    g = {(i, j): manifold.metric(p + np.array([i * h, j * h])) for i in (-1, 0, 1) for j in (-1, 0, 1)}
    E = lambda i, j: g[(i, j)][0, 0]
    F = lambda i, j: g[(i, j)][0, 1]
    G = lambda i, j: g[(i, j)][1, 1]
    E0, F0, G0 = E(0, 0), F(0, 0), G(0, 0)
    
    Eu = (E(1, 0) - E(-1, 0)) / (2 * h)
    Ev = (E(0, 1) - E(0, -1)) / (2 * h)
    Gu = (G(1, 0) - G(-1, 0)) / (2 * h)
    Gv = (G(0, 1) - G(0, -1)) / (2 * h)
    Fu = (F(1, 0) - F(-1, 0)) / (2 * h)
    Fv = (F(0, 1) - F(0, -1)) / (2 * h)
    Evv = (E(0, 1) - 2 * E0 + E(0, -1)) / (h * h)
    Guu = (G(1, 0) - 2 * G0 + G(-1, 0)) / (h * h)
    Euv = (E(1, 1) - E(1, -1) - E(-1, 1) + E(-1, -1)) / (4 * h * h)
    Fuv = (F(1, 1) - F(1, -1) - F(-1, 1) + F(-1, -1)) / (4 * h * h)
    
    A = np.array([[-0.5 * Evv + Fuv - 0.5 * Guu, 0.5 * Eu, Fu - 0.5 * Ev],
                  [Fv - 0.5 * Gu, E0, F0],
                  [0.5 * Gv, F0, G0]])
    B = np.array([[0.0, 0.5 * Ev, 0.5 * Gu],
                  [0.5 * Ev,  E0, F0],
                  [0.5 * Gu,  F0, G0]])
    
    denom = (E0 * G0 - F0 * F0) ** 2
    if denom <= 0:
        return 0.0
    return float((np.linalg.det(A) - np.linalg.det(B)) / denom)


def norm(field, mask):
    field_masked = field[mask]
    if field_masked.size == 0 or field_masked.max() <= 0:
        return np.zeros_like(field)
    scale = np.percentile(field_masked, 95)
    if scale <= 0:
        scale = field_masked.max()
    return np.clip(field / scale, 0.0, 1.0)


def norm_obstacle(dist_field, mask, res, bounds, obstacles):
    if obstacles:
        x0, x1, y0, y1 = bounds
        sigma = decay_rate * np.hypot(x1 - x0, y1 - y0)
        field_norm = np.where(mask, np.exp(-dist_field / sigma), 0.0)
        return np.where(mask, np.exp(-dist_field / sigma), 0.0)
    else:
        return np.zeros((res, res))
    

def difficulty_field(manifold: Manifold, bounds, obstacles, res=200, weights=(1.0, 1.0, 1.0)):
    x0, x1, y0, y1 = bounds
    xs = np.linspace(x0, x1, res)
    ys = np.linspace(y0, y1, res)
    grad = np.zeros((res, res))
    curv = np.zeros((res, res))
    dist = np.full((res, res), np.inf)
    mask = np.zeros((res, res), dtype=bool)
    
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            p = np.array([x, y])
            if manifold.in_domain(p) and not blocked(p, obstacles):
                mask[i, j] = True
                grad[i, j] = metric_grad_norm(manifold, p)
                curv[i, j] = abs(gaussian_curvature(manifold, p))
                dist[i, j] = obstacle_distance(p, obstacles)
    
    grad = norm(grad, mask)
    curv = norm(curv, mask)
    dist = norm_obstacle(dist, mask, res, bounds, obstacles)
    diff_field = np.where(mask, weights[0] * grad + weights[1] * curv + weights[2] * dist, 0.0)
    if diff_field.max() != 0:
        diff_field = diff_field / diff_field.max()
    return xs, ys, diff_field


def cell_level(x, y, dx, dy, interp, tau, max_levels):
    if tau <= 0:
        return max_levels
    cx, cy = x + 0.5 * dx, y + 0.5 * dy
    control_points = [(x, y), (x + dx, y), (x, y + dy), (x + dx, y + dy), (cx, cy)]
    d = max(float(interp([[px, py]])[0]) for px, py in control_points)
    return min(max_levels, int(d / tau))


def sample_points_adaptively(manifold: Manifold, bounds, obstacles, n0=10, max_levels=2,
    tau=0.33, weights=(1.0, 1.0, 1.0), field=None, res=200):
    x0, x1, y0, y1 = bounds
    dx = (x1 - x0) / n0
    dy = (y1 - y0) / n0
    fine = 2 ** max_levels
    fine_dx = dx / fine
    fine_dy = dy / fine
    nodes = {}
    
    if field is None:
        field = difficulty_field(manifold, bounds, obstacles, res, weights)
    xs, ys, diff_field = field
    interp = RegularGridInterpolator((xs, ys), diff_field, bounds_error=False, fill_value=0.0)
    
    for a in range(n0):
        for b in range(n0):
            x, y = x0 + a * dx, y0 + b * dy
            m = 2 ** cell_level(x, y, dx, dy, interp, tau, max_levels)
            
            for u in range(m + 1):
                for v in range(m + 1):
                    px, py = x + u * dx / m, y + v * dy / m
                    index = (round((px - x0) / fine_dx), round((py - y0) / fine_dy))
                    nodes[index] = np.array([px, py])
                    
    return [p for p in nodes.values() if manifold.in_domain(p) and not blocked(p, obstacles)]


def samlpe_boundary_points(manifold: Manifold, obstacles, bounds, offset=None, spacing=None):
    x0, x1, y0, y1 = bounds
    diag = np.hypot(x1 - x0, y1 - y0)
    points = []
    if offset is None:
        offset = 0.02 * diag
    if spacing is None:
        spacing = 0.03 * diag
        
    for o in obstacles:
        if isinstance(o, Circle):
            r = o.r + offset
            m = max(8, int(2 * np.pi * r / spacing))
            for t in np.linspace(0, 2 * np.pi, m, endpoint=False):
                p = np.array([o.cx + r * np.cos(t), o.cy + r * np.sin(t)])
                if manifold.in_domain(p) and not blocked(p, obstacles):
                    points.append(p)
                    
        elif isinstance(o, Rectangle):
            sx, ex, sy, ey = o.sx - offset, o.ex + offset, o.sy - offset, o.ey + offset
            nx = max(2, int((ex - sx) / spacing))
            ny = max(2, int((ey - sy) / spacing))
            potential_points = []
            
            for x in np.linspace(sx, ex, nx):
                potential_points.append((x, sy))
                potential_points.append((x, ey))
            for y in np.linspace(sy, ey, ny):
                potential_points.append((sx, y))
                potential_points.append((ex, y))
            for x, y in potential_points:
                p = np.array([x, y])
                if manifold.in_domain(p) and not blocked(p, obstacles):
                    points.append(p)
            
    return points

def graph_from_points(manifold: Manifold, points, obstacles, k=8):
    points = np.asarray([np.asarray(p, dtype=float) for p in points])
    n = len(points)
    tree = cKDTree(points)
    dists, indexs = tree.query(points, k=min(k + 1, n))
    indexs = np.atleast_2d(indexs)
    adj = [dict() for _ in range(n)]
    
    for i in range(n):
        for j in indexs[i]:
            j = int(j)
            if i != j and j not in adj[i] and segment_free(manifold, points[i], points[j], obstacles):
                w = manifold.segment_lenght(points[i], points[j])
                adj[i][j] = w
                adj[j][i] = w
                
    return [points[i] for i in range(n)], adj

def build_adaptive_graph(manifold, bounds, obstacles, n0=10, max_levels=2, tau=0.33, k=8, weights=(1.0, 1.0, 1.0), res=100, field=None,
                         boundary_offset=None, boundary_spacing=None):
    if field is None:
        field = difficulty_field(manifold, bounds, obstacles, res, weights)
        
    V_int = sample_points_adaptively(manifold, bounds, obstacles, n0, max_levels, tau, weights, field, res)
    V_bound = samlpe_boundary_points(manifold, obstacles, bounds, boundary_offset, boundary_spacing)
    return graph_from_points(manifold, V_int + V_bound, obstacles, k)


def build_uniform_graph(manifold, bounds, nx, ny, obstacles, k=8):
    x0, x1, y0, y1 = bounds
    points = []
    for x in np.linspace(x0, x1, nx):
        for y in np.linspace(y0, y1, ny):
            p = np.array([x, y])
            if manifold.in_domain(p) and not blocked(p, obstacles):
                points.append(p)
    return graph_from_points(manifold, points, obstacles, k)
