import heapq
import numpy as np
from graph_build import blocked

ds8 = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]

def quad(px, py, qx, qy, e, f, g):
    return e * px * qx + f * (px * qy + py * qx) + g * py * qy


def two_point(u1, u2, ax, ay, bx, by, e, f, g):
    cx, cy = ax - bx, ay - by
    a = quad(cx, cy, cx, cy, e, f, g)
    b = quad(cx, cy, bx, by, e, f, g)
    c = quad(bx, by, bx, by, e, f, g)
    d = u1 - u2
    k = a - d * d
    best = min(u2 + np.sqrt(c), u1 + np.sqrt(a + 2 * b + c))
    
    if a > 0 and abs(k) > 1e-14:
        disc = b * b * k * k - a * k * (b * b - d * d * c)
        if disc >= 0:
            sq = np.sqrt(disc)
            for s in ((-b * k + sq) / (a * k), (-b * k - sq) / (a * k)):
                if 0.0 < s < 1.0:
                    q = c + 2 * s * b + s * s * a
                    if q > 0:
                        val = u2 + s * d + np.sqrt(q)
                        if val < best:
                            best = val
    return best


def passable(available, nx, ny, a, b, dx, dy):
    if dx != 0 and dy != 0:
        if not (0 <= a + dx < nx) or not available[a + dx, b]:
            return False
        if not (0 <= b + dy < ny) or not available[a, b + dy]:
            return False
    return True


def solve_eikonal(manifold, bounds, nx, ny, sources, obstacles=()):
    x0, x1, y0, y1 = bounds
    xs, ys = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)
    hx, hy = (x1 - x0) / (nx - 1), (y1 - y0) / (ny - 1)
    e, f, g  = np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny))
    available = np.zeros((nx, ny), dtype=bool)

    for i in range(nx):
        for j in range(ny):
            p = np.array([xs[i], ys[j]])
            if manifold.in_domain(p) and not blocked(p, obstacles):
                g = manifold.metric(p)
                e[i, j], f[i, j], g[i, j] = g[0, 0], g[0, 1], g[1, 1]
                available[i, j] = True

    u = np.full((nx, ny), np.inf)
    accepted = np.zeros((nx, ny), dtype=bool)
    heap = []

    for s in sources:
        i = min(max(int(round((s[0] - x0) / hx)), 0), nx - 1)
        j = min(max(int(round((s[1] - y0) / hy)), 0), ny - 1)
        if available[i, j]:
            u[i, j] = 0.0
            heapq.heappush(heap, (0.0, i, j))

    while heap:
        val, i, j = heapq.heappop(heap)
        if accepted[i, j] or val > u[i, j]:
            continue
        accepted[i, j] = True

        for (dx, dy) in ds8:
            a, b = i + dx, j + dy
            if not (0 <= a < nx and 0 <= b < ny) or not available[a, b] or accepted[a, b]:
                continue
            e, f, g = e[a, b], f[a, b], g[a, b]
            best = np.inf
            
            for m in range(8):
                dx1, dy1 = ds8[m]
                dx2, dy2 = ds8[(m + 1) % 8]
                i1, j1 = a + dx1, b + dy1
                i2, j2 = a + dx2, b + dy2
                b1 = 0 <= i1 < nx and 0 <= j1 < ny and accepted[i1, j1] and passable(available, nx, ny, a, b, dx1, dy1)
                b2 = 0 <= i2 < nx and 0 <= j2 < ny and accepted[i2, j2] and passable(available, nx, ny, a, b, dx2, dy2)
                if not b1 and not b2:
                    continue
                ax, ay = dx1 * hx, dy1 * hy
                bx, by = dx2 * hx, dy2 * hy
                if b1 and b2:
                    cand = two_point(u[i1, j1], u[i2, j2], ax, ay, bx, by, e, f, g)
                elif b1:
                    cand = u[i1, j1] + np.sqrt(quad(ax, ay, ax, ay, e, f, g))
                else:
                    cand = u[i2, j2] + np.sqrt(quad(bx, by, bx, by, e, f, g))
                best = min(best, cand)

            if best < u[a, b]:
                u[a, b] = best
                heapq.heappush(heap, (best, a, b))

    u[~available] = np.inf
    return xs, ys, u


def sample_field(xs, ys, u, p):
    x, y = float(p[0]), float(p[1])
    i = np.clip(np.searchsorted(xs, x) - 1, 0, len(xs) - 2)
    j = np.clip(np.searchsorted(ys, y) - 1, 0, len(ys) - 2)
    tx = (x - xs[i]) / (xs[i + 1] - xs[i])
    ty = (y - ys[j]) / (ys[j + 1] - ys[j])
    c = u[i:i + 2, j:j + 2]
    if not np.isfinite(c).all():
        return float(c[int(round(tx)), int(round(ty))])
    return float((1 - tx) * (1 - ty) * c[0, 0] + tx * (1 - ty) * c[1, 0] + (1 - tx) * ty * c[0, 1] + tx * ty * c[1, 1])


def sample_points(xs, ys, u, points):
    points = np.asarray([np.asarray(p, dtype=float) for p in points])
    i = np.clip(np.searchsorted(xs, points[:, 0]) - 1, 0, len(xs) - 2)
    j = np.clip(np.searchsorted(ys, points[:, 1]) - 1, 0, len(ys) - 2)
    tx = (points[:, 0] - xs[i]) / (xs[i + 1] - xs[i])
    ty = (points[:, 1] - ys[j]) / (ys[j + 1] - ys[j])
    c = np.stack([u[i, j], u[i + 1, j], u[i, j + 1], u[i + 1, j + 1]])
    w = np.stack([(1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty])
    
    good = np.isfinite(c).all(0)
    out = np.where(good, (w * np.where(np.isfinite(c), c, 0.0)).sum(0), np.inf)
    bad = ~good
    if bad.any():
        out[bad] = np.where(np.isfinite(c[:, bad]), c[:, bad], np.inf).min(0)
    return out


def anisotropy_ratio(manifold, bounds, res=60):
    x0, x1, y0, y1 = bounds
    worst = 1.0
    for x in np.linspace(x0, x1, res):
        for y in np.linspace(y0, y1, res):
            p = np.array([x, y])
            if manifold.in_domain(p):
                w = np.linalg.eigvalsh(manifold.metric(p))
                if w[0] > 0:
                    worst = max(worst, float(np.sqrt(w[1] / w[0])))
    return worst
