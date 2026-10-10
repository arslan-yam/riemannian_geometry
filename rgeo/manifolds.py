"""Riemannian metrics in a 2-d chart: the Manifold base (lengths, energy, gradients) and the
analytic test manifolds (Euclidean, Poincare disk, Terrain, GaussianHills)."""
import numpy as np
from scipy.optimize import minimize


gauss_cache = {}
def gauss_legendre(n):
    if n not in gauss_cache:
        xs, ws = np.polynomial.legendre.leggauss(n)
        gauss_cache[n] = (0.5 * (xs + 1.0), 0.5 * ws)
    return gauss_cache[n]


def build_metric(fx, fy):
    g = np.empty((len(fx), 2, 2))
    g[:, 0, 0] = 1.0 + fx * fx
    g[:, 0, 1] = g[:, 1, 0] = fx * fy
    g[:, 1, 1] = 1.0 + fy * fy
    return g


def build_metric_grad(fx, fy, fxx, fxy, fyy):
    #g = I + grad f grad f^T, so d g / d x_k = h_k grad f^T + grad f h_k^T with h_k the k-th
    #column of the hessian. one call replaces four metric_batch calls of central differences
    gx = np.empty((len(fx), 2, 2))
    gx[:, 0, 0] = 2.0 * fx * fxx
    gx[:, 0, 1] = gx[:, 1, 0] = fxx * fy + fx * fxy
    gx[:, 1, 1] = 2.0 * fy * fxy
    gy = np.empty((len(fx), 2, 2))
    gy[:, 0, 0] = 2.0 * fx * fxy
    gy[:, 0, 1] = gy[:, 1, 0] = fxy * fy + fx * fyy
    gy[:, 1, 1] = 2.0 * fy * fyy
    return gx, gy


class Manifold:
    #a global lower bound of the metric speed, |v|_g >= lambda_min |v|. it turns the padding of
    #a separator window into a proof (notes/error_bound.md); None means "unknown"
    lambda_min = None

    def metric(self, p):
        raise NotImplementedError

    def in_domain(self, p):
        return True

    def in_domain_batch(self, points):
        return np.ones(len(points), dtype=bool)

    def segment_lenght(self, p, q, n=4):
        delta = q - p
        ts, ws = gauss_legendre(n)
        length = 0.0
        for t, w in zip(ts, ws):
            point = p + t * delta
            length += w * np.sqrt(delta.T @ self.metric(point) @ delta)
        return length

    def metric_batch(self, points): #override per manifold, this fallback is the slow path
        return np.array([self.metric(p) for p in points])

    def energy(self, points, m):
        delta = points[1:] - points[:-1]
        mid_point = 0.5 * (points[1:] + points[:-1])
        return m * float(np.einsum("ij,ijk,ik->", delta, self.metric_batch(mid_point), delta))

    def metric_grad(self, p, h=1e-4):
        dx, dy = np.array([h, 0.0]), np.array([0.0, h])
        gx = (self.metric(p + dx) - self.metric(p - dx)) / (2 * h)
        gy = (self.metric(p + dy) - self.metric(p - dy)) / (2 * h)
        return gx, gy

    def metric_grad_batch(self, points, h=1e-4):
        dx, dy = np.array([h, 0.0]), np.array([0.0, h])
        gx = (self.metric_batch(points + dx) - self.metric_batch(points - dx)) / (2 * h)
        gy = (self.metric_batch(points + dy) - self.metric_batch(points - dy)) / (2 * h)
        return gx, gy

    def energy_grad(self, points, m):
        delta = points[1:] - points[:-1]
        mid_point = 0.5 * (points[1:] + points[:-1])
        gd = 2.0 * np.einsum("ijk,ik->ij", self.metric_batch(mid_point), delta)
        gx, gy = self.metric_grad_batch(mid_point)
        gmid = 0.5 * np.stack([np.einsum("ij,ijk,ik->i", delta, gx, delta),
                               np.einsum("ij,ijk,ik->i", delta, gy, delta)], axis=1)
        grad = np.zeros((len(points), 2))
        grad[1:] += gd + gmid
        grad[:-1] += -gd + gmid
        return m * grad[1:-1].ravel()

    def energy_and_grad_points(self, points, m):
        #energy and its gradient from ONE metric evaluation at the midpoints. calling energy()
        #and energy_grad() separately evaluated the metric twice on the same points
        delta = points[1:] - points[:-1]
        mid_point = 0.5 * (points[1:] + points[:-1])
        g = self.metric_batch(mid_point)
        gd = np.einsum("ijk,ik->ij", g, delta)
        e = m * float(np.einsum("ij,ij->", delta, gd))
        gx, gy = self.metric_grad_batch(mid_point)
        gmid = 0.5 * np.stack([np.einsum("ij,ijk,ik->i", delta, gx, delta),
                               np.einsum("ij,ijk,ik->i", delta, gy, delta)], axis=1)
        grad = np.zeros((len(points), 2))
        grad[1:] += 2.0 * gd + gmid
        grad[:-1] += -2.0 * gd + gmid
        return e, m * grad[1:-1].ravel()

    def energy_and_grad(self, p, q, z, m):
        points = np.vstack([p, z.reshape(-1, 2), q])
        return self.energy_and_grad_points(points, m)

    def geodesic_path(self, p, q, m=8): #using minimize energy
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float)
        if m < 2:
            return np.vstack([p, q])
        ts = np.linspace(0.0, 1.0, m + 1)[1:-1]
        start = np.array([p + t * (q - p) for t in ts]).ravel()
        res = minimize(lambda z: self.energy_and_grad(p, q, z, m), start, method="L-BFGS-B", jac=True)
        return np.vstack([p, res.x.reshape(-1, 2), q])

    def geodesic_lenght(self, p, q, m=8):
        points = self.geodesic_path(p, q, m)
        return sum(self.segment_lenght(points[i], points[i + 1]) for i in range(len(points) - 1))

class Euclidean(Manifold):
    lambda_min = 1.0

    def metric(self, p):
        return np.eye(2)

    def metric_batch(self, points):
        return np.broadcast_to(np.eye(2), (len(points), 2, 2))
    
    def segment_lenght(self, p, q, n=10):
        return float(np.sqrt((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2))
    
class PoincareDisk(Manifold):
    lambda_min = 2.0 #lambda = 2/(1-|p|^2) >= 2, equality at the origin

    def __init__(self, rmax=0.95):
        self.rmax = rmax
        
    def metric(self, p):
        lmbd = 2.0 / (1.0 - float(p[0] ** 2 + p[1] ** 2))
        return (lmbd ** 2) * np.eye(2)
    
    def metric_batch(self, points):
        lmbd = 2.0 / (1.0 - (points ** 2).sum(1))
        return (lmbd ** 2)[:, None, None] * np.eye(2)

    def metric_grad_batch(self, points, h=None):
        #d(lambda^2)/dx = 2 lambda dlambda/dx and dlambda/dx = lambda^2 x, so 2 lambda^3 x
        lmbd = 2.0 / (1.0 - (points ** 2).sum(1))
        c = 2.0 * lmbd ** 3
        eye = np.eye(2)
        return (c * points[:, 0])[:, None, None] * eye, (c * points[:, 1])[:, None, None] * eye

    def in_domain(self, p):
        return p[0] ** 2 + p[1] ** 2 < self.rmax ** 2

    def in_domain_batch(self, points):
        return (points ** 2).sum(1) < self.rmax ** 2
    
    @staticmethod
    def exact_dist(u, v):
        u = complex(u[0], u[1])
        v = complex(v[0], v[1])
        num = 2.0 * abs(u - v) ** 2
        den = (1.0 - abs(u) ** 2) * (1.0 - abs(v) ** 2)
        return float(np.arccosh(1.0 + num / den))
    
    @staticmethod
    def exact_geodesic(u, v, n=100):
        u = complex(u[0], u[1])
        v = complex(v[0], v[1])
        T = lambda z: (z - u) / (1.0 - np.conj(u) * z)
        Ti = lambda w: (w + u) / (1.0 + np.conj(u) * w)
        w1 = T(v)
        ts = np.linspace(0.0, 1.0, n)
        zs = np.array([Ti(t * w1) for t in ts])
        return np.c_[zs.real, zs.imag]
    

class Terrain(Manifold):
    lambda_min = 1.0 #I + grad f grad f^T has smallest eigenvalue 1

    def __init__(self, amp=0.9, freq=1.4):
        self.amp = amp
        self.freq = freq
        
    def f(self, x, y):
        return self.amp * np.sin(self.freq * x) * np.cos(self.freq * y)
    
    def grad(self, p):
        x, y = p[0], p[1]
        fx = self.amp * self.freq * np.cos(self.freq * x) * np.cos(self.freq * y)
        fy = -self.amp * self.freq * np.sin(self.freq * x) * np.sin(self.freq * y)
        return fx, fy
    
    def metric(self, p):
        fx, fy = self.grad(p)
        return np.array([[1.0 + fx * fx, fx * fy],
                         [fx * fy, 1.0 + fy * fy]])

    def metric_batch(self, points):
        x, y = points[:, 0], points[:, 1]
        fx = self.amp * self.freq * np.cos(self.freq * x) * np.cos(self.freq * y)
        fy = -self.amp * self.freq * np.sin(self.freq * x) * np.sin(self.freq * y)
        return build_metric(fx, fy)

    def metric_grad_batch(self, points, h=None):
        x, y = points[:, 0], points[:, 1]
        a, w = self.amp, self.freq
        sx, cx, sy, cy = np.sin(w * x), np.cos(w * x), np.sin(w * y), np.cos(w * y)
        fx, fy = a * w * cx * cy, -a * w * sx * sy
        fxx, fxy, fyy = -a * w * w * sx * cy, -a * w * w * cx * sy, -a * w * w * sx * cy
        return build_metric_grad(fx, fy, fxx, fxy, fyy)


class GaussianHills(Manifold):
    lambda_min = 1.0

    def __init__(self, hills):
        self.hills = hills

    def f(self, x, y):
        z = 0.0
        for cx, cy, a, s in self.hills:
            z += a * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2.0 * s * s))
        return z

    def grad(self, p):
        x, y = p[0], p[1]
        fx = fy = 0.0
        for cx, cy, a, s in self.hills:
            e = a * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2.0 * s * s))
            fx += -e * (x - cx) / (s * s)
            fy += -e * (y - cy) / (s * s)
        return fx, fy

    def metric(self, p):
        fx, fy = self.grad(p)
        return np.array([[1.0 + fx * fx, fx * fy],
                         [fx * fy, 1.0 + fy * fy]])

    def metric_batch(self, points):
        x, y = points[:, 0], points[:, 1]
        fx = np.zeros(len(points))
        fy = np.zeros(len(points))
        for cx, cy, a, sg in self.hills:
            e = a * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2.0 * sg * sg))
            fx += -e * (x - cx) / (sg * sg)
            fy += -e * (y - cy) / (sg * sg)
        return build_metric(fx, fy)

    def metric_grad_batch(self, points, h=None):
        x, y = points[:, 0], points[:, 1]
        n = len(points)
        fx, fy, fxx, fxy, fyy = (np.zeros(n) for _ in range(5))
        for cx, cy, a, sg in self.hills:
            s2 = sg * sg
            u, v = x - cx, y - cy
            e = a * np.exp(-(u * u + v * v) / (2.0 * s2))
            fx += -e * u / s2
            fy += -e * v / s2
            fxx += e * (u * u / s2 - 1.0) / s2
            fyy += e * (v * v / s2 - 1.0) / s2
            fxy += e * u * v / (s2 * s2)
        return build_metric_grad(fx, fy, fxx, fxy, fyy)
