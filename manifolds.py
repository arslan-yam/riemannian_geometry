import numpy as np
from scipy.optimize import minimize


gauss_cache = {}
def gauss_legendre(n):
    if n not in gauss_cache:
        xs, ws = np.polynomial.legendre.leggauss(n)
        gauss_cache[n] = (0.5 * (xs + 1.0), 0.5 * ws)
    return gauss_cache[n]


class Manifold:
    def metric(self, p):
        raise NotImplementedError

    def in_domain(self, p):
        return True

    def segment_lenght(self, p, q, n=4):
        delta = q - p
        ts, ws = gauss_legendre(n)
        length = 0.0
        for t, w in zip(ts, ws):
            point = p + t * delta
            length += w * np.sqrt(delta.T @ self.metric(point) @ delta)
        return length

    def energy(self, points, m): 
        e = 0.0
        for i in range(m):
            delta = points[i + 1] - points[i]
            mid_point = 0.5 * (points[i] + points[i + 1])
            e += delta.T @ self.metric(mid_point) @ delta
        return m * e

    def metric_grad(self, p, h=1e-4):
        dx, dy = np.array([h, 0.0]), np.array([0.0, h])
        gx = (self.metric(p + dx) - self.metric(p - dx)) / (2 * h)
        gy = (self.metric(p + dy) - self.metric(p - dy)) / (2 * h)
        return gx, gy

    def energy_grad(self, points, m):
        grad = np.zeros((len(points), 2))
        for i in range(m):
            delta = points[i + 1] - points[i]
            mid_point = 0.5 * (points[i] + points[i + 1])
            gd = 2.0 * (self.metric(mid_point) @ delta)
            gx, gy = self.metric_grad(mid_point)
            gmid = 0.5 * np.array([delta @ gx @ delta, delta @ gy @ delta])
            grad[i + 1] += gd + gmid
            grad[i] += -gd + gmid
        return m * grad[1:-1].ravel()

    def energy_and_grad(self, p, q, z, m):
        points = np.vstack([p, z.reshape(-1, 2), q])
        return self.energy(points, m), self.energy_grad(points, m)

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
    def metric(self, p):
        return np.eye(2)
    
    def segment_lenght(self, p, q, n=10):
        return float(np.sqrt((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2))
    
class PoincareDisk(Manifold):
    def __init__(self, rmax=0.95):
        self.rmax = rmax
        
    def metric(self, p):
        lmbd = 2.0 / (1.0 - float(p[0] ** 2 + p[1] ** 2))
        return (lmbd ** 2) * np.eye(2)
    
    def in_domain(self, p):
        return p[0] ** 2 + p[1] ** 2 < self.rmax ** 2
    
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


class GaussianHills(Manifold):
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