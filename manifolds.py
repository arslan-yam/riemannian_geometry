import numpy as np

class Manifold:
    def metric(self, p):
        raise NotImplementedError
    
    def in_domain(self, p):
        return True
    
    def segment_lenght(self, p, q, n=10):
        delta = (q - p) / n
        length = 0.0
        for k in range(n):
            mid_point = p + (k + 0.5) * delta
            length += np.sqrt(delta.T @ self.metric(mid_point) @ delta)
        return length

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