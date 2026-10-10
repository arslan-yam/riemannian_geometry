"""Obstacles in the chart: point tests, batched tests, signed clearance for the refinement barrier."""
import numpy as np
from abc import ABC, abstractmethod


def blocked(p, obstacles):
    #True when the chart point p lies inside any obstacle
    for obstacle in obstacles:
        if obstacle.blocked(p):
            return True
    return False


class Obstacle(ABC):
    @abstractmethod
    def blocked(self, p):
        pass

    def blocked_batch(self, points): #override for speed, this fallback is the slow path
        return np.array([self.blocked(q) for q in points], dtype=bool)

    def clearance_batch(self, points):
        #signed distance to the boundary (negative inside) and its gradient, for all points at
        #once. the barrier in refine.py needs it, so every obstacle type has to provide one
        raise NotImplementedError(f"{type(self).__name__} needs clearance_batch for path refinement")
    
    
class Circle(Obstacle):
    def __init__(self, cx, cy, r):
        self.cx, self.cy, self.r = cx, cy, r
        
    def blocked(self, p):
        if (p[0] - self.cx) ** 2 + (p[1] - self.cy) ** 2 < self.r * self.r:
            return True
        return False

    def blocked_batch(self, points):
        return (points[:, 0] - self.cx) ** 2 + (points[:, 1] - self.cy) ** 2 < self.r * self.r

    def clearance_batch(self, points):
        v = points - np.array([self.cx, self.cy])
        d = np.hypot(v[:, 0], v[:, 1])
        g = np.where(d[:, None] > 1e-12, v / np.maximum(d, 1e-12)[:, None], np.array([1.0, 0.0]))
        return d - self.r, g

    
class Rectangle(Obstacle):
    def __init__(self, sx, ex, sy, ey):
        self.sx, self.ex, self.sy, self.ey = sx, ex, sy, ey

    def blocked(self, p):
        return self.sx <= p[0] <= self.ex and self.sy <= p[1] <= self.ey

    def blocked_batch(self, points):
        x, y = points[:, 0], points[:, 1]
        return (self.sx <= x) & (x <= self.ex) & (self.sy <= y) & (y <= self.ey)

    def clearance_batch(self, points):
        x, y = points[:, 0], points[:, 1]
        dx = np.maximum(np.maximum(self.sx - x, 0.0), x - self.ex)
        dy = np.maximum(np.maximum(self.sy - y, 0.0), y - self.ey)
        outside = (dx > 0.0) | (dy > 0.0)
        c = np.hypot(dx, dy)
        sx = np.where(x < self.sx, -1.0, np.where(x > self.ex, 1.0, 0.0))
        sy = np.where(y < self.sy, -1.0, np.where(y > self.ey, 1.0, 0.0))
        g = np.stack([sx * dx, sy * dy], 1) / np.maximum(c, 1e-12)[:, None]
        #inside: distance to the nearest side, gradient points out through that side
        sides = np.stack([x - self.sx, self.ex - x, y - self.sy, self.ey - y], 1)
        k = sides.argmin(1)
        cin = -sides[np.arange(len(x)), k]
        gin = np.array([[-1.0, 0.0], [1.0, 0.0], [0.0, -1.0], [0.0, 1.0]])[k]
        return np.where(outside, c, cin), np.where(outside[:, None], g, gin)
