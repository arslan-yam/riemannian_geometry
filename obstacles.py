import numpy as np
from abc import ABC, abstractmethod

class Obstacle(ABC):
    @abstractmethod
    def blocked(self, p):
        pass
    
    
class Circle(Obstacle):
    def __init__(self, cx, cy, r):
        self.cx, self.cy, self.r = cx, cy, r
        
    def blocked(self, p):
        if (p[0] - self.cx) ** 2 + (p[1] - self.cy) ** 2 < self.r * self.r:
            return True
        return False

    
class Rectangle(Obstacle):
    def __init__(self, sx, ex, sy, ey):
        self.sx, self.ex, self.sy, self.ey = sx, ex, sy, ey

    def blocked(self, p):
        return self.sx <= p[0] <= self.ex and self.sy <= p[1] <= self.ey
