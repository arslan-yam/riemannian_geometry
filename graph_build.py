import numpy as np
from collections import deque
from manifolds import Manifold
from typing import List
from obstacles import Obstacle

ds = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (2, -1), (1, -2)]
edge_sampling_params = [0.2, 0.4, 0.6, 0.8]

def blocked(p, obstacles: List[Obstacle]):
    for obstacle in obstacles:
        if obstacle.blocked(p):
            return True
    return False


def build_graph(manifold: Manifold, bounds, nx, ny, obstacles: List[Obstacle]):
    x0, x1, y0, y1 = bounds
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)
    idxs = -np.ones((nx, ny), dtype=int) #store the index of the point
    points = [] #storing all available points (not blocked by obstacles)
    
    for x in range(nx):
        for y in range(ny):
            p = np.array([xs[x], ys[y]])
            if manifold.in_domain(p) and not blocked(p, obstacles):
                idxs[x, y] = len(points)
                points.append(p)
                
    adj = [dict() for _ in range(len(points))] #for each point we store neiboring points connected by edge and edge weight
    
    for x in range(nx):
        for y in range(ny):
            idx = idxs[x, y]
            if idx == -1:
                continue
            
            for (dx, dy) in ds:
                if (0 <= x + dx < nx and 0 <= y + dy < ny):
                    neig_idx = idxs[x + dx, y + dy]
                    if neig_idx == -1:
                        continue
                    
                    can_connect = True
                    for param in edge_sampling_params:
                        q = (1 - param) * points[idx] + param * points[neig_idx]
                        if not manifold.in_domain(q) or blocked(q, obstacles):
                            can_connect = False
                            break
                    if not can_connect:
                        continue
                    
                    w = manifold.segment_lenght(points[idx], points[neig_idx])
                    adj[idx][neig_idx] = w
                    adj[neig_idx][idx] = w
                    
    return points, adj


def connected_comps(adj):
    n = len(adj)
    labels = np.full(n, -1)
    members = []
    for s in range(n):
        if labels[s] != -1:
            continue
        labels[s] = len(members)
        stack = [s]
        group = []
        while stack:
            v = stack.pop()
            group.append(v)
            for u in adj[v]:
                if labels[u] == -1:
                    labels[u] = len(members)
                    stack.append(u)
        members.append(group)
    return labels, members
    