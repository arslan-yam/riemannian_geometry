import numpy as np
from manifolds import Manifold
from typing import List
from obstacles import Obstacle

ds = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (2, -1), (1, -2)]
ts = [0.2, 0.4, 0.6, 0.8]

def blocked(p, obstacles: List[Obstacle]):
    for obstacle in obstacles:
        if obstacle.blocked(p):
            return True
    return False

def segment_free(manifold: Manifold, x, y, obstacles):
    for t in ts:
        q = (1.0 - t) * x + t * y
        if not manifold.in_domain(q) or blocked(q, obstacles):
            return False
    return True 


def sample_points_grid(manifold: Manifold, nx, ny, bounds, obstacles):
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
                
    return points, idxs
    
def sample_points_tiangular(manifold: Manifold, nx, ny, bounds, obstacles):
    x0, x1, y0, y1 = bounds
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    idxs = -np.ones((nx, ny), dtype=int) #store the index of the point
    points = [] #storing all available points (not blocked by obstacles)

    for j in range(ny):
        for i in range(nx):
            offset = 0.5 * dx if j % 2 else 0.0 #shift odd rows to make triangles
            p = np.array([x0 + i * dx + offset, y0 + j * dy])
            if manifold.in_domain(p) and not blocked(p, obstacles):
                idxs[i, j] = len(points)
                points.append(p)
                
    return points, idxs

def build_graph(manifold: Manifold, bounds, nx, ny, obstacles: List[Obstacle] = [], type = "rectangle"):
    points, idx = None, None
    if type == "rectangle":
        points, idxs = sample_points_grid(manifold, nx, ny, bounds, obstacles)                
    elif type == "triangular":
        points, idxs = sample_points_tiangular(manifold, nx, ny, bounds, obstacles)
    adj = [dict() for _ in range(len(points))] #for each point we store neiboring points connected by edge and edge weight

    if type == "rectangle":
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
                        if segment_free(manifold, points[idx], points[neig_idx], obstacles):
                            w = manifold.segment_lenght(points[idx], points[neig_idx])
                            adj[idx][neig_idx] = w
                            adj[neig_idx][idx] = w
                            
    elif type == "triangular":
        for x in range(nx):
            for y in range(ny):
                idx = idxs[x, y]
                if idx == -1:
                    continue
                if y % 2:
                    neigbs = [(x + 1, y), (x, y + 1), (x + 1, y + 1)]
                else:
                    neigbs = [(x + 1, y), (x - 1, y + 1), (x, y + 1)]

                for (dx, dy) in neigbs:
                    if 0 <= dx < nx and 0 <= dy < ny:
                        neig_idx = idxs[dx, dy]
                        if neig_idx == -1:
                            continue
                        if segment_free(manifold, points[idx], points[neig_idx], obstacles):
                            w = manifold.segment_lenght(points[idx], points[neig_idx])
                            adj[idx][neig_idx] = w
                            adj[neig_idx][idx] = w
                            
    return points, adj


def connected_comps(adj):
    n = len(adj)
    labels = np.full(n, -1, dtype=int)
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
