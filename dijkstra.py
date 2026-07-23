import heapq
import numpy as np

def dijkstra(adj, s, t=None):
    n = len(adj)
    dist = np.full(n, np.inf)
    pred = np.full(n, -1)
    dist[s] = 0.0
    q = [(0.0, s)]

    while q:
        curr_dist, v = heapq.heappop(q)
        if curr_dist > dist[v]:
            continue
        if t is not None and v == t:
            break
        
        for to, weigth in adj[v].items():
            if curr_dist + weigth < dist[to]:
                dist[to] = curr_dist + weigth
                pred[to] = v
                heapq.heappush(q, (curr_dist + weigth, to))
                
    return dist, pred


def get_path(pred, s, t):
    path = [t] 
    while path[-1] != s:
        p = int(pred[path[-1]])
        if p < 0:
            return None
        path.append(p)
    return path[::-1]