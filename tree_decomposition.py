import heapq

def tree_decomposition(adj):
    n = len(adj)
    H = [dict(to) for to in adj]
    heap = [(len(H[v]), v) for v in range(n)]
    heapq.heapify(heap)
    eliminated = [False] * n
    pi = [0] * n #elimination rank
    bags = [None] * n #X(v) for v in V
    iter = 0
    middle = {}
    
    while heap:
        deg, v = heap.heappop(heap)
        if eliminated[v]:
            continue
        if deg != len(H[v]):
            heapq.heappush(heap, (len(H[v]), v))
            continue
        
        iter += 1
        pi[v] = iter
        neigbs = list(H[v].items())
        bags[v] = neigbs
        
        # connecting every pair of remaining nodes
        for i in range(len(neigbs)):
            u, w_u = neigbs[i]
            G_u = H[u]
            for j in range(i + 1, len(neigbs)):
                x, w_x = neigbs[j]
                w = w_u + w_x
                curr_w = G_u.get(x)
                if curr_w is None or w < curr_w:
                    G_u[x] = w
                    H[x][u] = w
                    middle[frozenset((u, x))] = v
        
        #remove v from its neigbors            
        for u, w_u in neigbs:
            del H[u][v]
            heapq.heappush(heap, (len(H[u]), u))
            
        eliminated[v] = True
        H[v] = {}
        
    parent = [None] * n
    root = None
    
    for v in range(n):
        if bags[v]:
            parent[v] = min(bags[v], key=lambda x: pi[x[0]])[0]
        else:
            root = v
            
    return pi, bags, parent, root, middle
