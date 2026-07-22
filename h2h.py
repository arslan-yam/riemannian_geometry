import numpy as np
from tree_decomposition import tree_decomposition
from graph_build import connected_comps

class H2HIndex:
    def __init__(self, adj):
        pi, bags, parent, root, middle = tree_decomposition(adj)
        self.adj = self.adj
        self.n = len(adj)
        self.middle = middle
        self.root = root
        self.parent = parent
        
        children = [[] for _ in range(self.n)]
        for v in range(self.n):
            if parent[v] is not None:
                children[parent[v]].append(v)
        self.children = children
        
        anc = [None] * self.n
        depth = [0] * self.n
        anc[root] = np.array([root])
        order = [root]
        stack = [root]
        
        while stack:
            v = stack.pop()
            for c in children[v]:
                anc[c] = np.append(anc[v], c)
                depth[c] = depth[v] + 1
                order.append(c)
                stack.append(c)
        
        self.anc = anc
        self.depth = depth
        self.h = int(depth.max()) + 1
        self.width = max(len(b) for b in bags)
        self.Xv = [np.array([x for x, _ in bags[v]], dtype=np.int32) for v in range(self.n)]
        self.Xw = [np.array([w for _, w in bags[v]], dtype=float) for v in range(self.n)]
        
        dis = np.full((self.n, self.h), np.inf)
        mid = np.full((self.n, self.h), -1, dtype=np.int32)
        dis[root, 0] = 0.0
        
        for v in order:
            if v != root:
                l = len(anc[v])
                xv, xw = self.Xv[v], self.Xw[v]
                cand = np.empty((len(xv), l - 1))
                
                for j in range(len(xv)):
                    x = int(xv[j])
                    p = int(depth[x])
                    row = cand[j]
                    row[:p] = dis[x, :p]
                    row[p:] = dis[anc[v][p:l - 1], p]
                    row += xw[j]
                am = cand.argmin(axis=0)
                dis[v, :l - 1] = cand[am, np.arange(l - 1)]
                dis[v, l - 1] = 0.0
                mid[v, :l - 1] = xv[am]
                
        self.dis = dis
        self.mid = mid
        self.pos = [np.sort(np.concatenate([depth[self.Xv[v]], [depth[v]]])).astype(np.int64) for v in range(self.n)]
        self.build_lca()
        self.label_entries = int(depth.sum()) + self.n
        
    def build_lca(self):
        euler = [self.root]
        children =  self.children
        first = np.full(self.n, -1, dtype=np.int64)
        first[self.root] = 0
        stack = [(self.root, iter(children[self.root]))]
        
        while stack:
            v, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                stack.pop()
                if stack:
                    euler.append(stack[-1][0])
            else:
                if first[nxt] < 0:
                    first[nxt] = len(euler)
                euler.append(nxt)
                stack.append((nxt, iter(children[nxt])))
                
        euler = np.array(euler, dtype=np.int64)
        dep = self.depth[euler].astype(np.int64)
        m = len(euler)
        j = 1
        sp = np.zeros((m.bit_length(), m), dtype=np.int64)
        sp[0] = np.arange(m)
        
        while (1 << j) <= m:
            span, half = 1 << j, 1 << (j - 1)
            l = sp[j - 1, :m - span + 1]
            r = sp[j - 1, half:m - span + 1 + half]
            sp[j, :m - span + 1] = np.where(dep[l] <= dep[r], l, r)
            j += 1
            
        self.euler = euler
        self.first = first
        self.dep = dep
        self.sp = sp
        
    def lca(self, u, v):
        l, r = self.first[u], self.first[v]
        if l > r:
            l, r = r, l
        j = int(r - l + 1).bit_length() - 1
        a = self.sp[j, l]
        b = self.sp[j, r - (1 << j) + 1]
        i = a
        if self.dep[a] > self.dep[b]:
            i = b
        return int(self.euler[i])
    
    def query(self, s, t):
        if s == t:
            return 0.0
        lca = self.lca(s, t)
        p = self.pos[lca]
        dists = self.dis[s, p] + self.dis[t, p]
        return float(dists.min())
    
    def unpack(self, u, v):
        m = self.middle.get(frozenset((int(u), int(v))))
        if m is None:
            return [int(u), int(v)]
        return self.unpack(u, m)[:-1] + self.unpack(m, v)
        
    def path_up(self, v, i):
        a = int(self.anc[v][i])
        if a == v:
            return [v]
        x = int(self.mid[v, i])
        px = int(self.depth[x])
        head = self.unpack(v, x)
        if px > i:
            return head[:-1] + self.path_up(x, i)
        else:
            tail = self.path_up(a, px)
            return head[:-1] + tail[::-1] 
        
    def query_path(self, s, t):
        if s == t:
            return 0.0
        lca = self.lca(s, t)
        p = self.pos[lca]
        dists = self.dis[s, p] + self.dis[t, p]
        i = int(dists.argmin())
        c = int(self.pos[lca][i])
        l = self.path_up(s, c)
        r = self.path_up(t, c)
        return float(dists[i]), l + r[::-1][1:]
        
        
class MultiH2HIndex:
    def __init__(self, adj):
        self.adj = adj
        self.n = len(adj)
        self.comp, self.members = connected_comps(adj)
        self.local = np.full(self.n, -1)
        self.indices = []

        for group in self.members:
            g2l = {}
            for lid, g in enumerate(group):
                g2l[g] = lid
                self.local[g] = lid
            sub_adj = [dict() for _ in group]
            for g in group:
                lid = g2l[g]
                for u, w in adj[g].items():
                    sub_adj[lid][g2l[u]] = w
            self.indices.append(H2HIndex(sub_adj))

    def same_component(self, s, t):
        return self.comp[s] == self.comp[t]

    def query(self, s, t):
        if self.comp[s] != self.comp[t]:
            return float("inf")
        c = self.comp[s]
        return self.indices[c].query(self.local[s], self.local[t])

    def query_path(self, s, t):
        if self.comp[s] != self.comp[t]:
            return float("inf"), None
        c = self.comp[s]
        d, local_path = self.indices[c].query_path(self.local[s], self.local[t])
        group = self.members[c]
        return d, [group[i] for i in local_path]
