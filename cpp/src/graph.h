// The mesh graph (edges charged at 3-d length x the cheaper adjacent face weight, so graph paths are
// surface curves and graph distances are upper bounds) and a reusable single-source Dijkstra.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <unordered_map>
#include <vector>
#include "terrain.h"

namespace cci {

struct Graph {
    int n = 0;
    std::vector<int64_t> ptr;
    std::vector<int> adj;
    std::vector<double> wt;
    void from_mesh(const Terrain& T) {
        n = T.n;
        std::unordered_map<uint64_t, double> wmin; // undirected edge -> cheaper adjacent weight
        wmin.reserve(3 * T.m);
        for (int f = 0; f < T.m; ++f)
            for (int k = 0; k < 3; ++k) {
                uint64_t u = T.F[3 * f + k], v = T.F[3 * f + (k + 1) % 3];
                uint64_t key = std::min(u, v) * (uint64_t)T.n + std::max(u, v);
                auto it = wmin.find(key);
                if (it == wmin.end()) wmin[key] = T.w[f]; else it->second = std::min(it->second, T.w[f]);
            }
        std::vector<int64_t> deg(n + 1, 0);
        for (auto& kv : wmin) { deg[kv.first / T.n + 1]++; deg[kv.first % T.n + 1]++; }
        ptr.assign(n + 1, 0);
        for (int i = 0; i < n; ++i) ptr[i + 1] = ptr[i] + deg[i + 1];
        adj.resize(ptr[n]); wt.resize(ptr[n]);
        std::vector<int64_t> fill(ptr.begin(), ptr.end() - 1);
        for (auto& kv : wmin) {
            int u = (int)(kv.first / T.n), v = (int)(kv.first % T.n);
            double L = std::sqrt(std::pow(T.V[3 * u] - T.V[3 * v], 2) + std::pow(T.V[3 * u + 1] - T.V[3 * v + 1], 2) +
                                 std::pow(T.V[3 * u + 2] - T.V[3 * v + 2], 2)) * kv.second;
            adj[fill[u]] = v; wt[fill[u]++] = L;
            adj[fill[v]] = u; wt[fill[v]++] = L;
        }
    }
};

// Dijkstra with a distance limit (nodes beyond it stay INF); dist is valid until the next run
struct Dijkstra {
    const Graph* G;
    std::vector<double> dist;
    std::vector<int> touched;
    explicit Dijkstra(const Graph& g) : G(&g), dist(g.n, INF) {}
    void run(int src, double limit = INF) {
        for (int v : touched) dist[v] = INF;
        touched.clear();
        using item = std::pair<double, int>;
        std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
        dist[src] = 0; touched.push_back(src); pq.push({0.0, src});
        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d > dist[u]) continue;
            for (int64_t e = G->ptr[u]; e < G->ptr[u + 1]; ++e) {
                int v = G->adj[e];
                double nd = d + G->wt[e];
                if (nd > limit) continue;
                if (nd < dist[v]) {
                    if (dist[v] == INF) touched.push_back(v);
                    dist[v] = nd;
                    pq.push({nd, v});
                }
            }
        }
    }
};

}  // namespace cci
