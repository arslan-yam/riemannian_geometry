// Steiner graph solver for weighted terrains (src/steiner.h): the graph, its searches, and the Solver
// interface over it.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>
#include "solver.h"
#include "steiner.h"

namespace cci {

void SteinerGraph::build(const std::vector<double>& V, const std::vector<unsigned>& F, const std::vector<double>& W, int m_) {
    m = m_;
    n = (int)(V.size() / 3); nf = (int)(F.size() / 3);
    fw = W;
    // unique edges, sorted by (min, max) vertex
    std::vector<std::pair<uint64_t, int>> half(3 * (size_t)nf); // (key, 3 f + k)
    for (int f = 0; f < nf; ++f)
        for (int k = 0; k < 3; ++k) {
            uint64_t u = F[3 * f + k], v = F[3 * f + (k + 1) % 3];
            half[3 * (size_t)f + k] = {std::min(u, v) * (uint64_t)n + std::max(u, v), 3 * f + k};
        }
    std::sort(half.begin(), half.end());
    std::vector<int32_t> fe(3 * (size_t)nf); // edge id of (face, k)
    std::vector<uint64_t> ekey;
    for (size_t i = 0; i < half.size(); ++i) {
        if (i == 0 || half[i].first != half[i - 1].first) ekey.push_back(half[i].first);
        fe[half[i].second] = (int32_t)(ekey.size() - 1);
    }
    std::vector<std::pair<uint64_t, int>>().swap(half);
    const int64_t ne = (int64_t)ekey.size();
    std::vector<double> wmin(ne, INFINITY);
    for (int f = 0; f < nf; ++f) for (int k = 0; k < 3; ++k) wmin[fe[3 * (size_t)f + k]] = std::min(wmin[fe[3 * (size_t)f + k]], fw[f]);
    N = n + ne * m;
    if (N >= (int64_t)INT32_MAX) throw std::length_error("Steiner graph: more than 2^31 nodes");
    P.resize(3 * N);
    std::copy(V.begin(), V.end(), P.begin());
    for (int64_t e = 0; e < ne; ++e) {
        int64_t a = (int64_t)(ekey[e] / n), b = (int64_t)(ekey[e] % n);
        for (int i = 0; i < m; ++i) {
            double t = (i + 1.0) / (m + 1.0);
            for (int c = 0; c < 3; ++c) P[3 * (n + e * m + i) + c] = V[3 * a + c] + t * (V[3 * b + c] - V[3 * a + c]);
        }
    }
    // boundary nodes of each face: v0, v1, v2, then the points of edges (v0 v1), (v1 v2), (v2 v0) in that
    // direction; mask = the face edges a node lies on
    const int Bn = B();
    fpts.resize((size_t)nf * Bn);
    std::vector<int> mask(Bn, 0);
    mask[0] = 1 | 4; mask[1] = 1 | 2; mask[2] = 2 | 4;
    for (int k = 0; k < 3; ++k) for (int i = 0; i < m; ++i) mask[3 + k * m + i] = 1 << k;
    std::vector<std::pair<int, int>> pairs; // index pairs joined inside a face
    for (int i = 0; i < Bn; ++i) for (int j = i + 1; j < Bn; ++j) if (!(mask[i] & mask[j])) pairs.push_back({i, j});
    for (int f = 0; f < nf; ++f) {
        int32_t* q = &fpts[(size_t)f * Bn];
        for (int k = 0; k < 3; ++k) q[k] = (int32_t)F[3 * f + k];
        for (int k = 0; k < 3; ++k) {
            int64_t e = fe[3 * (size_t)f + k];
            bool fwd = (uint64_t)F[3 * f + k] == ekey[e] / n; // stored as (min -> max)
            for (int i = 0; i < m; ++i) q[3 + k * m + i] = (int32_t)(n + e * m + (fwd ? i : m - 1 - i));
        }
    }
    // CSR: degrees, then fill
    std::vector<int64_t> deg(N + 1, 0);
    auto chain_node = [&](int64_t e, int i) -> int64_t { // i = 0 .. m + 1 along the edge (min -> max)
        if (i == 0) return (int64_t)(ekey[e] / n);
        if (i == m + 1) return (int64_t)(ekey[e] % n);
        return n + e * m + (i - 1);
    };
    for (int64_t e = 0; e < ne; ++e) for (int i = 0; i <= m; ++i) { deg[chain_node(e, i) + 1]++; deg[chain_node(e, i + 1) + 1]++; }
    for (int f = 0; f < nf; ++f) {
        const int32_t* q = &fpts[(size_t)f * Bn];
        for (auto& pr : pairs) { deg[q[pr.first] + 1]++; deg[q[pr.second] + 1]++; }
    }
    ptr.assign(N + 1, 0);
    for (int64_t i = 0; i < N; ++i) ptr[i + 1] = ptr[i] + deg[i + 1];
    std::vector<int64_t>().swap(deg);
    adj.resize(ptr[N]); wt.resize(ptr[N]);
    std::vector<int64_t> fill(ptr.begin(), ptr.end() - 1);
    auto add = [&](int64_t a, int64_t b, double w) {
        double L = w * std::sqrt((P[3 * a] - P[3 * b]) * (P[3 * a] - P[3 * b]) + (P[3 * a + 1] - P[3 * b + 1]) * (P[3 * a + 1] - P[3 * b + 1]) +
                                 (P[3 * a + 2] - P[3 * b + 2]) * (P[3 * a + 2] - P[3 * b + 2]));
        adj[fill[a]] = (int32_t)b; wt[fill[a]++] = L;
        adj[fill[b]] = (int32_t)a; wt[fill[b]++] = L;
    };
    for (int64_t e = 0; e < ne; ++e) for (int i = 0; i <= m; ++i) add(chain_node(e, i), chain_node(e, i + 1), wmin[e]);
    for (int f = 0; f < nf; ++f) {
        const int32_t* q = &fpts[(size_t)f * Bn];
        for (auto& pr : pairs) add(q[pr.first], q[pr.second], fw[f]);
    }
}

void SteinerGraph::init(State& S) const { S.D.assign(N, INFINITY); S.done.assign(N, 0); S.want.assign(N, 0); S.touched.clear(); }

void SteinerGraph::search(State& S, const Seeds& seeds, const std::vector<int32_t>& stop) const {
    if ((int64_t)S.D.size() != N) init(S);
    for (int32_t v : S.touched) { S.D[v] = INFINITY; S.done[v] = 0; }
    S.touched.clear();
    int64_t remaining = 0;
    for (int32_t v : stop) if (!S.want[v]) { S.want[v] = 1; ++remaining; }
    using item = std::pair<double, int32_t>;
    std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
    for (auto& sd : seeds)
        if (sd.second < S.D[sd.first]) { if (S.D[sd.first] == INFINITY) S.touched.push_back(sd.first); S.D[sd.first] = sd.second; pq.push({sd.second, sd.first}); }
    while (!pq.empty() && remaining > 0) {
        auto [d, u] = pq.top(); pq.pop();
        if (S.done[u] || d > S.D[u]) continue;
        S.done[u] = 1;
        if (S.want[u]) --remaining;
        for (int64_t k = ptr[u]; k < ptr[u + 1]; ++k) {
            int32_t v = adj[k];
            double nd = d + wt[k];
            if (nd < S.D[v]) { if (S.D[v] == INFINITY) S.touched.push_back(v); S.D[v] = nd; pq.push({nd, v}); }
        }
    }
    for (int32_t v : stop) S.want[v] = 0;
}

void SteinerGraph::point_seeds(int f, double x, double y, double z, Seeds& out) const {
    out.clear();
    for (int i = 0; i < B(); ++i) {
        int32_t v = fpts[(size_t)f * B() + i];
        double dx = P[3 * v] - x, dy = P[3 * v + 1] - y, dz = P[3 * v + 2] - z;
        out.push_back({v, fw[f] * std::sqrt(dx * dx + dy * dy + dz * dz)});
    }
}

double SteinerGraph::point_dist(const State& S, int f, double x, double y, double z) const {
    double best = INFINITY;
    for (int i = 0; i < B(); ++i) {
        int32_t v = fpts[(size_t)f * B() + i];
        double dx = P[3 * v] - x, dy = P[3 * v + 1] - y, dz = P[3 * v + 2] - z;
        best = std::min(best, S.D[v] + fw[f] * std::sqrt(dx * dx + dy * dy + dz * dz));
    }
    return best;
}

// the Solver interface over a graph of its own (sub-meshes, the truth on a whole mesh)
struct SteinerSolver : Solver {
    SteinerGraph G;
    SteinerGraph::State S;
    SteinerSolver(const std::vector<double>& V, const std::vector<unsigned>& F, const std::vector<double>& W, int m) { G.build(V, F, W, m); G.init(S); }
    void dist(unsigned src, const std::vector<unsigned>& targets, double* out) override {
        std::vector<int32_t> stop(targets.begin(), targets.end());
        G.search(S, {{(int32_t)src, 0.0}}, stop);
        for (size_t i = 0; i < targets.size(); ++i) out[i] = S.D[targets[i]];
    }
    // surface points: seeded at the boundary nodes of the source's face; a target is reached through a
    // boundary node of its face, or directly when both lie in the same face
    bool point_dists(const SurfPt& s, const std::vector<SurfPt>& targets, double* out) override {
        SteinerGraph::Seeds seeds;
        G.point_seeds((int)s.face, s.x, s.y, s.z, seeds);
        std::vector<int32_t> stop;
        for (auto& t : targets) for (int i = 0; i < G.B(); ++i) stop.push_back(G.fpts[(size_t)t.face * G.B() + i]);
        G.search(S, seeds, stop);
        for (size_t k = 0; k < targets.size(); ++k) {
            const SurfPt& t = targets[k];
            double best = G.point_dist(S, (int)t.face, t.x, t.y, t.z);
            if (t.face == s.face) best = std::min(best, G.fw[s.face] * std::sqrt((s.x - t.x) * (s.x - t.x) + (s.y - t.y) * (s.y - t.y) + (s.z - t.z) * (s.z - t.z)));
            out[k] = best;
        }
        return true;
    }
    const char* name() const override { return "Steiner graph"; }
};

Solver* make_steiner(const std::vector<double>& V, const std::vector<unsigned>& F, const std::vector<double>& W, int m) {
    return new SteinerSolver(V, F, W, m);
}

// bytes of a Steiner solver per face (for the build's memory budget)
double steiner_bytes_per_face(int m) {
    const int B = 3 + 3 * m;
    double pairs = B * (B - 1) / 2.0 - 3.0 * (m + 2) * (m + 1) / 2.0;
    return pairs * 2 * (4 + 8) + 1.5 * m * (24 + 8 + 8 + 2) + B * 4 + 64;
}

}  // namespace cci
