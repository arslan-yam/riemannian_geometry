// SEP baseline: our reimplementation of Gotsman and Hormann, "Compressing geodesic information for fast
// point-to-point geodesic distance queries" (SIGGRAPH Asia 2022), SEPfun variant, with exact VTP distances.
// Their code is MATLAB and not public; this follows Sections 3-4 of the paper:
//   - nested bisection by a Fiedler-vector graph cut of each patch's vertex graph, down to Vmax = 30
//     vertices per leaf; disconnected patches are split by components with an empty separator;
//   - separator curves through the midpoints of the cut edges, sub-sampled by three and smoothed by a few
//     filter-and-project steps, then sampled at a regular arc spacing; the faces the cut crosses belong to
//     both children (their other vertices are "ghost" vertices of the child);
//   - for every vertex of a face holding a sample, one exact sweep over the patch (d restricted to P(N));
//     the distance from a vertex x to a sample s in face [v1 v2 v3] by Rustamov's formula
//       d(x,s)^2 ~ sum_i l_i d(x,v_i)^2 - sum_{i<j} l_i l_j |v_i - v_j|^2;
//   - per vertex and ancestor curve, the better of a weighted least-squares polynomial fit to d or to d^2,
//     weight 1/d, degree raised from 2 until the mean relative error is below 0.5% or the degree is 12;
//     loops are shifted so the largest distance sits at the interval end; the minimum of each fit is kept
//     for the quick reject;
//   - a query scans the curves of all common ancestors bottom-up, skips a curve when m_x + m_y >= best,
//     and minimises fit_x + fit_y on a parameter grid of spacing 0.005 (mesh scaled to the unit cube);
//     same-leaf pairs use the Euclidean distance (their setting) or, with --leaf-exact, a stored exact one.
// Choices the paper leaves open, ours: the Fiedler vector of the random-walk Laplacian by multilevel
// coarsening and lazy power iteration (any graph cut is allowed there); smoothing by five Laplacian steps in
// the chart, each point lifted onto the surface and kept unsmoothed if it would leave the patch; samples a
// third of a median edge apart (three points per cut edge), at least two samples per coefficient; a Chebyshev basis; weight 1/d^3 for the d^2 fit (the d fit's 1/d in relative terms).
//
// usage: sep_bench mesh.off [pairs truth] [--json out.json] [--vmax 30] [--gap 0.333] [--tol 0.005]
//        [--maxdeg 12] [--leaf-exact] [--random K] (K random pairs with exact MMP truth, for testing)
#include "../src/mesh_io.h"
#include "../src/solver.h"
#include "../src/terrain.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cci;
typedef std::chrono::steady_clock Clock;
static double secs(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// ---------------------------------------------------------------- Fiedler vector (multilevel)
struct WGraph {  // CSR with edge and node weights
    std::vector<int> ptr, adj;
    std::vector<double> w, mass;
    int n() const { return (int)mass.size(); }
};

// heavy-edge matching; returns the coarse graph and the fine -> coarse map
static WGraph coarsen(const WGraph& G, std::vector<int>& cmap, std::mt19937& rng) {
    int n = G.n();
    std::vector<int> order(n), match(n, -1);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    for (int u : order) {
        if (match[u] >= 0) continue;
        int best = -1; double bw = -1;
        for (int k = G.ptr[u]; k < G.ptr[u + 1]; ++k) {
            int v = G.adj[k];
            if (v != u && match[v] < 0 && G.w[k] > bw) { bw = G.w[k]; best = v; }
        }
        match[u] = best >= 0 ? best : u;
        if (best >= 0) match[best] = u;
    }
    cmap.assign(n, -1);
    int nc = 0;
    for (int u = 0; u < n; ++u)
        if (cmap[u] < 0) { cmap[u] = nc; if (match[u] != u) cmap[match[u]] = nc; ++nc; }
    WGraph C;
    C.mass.assign(nc, 0);
    for (int u = 0; u < n; ++u) C.mass[cmap[u]] += G.mass[u];
    std::vector<std::vector<std::pair<int, double>>> rows(nc);
    for (int u = 0; u < n; ++u)
        for (int k = G.ptr[u]; k < G.ptr[u + 1]; ++k) {
            int a = cmap[u], b = cmap[G.adj[k]];
            if (a != b) rows[a].push_back({b, G.w[k]});
        }
    C.ptr.assign(nc + 1, 0);
    for (int a = 0; a < nc; ++a) {
        auto& r = rows[a];
        std::sort(r.begin(), r.end());
        int k = 0;
        for (size_t i = 0; i < r.size(); ++i) {
            if (k > 0 && r[k - 1].first == r[i].first) r[k - 1].second += r[i].second;
            else r[k++] = r[i];
        }
        r.resize(k);
        C.ptr[a + 1] = C.ptr[a] + k;
    }
    for (auto& r : rows) for (auto& e : r) { C.adj.push_back(e.first); C.w.push_back(e.second); }
    return C;
}

// lazy random-walk power iteration with the stationary (constant) component removed
static void smooth_fiedler(const WGraph& G, std::vector<double>& x, int iters) {
    int n = G.n();
    std::vector<double> deg(n, 0), y(n);
    double dsum = 0;
    for (int u = 0; u < n; ++u) { for (int k = G.ptr[u]; k < G.ptr[u + 1]; ++k) deg[u] += G.w[k]; dsum += deg[u]; }
    for (int it = 0; it < iters; ++it) {
        for (int u = 0; u < n; ++u) {
            double s = 0;
            for (int k = G.ptr[u]; k < G.ptr[u + 1]; ++k) s += G.w[k] * x[G.adj[k]];
            y[u] = deg[u] > 0 ? 0.5 * x[u] + 0.5 * s / deg[u] : x[u];
        }
        double mean = 0;  // remove the degree-weighted mean (the stationary direction)
        for (int u = 0; u < n; ++u) mean += deg[u] * y[u];
        mean /= std::max(dsum, 1e-300);
        double nrm = 0;
        for (int u = 0; u < n; ++u) { y[u] -= mean; nrm += deg[u] * y[u] * y[u]; }
        nrm = std::sqrt(std::max(nrm, 1e-300));
        for (int u = 0; u < n; ++u) x[u] = y[u] / nrm;
    }
}

// dense symmetric eigen-decomposition (cyclic Jacobi), eigenvalues ascending with vectors in columns
static void jacobi_eig(std::vector<double> A, int n, std::vector<double>& val, std::vector<double>& vec) {
    vec.assign(n * n, 0);
    for (int i = 0; i < n; ++i) vec[i * n + i] = 1;
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0;
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) off += A[i * n + j] * A[i * n + j];
        if (off < 1e-22) break;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) {
                double apq = A[p * n + q];
                if (std::fabs(apq) < 1e-300) continue;
                double th = 0.5 * (A[q * n + q] - A[p * n + p]) / apq;
                double t = (th >= 0 ? 1 : -1) / (std::fabs(th) + std::sqrt(th * th + 1));
                double c = 1 / std::sqrt(t * t + 1), s = t * c;
                for (int k = 0; k < n; ++k) {
                    double akp = A[k * n + p], akq = A[k * n + q];
                    A[k * n + p] = c * akp - s * akq; A[k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    double apk = A[p * n + k], aqk = A[q * n + k];
                    A[p * n + k] = c * apk - s * aqk; A[q * n + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) {
                    double vkp = vec[k * n + p], vkq = vec[k * n + q];
                    vec[k * n + p] = c * vkp - s * vkq; vec[k * n + q] = s * vkp + c * vkq;
                }
            }
    }
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return A[a * n + a] < A[b * n + b]; });
    std::vector<double> v2(n * n);
    val.resize(n);
    for (int j = 0; j < n; ++j) { val[j] = A[idx[j] * n + idx[j]]; for (int i = 0; i < n; ++i) v2[i * n + j] = vec[i * n + idx[j]]; }
    vec = v2;
}

// Fiedler vector of the random-walk Laplacian of a connected graph
static std::vector<double> fiedler(const WGraph& G0, std::mt19937& rng) {
    std::vector<WGraph> levels{G0};
    std::vector<std::vector<int>> maps;
    while (levels.back().n() > 200) {
        std::vector<int> cm;
        WGraph C = coarsen(levels.back(), cm, rng);
        if (C.n() > 0.9 * levels.back().n()) break;
        maps.push_back(cm);
        levels.push_back(std::move(C));
    }
    const WGraph& C = levels.back();
    int nc = C.n();
    // normalised Laplacian I - D^-1/2 W D^-1/2 on the coarsest graph; Fiedler = D^-1/2 times its 2nd vector
    std::vector<double> deg(nc, 0), A(nc * nc, 0);
    for (int u = 0; u < nc; ++u) for (int k = C.ptr[u]; k < C.ptr[u + 1]; ++k) deg[u] += C.w[k];
    for (int u = 0; u < nc; ++u) {
        A[u * nc + u] = 1;
        for (int k = C.ptr[u]; k < C.ptr[u + 1]; ++k) {
            int v = C.adj[k];
            A[u * nc + v] -= C.w[k] / std::sqrt(std::max(deg[u] * deg[v], 1e-300));
        }
    }
    std::vector<double> val, vec, x(nc);
    jacobi_eig(A, nc, val, vec);
    for (int u = 0; u < nc; ++u) x[u] = vec[u * nc + std::min(1, nc - 1)] / std::sqrt(std::max(deg[u], 1e-300));
    smooth_fiedler(C, x, 20);
    for (int l = (int)levels.size() - 2; l >= 0; --l) {
        std::vector<double> xf(levels[l].n());
        for (int u = 0; u < levels[l].n(); ++u) xf[u] = x[maps[l][u]];
        x.swap(xf);
        smooth_fiedler(levels[l], x, 30);
    }
    return x;
}

// ---------------------------------------------------------------- polynomial fits (Chebyshev basis)
static inline double cheb(const double* c, int deg, double u) {  // Clenshaw
    double b1 = 0, b2 = 0;
    for (int k = deg; k >= 1; --k) { double b0 = 2 * u * b1 - b2 + c[k]; b2 = b1; b1 = b0; }
    return u * b1 - b2 + c[0];
}

// weighted least squares fit of degree deg; false if singular
static bool lsq(const std::vector<double>& u, const std::vector<double>& f, const std::vector<double>& w, int deg, double* c) {
    int m = deg + 1;
    double M[13 * 13] = {0}, r[13] = {0}, T[13];
    for (size_t i = 0; i < u.size(); ++i) {
        T[0] = 1; if (m > 1) T[1] = u[i];
        for (int k = 2; k < m; ++k) T[k] = 2 * u[i] * T[k - 1] - T[k - 2];
        for (int a = 0; a < m; ++a) { r[a] += w[i] * T[a] * f[i]; for (int b = 0; b <= a; ++b) M[a * m + b] += w[i] * T[a] * T[b]; }
    }
    for (int a = 0; a < m; ++a) for (int b = a + 1; b < m; ++b) M[a * m + b] = M[b * m + a];
    // Cholesky
    double L[13 * 13] = {0};
    for (int j = 0; j < m; ++j) {
        double s = M[j * m + j];
        for (int k = 0; k < j; ++k) s -= L[j * m + k] * L[j * m + k];
        if (s <= 1e-300) return false;
        L[j * m + j] = std::sqrt(s);
        for (int i = j + 1; i < m; ++i) {
            double t = M[i * m + j];
            for (int k = 0; k < j; ++k) t -= L[i * m + k] * L[j * m + k];
            L[i * m + j] = t / L[j * m + j];
        }
    }
    double y[13];
    for (int i = 0; i < m; ++i) { double t = r[i]; for (int k = 0; k < i; ++k) t -= L[i * m + k] * y[k]; y[i] = t / L[i * m + i]; }
    for (int i = m - 1; i >= 0; --i) { double t = y[i]; for (int k = i + 1; k < m; ++k) t -= L[k * m + i] * c[k]; c[i] = t / L[i * m + i]; }
    return true;
}

// ---------------------------------------------------------------- index
struct Curve { double len = 0; bool loop = false; int grid = 1; std::vector<int> sf; std::vector<double> sl, st; };  // samples kept for debugging  // grid: query evaluations along the curve
struct Fit { int deg = 0; bool sq = false; double shift = 0, mn = 0; uint32_t coef = 0; };  // coef: offset into coeffs
struct Node { int parent = -1, depth = 0; std::vector<int> curves; };

struct Sep {
    Terrain T; FaceGrid G;
    int vmax = 30, maxdeg = 12; double gap = 1.0 / 3, tol = 0.005, scale = 1, edge3 = 1;
    std::vector<Node> nodes; std::vector<Curve> curves;
    std::vector<int> leaf;                       // vertex -> leaf node
    std::vector<std::vector<Fit>> fits;          // vertex -> fits of its ancestor curves, root first
    std::vector<std::vector<int>> fit_off;       // vertex -> index of the first fit of the ancestor at each depth
    std::vector<double> coeffs;
    std::unordered_map<uint64_t, float> leafD;   // exact same-leaf distances (--leaf-exact)
    double t_tree = 0, t_sweep = 0, t_fit = 0; long sweeps = 0;

    double P3(int v, int k) const { return T.V[3 * v + k]; }

    // node of depth d on the path from leaf l to the root
    int anc(int l, int d) const { while (nodes[l].depth > d) l = nodes[l].parent; return l; }
};

// one patch: owned vertices (those whose leaf lies below) and the faces of P(N)
struct Patch { int node; std::vector<int> own, faces; };

static void split(Sep& S, Patch& P, std::vector<Patch>& todo, std::mt19937& rng) {
    const Terrain& T = S.T;
    int node = P.node;
    if ((int)P.own.size() <= S.vmax) {
        for (int v : P.own) S.leaf[v] = node;
        return;
    }
    auto t0 = Clock::now();
    // vertex graph on the owned vertices
    std::unordered_map<int, int> loc;
    for (size_t i = 0; i < P.own.size(); ++i) loc[P.own[i]] = (int)i;
    int n = (int)P.own.size();
    std::vector<std::vector<int>> nb(n);
    for (int f : P.faces)
        for (int k = 0; k < 3; ++k) {
            auto a = loc.find(T.F[3 * f + k]), b = loc.find(T.F[3 * f + (k + 1) % 3]);
            if (a != loc.end() && b != loc.end()) { nb[a->second].push_back(b->second); nb[b->second].push_back(a->second); }
        }
    for (auto& r : nb) { std::sort(r.begin(), r.end()); r.erase(std::unique(r.begin(), r.end()), r.end()); }
    // components
    std::vector<int> comp(n, -1); int nc = 0;
    for (int s = 0; s < n; ++s) {
        if (comp[s] >= 0) continue;
        std::vector<int> st{s}; comp[s] = nc;
        while (!st.empty()) { int u = st.back(); st.pop_back(); for (int v : nb[u]) if (comp[v] < 0) { comp[v] = nc; st.push_back(v); } }
        ++nc;
    }
    std::vector<char> side(n, 0);
    bool has_sep = false;
    if (nc > 1) {  // assign whole components, largest first, to the lighter side: empty separator
        std::vector<int> sz(nc, 0); for (int c : comp) sz[c]++;
        std::vector<int> ord(nc); std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(), [&](int a, int b) { return sz[a] > sz[b]; });
        std::vector<char> cs(nc); long w0 = 0, w1 = 0;
        for (int c : ord) { if (w0 <= w1) { cs[c] = 0; w0 += sz[c]; } else { cs[c] = 1; w1 += sz[c]; } }
        for (int i = 0; i < n; ++i) side[i] = cs[comp[i]];
    } else {
        WGraph W; W.mass.assign(n, 1); W.ptr.assign(n + 1, 0);
        for (int i = 0; i < n; ++i) { W.ptr[i + 1] = W.ptr[i] + (int)nb[i].size(); for (int j : nb[i]) { W.adj.push_back(j); W.w.push_back(1); } }
        std::vector<double> x = fiedler(W, rng);
        std::vector<int> ord(n); std::iota(ord.begin(), ord.end(), 0);
        std::nth_element(ord.begin(), ord.begin() + n / 2, ord.end(), [&](int a, int b) { return x[a] < x[b]; });
        for (int i = n / 2; i < n; ++i) side[ord[i]] = 1;
        has_sep = true;
    }
    auto sideof = [&](unsigned v) -> int { auto it = loc.find(v); return it == loc.end() ? -1 : side[it->second]; };
    // children: faces with an owned vertex on that side (faces the cut crosses go to both)
    Patch C[2];
    for (int k = 0; k < 2; ++k) {
        C[k].node = (int)S.nodes.size();
        Node nd; nd.parent = node; nd.depth = S.nodes[node].depth + 1; S.nodes.push_back(nd);
    }
    for (int i = 0; i < n; ++i) C[(int)side[i]].own.push_back(P.own[i]);
    for (int f : P.faces) {
        bool h[2] = {false, false};
        for (int k = 0; k < 3; ++k) { int s = sideof(T.F[3 * f + k]); if (s >= 0) h[s] = true; }
        for (int k = 0; k < 2; ++k) if (h[k]) C[k].faces.push_back(f);
    }
    S.t_tree += secs(t0);
    if (has_sep) {
        // ---- separator curves through the midpoints of the cut edges
        t0 = Clock::now();
        std::vector<char> inP(T.m, 0);
        for (int f : P.faces) inP[f] = 1;
        std::unordered_map<uint64_t, int> mid;  // cut edge -> midpoint id
        std::vector<std::pair<int, int>> mids;  // the edge's vertices
        std::vector<int> mface;                 // a face containing the midpoint
        std::vector<std::vector<int>> link;
        auto key = [](unsigned a, unsigned b) { if (a > b) std::swap(a, b); return ((uint64_t)a << 32) | b; };
        for (int f : P.faces) {
            int ids[3], c = 0;
            for (int k = 0; k < 3; ++k) {
                unsigned a = T.F[3 * f + k], b = T.F[3 * f + (k + 1) % 3];
                int sa = sideof(a), sb = sideof(b);
                if (sa < 0 || sb < 0 || sa == sb) continue;
                auto it = mid.find(key(a, b));
                int id;
                if (it == mid.end()) { id = (int)mids.size(); mid[key(a, b)] = id; mids.push_back({(int)a, (int)b}); mface.push_back(f); link.emplace_back(); }
                else id = it->second;
                ids[c++] = id;
            }
            if (c == 2) { link[ids[0]].push_back(ids[1]); link[ids[1]].push_back(ids[0]); }
        }
        // chains: paths from the ends first, then loops
        std::vector<char> used(mids.size(), 0);
        std::vector<std::pair<std::vector<int>, bool>> chains;
        for (int pass = 0; pass < 2; ++pass)
            for (size_t s = 0; s < mids.size(); ++s) {
                if (used[s] || (pass == 0 && link[s].size() == 2)) continue;
                std::vector<int> ch{(int)s}; used[s] = 1;
                int prev = -1, cur = (int)s;
                while (true) {
                    int nx = -1;
                    for (int v : link[cur]) if (v != prev && !used[v]) { nx = v; break; }
                    if (nx < 0) break;
                    used[nx] = 1; ch.push_back(nx); prev = cur; cur = nx;
                }
                chains.push_back({ch, pass == 1});
            }
        // per curve: points (chart), smoothed, lifted, sampled
        for (auto& chp : chains) {
            auto& ch = chp.first; bool loop = chp.second && ch.size() > 2;
            std::vector<double> X, Y; std::vector<int> FF;
            auto midpt = [&](int id, double& x, double& y) {
                int a = mids[id].first, b = mids[id].second;
                x = 0.5 * (T.xy[2 * a] + T.xy[2 * b]); y = 0.5 * (T.xy[2 * a + 1] + T.xy[2 * b + 1]);
            };
            size_t K = ch.size();
            for (size_t i = 0; i < K; ++i) {  // sub-sample by three
                double x0, y0; midpt(ch[i], x0, y0);
                X.push_back(x0); Y.push_back(y0); FF.push_back(mface[ch[i]]);
                if (i + 1 < K || loop) {
                    double x1, y1; midpt(ch[(i + 1) % K], x1, y1);
                    for (int j = 1; j < 3; ++j) { X.push_back(x0 + (x1 - x0) * j / 3.0); Y.push_back(y0 + (y1 - y0) * j / 3.0); FF.push_back(-1); }
                }
            }
            size_t N = X.size();
            std::vector<double> X0 = X, Y0 = Y;
            for (int it = 0; it < 5; ++it) {  // filter and project back (onto the patch)
                std::vector<double> X1 = X, Y1 = Y;
                for (size_t i = 0; i < N; ++i) {
                    if (!loop && (i == 0 || i + 1 == N)) continue;
                    size_t a = (i + N - 1) % N, b = (i + 1) % N;
                    X1[i] = 0.5 * X[i] + 0.25 * (X[a] + X[b]); Y1[i] = 0.5 * Y[i] + 0.25 * (Y[a] + Y[b]);
                }
                X = X1; Y = Y1;
            }
            std::vector<double> PX, PY, PZ; std::vector<int> PF;
            for (size_t i = 0; i < N; ++i) {
                int f = locate(T, S.G, X[i], Y[i]);
                if (f < 0 || !inP[f]) { X[i] = X0[i]; Y[i] = Y0[i]; f = locate(T, S.G, X[i], Y[i]); }
                if (f < 0 || !inP[f]) continue;
                PX.push_back(X[i]); PY.push_back(Y[i]); PZ.push_back(T.z_at(f, X[i], Y[i])); PF.push_back(f);
            }
            if (PX.size() < 2) continue;
            if (loop) { PX.push_back(PX[0]); PY.push_back(PY[0]); PZ.push_back(PZ[0]); PF.push_back(PF[0]); }
            // regular samples along the lifted polyline
            std::vector<double> cum(PX.size(), 0);
            for (size_t i = 1; i < PX.size(); ++i) cum[i] = cum[i - 1] + std::sqrt((PX[i] - PX[i - 1]) * (PX[i] - PX[i - 1]) + (PY[i] - PY[i - 1]) * (PY[i] - PY[i - 1]) + (PZ[i] - PZ[i - 1]) * (PZ[i] - PZ[i - 1]));
            double L = cum.back(), h = S.gap * S.edge3;
            int ns = std::max(2, (int)std::ceil(L / h) + (loop ? 0 : 1));
            std::vector<int> sf; std::vector<double> sl, st;  // sample face, barycentrics, parameter
            size_t seg = 0;
            for (int j = 0; j < ns; ++j) {
                double t = loop ? L * j / ns : L * j / (ns - 1);
                while (seg + 2 < cum.size() && cum[seg + 1] < t) ++seg;
                double a = cum[seg + 1] > cum[seg] ? (t - cum[seg]) / (cum[seg + 1] - cum[seg]) : 0;
                double x = PX[seg] + a * (PX[seg + 1] - PX[seg]), y = PY[seg] + a * (PY[seg + 1] - PY[seg]);
                int f = locate(T, S.G, x, y);
                if (f < 0 || !inP[f]) { f = PF[a < 0.5 ? seg : seg + 1]; x = a < 0.5 ? PX[seg] : PX[seg + 1]; y = a < 0.5 ? PY[seg] : PY[seg + 1]; }
                double l[3]; bary(T, f, x, y, l);
                for (int k = 0; k < 3; ++k) l[k] = std::max(0.0, l[k]);
                double ls = l[0] + l[1] + l[2];
                sf.push_back(f); for (int k = 0; k < 3; ++k) sl.push_back(l[k] / ls); st.push_back(t);
            }
            int ci = (int)S.curves.size();
            Curve cv; cv.len = L; cv.loop = loop;
            cv.grid = std::max(2, (int)std::ceil(L / S.scale / 0.005) + 1);
            if (getenv("SEP_DEBUG")) { cv.sf = sf; cv.sl = sl; cv.st = st; }
            S.curves.push_back(cv);
            S.nodes[node].curves.push_back(ci);
            S.t_tree += secs(t0);
            // ---- exact sweeps over the patch from the vertices of the sample faces
            t0 = Clock::now();
            std::vector<int> src;
            for (int f : sf) for (int k = 0; k < 3; ++k) src.push_back((int)T.F[3 * f + k]);
            std::sort(src.begin(), src.end()); src.erase(std::unique(src.begin(), src.end()), src.end());
            std::vector<int> pv;  // patch vertices, local numbering
            std::unordered_map<int, int> pl;
            std::vector<double> V; std::vector<unsigned> F;
            for (int f : P.faces)
                for (int k = 0; k < 3; ++k) {
                    int v = (int)T.F[3 * f + k];
                    auto it = pl.find(v);
                    if (it == pl.end()) { it = pl.emplace(v, (int)pv.size()).first; pv.push_back(v); for (int d = 0; d < 3; ++d) V.push_back(T.V[3 * v + d]); }
                    F.push_back((unsigned)it->second);
                }
            std::unique_ptr<Solver> sol(make_vtp(V, F, true));
            if (!sol) sol.reset(make_mmp(V, F, false));
            std::vector<unsigned> tg(P.own.size());
            for (size_t i = 0; i < P.own.size(); ++i) tg[i] = (unsigned)pl.at(P.own[i]);
            std::unordered_map<int, int> si;
            std::vector<float> D(src.size() * P.own.size());
            std::vector<double> out(P.own.size());
            for (size_t s = 0; s < src.size(); ++s) {
                si[src[s]] = (int)s;
                sol->dist((unsigned)pl.at(src[s]), tg, out.data());
                for (size_t i = 0; i < out.size(); ++i) D[s * P.own.size() + i] = (float)out[i];
                ++S.sweeps;
            }
            S.t_sweep += secs(t0);
            // ---- per owned vertex: Rustamov distances to the samples, polynomial fit
            t0 = Clock::now();
            int ns2 = (int)sf.size();
            std::vector<double> u(ns2), d(ns2), f2(ns2), w(ns2), wq(ns2);
            double e3[3];  // squared 3D edge lengths of each sample face, recomputed per sample below
            for (size_t i = 0; i < P.own.size(); ++i) {
                int x = P.own[i];
                for (int j = 0; j < ns2; ++j) {
                    int f = sf[j];
                    double dv[3];
                    for (int k = 0; k < 3; ++k) dv[k] = D[si[(int)T.F[3 * f + k]] * P.own.size() + i];
                    for (int k = 0; k < 3; ++k) {
                        int a = (int)T.F[3 * f + k], b = (int)T.F[3 * f + (k + 1) % 3];
                        e3[k] = 0; for (int q = 0; q < 3; ++q) e3[k] += (S.P3(a, q) - S.P3(b, q)) * (S.P3(a, q) - S.P3(b, q));
                    }
                    const double* l = &sl[3 * j];
                    double v2 = l[0] * dv[0] * dv[0] + l[1] * dv[1] * dv[1] + l[2] * dv[2] * dv[2]
                              - l[0] * l[1] * e3[0] - l[1] * l[2] * e3[1] - l[2] * l[0] * e3[2];
                    d[j] = std::sqrt(std::max(v2, 0.0));
                }
                // loops: shift so the largest distance sits at the interval end
                double shift = 0; std::vector<double> dd = d, tt = st;
                if (loop) {
                    int am = (int)(std::max_element(d.begin(), d.end()) - d.begin());
                    shift = st[am];
                    for (int j = 0; j < ns2; ++j) { int jj = (am + j) % ns2; dd[j] = d[jj]; tt[j] = std::fmod(st[jj] - shift + L, L); }
                }
                double dmin = 1e-9 * S.scale;
                for (int j = 0; j < ns2; ++j) {
                    double dj = std::max(dd[j], dmin);
                    u[j] = 2 * tt[j] / L - 1; f2[j] = dj * dj; w[j] = 1 / dj; wq[j] = 1 / (dj * dj * dj);
                    dd[j] = dj;
                }
                double best_err = 1e300; Fit bf; std::vector<double> bc;
                int top = std::min(S.maxdeg, (ns2 - 1) / 2);  // at least two samples per coefficient: no interpolation
                for (int deg = std::min(2, top); deg <= top; ++deg) {
                    for (int sq = 0; sq < 2; ++sq) {
                        double c[13];
                        if (!lsq(u, sq ? f2 : dd, sq ? wq : w, deg, c)) continue;
                        double err = 0;
                        for (int j = 0; j < ns2; ++j) {
                            double p = cheb(c, deg, u[j]);
                            if (sq) p = std::sqrt(std::max(p, 0.0));
                            err += std::fabs(p - dd[j]) / dd[j];
                        }
                        err /= ns2;
                        if (err < best_err) { best_err = err; bf.deg = deg; bf.sq = sq; bc.assign(c, c + deg + 1); }
                    }
                    if (best_err < S.tol) break;
                }
                bf.shift = shift; bf.coef = (uint32_t)S.coeffs.size();
                S.coeffs.insert(S.coeffs.end(), bc.begin(), bc.end());
                // minimum on the query grid, for the quick reject
                double mn = 1e300;
                for (int g = 0; g < cv.grid; ++g) {
                    double p = cheb(&S.coeffs[bf.coef], bf.deg, 2.0 * g / (cv.grid - 1) - 1);
                    if (bf.sq) p = std::sqrt(std::max(p, 0.0));
                    mn = std::min(mn, p);
                }
                bf.mn = mn;
                S.fits[x].push_back(bf);
            }
            S.t_fit += secs(t0);
            t0 = Clock::now();
        }
    }
    // fit offsets: this node's curves start at the current end of each owned vertex's list
    todo.push_back(std::move(C[0])); todo.push_back(std::move(C[1]));
}

// value of x's fit on curve c at curve parameter t in [0, L]
static inline double eval(const Sep& S, const Fit& f, const Curve& c, double t) {
    double tt = c.loop ? std::fmod(t - f.shift + 2 * c.len, c.len) : t;
    double p = cheb(&S.coeffs[f.coef], f.deg, 2 * tt / c.len - 1);
    return f.sq ? std::sqrt(std::max(p, 0.0)) : p;
}

static double query(const Sep& S, int x, int y, int* same_leaf = nullptr) {
    if (x == y) return 0;
    int lx = S.leaf[x], ly = S.leaf[y];
    if (lx == ly) {
        if (same_leaf) *same_leaf = 1;
        if (!S.leafD.empty()) {
            auto it = S.leafD.find(((uint64_t)std::min(x, y) << 32) | (uint64_t)std::max(x, y));
            if (it != S.leafD.end()) return it->second;
        }
        double s = 0; for (int k = 0; k < 3; ++k) s += (S.P3(x, k) - S.P3(y, k)) * (S.P3(x, k) - S.P3(y, k));
        return std::sqrt(s);
    }
    // lowest common ancestor
    int a = lx, b = ly;
    while (S.nodes[a].depth > S.nodes[b].depth) a = S.nodes[a].parent;
    while (S.nodes[b].depth > S.nodes[a].depth) b = S.nodes[b].parent;
    while (a != b) { a = S.nodes[a].parent; b = S.nodes[b].parent; }
    double best = 1e300;
    for (int d = S.nodes[a].depth; d >= 0; --d) {  // bottom-up over the common ancestors
        int nd = S.anc(a, d);
        const auto& cs = S.nodes[nd].curves;
        int ox = S.fit_off[x][d], oy = S.fit_off[y][d];
        for (size_t k = 0; k < cs.size(); ++k) {
            const Fit &fx = S.fits[x][ox + k], &fy = S.fits[y][oy + k];
            if (fx.mn + fy.mn >= best) continue;  // quick reject
            const Curve& c = S.curves[cs[k]];
            for (int g = 0; g < c.grid; ++g) {
                double t = c.len * g / (c.grid - 1);
                best = std::min(best, eval(S, fx, c, t) + eval(S, fy, c, t));
            }
        }
    }
    return best;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: sep_bench mesh.off [pairs truth] [--json out] [--vmax 30] [--gap 0.333] [--tol 0.005] [--maxdeg 12] [--leaf-exact] [--random K]\n"); return 1; }
    std::string mesh = argv[1], pairs_f, truth_f, json;
    Sep S; bool leaf_exact = false; int nrand = 0;
    int i = 2;
    if (argc > 3 && argv[2][0] != '-') { pairs_f = argv[2]; truth_f = argv[3]; i = 4; }
    for (; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--json") json = nx();
        else if (a == "--vmax") S.vmax = std::stoi(nx());
        else if (a == "--gap") S.gap = std::stod(nx());
        else if (a == "--tol") S.tol = std::stod(nx());
        else if (a == "--maxdeg") S.maxdeg = std::min(12, std::stoi(nx()));
        else if (a == "--leaf-exact") leaf_exact = true;
        else if (a == "--random") nrand = std::stoi(nx());
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 1; }
    }
    std::vector<double> V; std::vector<unsigned> F;
    read_off(mesh, V, F);
    S.T.init(V, F); S.G.init(S.T);
    const Terrain& T = S.T;
    {
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (int v = 0; v < T.n; ++v) for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], T.V[3 * v + k]); hi[k] = std::max(hi[k], T.V[3 * v + k]); }
        S.scale = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
        std::vector<double> el;
        for (int f = 0; f < T.m; ++f) {
            double s = 0; for (int k = 0; k < 3; ++k) s += (V[3 * F[3 * f] + k] - V[3 * F[3 * f + 1] + k]) * (V[3 * F[3 * f] + k] - V[3 * F[3 * f + 1] + k]);
            el.push_back(std::sqrt(s));
        }
        std::nth_element(el.begin(), el.begin() + el.size() / 2, el.end());
        S.edge3 = el[el.size() / 2];
    }
    // ---- build
    auto tb = Clock::now();
    std::mt19937 rng(1);
    S.leaf.assign(T.n, -1); S.fits.resize(T.n);
    S.nodes.push_back(Node());
    std::vector<Patch> todo(1);
    todo[0].node = 0; todo[0].own.resize(T.n); std::iota(todo[0].own.begin(), todo[0].own.end(), 0);
    todo[0].faces.resize(T.m); std::iota(todo[0].faces.begin(), todo[0].faces.end(), 0);
    // depth-first so that each vertex receives its fits root first; a FIFO would do as well, since the fits
    // of a vertex come only from its own ancestors, which are processed in depth order along any order
    while (!todo.empty()) {
        Patch P = std::move(todo.back()); todo.pop_back();
        split(S, P, todo, rng);
    }
    // per vertex: offset of the first fit of the ancestor at each depth
    S.fit_off.resize(T.n);
    int maxd = 0; double sumd = 0;
    for (int v = 0; v < T.n; ++v) {
        int l = S.leaf[v], dl = S.nodes[l].depth;
        maxd = std::max(maxd, dl); sumd += dl;
        S.fit_off[v].assign(dl + 1, 0);
        int off = 0;
        for (int d = 0; d <= dl; ++d) { S.fit_off[v][d] = off; off += (int)S.nodes[S.anc(l, d)].curves.size(); }
        if (off != (int)S.fits[v].size()) { std::fprintf(stderr, "fit bookkeeping mismatch at vertex %d (%d vs %zu)\n", v, off, S.fits[v].size()); return 2; }
    }
    double t_leaf = 0;
    if (leaf_exact) {  // exact same-leaf distances, MMP with early stop on the whole mesh
        auto t0 = Clock::now();
        std::unique_ptr<Solver> mmp(make_mmp(V, F, true));
        std::vector<std::vector<int>> members(S.nodes.size());
        for (int v = 0; v < T.n; ++v) members[S.leaf[v]].push_back(v);
        for (auto& mb : members)
            for (size_t a = 0; a + 1 < mb.size(); ++a) {
                std::vector<unsigned> tg(mb.begin() + a + 1, mb.end());
                std::vector<double> out(tg.size());
                mmp->dist((unsigned)mb[a], tg, out.data());
                for (size_t k = 0; k < tg.size(); ++k) S.leafD[((uint64_t)std::min(mb[a], (int)tg[k]) << 32) | (uint64_t)std::max(mb[a], (int)tg[k])] = (float)out[k];
            }
        t_leaf = secs(t0);
    }
    double t_build = secs(tb);
    long nfit = 0; for (auto& fv : S.fits) nfit += (long)fv.size();
    // storage as counted in the paper: doubles per vertex (coefficients; plus the shift of loop fits)
    double doubles = (double)S.coeffs.size();
    for (int v = 0; v < T.n; ++v) { int l = S.leaf[v]; for (int d = 0, o = 0; d <= S.nodes[l].depth; ++d) for (int c : S.nodes[S.anc(l, d)].curves) { if (S.curves[c].loop) doubles += 1; ++o; } }
    double bytes = 8 * doubles + 2.0 * nfit + 4.0 * T.n;  // + one byte each for degree and d/d^2, the leaf id
    if (leaf_exact) bytes += 4.0 * S.leafD.size();
    int leaves = 0; for (size_t k = 0; k < S.nodes.size(); ++k) if (S.nodes[k].curves.empty()) ++leaves;
    std::printf("SEP %s: n %d, %zu nodes, depth avg %.1f max %d, %zu curves, %ld sweeps; build %.1f s (tree %.1f, sweeps %.1f, fits %.1f, leaf %.1f); %.1f doubles/vertex, %.1f MB\n",
                mesh.c_str(), T.n, S.nodes.size(), sumd / T.n, maxd, S.curves.size(), S.sweeps, t_build, S.t_tree, S.t_sweep, S.t_fit, t_leaf, doubles / T.n, bytes / 1e6);
    // ---- queries
    std::vector<std::pair<int, int>> qp; std::vector<double> tr;
    if (!pairs_f.empty()) {
        std::ifstream pf(pairs_f), tf(truth_f); int s, t;
        while (pf >> s >> t) { qp.push_back({s, t}); std::string rest; std::getline(pf, rest); }
        tr.resize(qp.size()); for (auto& d : tr) tf >> d;
    } else if (nrand > 0) {
        std::uniform_int_distribution<int> U(0, T.n - 1);
        std::unique_ptr<Solver> mmp(make_mmp(V, F, true));
        for (int k = 0; k < nrand; ++k) {
            int s = U(rng), t = U(rng); if (s == t) continue;
            std::vector<unsigned> tg{(unsigned)t}; double d;
            mmp->dist((unsigned)s, tg, &d);
            qp.push_back({s, t}); tr.push_back(d);
        }
    }
    if (qp.empty()) return 0;
    std::vector<double> est(qp.size()); int same = 0;
    auto tq = Clock::now();
    for (size_t k = 0; k < qp.size(); ++k) est[k] = query(S, qp[k].first, qp[k].second);
    double q_us = 1e6 * secs(tq) / qp.size();
    std::vector<double> ae; double esum = 0, emax = -1e300, emin = 1e300;
    for (size_t k = 0; k < qp.size(); ++k) {
        int sl = 0; query(S, qp[k].first, qp[k].second, &sl); same += sl;
        double e = est[k] / tr[k] - 1;
        esum += e; emax = std::max(emax, e); emin = std::min(emin, e); ae.push_back(std::fabs(e));
    }
    if (getenv("SEP_DEBUG")) {
        size_t w = 0; for (size_t k = 0; k < qp.size(); ++k) if (std::fabs(est[k] / tr[k] - 1) > std::fabs(est[w] / tr[w] - 1)) w = k;
        int x = qp[w].first, y = qp[w].second;
        std::printf("worst %d %d est %.6g truth %.6g leaves %d %d depth %d %d\n", x, y, est[w], tr[w], S.leaf[x], S.leaf[y], S.nodes[S.leaf[x]].depth, S.nodes[S.leaf[y]].depth);
        for (int d = 0; d <= S.nodes[S.leaf[x]].depth; ++d) {
            int nd = S.anc(S.leaf[x], d); if (S.anc(S.leaf[y], std::min(d, S.nodes[S.leaf[y]].depth)) != nd) break;
            for (size_t k = 0; k < S.nodes[nd].curves.size(); ++k) {
                const Fit &fx = S.fits[x][S.fit_off[x][d] + k], &fy = S.fits[y][S.fit_off[y][d] + k]; const Curve& c = S.curves[S.nodes[nd].curves[k]];
                double b = 1e300; for (int g = 0; g < c.grid; ++g) { double t = c.len * g / (c.grid - 1); b = std::min(b, eval(S, fx, c, t) + eval(S, fy, c, t)); }
                if (getenv("SEP_DEBUG2") && atoi(getenv("SEP_DEBUG2")) == S.nodes[nd].curves[k]) {
                    std::unique_ptr<Solver> mm(make_mmp(V, F, false));
                    std::vector<unsigned> all(T.n); std::iota(all.begin(), all.end(), 0);
                    std::vector<double> dx(T.n), dy(T.n); mm->dist(x, all, dx.data()); mm->dist(y, all, dy.data());
                    for (size_t j = 0; j < c.sf.size(); ++j) {
                        int f = c.sf[j]; double rx = 0, ry = 0, e = 0;
                        for (int q = 0; q < 3; ++q) { int a = T.F[3*f+q], b = T.F[3*f+(q+1)%3]; double ee = 0; for (int z = 0; z < 3; ++z) ee += (S.P3(a,z)-S.P3(b,z))*(S.P3(a,z)-S.P3(b,z));
                            rx += c.sl[3*j+q]*dx[a]*dx[a]; ry += c.sl[3*j+q]*dy[a]*dy[a]; e += c.sl[3*j+q]*c.sl[3*j+(q+1)%3]*ee; }
                        std::printf("    t %.1f face %d  x: rust %.2f fit %.2f   y: rust %.2f fit %.2f\n", c.st[j], f, std::sqrt(std::max(0.0, rx-e)), eval(S, fx, c, c.st[j]), std::sqrt(std::max(0.0, ry-e)), eval(S, fy, c, c.st[j]));
                    }
                }
                std::printf("  depth %d curve %d len %.4g loop %d grid %d: degx %d sq %d mn %.4g | degy %d sq %d mn %.4g -> %.6g\n", d, S.nodes[nd].curves[k], c.len, c.loop, c.grid, fx.deg, fx.sq, fx.mn, fy.deg, fy.sq, fy.mn, b);
            }
        }
    }
    std::sort(ae.begin(), ae.end());
    double am = 0; for (double x : ae) am += x / ae.size();
    size_t m = qp.size();
    std::printf("SEP queries: %zu (%d same-leaf), %.2f us, err mean |.| %.4f%%, p99 %.4f%%, max |.| %.4f%% (signed max %.4f%%, min %.4f%%)\n",
                m, same, q_us, 100 * am, 100 * ae[(size_t)(0.99 * (m - 1))], 100 * ae.back(), 100 * emax, 100 * emin);
    if (!json.empty()) {
        std::ofstream o(json);
        o << "{\"tool\":\"sep_bench\",\"mesh\":\"" << mesh << "\",\"n\":" << T.n << ",\"faces\":" << T.m << ",\"vmax\":" << S.vmax
          << ",\"gap_edges\":" << S.gap << ",\"tol\":" << S.tol << ",\"maxdeg\":" << S.maxdeg << ",\"leaf_exact\":" << (leaf_exact ? "true" : "false")
          << ",\"nodes\":" << S.nodes.size() << ",\"leaves\":" << leaves << ",\"curves\":" << S.curves.size() << ",\"depth_avg\":" << sumd / T.n << ",\"depth_max\":" << maxd
          << ",\"sweeps\":" << S.sweeps << ",\"build_s\":" << t_build << ",\"tree_s\":" << S.t_tree << ",\"sweep_s\":" << S.t_sweep << ",\"fit_s\":" << S.t_fit << ",\"leaf_s\":" << t_leaf
          << ",\"doubles_per_vertex\":" << doubles / T.n << ",\"index_bytes\":" << (long long)bytes
          << ",\"queries\":" << m << ",\"same_leaf\":" << same << ",\"query_us\":" << q_us
          << ",\"err_mean\":" << am << ",\"err_p99\":" << ae[(size_t)(0.99 * (m - 1))] << ",\"err_max\":" << ae.back()
          << ",\"err_mean_signed\":" << esum / m << ",\"err_max_signed\":" << emax << ",\"err_min_signed\":" << emin << "}\n";
    }
}
