// Terrain geometry for the C++ Certified Cut Index: the mesh in its (x, y) chart, face planes and
// weights, a face grid (faces meeting a box), a vertex grid (radius / nearest queries), tolerant
// point location and the exact lifted length of a chart segment (face walk).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

namespace cci {

constexpr double INF = std::numeric_limits<double>::infinity();

struct Terrain {
    int n = 0, m = 0;
    std::vector<double> V;       // n x 3, translated (chart minimum at the origin)
    std::vector<unsigned> F;     // m x 3
    std::vector<double> xy;      // n x 2
    std::vector<double> a, b, c; // face planes z = a x + b y + c
    std::vector<double> w;       // face weights (1 when unweighted)
    std::vector<int> nbr;        // m x 3: face across edge k (vertex k -> k+1), -1 on the boundary
    double edge = 0;             // median chart length of the edge F[:,0] -> F[:,1]
    double lam = 1;              // min face weight
    double lo[2], hi[2];         // chart bounding box

    void init(const std::vector<double>& V_, const std::vector<unsigned>& F_, const std::vector<double>* weights = nullptr) {
        V = V_; F = F_;
        n = (int)(V.size() / 3); m = (int)(F.size() / 3);
        xy.resize(2 * n);
        lo[0] = lo[1] = INF; hi[0] = hi[1] = -INF;
        for (int i = 0; i < n; ++i) {
            xy[2 * i] = V[3 * i]; xy[2 * i + 1] = V[3 * i + 1];
            for (int k = 0; k < 2; ++k) { lo[k] = std::min(lo[k], xy[2 * i + k]); hi[k] = std::max(hi[k], xy[2 * i + k]); }
        }
        a.resize(m); b.resize(m); c.resize(m);
        std::vector<double> el(m);
        for (int f = 0; f < m; ++f) {
            const double *P0 = &V[3 * F[3 * f]], *P1 = &V[3 * F[3 * f + 1]], *P2 = &V[3 * F[3 * f + 2]];
            double x0 = P0[0], y0 = P0[1], z0 = P0[2];
            double det = (P1[0] - x0) * (P2[1] - y0) - (P2[0] - x0) * (P1[1] - y0);
            if (std::fabs(det) < 1e-300) det = 1e-300;
            a[f] = ((P1[2] - z0) * (P2[1] - y0) - (P2[2] - z0) * (P1[1] - y0)) / det;
            b[f] = ((P1[0] - x0) * (P2[2] - z0) - (P2[0] - x0) * (P1[2] - z0)) / det;
            c[f] = z0 - a[f] * x0 - b[f] * y0;
            el[f] = std::hypot(P1[0] - x0, P1[1] - y0);
        }
        std::vector<double> s(el);
        std::nth_element(s.begin(), s.begin() + m / 2, s.end());
        edge = s[m / 2];
        if (m % 2 == 0) { // numpy median of an even count: mean of the two middle values
            double lo2 = *std::max_element(s.begin(), s.begin() + m / 2);
            edge = 0.5 * (edge + lo2);
        }
        w.assign(m, 1.0);
        if (weights) w = *weights;
        lam = *std::min_element(w.begin(), w.end());
        // face adjacency through shared edges
        nbr.assign(3 * m, -1);
        std::map<std::pair<unsigned, unsigned>, int> half;
        for (int f = 0; f < m; ++f)
            for (int k = 0; k < 3; ++k) {
                unsigned u = F[3 * f + k], v = F[3 * f + (k + 1) % 3];
                auto key = std::make_pair(std::min(u, v), std::max(u, v));
                auto it = half.find(key);
                if (it == half.end()) half[key] = 3 * f + k;
                else { nbr[3 * f + k] = it->second / 3; nbr[it->second] = f; }
            }
    }
    double z_at(int f, double x, double y) const { return a[f] * x + b[f] * y + c[f]; }
};

// barycentric coordinates of (x, y) in face f
inline void bary(const Terrain& T, int f, double x, double y, double* l) {
    const double *A = &T.xy[2 * T.F[3 * f]], *B = &T.xy[2 * T.F[3 * f + 1]], *C = &T.xy[2 * T.F[3 * f + 2]];
    double det = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);
    double l1 = ((x - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (y - A[1])) / det;
    double l2 = ((B[0] - A[0]) * (y - A[1]) - (x - A[0]) * (B[1] - A[1])) / det;
    l[0] = 1.0 - l1 - l2; l[1] = l1; l[2] = l2;
}

struct FaceGrid {
    double cell = 1, org[2] = {0, 0};
    int nc[2] = {1, 1};
    std::vector<double> flo, fhi;     // m x 2 face bounding boxes
    std::vector<int64_t> ptr;         // CSR over cells
    std::vector<int> faces;
    mutable std::vector<int> stamp;   // de-duplication marks
    mutable int cur = 0;

    void init(const Terrain& T) {
        flo.resize(2 * T.m); fhi.resize(2 * T.m);
        for (int f = 0; f < T.m; ++f)
            for (int k = 0; k < 2; ++k) {
                double v0 = T.xy[2 * T.F[3 * f] + k], v1 = T.xy[2 * T.F[3 * f + 1] + k], v2 = T.xy[2 * T.F[3 * f + 2] + k];
                flo[2 * f + k] = std::min({v0, v1, v2}); fhi[2 * f + k] = std::max({v0, v1, v2});
            }
        cell = 2.0 * T.edge;
        for (int k = 0; k < 2; ++k) org[k] = T.lo[k] - 1e-9 * std::max(1.0, cell);
        int mx[2] = {0, 0};
        for (int f = 0; f < T.m; ++f)
            for (int k = 0; k < 2; ++k) mx[k] = std::max(mx[k], (int)std::floor((fhi[2 * f + k] - org[k]) / cell));
        nc[0] = mx[0] + 1; nc[1] = mx[1] + 1;
        std::vector<int64_t> cnt((size_t)nc[0] * nc[1] + 1, 0);
        auto range = [&](int f, int k, int& i0, int& i1) {
            i0 = (int)std::floor((flo[2 * f + k] - org[k]) / cell); i1 = (int)std::floor((fhi[2 * f + k] - org[k]) / cell);
        };
        for (int f = 0; f < T.m; ++f) {
            int x0, x1, y0, y1; range(f, 0, x0, x1); range(f, 1, y0, y1);
            for (int cx = x0; cx <= x1; ++cx) for (int cy = y0; cy <= y1; ++cy) cnt[(size_t)cx * nc[1] + cy + 1]++;
        }
        ptr.assign(cnt.size(), 0);
        for (size_t i = 1; i < cnt.size(); ++i) ptr[i] = ptr[i - 1] + cnt[i];
        faces.resize(ptr.back());
        std::vector<int64_t> fill(ptr.begin(), ptr.end() - 1);
        for (int f = 0; f < T.m; ++f) {
            int x0, x1, y0, y1; range(f, 0, x0, x1); range(f, 1, y0, y1);
            for (int cx = x0; cx <= x1; ++cx) for (int cy = y0; cy <= y1; ++cy) faces[fill[(size_t)cx * nc[1] + cy]++] = f;
        }
        stamp.assign(T.m, 0);
    }
    // faces whose bounding box meets [x0, x1] x [y0, y1], ascending
    std::vector<int> in_box(double x0, double x1, double y0, double y1) const {
        auto clampi = [](long v, int lo_, int hi_) { return (int)std::max<long>(lo_, std::min<long>(hi_, v)); };
        int a0 = clampi((long)std::floor((x0 - org[0]) / cell), 0, nc[0] - 1), a1 = clampi((long)std::floor((x1 - org[0]) / cell), 0, nc[0] - 1);
        int b0 = clampi((long)std::floor((y0 - org[1]) / cell), 0, nc[1] - 1), b1 = clampi((long)std::floor((y1 - org[1]) / cell), 0, nc[1] - 1);
        if (++cur == std::numeric_limits<int>::max()) { std::fill(stamp.begin(), stamp.end(), 0); cur = 1; }
        std::vector<int> out;
        for (int cx = a0; cx <= a1; ++cx) {
            int64_t s = ptr[(size_t)cx * nc[1] + b0], e = ptr[(size_t)cx * nc[1] + b1 + 1];
            for (int64_t i = s; i < e; ++i) {
                int f = faces[i];
                if (stamp[f] == cur) continue;
                stamp[f] = cur;
                if (flo[2 * f] <= x1 && fhi[2 * f] >= x0 && flo[2 * f + 1] <= y1 && fhi[2 * f + 1] >= y0) out.push_back(f);
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }
};

// tolerant point location: the face of maximum min-barycentric among the candidates of the cell,
// accepted if that minimum is >= -tol
inline int locate(const Terrain& T, const FaceGrid& G, double x, double y, double tol = 1e-7) {
    double r = 1e-9 * std::max(1.0, T.edge);
    std::vector<int> cand = G.in_box(x - r, x + r, y - r, y + r);
    int best = -1; double bm = -INF, l[3];
    for (int f : cand) {
        bary(T, f, x, y, l);
        double mn = std::min({l[0], l[1], l[2]});
        if (mn > bm) { bm = mn; best = f; }
    }
    return bm >= -tol ? best : -1;
}

// vertex grid: radius queries and nearest vertex
struct PointGrid {
    double cell = 1, org[2] = {0, 0};
    int nc[2] = {1, 1};
    std::vector<int64_t> ptr;
    std::vector<int> ids;
    const std::vector<double>* P = nullptr;
    void init(const std::vector<double>& xy, double cell_) {
        P = &xy;
        size_t n = xy.size() / 2;
        double lo_[2] = {INF, INF}, hi_[2] = {-INF, -INF};
        for (size_t i = 0; i < n; ++i) for (int k = 0; k < 2; ++k) { lo_[k] = std::min(lo_[k], xy[2 * i + k]); hi_[k] = std::max(hi_[k], xy[2 * i + k]); }
        cell = cell_;
        for (int k = 0; k < 2; ++k) { org[k] = lo_[k]; nc[k] = (int)std::floor((hi_[k] - lo_[k]) / cell) + 1; }
        std::vector<int64_t> cnt((size_t)nc[0] * nc[1] + 1, 0);
        for (size_t i = 0; i < n; ++i) cnt[cell_of(xy[2 * i], xy[2 * i + 1]) + 1]++;
        ptr.assign(cnt.size(), 0);
        for (size_t i = 1; i < cnt.size(); ++i) ptr[i] = ptr[i - 1] + cnt[i];
        ids.resize(n);
        std::vector<int64_t> fill(ptr.begin(), ptr.end() - 1);
        for (size_t i = 0; i < n; ++i) ids[fill[cell_of(xy[2 * i], xy[2 * i + 1])]++] = (int)i;
    }
    size_t cell_of(double x, double y) const {
        int cx = std::max(0, std::min(nc[0] - 1, (int)std::floor((x - org[0]) / cell)));
        int cy = std::max(0, std::min(nc[1] - 1, (int)std::floor((y - org[1]) / cell)));
        return (size_t)cx * nc[1] + cy;
    }
    // ids with |p - q| <= r (ascending)
    std::vector<int> radius(double x, double y, double r) const {
        std::vector<int> out;
        int x0 = std::max(0, (int)std::floor((x - r - org[0]) / cell)), x1 = std::min(nc[0] - 1, (int)std::floor((x + r - org[0]) / cell));
        int y0 = std::max(0, (int)std::floor((y - r - org[1]) / cell)), y1 = std::min(nc[1] - 1, (int)std::floor((y + r - org[1]) / cell));
        double r2 = r * r;
        for (int cx = x0; cx <= x1; ++cx)
            for (int cy = y0; cy <= y1; ++cy)
                for (int64_t i = ptr[(size_t)cx * nc[1] + cy]; i < ptr[(size_t)cx * nc[1] + cy + 1]; ++i) {
                    int id = ids[i];
                    double dx = (*P)[2 * id] - x, dy = (*P)[2 * id + 1] - y;
                    if (dx * dx + dy * dy <= r2) out.push_back(id);
                }
        std::sort(out.begin(), out.end());
        return out;
    }
    int nearest(double x, double y) const {
        int best = -1; double bd = INF;
        int cx = (int)std::floor((x - org[0]) / cell), cy = (int)std::floor((y - org[1]) / cell);
        for (int ring = 0; ring < std::max(nc[0], nc[1]) + 2; ++ring) {
            for (int i = cx - ring; i <= cx + ring; ++i)
                for (int j = cy - ring; j <= cy + ring; ++j) {
                    if (std::max(std::abs(i - cx), std::abs(j - cy)) != ring) continue;
                    if (i < 0 || j < 0 || i >= nc[0] || j >= nc[1]) continue;
                    for (int64_t k = ptr[(size_t)i * nc[1] + j]; k < ptr[(size_t)i * nc[1] + j + 1]; ++k) {
                        int id = ids[k];
                        double dx = (*P)[2 * id] - x, dy = (*P)[2 * id + 1] - y, d = dx * dx + dy * dy;
                        if (d < bd) { bd = d; best = id; }
                    }
                }
            if (best >= 0 && std::sqrt(bd) <= ring * cell) break; //every unvisited cell is farther than ring*cell
        }
        return best;
    }
};

// exact lifted (weighted) length of the chart segment A -> B, walking the faces it crosses;
// INF if it leaves the mesh. Endpoints at vertices are handled by starting a hair inside.
inline double walk_length(const Terrain& T, const FaceGrid& G, const double* A, const double* B, int max_hops = 4096) {
    double d[2] = {B[0] - A[0], B[1] - A[1]};
    double L = std::hypot(d[0], d[1]);
    if (L == 0) return 0;
    double u[2] = {d[0] / L, d[1] / L};
    const double shrink = 1e-7;
    double s = shrink, s_end = 1.0 - shrink; // interior part, the two hair-ends are added back below
    double P[2] = {A[0] + s * d[0], A[1] + s * d[1]};
    int f = locate(T, G, P[0], P[1], 1e-7);
    if (f < 0) return INF;
    double total = 0;
    for (int hop = 0; hop < max_hops; ++hop) {
        // barycentric coordinates along the segment are affine in the parameter
        double l0[3], l1[3];
        bary(T, f, A[0] + s * d[0], A[1] + s * d[1], l0);
        bary(T, f, A[0] + (s + 1.0) * d[0], A[1] + (s + 1.0) * d[1], l1);
        double s_exit = INF; int k_exit = -1;
        for (int i = 0; i < 3; ++i) {
            double dl = l1[i] - l0[i];
            if (dl < 0) {
                double t = s - l0[i] / dl;
                if (t < s_exit) { s_exit = t; k_exit = i; }
            }
        }
        double speed = T.w[f] * std::sqrt(1.0 + std::pow(T.a[f] * u[0] + T.b[f] * u[1], 2));
        if (s_exit >= s_end - 1e-12) {
            total += (s_end - s) * L * speed;
            // hair-ends: the first and last 1e-7 of the segment at the speeds of the end faces
            int fa = locate(T, G, A[0] + shrink * d[0], A[1] + shrink * d[1], 1e-7);
            double sa = fa >= 0 ? T.w[fa] * std::sqrt(1.0 + std::pow(T.a[fa] * u[0] + T.b[fa] * u[1], 2)) : speed;
            return total + shrink * L * (sa + speed);
        }
        if (s_exit < s) s_exit = s;
        total += (s_exit - s) * L * speed;
        s = s_exit;
        // vertex i reached zero weight: leave through the edge opposite vertex i (vertices k_exit+1, k_exit+2)
        int edge_k = (k_exit + 1) % 3; // edge from vertex k_exit+1 to k_exit+2
        int g = T.nbr[3 * f + edge_k];
        if (g < 0) return INF;
        f = g;
    }
    return INF;
}

}  // namespace cci
