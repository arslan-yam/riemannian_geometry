#include "certified.h"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>

namespace cci {

using clk = std::chrono::steady_clock;
static double secs(clk::time_point a) { return std::chrono::duration<double>(clk::now() - a).count(); }

// body(tid) on nt threads (inline when nt == 1)
template <class F> static void run_threads(int nt, F&& body) {
    if (nt <= 1) { body(0); return; }
    std::vector<std::thread> th;
    for (int t = 0; t < nt; ++t) th.emplace_back([&body, t]() { body(t); });
    for (auto& x : th) x.join();
}
// body(i, tid) for i < n, handed out one at a time (dynamic scheduling); in order when nt == 1
template <class F> static void parallel_for(int nt, size_t n, F&& body) {
    std::atomic<size_t> next{0};
    run_threads((int)std::min<size_t>((size_t)std::max(nt, 1), std::max<size_t>(n, 1)), [&](int tid) {
        for (size_t i; (i = next++) < n;) body(i, tid);
    });
}

// numpy.interp for increasing xp
static double interp1(double x, const std::vector<double>& xp, const std::vector<double>& fp) {
    size_t n = xp.size();
    if (x <= xp[0]) return fp[0];
    if (x >= xp[n - 1]) return fp[n - 1];
    size_t j = std::upper_bound(xp.begin(), xp.end(), x) - xp.begin();
    double t = (x - xp[j - 1]) / (xp[j] - xp[j - 1]);
    return fp[j - 1] + t * (fp[j] - fp[j - 1]);
}

CertifiedIndex::CertifiedIndex(const Terrain& t, const Options& o) : T(t), opt(o) {
    auto t0 = clk::now();
    eps = opt.eps;
    lam = T.lam;
    g = opt.gap > 0 ? opt.gap : opt.gap_factor * T.edge;
    Rn = g / (eps * lam);
    Rq = Rn * (1.0 + 1e-9);
    slot_of.assign(T.n, -1);
    if (opt.poi.empty()) {
        for (int i = 0; i < T.n; ++i) slot_of[i] = i;
        NP = T.n;
    } else {
        if (opt.hierarchy == "shifted") throw std::invalid_argument("POI mode needs the kd hierarchy");
        poi_list = opt.poi;
        std::sort(poi_list.begin(), poi_list.end());
        poi_list.erase(std::unique(poi_list.begin(), poi_list.end()), poi_list.end());
        if (poi_list.front() < 0 || poi_list.back() >= T.n) throw std::invalid_argument("POI out of range");
        for (size_t i = 0; i < poi_list.size(); ++i) slot_of[poi_list[i]] = (int32_t)i;
        NP = (int)poi_list.size();
    }
    if (opt.store != "dense" && opt.store != "code16" && opt.store != "code8") throw std::invalid_argument("store: dense | code16 | code8");
    grid.init(T);
    pgrid.init(T.xy, std::max(Rq, T.edge));
    G.from_mesh(T);
    for (int k = 0; k < 2; ++k) {
        ucoord[k].resize(T.n);
        for (int i = 0; i < T.n; ++i) ucoord[k][i] = T.xy[2 * i + k];
        std::sort(ucoord[k].begin(), ucoord[k].end());
        ucoord[k].erase(std::unique(ucoord[k].begin(), ucoord[k].end()), ucoord[k].end());
    }
    timing["setup"] = secs(t0);
    t0 = clk::now(); if (opt.hierarchy == "shifted") build_shifted(); else build_tree(); timing["tree"] = secs(t0);
    t0 = clk::now(); cut_geometry(); timing["cut_geometry"] = secs(t0);
    t0 = clk::now(); labels_phase(); timing["labels"] = secs(t0);
    t0 = clk::now(); near_phase(); timing["near_tables"] = secs(t0);
    t0 = clk::now(); feet_phase(); timing["feet"] = secs(t0);
    if (opt.store == "code16") storage = CODE16;
    if (opt.store == "code8") storage = CODE8;
}

// ---------------------------------------------------------------- separator tree
double CertifiedIndex::split_value(std::vector<double> cs, int axis) const {
    std::sort(cs.begin(), cs.end());
    size_t k = cs.size() / 2;
    double lo = cs[k - 1], hi = cs[k];
    if (lo == hi) { // ties at the median: the nearest distinct split
        long best = -1, bd = -1;
        for (size_t j = 0; j < cs.size(); ++j)
            if (cs[j] != cs[k]) {
                long dj = std::labs((long)j - (long)k);
                if (best < 0 || dj < bd) { best = (long)j; bd = dj; }
            }
        if (best < 0) return NAN;
        if ((size_t)best < k) { lo = cs[best]; hi = cs[best + 1]; }
        else { lo = cs[best - 1]; hi = cs[best]; }
    }
    const auto& U = ucoord[axis];
    std::vector<double> pts{lo};
    for (auto it = std::upper_bound(U.begin(), U.end(), lo); it != U.end() && *it < hi; ++it) pts.push_back(*it);
    pts.push_back(hi);
    size_t i = 0; double gbest = -1;
    for (size_t j = 0; j + 1 < pts.size(); ++j)
        if (pts[j + 1] - pts[j] > gbest) { gbest = pts[j + 1] - pts[j]; i = j; }
    return 0.5 * (pts[i] + pts[i + 1]);
}

void CertifiedIndex::build_tree() {
    const int n = NP, maxd = 64; // rows: indexed vertices (all of them, or the POIs)
    std::vector<int32_t> ancF((size_t)n * maxd, -1), rowF((size_t)n * maxd, -1);
    depth_of.assign(n, 0);
    struct Item { std::vector<int> M; int depth; };
    std::vector<Item> stack;
    std::vector<int> all(poi_list);
    if (all.empty()) { all.resize(T.n); std::iota(all.begin(), all.end(), 0); }
    stack.push_back({all, 0});
    while (!stack.empty()) {
        Item it = std::move(stack.back());
        stack.pop_back();
        const auto& M = it.M;
        int depth = it.depth;
        double lo[2] = {INF, INF}, hi[2] = {-INF, -INF};
        for (int u : M) for (int k = 0; k < 2; ++k) { lo[k] = std::min(lo[k], T.xy[2 * u + k]); hi[k] = std::max(hi[k], T.xy[2 * u + k]); }
        double v = NAN; int a = 0;
        if (M.size() > 1 && std::hypot(hi[0] - lo[0], hi[1] - lo[1]) > Rn && depth < maxd) {
            a = (hi[1] - lo[1] > hi[0] - lo[0]) ? 1 : 0;
            std::vector<double> c(M.size());
            for (size_t i = 0; i < M.size(); ++i) c[i] = T.xy[2 * M[i] + a];
            v = split_value(c, a);
            if (std::isnan(v)) {
                a = 1 - a;
                for (size_t i = 0; i < M.size(); ++i) c[i] = T.xy[2 * M[i] + a];
                v = split_value(c, a);
            }
        }
        if (std::isnan(v)) {
            for (int u : M) depth_of[slot_of[u]] = depth;
            leaves.push_back(M);
            continue;
        }
        int cid = (int)cuts.size();
        Cut C;
        C.axis = a; C.v = v; C.depth = depth; C.M = M;
        for (int k = 0; k < 2; ++k) { C.lo[k] = lo[k]; C.hi[k] = hi[k]; }
        std::vector<int> Mg, Ml;
        for (size_t r = 0; r < M.size(); ++r) {
            int u = M[r];
            ancF[(size_t)slot_of[u] * maxd + depth] = cid; rowF[(size_t)slot_of[u] * maxd + depth] = (int32_t)r;
            double x = T.xy[2 * u + a];
            if (x > v) Mg.push_back(u); else if (x < v) Ml.push_back(u);
        }
        cuts.push_back(std::move(C));
        stack.push_back({Mg, depth + 1});
        stack.push_back({Ml, depth + 1});
    }
    D = *std::max_element(depth_of.begin(), depth_of.end()) + 1;
    anc.assign((size_t)n * D, -1); row.assign((size_t)n * D, -1);
    for (int u = 0; u < n; ++u)
        for (int k = 0; k < D; ++k) { anc[(size_t)u * D + k] = ancF[(size_t)u * maxd + k]; row[(size_t)u * D + k] = rowF[(size_t)u * maxd + k]; }
}

// ---------------------------------------------------------------- shifted hierarchy
// copies regular quadtrees, copy j shifted by j/3 of the top cell (2-d shifting lemma: any two points at
// L-inf distance r share a cell of side <= 6r in one copy). A cell is a leaf when its diagonal is <= R
// or it holds one vertex. Each cell gets its vertical / horizontal split line as a cut, but only if it
// has members on both sides of it. The grid is offset by a tiny irrational amount so no vertex lies on
// any line.
void CertifiedIndex::build_shifted() {
    const int n = T.n, maxd = 48, ncopy = std::max(1, opt.copies);
    double E = std::max(T.hi[0] - T.lo[0], T.hi[1] - T.lo[1]);
    double S0 = E * ncopy * (1 + 1e-9) + 1e-6;   // copy j covers [o_j, o_j + S0] with o_j = lo - j S0 / ncopy
    std::vector<int> all(n);
    std::iota(all.begin(), all.end(), 0);
    copies.assign(ncopy, ShiftCopy());
    for (int j = 0; j < ncopy; ++j) {
        ShiftCopy& Cp = copies[j];
        // origin with a tiny irrational offset; retried until no vertex is within 1e-9 edges of a line
        double off = 0.6180339887 * 1e-4 * T.edge * (j + 1);
        double org[2];
        for (int attempt = 0; attempt < 20; ++attempt, off *= 1.37) {
            for (int k = 0; k < 2; ++k) org[k] = T.lo[k] - j * S0 / ncopy - off;
            bool clash = false;
            for (int lev = 0; lev < maxd && !clash; ++lev) {
                double Sl = S0 / std::ldexp(1.0, lev);
                if (Sl < 0.5 * Rn) break;
                for (int i = 0; i < n && !clash; ++i)
                    for (int k = 0; k < 2; ++k) {
                        double f = (T.xy[2 * i + k] - org[k]) / Sl;
                        if (std::fabs(f - std::round(f)) * Sl < 1e-9 * std::max(1.0, T.edge)) { clash = true; break; }
                    }
            }
            if (!clash) break;
        }
        std::vector<int32_t> ancF((size_t)n * maxd, -1), rowF((size_t)n * maxd, -1);
        Cp.depth_of.assign(n, 0);
        struct Item { std::vector<int> M; int depth; double x0, y0, S; };
        std::vector<Item> stack{{all, 0, org[0], org[1], S0}};
        int maxdepth = 0;
        while (!stack.empty()) {
            Item it = std::move(stack.back());
            stack.pop_back();
            const auto& M = it.M;
            int depth = it.depth;
            int node = (int)Cp.nodes.size();
            ShiftNode N;
            N.side = it.S; N.mid[0] = it.x0 + it.S / 2; N.mid[1] = it.y0 + it.S / 2;
            for (size_t r = 0; r < M.size(); ++r) { ancF[(size_t)M[r] * maxd + depth] = node; rowF[(size_t)M[r] * maxd + depth] = (int32_t)r; }
            for (int u : M) Cp.depth_of[u] = depth;
            maxdepth = std::max(maxdepth, depth);
            bool leaf = M.size() <= 1 || it.S * std::sqrt(2.0) <= Rn || depth + 1 >= maxd;
            if (!leaf) {
                double lo[2] = {INF, INF}, hi[2] = {-INF, -INF};
                for (int u : M) for (int k = 0; k < 2; ++k) { lo[k] = std::min(lo[k], T.xy[2 * u + k]); hi[k] = std::max(hi[k], T.xy[2 * u + k]); }
                for (int a = 0; a < 2; ++a) {
                    if (!(lo[a] < N.mid[a] && N.mid[a] < hi[a])) continue; // nobody on one side: no pair split here
                    Cut C;
                    C.axis = a; C.v = N.mid[a]; C.depth = depth; C.M = M; C.cell_side = it.S;
                    for (int k = 0; k < 2; ++k) { C.lo[k] = lo[k]; C.hi[k] = hi[k]; }
                    N.cut[a] = (int)cuts.size();
                    cuts.push_back(std::move(C));
                }
                std::vector<int> Q[4];
                for (int u : M) Q[(T.xy[2 * u] > N.mid[0] ? 1 : 0) + (T.xy[2 * u + 1] > N.mid[1] ? 2 : 0)].push_back(u);
                double h = it.S / 2;
                for (int q = 3; q >= 0; --q)
                    if (!Q[q].empty()) stack.push_back({Q[q], depth + 1, it.x0 + (q & 1 ? h : 0), it.y0 + (q & 2 ? h : 0), h});
            } else if (j == 0) {
                leaves.push_back(M);
            }
            Cp.nodes.push_back(N);
        }
        Cp.D = maxdepth + 1;
        Cp.anc.assign((size_t)n * Cp.D, -1); Cp.row.assign((size_t)n * Cp.D, -1);
        for (int u = 0; u < n; ++u)
            for (int k = 0; k < Cp.D; ++k) { Cp.anc[(size_t)u * Cp.D + k] = ancF[(size_t)u * maxd + k]; Cp.row[(size_t)u * Cp.D + k] = rowF[(size_t)u * maxd + k]; }
    }
    depth_of = copies[0].depth_of;
}

// ---------------------------------------------------------------- cut lines
std::vector<Run> CertifiedIndex::line_pieces(int a, double v, double lo, double hi) const {
    int b = 1 - a;
    std::vector<int> faces = a == 0 ? grid.in_box(v, v, lo, hi) : grid.in_box(lo, hi, v, v);
    struct Piece { double u0, u1, speed; };
    std::vector<Piece> P;
    for (int f : faces) {
        double Q[3][2], da[3];
        for (int i = 0; i < 3; ++i) { Q[i][0] = T.xy[2 * T.F[3 * f + i]]; Q[i][1] = T.xy[2 * T.F[3 * f + i] + 1]; da[i] = Q[i][a] - v; }
        double us[3]; int cnt = 0;
        const int E[3][2] = {{0, 1}, {1, 2}, {2, 0}};
        for (int e = 0; e < 3; ++e) {
            int i = E[e][0], j = E[e][1];
            if (da[i] * da[j] < 0) {
                double t = da[i] / (da[i] - da[j]);
                us[cnt++] = Q[i][b] + t * (Q[j][b] - Q[i][b]);
            }
        }
        if (cnt != 2) continue; // no vertex lies on the line: a crossed face is crossed through two edges
        double u0 = std::min(us[0], us[1]), u1 = std::max(us[0], us[1]);
        double u0c = std::max(u0, lo), u1c = std::min(u1, hi);
        if (!(u1c > u0c)) continue;
        double grad = a == 0 ? T.b[f] : T.a[f];
        P.push_back({u0c, u1c, T.w[f] * std::sqrt(1.0 + grad * grad)});
    }
    std::vector<Run> runs;
    if (P.empty()) return runs;
    std::stable_sort(P.begin(), P.end(), [](const Piece& x, const Piece& y) { return x.u0 < y.u0; });
    double tol = 1e-9 * std::max(1.0, T.edge);
    size_t start = 0;
    for (size_t i = 1; i <= P.size(); ++i) {
        if (i == P.size() || P[i].u0 > P[i - 1].u1 + tol) {
            Run r;
            r.ub.push_back(P[start].u0); r.ab.push_back(0.0);
            double acc = 0;
            for (size_t j = start; j < i; ++j) {
                acc += (P[j].u1 - P[j].u0) * P[j].speed;
                r.ub.push_back(j + 1 < i ? P[j + 1].u0 : P[j].u1);
                r.ab.push_back(acc);
                r.smax = std::max(r.smax, P[j].speed);
            }
            runs.push_back(std::move(r));
            start = i;
        }
    }
    return runs;
}

std::vector<double> CertifiedIndex::place(const Run& R, double r_lo, double r_hi, double g) const {
    const auto &ub = R.ub, &ab = R.ab;
    double U0 = ub.front(), U1 = ub.back();
    auto arc = [&](double u) { return interp1(u, ub, ab); };
    auto inv = [&](double s) { return interp1(s, ab, ub); };
    std::vector<double> pos{U0, U1};
    double k2 = 2.0 * eps * lam;
    double rl = std::max(U0, r_lo), rh = std::min(U1, r_hi);
    if (rh > rl) {
        double A0 = arc(rl), A1 = arc(rh);
        int m = std::max(1, (int)std::ceil((A1 - A0) / g));
        for (int i = 0; i <= m; ++i) pos.push_back(inv(A0 + (A1 - A0) * i / m));
    }
    if (U0 < r_lo) {
        double cur = std::min(U1, r_lo);
        pos.push_back(cur);
        while (cur > U0) {
            double c = r_lo - cur;
            double nxt = arc(cur) - std::max(g, k2 * c);
            if (nxt <= 0.0) break;
            cur = inv(nxt);
            pos.push_back(cur);
        }
    }
    if (U1 > r_hi) {
        double cur = std::max(U0, r_hi);
        pos.push_back(cur);
        while (cur < U1) {
            double c = cur - r_hi;
            double nxt = arc(cur) + std::max(g, k2 * c);
            if (nxt >= ab.back()) break;
            cur = inv(nxt);
            pos.push_back(cur);
        }
    }
    for (auto& p : pos) p = std::min(std::max(p, U0), U1);
    std::sort(pos.begin(), pos.end());
    pos.erase(std::unique(pos.begin(), pos.end()), pos.end());
    std::vector<double> out;
    double thr = 1e-9 * std::max(1.0, T.edge);
    for (size_t i = 0; i < pos.size(); ++i)
        if (i == 0 || pos[i] - pos[i - 1] > thr) out.push_back(pos[i]);
    return out;
}

void CertifiedIndex::cut_geometry() {
    const int nt = std::max(1, opt.threads);
    std::vector<Dijkstra> dijs(nt, Dijkstra(G));
    parallel_for(nt, cuts.size(), [&](size_t cid, int tid) {
        Cut& C = cuts[cid];
        Dijkstra& dij = dijs[tid];
        int a = C.axis, b = 1 - a;
        double r_lo = C.lo[b], r_hi = C.hi[b];
        double ctr[2] = {0.5 * (C.lo[0] + C.hi[0]), 0.5 * (C.lo[1] + C.hi[1])};
        int hub = C.M[0]; double bd = INF;
        for (int u : C.M) {
            double d = std::hypot(T.xy[2 * u] - ctr[0], T.xy[2 * u + 1] - ctr[1]);
            if (d < bd) { bd = d; hub = u; }
        }
        double lim = opt.hub_slack * std::hypot(C.hi[0] - C.lo[0], C.hi[1] - C.lo[1]) * 2.0 + 10 * T.edge;
        dij.run(hub, lim);
        bool fin = true;
        for (int u : C.M) if (!std::isfinite(dij.dist[u])) { fin = false; break; }
        if (!fin) dij.run(hub);
        double mu = 0;
        for (int u : C.M) mu = std::max(mu, dij.dist[u]);
        double off = mu / lam;
        double wlo = std::max(r_lo - off, T.lo[b]), whi = std::min(r_hi + off, T.hi[b]);
        C.runs = line_pieces(a, C.v, wlo, whi);
        // the gap inside the box: fixed g; capped at k_max portals (version A: the arc of the box part over
        // k_max); or scale-dependent max(g, eps lam S / 6) on the shifted hierarchy (set with cell_side)
        double in_box = 0;
        for (auto& R : C.runs) {
            double rl = std::max(R.ub.front(), r_lo), rh = std::min(R.ub.back(), r_hi);
            if (rh > rl) in_box += interp1(rh, R.ub, R.ab) - interp1(rl, R.ub, R.ab);
        }
        C.gap_box = g;
        if (opt.k_max > 0) C.gap_box = std::max(g, in_box / opt.k_max);
        if (C.cell_side > 0) C.gap_box = std::max(g, eps * lam * C.cell_side / 6.0);
        C.pos.clear(); C.arc.clear(); C.run.clear();
        for (size_t r = 0; r < C.runs.size(); ++r) {
            auto p = place(C.runs[r], r_lo, r_hi, C.gap_box);
            for (double x : p) { C.pos.push_back(x); C.arc.push_back(interp1(x, C.runs[r].ub, C.runs[r].ab)); C.run.push_back((int64_t)r); }
        }
        int K = C.K();
        C.X.assign(2 * K, 0);
        for (int j = 0; j < K; ++j) { C.X[2 * j + a] = C.v; C.X[2 * j + b] = C.pos[j]; }
        C.mu = mu; C.hub = hub;
    });
}

// ---------------------------------------------------------------- labels
void CertifiedIndex::submesh(const std::vector<int>& faces, const std::vector<double>& X, const double* along,
                             std::vector<double>& V2, std::vector<unsigned>& F2, std::vector<int>& used,
                             std::vector<int>& pidx, std::vector<double>& Xu) const {
    // local vertices
    used.clear();
    for (int f : faces) for (int k = 0; k < 3; ++k) used.push_back((int)T.F[3 * f + k]);
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());
    auto loc = [&](int v) { return (int)(std::lower_bound(used.begin(), used.end(), v) - used.begin()); };
    V2.resize(3 * used.size());
    for (size_t i = 0; i < used.size(); ++i) for (int k = 0; k < 3; ++k) V2[3 * i + k] = T.V[3 * used[i] + k];
    int K = (int)X.size() / 2;
    Xu = X;
    pidx.assign(K, -1);
    std::vector<int> gf(K, -1), lf(K, -1);
    for (int j = 0; j < K; ++j) {
        gf[j] = locate(T, grid, X[2 * j], X[2 * j + 1]);
        if (gf[j] >= 0) {
            auto it = std::lower_bound(faces.begin(), faces.end(), gf[j]);
            if (it != faces.end() && *it == gf[j]) lf[j] = (int)(it - faces.begin());
        }
    }
    // keep points off the edges of their face: move along the cut by 1e-5 edge lengths
    for (int j = 0; j < K; ++j) {
        if (lf[j] < 0) continue;
        int f = gf[j];
        const double *A = &T.xy[2 * T.F[3 * f]], *B = &T.xy[2 * T.F[3 * f + 1]];
        double scale = std::hypot(B[0] - A[0], B[1] - A[1]);
        for (int it = 0; it < 8; ++it) {
            double l[3]; bary(T, f, Xu[2 * j], Xu[2 * j + 1], l);
            if (std::min({l[0], l[1], l[2]}) >= 1e-7) break;
            double p[2] = {Xu[2 * j] + 1e-5 * scale * along[0], Xu[2 * j + 1] + 1e-5 * scale * along[1]};
            double q[2] = {Xu[2 * j] - 1e-5 * scale * along[0], Xu[2 * j + 1] - 1e-5 * scale * along[1]};
            double lp[3], lq[3]; bary(T, f, p[0], p[1], lp); bary(T, f, q[0], q[1], lq);
            if (std::min({lp[0], lp[1], lp[2]}) >= std::min({lq[0], lq[1], lq[2]})) { Xu[2 * j] = p[0]; Xu[2 * j + 1] = p[1]; }
            else { Xu[2 * j] = q[0]; Xu[2 * j + 1] = q[1]; }
        }
    }
    // new vertices on their face planes, then each face with points is split point by point (1 -> 3)
    std::vector<std::vector<int>> per_face(faces.size());
    int nV0 = (int)used.size();
    for (int j = 0; j < K; ++j) {
        if (lf[j] < 0) continue;
        pidx[j] = (int)(V2.size() / 3);
        V2.push_back(Xu[2 * j]); V2.push_back(Xu[2 * j + 1]); V2.push_back(T.z_at(gf[j], Xu[2 * j], Xu[2 * j + 1]));
        per_face[lf[j]].push_back(j);
    }
    (void)nV0;
    F2.clear();
    F2.reserve(3 * (faces.size() + 2 * K));
    auto lbary = [&](unsigned ia, unsigned ib, unsigned ic, double x, double y, double* l) {
        const double *A = &V2[3 * ia], *B = &V2[3 * ib], *C = &V2[3 * ic];
        double det = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);
        double l1 = ((x - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (y - A[1])) / det;
        double l2 = ((B[0] - A[0]) * (y - A[1]) - (x - A[0]) * (B[1] - A[1])) / det;
        l[0] = 1 - l1 - l2; l[1] = l1; l[2] = l2;
    };
    for (size_t fi = 0; fi < faces.size(); ++fi) {
        int f = faces[fi];
        unsigned c0 = loc(T.F[3 * f]), c1 = loc(T.F[3 * f + 1]), c2 = loc(T.F[3 * f + 2]);
        if (per_face[fi].empty()) { F2.push_back(c0); F2.push_back(c1); F2.push_back(c2); continue; }
        std::vector<unsigned> tri{c0, c1, c2};
        for (int j : per_face[fi]) {
            unsigned p = (unsigned)pidx[j];
            size_t best = 0; double bm = -INF;
            for (size_t t = 0; t < tri.size(); t += 3) {
                double l[3]; lbary(tri[t], tri[t + 1], tri[t + 2], V2[3 * p], V2[3 * p + 1], l);
                double mn = std::min({l[0], l[1], l[2]});
                if (mn > bm) { bm = mn; best = t; }
            }
            unsigned A = tri[best], B = tri[best + 1], C = tri[best + 2];
            tri[best + 2] = p;                       // (A, B, p)
            tri.insert(tri.end(), {B, C, p, C, A, p}); // (B, C, p), (C, A, p)
        }
        F2.insert(F2.end(), tri.begin(), tri.end());
    }
    // one orientation for every triangle (counter-clockwise in the chart): VTP's window lists need
    // consistently oriented faces and abort otherwise; distances do not depend on it
    for (size_t t = 0; t < F2.size(); t += 3) {
        const double *A = &V2[3 * F2[t]], *B = &V2[3 * F2[t + 1]], *C = &V2[3 * F2[t + 2]];
        if ((B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]) < 0) std::swap(F2[t + 1], F2[t + 2]);
    }
}

void CertifiedIndex::finalize_cut(Cut& C, std::vector<float>& Lc, std::vector<double>& Xu, int K) {
    int a = C.axis, b = 1 - a;
    size_t nm = C.M.size();
    std::vector<int> keep;
    for (int j = 0; j < K; ++j) {
        bool ok = true;
        for (size_t r = 0; r < nm; ++r) if (!std::isfinite(Lc[r * K + j])) { ok = false; break; }
        if (ok) keep.push_back(j);
    }
    std::vector<std::pair<double, int>> ord;
    for (int j : keep) ord.push_back({Xu[2 * j + b], j});
    std::stable_sort(ord.begin(), ord.end(), [](auto& x, auto& y) { return x.first < y.first; });
    int K2 = (int)ord.size();
    std::vector<float> L2(nm * K2);
    std::vector<double> X2(2 * K2), pos(K2);
    for (int jj = 0; jj < K2; ++jj) {
        int j = ord[jj].second;
        for (size_t r = 0; r < nm; ++r) L2[r * K2 + jj] = Lc[r * K + j];
        X2[2 * jj] = Xu[2 * j]; X2[2 * jj + 1] = Xu[2 * j + 1]; pos[jj] = ord[jj].first;
    }
    std::vector<int64_t> rid(K2, 0);
    for (int j = 0; j < K2; ++j) {
        int r = 0;
        for (size_t q = 0; q < C.runs.size(); ++q) if (C.runs[q].ub[0] <= pos[j]) r = (int)q;
        rid[j] = r;
    }
    std::vector<double> arcs(K2, 0.0);
    C.end_idx.clear(); C.end_gap.clear(); C.lb_ok = true;
    for (size_t r = 0; r < C.runs.size(); ++r) {
        std::vector<int> m;
        for (int j = 0; j < K2; ++j) if (rid[j] == (int64_t)r) m.push_back(j);
        if (m.empty()) { C.lb_ok = false; continue; }
        for (int j : m) arcs[j] = interp1(pos[j], C.runs[r].ub, C.runs[r].ab);
        double lead = arcs[m.front()] - C.runs[r].ab.front(), trail = C.runs[r].ab.back() - arcs[m.back()];
        if (lead > 0) { C.end_idx.push_back(m.front()); C.end_gap.push_back(2.0 * lead); }
        if (trail > 0) { C.end_idx.push_back(m.back()); C.end_gap.push_back(2.0 * trail); }
    }
    C.gap.assign(K2 > 1 ? K2 - 1 : 0, INF);
    for (int j = 0; j + 1 < K2; ++j) if (rid[j + 1] == rid[j]) C.gap[j] = arcs[j + 1] - arcs[j];
    C.L = std::move(L2); C.X = std::move(X2); C.pos = std::move(pos); C.arc = std::move(arcs); C.run = std::move(rid);
}

// the faces meeting the members' and portals' box dilated by R = (max UB(portal, hub) + mu) / (2 lambda)
// (ellipse lemma), with the portals inserted as vertices; the members' local ids; VTP or MMP
void CertifiedIndex::prepare_cut(const Cut& C, Dijkstra& dij, CutWork& W) const {
    int K = C.K(), a = C.axis;
    size_t nm = C.M.size();
    double lim = opt.hub_slack * std::hypot(C.hi[0] - C.lo[0], C.hi[1] - C.lo[1]) * 2.0 + 10 * T.edge + 2 * C.mu;
    dij.run(C.hub, lim);
    // upper bounds of d(portal, hub): straight walk, or a corner of the portal's face + its graph distance
    auto ubounds = [&](std::vector<double>& U) {
        U.assign(K, INF);
        const double* H = &T.xy[2 * C.hub];
        for (int j = 0; j < K; ++j) {
            const double* x = &C.X[2 * j];
            U[j] = walk_length(T, grid, x, H);
            int f = locate(T, grid, x[0], x[1]);
            if (f >= 0) {
                double z = T.z_at(f, x[0], x[1]);
                for (int k = 0; k < 3; ++k) {
                    int vtx = (int)T.F[3 * f + k];
                    double d = T.w[f] * std::sqrt(std::pow(T.V[3 * vtx] - x[0], 2) + std::pow(T.V[3 * vtx + 1] - x[1], 2) + std::pow(T.V[3 * vtx + 2] - z, 2));
                    U[j] = std::min(U[j], d + dij.dist[vtx]);
                }
            }
        }
    };
    std::vector<double> U;
    ubounds(U);
    if (!std::all_of(U.begin(), U.end(), [](double x) { return std::isfinite(x); })) { dij.run(C.hub); ubounds(U); }
    double R = (*std::max_element(U.begin(), U.end()) + C.mu) / (2.0 * lam);
    double blo[2] = {INF, INF}, bhi[2] = {-INF, -INF};
    for (int u : C.M) for (int k = 0; k < 2; ++k) { blo[k] = std::min(blo[k], T.xy[2 * u + k]); bhi[k] = std::max(bhi[k], T.xy[2 * u + k]); }
    for (int j = 0; j < K; ++j) for (int k = 0; k < 2; ++k) { blo[k] = std::min(blo[k], C.X[2 * j + k]); bhi[k] = std::max(bhi[k], C.X[2 * j + k]); }
    std::vector<int> faces = grid.in_box(blo[0] - R, bhi[0] + R, blo[1] - R, bhi[1] + R);
    double along[2] = {0, 0}; along[1 - a] = 1.0;
    submesh(faces, C.X, along, W.V2, W.F2, W.used, W.pidx, W.Xu);
    W.tgt.resize(nm);
    for (size_t r = 0; r < nm; ++r) W.tgt[r] = (unsigned)(std::lower_bound(W.used.begin(), W.used.end(), C.M[r]) - W.used.begin());
    W.vtp = opt.solver == "vtp";
    if (opt.solver == "hybrid") {
        size_t nbox = grid.in_box(blo[0], bhi[0], blo[1], bhi[1]).size();
        W.vtp = faces.size() <= opt.hybrid_ratio * nbox;
    }
    W.nfaces = faces.size();
}

// VTP on the sub-mesh, or MMP with early stop (also when VTP rejects the sub-mesh)
Solver* CertifiedIndex::cut_solver(const CutWork& W, bool& vtp) const {
    Solver* s = vtp ? make_vtp(W.V2, W.F2) : nullptr;
    if (!s) { vtp = false; s = make_mmp(W.V2, W.F2, true); }
    return s;
}

void CertifiedIndex::finish_cut(Cut& C, std::vector<float>& Lc, CutWork& W, int K) {
    C.submesh_faces = (int)W.nfaces;
    finalize_cut(C, Lc, W.Xu, K);
    int bits = opt.store == "code8" ? 8 : opt.store == "code16" ? 16 : 0;
    if (bits && encode_cut(C, bits)) std::vector<float>().swap(C.L); // coded: the float32 labels go
}

// One sweep per portal on the cut's sub-mesh. With several threads the top cuts (largest sub-meshes)
// are done one at a time with their portals spread over the threads (a solver per thread), the rest
// one cut per thread, largest first. Every sweep is the same as on one thread, so are the labels.
void CertifiedIndex::labels_phase() {
    const int nt = std::max(1, opt.threads);
    std::atomic<int64_t> n_mmp{0}, n_vtp{0};
    std::mutex log_mu;
    std::vector<Dijkstra> dijs(nt, Dijkstra(G));
    auto report = [&](size_t cid, const Cut& C, int K, bool vtp) {
        if (!opt.verbose || cid % 100 != 0) return;
        std::lock_guard<std::mutex> lk(log_mu);
        std::fprintf(stderr, "  cut %zu/%zu: %zu members, %d portals, sub-mesh %d faces, %s\n", cid, cuts.size(), C.M.size(), K, C.submesh_faces, vtp ? "VTP" : "MMP");
    };
    // portals of cut cid swept by `workers` threads (the calling thread alone when workers == 1)
    auto do_cut = [&](size_t cid, int tid, int workers) {
        Cut& C = cuts[cid];
        int K = C.K();
        if (K == 0) { std::vector<float> L0; std::vector<double> X0; finalize_cut(C, L0, X0, 0); return; }
        CutWork W;
        prepare_cut(C, dijs[tid], W);
        size_t nm = C.M.size();
        std::vector<float> Lc(nm * K, std::numeric_limits<float>::infinity());
        std::atomic<int> next{0};
        std::vector<char> vtp_used(workers, 0);
        run_threads(workers, [&](int w) {
            bool vtp = W.vtp;
            std::unique_ptr<Solver> sol(cut_solver(W, vtp));
            vtp_used[w] = vtp;
            std::vector<double> out(nm);
            for (int j; (j = next++) < K;) {
                if (W.pidx[j] < 0) continue;
                sol->dist((unsigned)W.pidx[j], W.tgt, out.data());
                for (size_t r = 0; r < nm; ++r) Lc[r * K + j] = f32_up(out[r]);
                (vtp ? n_vtp : n_mmp)++;
            }
        });
        C.used_vtp = vtp_used[0];
        finish_cut(C, Lc, W, K);
        report(cid, C, K, C.used_vtp);
    };
    if (nt == 1) {
        for (size_t cid = 0; cid < cuts.size(); ++cid) do_cut(cid, 0, 1);
    } else {
        int split_depth = 1;
        while ((1 << split_depth) < nt) ++split_depth;
        std::vector<size_t> rest;
        for (size_t cid = 0; cid < cuts.size(); ++cid) {
            if (cuts[cid].depth <= split_depth) do_cut(cid, 0, nt);
            else rest.push_back(cid);
        }
        std::stable_sort(rest.begin(), rest.end(), [&](size_t x, size_t y) { return cuts[x].depth < cuts[y].depth; });
        parallel_for(nt, rest.size(), [&](size_t i, int tid) { do_cut(rest[i], tid, 1); });
    }
    sweeps_mmp += n_mmp; sweeps_vtp += n_vtp;
}

// ---------------------------------------------------------------- near tables
void CertifiedIndex::near_phase() {
    std::vector<std::vector<int32_t>> ids(NP);
    std::vector<std::vector<float>> dd(NP);
    const int nt = std::max(1, opt.threads);
    std::vector<Dijkstra> dijs(nt, Dijkstra(G));
    parallel_for(nt, leaves.size(), [&](size_t li, int tid) {
        const auto& M = leaves[li];
        Dijkstra& dij = dijs[tid];
        std::vector<double> V2, Xu, out;
        std::vector<unsigned> F2;
        std::vector<int> used, pidx;
        std::vector<std::vector<int>> Tg(M.size());
        std::vector<int> allt(M.begin(), M.end());
        for (size_t i = 0; i < M.size(); ++i) {
            int s = M[i];
            for (int t : pgrid.radius(T.xy[2 * s], T.xy[2 * s + 1], Rq)) if (t != s && slot_of[t] >= 0) Tg[i].push_back(t);
            allt.insert(allt.end(), Tg[i].begin(), Tg[i].end());
        }
        std::sort(allt.begin(), allt.end());
        allt.erase(std::unique(allt.begin(), allt.end()), allt.end());
        double Umax = 0;
        for (size_t i = 0; i < M.size(); ++i) {
            if (Tg[i].empty()) continue;
            dij.run(M[i], opt.hub_slack * Rq + 10 * T.edge);
            bool fin = true;
            for (int t : Tg[i]) if (!std::isfinite(dij.dist[t])) { fin = false; break; }
            if (!fin) dij.run(M[i]);
            for (int t : Tg[i]) Umax = std::max(Umax, dij.dist[t]);
        }
        if (Umax == 0) return; // no near pairs in this leaf
        double R = Umax / (2.0 * lam);
        double blo[2] = {INF, INF}, bhi[2] = {-INF, -INF};
        for (int u : allt) for (int k = 0; k < 2; ++k) { blo[k] = std::min(blo[k], T.xy[2 * u + k]); bhi[k] = std::max(bhi[k], T.xy[2 * u + k]); }
        std::vector<int> faces = grid.in_box(blo[0] - R, bhi[0] + R, blo[1] - R, bhi[1] + R);
        double along[2] = {1, 0};
        submesh(faces, {}, along, V2, F2, used, pidx, Xu);
        std::unique_ptr<Solver> sol(make_mmp(V2, F2, true));
        auto loc = [&](int v) { return (unsigned)(std::lower_bound(used.begin(), used.end(), v) - used.begin()); };
        for (size_t i = 0; i < M.size(); ++i) {
            int s = M[i];
            if (Tg[i].empty()) continue;
            std::vector<unsigned> tl(Tg[i].size());
            for (size_t k = 0; k < tl.size(); ++k) tl[k] = loc(Tg[i][k]);
            out.resize(tl.size());
            sol->dist(loc(s), tl, out.data());
            int32_t r = slot_of[s];
            ids[r].assign(Tg[i].begin(), Tg[i].end());
            dd[r].resize(tl.size());
            for (size_t k = 0; k < tl.size(); ++k) dd[r][k] = f32_up(out[k]);
        }
    });
    near_ptr.assign(NP + 1, 0);
    for (int i = 0; i < NP; ++i) near_ptr[i + 1] = near_ptr[i] + (int64_t)ids[i].size();
    near_ids.resize(near_ptr[NP]); near_d.resize(near_ptr[NP]);
    for (int i = 0; i < NP; ++i) {
        std::copy(ids[i].begin(), ids[i].end(), near_ids.begin() + near_ptr[i]);
        std::copy(dd[i].begin(), dd[i].end(), near_d.begin() + near_ptr[i]);
    }
}

// ---------------------------------------------------------------- feet (model interpolation)
void CertifiedIndex::feet_phase() {
    for (auto& C : cuts) {
        int a = C.axis, b = 1 - a;
        size_t nm = C.M.size();
        C.foot_r.assign(nm, 0); C.foot_arc.assign(nm, NAN); C.foot_run.assign(nm, -1);
        for (size_t i = 0; i < nm; ++i) {
            const double* P = &T.xy[2 * C.M[i]];
            double Q[2] = {P[0], P[1]};
            Q[a] = C.v;
            double r = walk_length(T, grid, P, Q);
            if (!std::isfinite(r)) r = lam * std::fabs(P[a] - C.v);
            C.foot_r[i] = (float)r;
            for (size_t k = 0; k < C.runs.size(); ++k) {
                const auto& ub = C.runs[k].ub;
                if (P[b] >= ub.front() && P[b] <= ub.back()) { C.foot_arc[i] = (float)interp1(P[b], ub, C.runs[k].ab); C.foot_run[i] = (int64_t)k; }
            }
        }
    }
}

// ---------------------------------------------------------------- query
static inline double hyp(double u, double xc, double r) { double d = u - xc; return std::sqrt(r * r + d * d); }

static double model_min(const double* x, const double* y, double xs, double r1, double xt, double r2) {
    double g0 = y[0] - hyp(x[0], xs, r1) - hyp(x[0], xt, r2);
    double g1 = y[1] - hyp(x[1], xs, r1) - hyp(x[1], xt, r2);
    double g2 = y[2] - hyp(x[2], xs, r1) - hyp(x[2], xt, r2);
    double d1 = (g1 - g0) / (x[1] - x[0]), d2 = (g2 - g1) / (x[2] - x[1]);
    double a = (d2 - d1) / (x[2] - x[0]), b = d1 - a * (x[0] + x[1]), c = g0 - a * x[0] * x[0] - b * x[0];
    double lo = x[0], hi = x[2], u, D, DD;
    auto der = [&](double uu) {
        double M1 = hyp(uu, xs, r1), M2 = hyp(uu, xt, r2);
        D = (uu - xs) / M1 + (uu - xt) / M2 + 2 * a * uu + b;
        DD = r1 * r1 / (M1 * M1 * M1) + r2 * r2 / (M2 * M2 * M2) + 2 * a;
    };
    der(lo);
    if (D >= 0) u = lo;
    else {
        der(hi);
        if (D <= 0) u = hi;
        else {
            u = x[1];
            for (int it = 0; it < 40; ++it) {
                der(u);
                if (D < 0) lo = u; else hi = u;
                double un = DD > 0 ? u - D / DD : 0.5 * (lo + hi);
                if (!(lo < un && un < hi)) un = 0.5 * (lo + hi);
                if (std::fabs(un - u) < 1e-9 * std::max(1.0, std::fabs(u))) { u = un; break; }
                u = un;
            }
        }
    }
    double base = hyp(u, xs, r1) + hyp(u, xt, r2);
    double val = std::max(base + a * u * u + b * u + c, base);
    return std::min(y[1], val);
}

double CertifiedIndex::interp(const Cut& C, const double* f, int K, int i, int rs, int rt) const {
    const auto& run = C.run;
    if (!(0 < i && i < K - 1) || !(run[i - 1] == run[i] && run[i] == run[i + 1])) return f[i];
    double x[3] = {C.arc[i - 1], C.arc[i], C.arc[i + 1]}, y[3] = {f[i - 1], f[i], f[i + 1]};
    if (opt.interp == 2 && C.foot_run[rs] == run[i] && C.foot_run[rt] == run[i]) {
        double xs = C.foot_arc[rs], r1 = C.foot_r[rs], xt = C.foot_arc[rt], r2 = C.foot_r[rt];
        if (std::min(r1, r2) < 2.0 * (x[2] - x[0])) return model_min(x, y, xs, r1, xt, r2);
    }
    double d1 = (y[1] - y[0]) / (x[1] - x[0]), d2 = (y[2] - y[1]) / (x[2] - x[1]);
    double a2 = (d2 - d1) / (x[2] - x[0]);
    if (a2 <= 0) return y[1];
    double xv = 0.5 * (x[0] + x[1]) - d1 / (2.0 * a2);
    if (!(x[0] <= xv && xv <= x[2])) return y[1];
    double fv = y[0] + d1 * (xv - x[0]) + a2 * (xv - x[0]) * (xv - x[1]);
    double cap = 0.25 * std::min(y[0] - y[1], y[2] - y[1]);
    return std::max(fv, y[1] - cap);
}

// the cut that answers (s, t): the deepest common ancestor cut (kd), or on the shifted hierarchy the
// copy whose common cell is smallest and that cell's line separating s and t
bool CertifiedIndex::select_cut(int s, int t, int& cid, int& rs, int& rt) const {
    if (copies.empty()) {
        const size_t ss = (size_t)slot_of[s], st = (size_t)slot_of[t];
        const int32_t *as = &anc[ss * D], *at = &anc[st * D];
        int k = -1;
        for (int j = 0; j < D; ++j) if (as[j] != at[j]) { k = j; break; }
        if (k < 0) k = std::min(depth_of[ss], depth_of[st]);
        if (k <= 0) return false;
        cid = as[k - 1]; rs = row[ss * D + k - 1]; rt = row[st * D + k - 1];
        return cid >= 0;
    }
    double best = INF;
    cid = -1;
    for (const auto& Cp : copies) {
        const int32_t *as = &Cp.anc[(size_t)s * Cp.D], *at = &Cp.anc[(size_t)t * Cp.D];
        int k = -1;
        for (int j = 0; j < Cp.D; ++j) if (as[j] != at[j]) { k = j; break; }
        if (k <= 0) continue; // same leaf in this copy (a near pair) or different roots (cannot happen)
        const ShiftNode& N = Cp.nodes[as[k - 1]];
        if (N.side >= best) continue;
        int a = ((T.xy[2 * s] > N.mid[0]) != (T.xy[2 * t] > N.mid[0])) ? 0 : 1;
        if (N.cut[a] < 0) continue;
        best = N.side; cid = N.cut[a];
        rs = Cp.row[(size_t)s * Cp.D + k - 1]; rt = Cp.row[(size_t)t * Cp.D + k - 1];
    }
    return cid >= 0;
}

// LB, UB and estimate from one cut. f = s-row + t-row (dense, decoded 16-bit, or reconstructed from
// covers: upper reconstruction for UB and the estimate, lower reconstruction for the interval LB);
// each LB term subtracts the per-label error bound of the storage
void CertifiedIndex::eval_cut(const Cut& C, int s, int t, int rs, int rt, double* out) const {
    auto chord = [&]() {
        const double *A = &T.V[3 * s], *B = &T.V[3 * t];
        return lam * std::sqrt(std::pow(A[0] - B[0], 2) + std::pow(A[1] - B[1], 2) + std::pow(A[2] - B[2], 2));
    };
    int K = C.K();
    thread_local std::vector<double> f, fl, err, a1, a2;
    f.resize(std::max(K, 1)); fl.resize(std::max(K, 1)); err.assign(std::max(K, 1), 0.0);
    double* flo = f.data(); // lower f for the LB terms (same array unless covers)
    if (storage == DENSE || (storage == CODE16 && C.Lq.empty()) || (storage == CODE8 && C.Lq8.empty()) || (storage == COVER && C.cov_ptr.empty())) {
        const float *Ls = &C.L[(size_t)rs * K], *Lt = &C.L[(size_t)rt * K];
        for (int j = 0; j < K; ++j) f[j] = (double)Ls[j] + (double)Lt[j];
    } else if (storage == CODE8) {
        const int8_t *Qs = &C.Lq8[(size_t)rs * K], *Qt = &C.Lq8[(size_t)rt * K];
        const float *As = &C.La[(size_t)rs * C.n_anchor], *At = &C.La[(size_t)rt * C.n_anchor];
        double ds = 0, dt = 0;
        int slot = 0;
        for (int j = 0; j < K; ++j) {
            if (j == 0 || C.run[j] != C.run[j - 1]) { ds = As[slot]; dt = At[slot]; ++slot; }
            else { ds += Qs[j] * C.qstep[j]; dt += Qt[j] * C.qstep[j]; }
            f[j] = ds + dt;
            err[j] = 2 * C.qerr[j];
        }
    } else if (storage == CODE16) {
        const int16_t *Qs = &C.Lq[(size_t)rs * K], *Qt = &C.Lq[(size_t)rt * K];
        const float *As = &C.La[(size_t)rs * C.n_anchor], *At = &C.La[(size_t)rt * C.n_anchor];
        double ds = 0, dt = 0;
        int slot = 0;
        for (int j = 0; j < K; ++j) {
            if (j == 0 || C.run[j] != C.run[j - 1]) { ds = As[slot]; dt = At[slot]; ++slot; }
            else { ds += Qs[j] * C.qstep[j]; dt += Qt[j] * C.qstep[j]; }
            f[j] = ds + dt;
            err[j] = 2 * C.qerr[j];
        }
    } else { // covers
        auto recon = [&](int r, std::vector<double>& up, std::vector<double>& lo) {
            up.assign(K, INF); lo.assign(K, -INF);
            int b = C.cov_ptr[r], e = C.cov_ptr[r + 1];
            // forward and backward passes along each run: up = min_q l_q + |A - A_q|, lo = max_q l_q - |A - A_q|
            double cu = INF, cl = -INF; int c = b;
            for (int j = 0; j < K; ++j) {
                if (j > 0 && C.run[j] == C.run[j - 1]) { double dA = C.arc[j] - C.arc[j - 1]; cu += dA; cl -= dA; }
                else { cu = INF; cl = -INF; }
                while (c < e && C.cov_idx[c] < j) ++c;
                if (c < e && C.cov_idx[c] == j) { cu = std::min(cu, (double)C.cov_val[c]); cl = std::max(cl, (double)C.cov_val[c]); }
                up[j] = std::min(up[j], cu); lo[j] = std::max(lo[j], cl);
            }
            cu = INF; cl = -INF; c = e - 1;
            for (int j = K - 1; j >= 0; --j) {
                if (j < K - 1 && C.run[j] == C.run[j + 1]) { double dA = C.arc[j + 1] - C.arc[j]; cu += dA; cl -= dA; }
                else { cu = INF; cl = -INF; }
                while (c >= b && C.cov_idx[c] > j) --c;
                if (c >= b && C.cov_idx[c] == j) { cu = std::min(cu, (double)C.cov_val[c]); cl = std::max(cl, (double)C.cov_val[c]); }
                up[j] = std::min(up[j], cu); lo[j] = std::max(lo[j], cl);
            }
        };
        thread_local std::vector<double> us, ls_, ut, lt;
        recon(rs, us, ls_); recon(rt, ut, lt);
        for (int j = 0; j < K; ++j) { f[j] = us[j] + ut[j]; fl[j] = ls_[j] + lt[j]; }
        flo = fl.data();
    }
    double ub = INF, lb = INF;
    int imin = 0;
    bool use_lb = exact && C.lb_ok;
    for (int j = 0; j < K; ++j) {
        if (f[j] < ub) { ub = f[j]; imin = j; }
        if (use_lb && j > 0 && std::isfinite(C.gap[j - 1]))
            lb = std::min(lb, 0.5 * (flo[j - 1] + flo[j]) - C.gap[j - 1] - 0.5 * (err[j - 1] + err[j]));
    }
    if (use_lb) for (size_t e = 0; e < C.end_idx.size(); ++e) lb = std::min(lb, flo[C.end_idx[e]] - C.end_gap[e] - err[C.end_idx[e]]);
    if (!use_lb || !std::isfinite(lb)) lb = 0;
    lb *= 1.0 - 2 * F32_REL;
    lb = std::max(lb, chord());
    double est = ub;
    if (opt.interp) {
        if (storage == CODE8 || storage == CODE16) {
            // coded labels sit above the stored ones by half a quantum on average: the estimate (not the
            // bounds) is taken from the de-biased values
            thread_local std::vector<double> fe;
            fe.resize(std::max(K, 1));
            int ie = 0;
            for (int j = 0; j < K; ++j) { fe[j] = f[j] - 0.5 * err[j]; if (fe[j] < fe[ie]) ie = j; }
            est = interp(C, fe.data(), K, ie, rs, rt);
        } else {
            est = interp(C, f.data(), K, imin, rs, rt);
        }
    }
    out[0] = lb; out[1] = ub; out[2] = std::min(std::max(est, lb), ub);
}

void CertifiedIndex::bounds(int s, int t, double* out) const {
    if (s == t) { out[0] = out[1] = out[2] = 0; return; }
    const int32_t ss = slot_of[s];
    if (ss < 0 || slot_of[t] < 0) { out[0] = 0; out[1] = out[2] = INF; return; } // POI mode: not indexed
    double dx = T.xy[2 * s] - T.xy[2 * t], dy = T.xy[2 * s + 1] - T.xy[2 * t + 1];
    if (dx * dx + dy * dy <= Rq * Rq) {
        auto b0 = near_ids.begin() + near_ptr[ss], b1 = near_ids.begin() + near_ptr[ss + 1];
        auto it = std::lower_bound(b0, b1, t);
        if (it != b1 && *it == t) {
            double d = near_d[it - near_ids.begin()];
            const double *A = &T.V[3 * s], *B = &T.V[3 * t];
            out[0] = exact ? d * (1.0 - F32_REL) : lam * std::sqrt(std::pow(A[0] - B[0], 2) + std::pow(A[1] - B[1], 2) + std::pow(A[2] - B[2], 2));
            out[1] = out[2] = d;
            return;
        }
    }
    int cid, rs, rt;
    if (!select_cut(s, t, cid, rs, rt)) { out[0] = 0; out[1] = out[2] = INF; return; }
    eval_cut(cuts[cid], s, t, rs, rt, out);
}

// exact distance of a pair: one early-stopping MMP sweep on the faces meeting bbox(s, t) dilated by
// ub / (2 lambda) (the ellipse lemma: every path of length <= ub stays inside)
double CertifiedIndex::exact_pair(int s, int t, double ub) const {
    double R = ub / (2.0 * lam) + 1e-6 * T.edge;
    double lo[2], hi[2];
    for (int k = 0; k < 2; ++k) { lo[k] = std::min(T.xy[2 * s + k], T.xy[2 * t + k]) - R; hi[k] = std::max(T.xy[2 * s + k], T.xy[2 * t + k]) + R; }
    std::vector<int> faces = grid.in_box(lo[0], hi[0], lo[1], hi[1]);
    std::vector<double> V2, Xu; std::vector<unsigned> F2; std::vector<int> used, pidx;
    double along[2] = {1, 0};
    submesh(faces, {}, along, V2, F2, used, pidx, Xu);
    auto loc = [&](int v) { return (unsigned)(std::lower_bound(used.begin(), used.end(), v) - used.begin()); };
    std::unique_ptr<Solver> sol(make_mmp(V2, F2, true));
    double d;
    std::vector<unsigned> tg{loc(t)};
    sol->dist(loc(s), tg, &d);
    return d;
}

bool CertifiedIndex::bounds_refined(int s, int t, double* out) const {
    bounds(s, t, out);
    if (s == t || out[1] <= (1.0 + eps) * out[0]) return false;
    if (!exact) return false; // weighted terrains: no exact solver
    double d = exact_pair(s, t, std::isfinite(out[1]) ? out[1] : 4 * lam * std::hypot(T.xy[2 * s] - T.xy[2 * t], T.xy[2 * s + 1] - T.xy[2 * t + 1]) + 10 * T.edge);
    out[0] = d * (1 - 1e-12); out[1] = d * (1 + 1e-12); out[2] = d;
    return true;
}

// 16-bit coding: along a run, labels are 1-Lipschitz in arc length (the cut is a path), so the step to the
// next label lies in [-(A_j + e_{j-1}), A_j] with e the previous decoding error. Each step is stored as
// ceil(step / q_j) with q_j = (A_j + e_{j-1}) / 32767 (shared by all rows of the cut), decoded values are
// re-used (closed loop), so every decoded label is >= the stored one and exceeds it by < q_j
bool CertifiedIndex::encode_cut(Cut& C, int bits) {
    const double qmax = bits == 8 ? 127.0 : 32767.0, qmin = bits == 8 ? -128.0 : -32768.0;
    int K = C.K();
    size_t nm = C.M.size();
    if (K == 0) { C.Lq.clear(); C.Lq8.clear(); C.La.clear(); return false; }
    float Lmax = 0;
    for (float x : C.L) if (std::isfinite(x)) Lmax = std::max(Lmax, x);
    double ulp = std::ldexp((double)Lmax, -22); // a float32 anchor's error, generously
    C.qstep.assign(K, 0.0); C.qerr.assign(K, ulp);
    C.n_anchor = 0;
    for (int j = 0; j < K; ++j) {
        if (j == 0 || C.run[j] != C.run[j - 1]) { C.n_anchor++; continue; }
        double A = C.arc[j] - C.arc[j - 1];
        C.qstep[j] = (A * (1 + 1e-6) + C.qerr[j - 1] + 4 * ulp) / qmax;
        C.qerr[j] = C.qstep[j];
    }
    if (bits == 8) { C.Lq8.assign(nm * K, 0); C.Lq.clear(); } else { C.Lq.assign(nm * K, 0); C.Lq8.clear(); }
    C.La.assign(nm * C.n_anchor, 0.f);
    bool ok = true;
    for (size_t r = 0; r < nm && ok; ++r) {
        const float* Lr = &C.L[r * K];
        double dec = 0;
        int slot = 0;
        for (int j = 0; j < K; ++j) {
            if (j == 0 || C.run[j] != C.run[j - 1]) {
                C.La[r * C.n_anchor + slot++] = Lr[j];
                dec = Lr[j];
            } else {
                double q = std::ceil(((double)Lr[j] - dec) / C.qstep[j]);
                if (q < qmin || q > qmax || !std::isfinite(q)) { ok = false; break; } // cannot happen (1-Lipschitz rows)
                if (bits == 8) C.Lq8[r * K + j] = (int8_t)q; else C.Lq[r * K + j] = (int16_t)q;
                dec += q * C.qstep[j];
            }
        }
    }
    if (!ok) { C.Lq.clear(); C.Lq8.clear(); C.La.clear(); } // this cut stays dense
    return ok;
}

void CertifiedIndex::encode16(int bits) {
    for (auto& C : cuts) {
        if (C.L.empty() && C.K() > 0 && !C.M.empty()) { // coded during the build: the labels are gone
            bool have = bits == 8 ? !C.Lq8.empty() : !C.Lq.empty();
            if (!have) std::fprintf(stderr, "encode16(%d): cut coded at the other width during the build, left as is\n", bits);
            continue;
        }
        encode_cut(C, bits);
    }
}

// eps_c-covers: per row keep a subset Q of portals such that every portal p of the same run has some q
// in Q with l(q) + A(q, p) <= (1 + eps_c) l(p). Greedy outwards from each run's minimum (valid for any
// row; near-minimal for the V-shaped rows of a terrain)
void CertifiedIndex::build_covers(double eps_c) {
    for (auto& C : cuts)
        if (C.L.empty() && C.K() > 0 && !C.M.empty()) { std::fprintf(stderr, "build_covers: labels were coded during the build (store != dense)\n"); return; }
    cover_eps = eps_c;
    cover_entries = 0;
    for (auto& C : cuts) {
        int K = C.K();
        size_t nm = C.M.size();
        C.cov_ptr.assign(nm + 1, 0); C.cov_idx.clear(); C.cov_val.clear();
        std::vector<char> keep(K);
        for (size_t r = 0; r < nm; ++r) {
            const float* Lr = &C.L[r * K];
            std::fill(keep.begin(), keep.end(), 0);
            int b = 0;
            while (b < K) {
                int e = b;
                while (e + 1 < K && C.run[e + 1] == C.run[b]) ++e;
                int m = b;
                for (int j = b; j <= e; ++j) if (Lr[j] < Lr[m]) m = j;
                keep[m] = 1;
                double c = (double)Lr[m] - C.arc[m]; // best l_q - A_q of kept portals left of j (going right)
                for (int j = m + 1; j <= e; ++j) {
                    if ((1 + eps_c) * Lr[j] < c + C.arc[j]) { keep[j] = 1; c = std::min(c, (double)Lr[j] - C.arc[j]); }
                }
                c = (double)Lr[m] + C.arc[m];       // best l_q + A_q going left
                for (int j = m - 1; j >= b; --j) {
                    if ((1 + eps_c) * Lr[j] < c - C.arc[j]) { keep[j] = 1; c = std::min(c, (double)Lr[j] + C.arc[j]); }
                }
                b = e + 1;
            }
            for (int j = 0; j < K; ++j) if (keep[j]) { C.cov_idx.push_back((uint16_t)j); C.cov_val.push_back(Lr[j]); }
            C.cov_ptr[r + 1] = (int32_t)C.cov_idx.size();
        }
        cover_entries += (int64_t)C.cov_idx.size();
    }
}

size_t CertifiedIndex::nbytes() const {
    size_t b = 0;
    for (auto& C : cuts) {
        if (storage == CODE16 && !C.Lq.empty()) b += C.Lq.size() * 2 + C.La.size() * 4 + C.qstep.size() * 8;
        else if (storage == CODE8 && !C.Lq8.empty()) b += C.Lq8.size() + C.La.size() * 4 + C.qstep.size() * 8;
        else if (storage == COVER && !C.cov_ptr.empty()) b += C.cov_idx.size() * 6 + C.cov_ptr.size() * 4;
        else b += C.L.size() * 4;
        b += C.arc.size() * 8 + C.run.size() * 8 + C.gap.size() * 8 + C.X.size() * 8;
        b += C.end_idx.size() * 16 + 1;
        b += C.foot_r.size() * 4 + C.foot_arc.size() * 4 + C.foot_run.size() * 2;
    }
    b += near_ids.size() * 4 + near_d.size() * 4 + near_ptr.size() * 8;
    b += anc.size() * 4 + row.size() * 4;
    if (!poi_list.empty()) b += slot_of.size() * 4; // vertex -> POI row
    for (auto& Cp : copies) b += Cp.anc.size() * 4 + Cp.row.size() * 4 + Cp.nodes.size() * sizeof(ShiftNode);
    return b;
}

}  // namespace cci
