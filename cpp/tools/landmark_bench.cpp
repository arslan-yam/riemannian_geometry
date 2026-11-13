// Exact-landmark interval baseline (review item 3; ALT-style bounds, Goldberg and Harrelson 2005, on the exact
// surface metric, as in Burghard and Klein 2015): L landmarks chosen by farthest-point sampling, one exact sweep
// each over the whole terrain (VTP, MMP if VTP rejects the mesh), every vertex stores its exact distance to every
// landmark (float32 rounded up). A query (s, t) returns
//   UB = min_l d(l,s) + d(l,t)                      (the length of a path through l)
//   LB = max(|S - T|, max_l |d(l,s) - d(l,t)|)      (triangle inequality; chord: 3-D length dominates it)
// with the float32 margin subtracted from LB, so the interval is sound like ours. It isolates what the cut
// certificate adds over a plain landmark interval at the same memory. Reported per L (prefixes of one sampling):
// interval width, the fraction certified at 1 + eps (UB <= (1+eps) LB), lookup time, and certify-or-refine cost
// (an exact early-stopping MMP sweep for each pair that is not certified, timed on a sample). With --prox, kNN
// queries over the exported objects of cpp/tools/proximity (--export): candidates decided from the intervals,
// ambiguous ones refined by one early-stopping MMP sweep to them, and the answer checked against the truth.
// usage: landmark_bench mesh.off pairs truth [--L 16,64,256,1024] [--eps 0.2] [--prox PREFIX.pairs] [--k 1,10,100]
//                       [--refine-sample 100] [--seed 7] [--json out.json]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "../src/mesh_io.h"
#include "../src/solver.h"
using namespace cci;

typedef std::chrono::steady_clock clk;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }
static std::vector<double> parse_list(const std::string& s) {
    std::vector<double> v; std::stringstream ss(s); std::string x;
    while (std::getline(ss, x, ',')) if (!x.empty()) v.push_back(std::stod(x));
    return v;
}
static float f32_up(double x) { float f = (float)x; if ((double)f < x) f = std::nextafter(f, INFINITY); return f; }
static const double REL = 1.0 / (1 << 22);  // float32 rounding margin (two rounded-up values)

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: landmark_bench mesh.off pairs truth [--L 16,64,256,1024] [--eps 0.2] [--prox P.pairs] [--k 1,10,100] [--refine-sample 100] [--seed 7] [--json out]\n"); return 1; }
    std::string Ls = "16,64,256,1024", json, prox, ks_s = "1,10,100"; double eps = 0.2; int seed = 7, rsample = 100;
    for (int i = 4; i + 1 < argc; ++i) {
        std::string a = argv[i]; const char* v = argv[i + 1];
        if (a == "--L") Ls = v; else if (a == "--eps") eps = std::stod(v); else if (a == "--json") json = v;
        else if (a == "--prox") prox = v; else if (a == "--k") ks_s = v; else if (a == "--seed") seed = std::stoi(v);
        else if (a == "--refine-sample") rsample = std::stoi(v);
    }
    std::vector<int> Lv; for (double x : parse_list(Ls)) Lv.push_back((int)x);
    std::sort(Lv.begin(), Lv.end());
    const int Lmax = Lv.back();
    std::vector<int> ks; for (double x : parse_list(ks_s)) ks.push_back((int)x);

    std::vector<double> V; std::vector<unsigned> F;
    read_off(argv[1], V, F);
    const int n = (int)V.size() / 3;
    std::vector<std::pair<int, int>> qp; std::vector<double> tr;
    { std::ifstream pf(argv[2]); int s, t; std::string rest; while (pf >> s >> t) { qp.push_back({s, t}); std::getline(pf, rest); } }
    { std::ifstream tf(argv[3]); double d; while (tf >> d) tr.push_back(d); }
    if (tr.size() < qp.size()) qp.resize(tr.size());
    const size_t m = qp.size();

    // ---- landmarks: farthest-point sampling, one full exact sweep each
    std::vector<unsigned> F2(F);
    for (size_t f = 0; f < F2.size(); f += 3) { // counter-clockwise in the chart, as VTP needs
        const double *A = &V[3 * F2[f]], *B = &V[3 * F2[f + 1]], *C = &V[3 * F2[f + 2]];
        if ((B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]) < 0) std::swap(F2[f + 1], F2[f + 2]);
    }
    std::unique_ptr<Solver> sol(make_vtp(V, F2, true));
    if (!sol) sol.reset(make_mmp(V, F, false));
    std::vector<unsigned> all(n); for (int i = 0; i < n; ++i) all[i] = i;
    std::vector<float> D((size_t)n * Lmax);  // vertex-major: the row of a vertex is contiguous
    std::vector<double> mind(n, INFINITY), row(n), cum_s(Lmax + 1, 0.0);
    std::mt19937_64 rng(seed);
    int lm = (int)(rng() % n);
    std::vector<int> lms;
    auto t_build = clk::now();
    for (int l = 0; l < Lmax; ++l) {
        lms.push_back(lm);
        sol->dist(lm, all, row.data());
        int far = 0;
        for (int v = 0; v < n; ++v) {
            D[(size_t)v * Lmax + l] = f32_up(row[v]);
            mind[v] = std::min(mind[v], row[v]);
            if (mind[v] > mind[far]) far = v;
        }
        lm = far;
        cum_s[l + 1] = since(t_build);
        if ((l + 1) % 64 == 0) { std::printf("  %d landmarks, %.1f s\n", l + 1, cum_s[l + 1]); std::fflush(stdout); }
    }
    std::printf("landmark_bench %s: n %d, %d landmarks by %s in %.1f s (%.2f s per sweep)\n", argv[1], n, Lmax, sol->name(), cum_s[Lmax], cum_s[Lmax] / Lmax);
    std::fflush(stdout);

    auto chord = [&](int s, int t) { double dx = V[3 * s] - V[3 * t], dy = V[3 * s + 1] - V[3 * t + 1], dz = V[3 * s + 2] - V[3 * t + 2]; return std::sqrt(dx * dx + dy * dy + dz * dz); };
    auto bounds = [&](int s, int t, int L, double& lb, double& ub) {
        if (s == t) { lb = ub = 0; return; }
        const float *a = &D[(size_t)s * Lmax], *b = &D[(size_t)t * Lmax];
        double u = INFINITY, g = 0;  // in double: a float sum could round below the true sum
        for (int l = 0; l < L; ++l) { double x = a[l], y = b[l]; u = std::min(u, x + y); g = std::max(g, std::fabs(x - y) - REL * std::max(x, y)); }
        ub = u; lb = std::max((double)g, chord(s, t) * (1 - 1e-12));
    };
    std::unique_ptr<Solver> mmp(make_mmp(V, F, true));  // the refine: an early-stopping exact sweep

    std::ostringstream js;
    js << "{\"tool\":\"landmark_bench\",\"mesh\":\"" << argv[1] << "\",\"n\":" << n << ",\"eps\":" << eps << ",\"solver\":\"" << sol->name()
       << "\",\"pairs\":" << m << ",\"landmarks\":[";
    for (size_t li = 0; li < Lv.size(); ++li) {
        const int L = Lv[li];
        std::vector<double> lb(m), ub(m);
        for (size_t i = 0; i < m; ++i) bounds(qp[i].first, qp[i].second, L, lb[i], ub[i]);
        double best = 1e30; volatile double sink = 0;
        for (int rep = 0; rep < 4; ++rep) { // one warm-up pass, best of three
            auto t0 = clk::now(); double lo, hi;
            for (auto& q : qp) { bounds(q.first, q.second, L, lo, hi); sink = sink + lo + hi; }
            if (rep) best = std::min(best, since(t0) / m);
        }
        bool sound = true; int cert = 0; double ub_mean = 0, ub_max = 0, lb_mean = 0;
        std::vector<double> widths, ratios;
        std::vector<int> uncert;
        for (size_t i = 0; i < m; ++i) {
            double d = tr[i];
            if (!(lb[i] <= d * (1 + 1e-9) && d <= ub[i] * (1 + 1e-9))) {
                if (sound) std::fprintf(stderr, "UNSOUND L=%d pair %zu (%d, %d): LB %.9g d %.9g UB %.9g chord %.9g\n", L, i, qp[i].first, qp[i].second, lb[i], d, ub[i], chord(qp[i].first, qp[i].second));
                sound = false;
            }
            ub_mean += (ub[i] / d - 1) / m; ub_max = std::max(ub_max, ub[i] / d - 1); lb_mean += (1 - lb[i] / d) / m;
            widths.push_back(lb[i] > 0 ? ub[i] / lb[i] - 1 : INFINITY);
            if (ub[i] <= (1 + eps) * lb[i]) ++cert; else uncert.push_back((int)i);
        }
        std::sort(widths.begin(), widths.end());
        // refine cost: an exact early-stopping sweep per uncertified pair, timed on a sample
        double rsum = 0; int rn = 0;
        for (int j = 0; j < (int)uncert.size() && rn < rsample; ++j) {
            int i = uncert[(size_t)j * uncert.size() / std::min<size_t>(uncert.size(), rsample)];
            std::vector<unsigned> tg{(unsigned)qp[i].second}; double out;
            auto t0 = clk::now(); mmp->dist(qp[i].first, tg, &out); rsum += since(t0); ++rn;
            if (std::fabs(out / tr[i] - 1) > 1e-6) std::fprintf(stderr, "refine mismatch pair %d: %.9g vs %.9g\n", i, out, tr[i]);
        }
        double rmean = rn ? rsum / rn : 0, frac = (double)uncert.size() / m, e2e = best + frac * rmean;
        std::printf("L=%d (%.1f MB, build %.1f s): lookup %.3f us, sound %s, UB err mean %.4f%% max %.3f%%, LB err mean %.4f%%, width median %.3f%% p90 %.3f%%, "
                    "certified at 1+eps %d/%zu, refine %.1f ms per pair (sample %d), end-to-end mean %.2f ms\n",
                    L, (double)n * L * 4 / 1e6, cum_s[L], best * 1e6, sound ? "yes" : "NO", 100 * ub_mean, 100 * ub_max, 100 * lb_mean,
                    100 * widths[m / 2], 100 * widths[(size_t)(0.9 * (m - 1))], cert, m, 1e3 * rmean, rn, 1e3 * e2e);
        std::fflush(stdout);
        js << (li ? "," : "") << "{\"L\":" << L << ",\"index_bytes\":" << (long long)n * L * 4 << ",\"build_s\":" << cum_s[L] << ",\"query_us\":" << best * 1e6
           << ",\"sound\":" << (sound ? "true" : "false") << ",\"ub_err_mean\":" << ub_mean << ",\"ub_err_max\":" << ub_max << ",\"lb_err_mean\":" << lb_mean
           << ",\"width_median\":" << widths[m / 2] << ",\"width_p90\":" << widths[(size_t)(0.9 * (m - 1))] << ",\"certified\":" << cert
           << ",\"refine_ms\":" << 1e3 * rmean << ",\"refine_sample\":" << rn << ",\"end_to_end_ms\":" << 1e3 * e2e << ",\"knn\":[";

        // ---- certified kNN over the exported objects
        if (!prox.empty()) {
            std::vector<int> qv; std::vector<std::vector<std::pair<int, double>>> objs;
            { std::ifstream in(prox); long q, o; double d; std::unordered_map<long, int> at;
              while (in >> q >> o >> d) { auto it = at.find(q); if (it == at.end()) { at[q] = (int)qv.size(); qv.push_back((int)q); objs.emplace_back(); it = at.find(q); } objs[it->second].push_back({(int)o, d}); } }
            const int Q = (int)qv.size();
            for (size_t ki = 0; ki < ks.size(); ++ki) {
                const int k = ks[ki];
                double t_int = 0, t_ref = 0, amb_sum = 0; int refined = 0, exact = 0;
                for (int qi = 0; qi < Q; ++qi) {
                    auto& ob = objs[qi]; const int M = (int)ob.size(), kk = std::min(k, M);
                    auto t0 = clk::now();
                    std::vector<double> lo(M), hi(M);
                    for (int j = 0; j < M; ++j) bounds(qv[qi], ob[j].first, L, lo[j], hi[j]);
                    std::vector<double> slo(lo), shi(hi);
                    std::nth_element(slo.begin(), slo.begin() + (kk - 1), slo.end());
                    std::nth_element(shi.begin(), shi.begin() + (kk - 1), shi.end());
                    // out: at least k objects certainly closer (k-th smallest UB below its LB); in: fewer than k other
                    // objects can be closer (we count those whose LB is below its UB)
                    std::sort(slo.begin(), slo.end());
                    const double kub = shi[kk - 1];
                    std::vector<int> in_, amb;
                    for (int j = 0; j < M; ++j) {
                        if (lo[j] > kub) continue;  // out
                        int could = (int)(std::lower_bound(slo.begin(), slo.end(), hi[j]) - slo.begin()) - (lo[j] < hi[j] ? 1 : 0);  // others with LB < UB_j
                        if (could < kk) in_.push_back(j); else amb.push_back(j);
                    }
                    t_int += since(t0);
                    std::vector<std::pair<double, int>> fin;
                    for (int j : in_) fin.push_back({hi[j], j});
                    if ((int)in_.size() < kk) { // refine the ambiguous ones by one exact sweep to them
                        auto t1 = clk::now();
                        std::vector<unsigned> tg; for (int j : amb) tg.push_back((unsigned)ob[j].first);
                        std::vector<double> out(tg.size());
                        mmp->dist(qv[qi], tg, out.data());
                        std::vector<std::pair<double, int>> ex;
                        for (size_t a = 0; a < amb.size(); ++a) ex.push_back({out[a], amb[a]});
                        std::sort(ex.begin(), ex.end());
                        for (int a = 0; a < kk - (int)in_.size() && a < (int)ex.size(); ++a) fin.push_back(ex[a]);
                        t_ref += since(t1); ++refined; amb_sum += amb.size();
                    }
                    // check: the true k nearest (ties by distance)
                    std::vector<double> dd; for (auto& o : ob) dd.push_back(o.second);
                    std::vector<double> sd(dd); std::nth_element(sd.begin(), sd.begin() + (kk - 1), sd.end());
                    const double dk = sd[kk - 1]; bool ok = (int)fin.size() == kk;
                    for (auto& f : fin) ok = ok && dd[f.second] <= dk * (1 + 1e-6);
                    exact += ok;
                }
                std::printf("  kNN k=%d over %zu objects, %d queries: intervals %.3f ms, refined %d queries (%.1f ambiguous objects), refine %.2f ms per refined query, "
                            "end-to-end mean %.3f ms, exact answers %d/%d\n", k, objs.empty() ? 0 : objs[0].size(), Q, 1e3 * t_int / Q, refined,
                            refined ? amb_sum / refined : 0.0, refined ? 1e3 * t_ref / refined : 0.0, 1e3 * (t_int + t_ref) / Q, exact, Q);
                std::fflush(stdout);
                js << (ki ? "," : "") << "{\"k\":" << k << ",\"queries\":" << Q << ",\"interval_ms\":" << 1e3 * t_int / Q << ",\"refined_queries\":" << refined
                   << ",\"ambiguous_mean\":" << (refined ? amb_sum / refined : 0.0) << ",\"refine_ms\":" << (refined ? 1e3 * t_ref / refined : 0.0)
                   << ",\"end_to_end_ms\":" << 1e3 * (t_int + t_ref) / Q << ",\"exact_answers\":" << exact << "}";
            }
        }
        js << "]}";
        if (!json.empty()) { std::ofstream o(json); o << js.str() << "],\"complete\":" << (li + 1 == Lv.size() ? "true" : "false") << "}\n"; }
    }
    return 0;
}
