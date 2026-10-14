// Proximity queries with the certified index: k nearest objects and range queries from a query vertex,
// against exact on-the-fly search (one full VTP sweep per query vertex).
//
// Objects: M random vertices (or a POI file). For a query vertex q:
//  * kNN, certified: the objects within growing chart radius until k are found give tau = the k-th
//    smallest UB; every object with lambda |q - o| <= tau is a candidate (LB >= lambda |q - o|), and gets
//    its interval. An object is surely in if fewer than k other candidates have LB < its UB, surely out if
//    at least k candidates have UB < its LB; the rest are refined together by one exact MMP sweep on the
//    ellipse sub-mesh (CertifiedIndex::exact_from). The answer set is exact.
//  * kNN, estimate: the k candidates of smallest estimate (no refine).
//  * range r, certified: candidates within lambda |q - o| <= r; in if UB <= r, out if LB > r, the rest refined.
//  * range r, estimate: estimate <= r.
//  * certified approximate (no refine, the candidate pass only): kNN = the k candidates of smallest UB,
//    each provably within (1 + eps) of the true k-th distance (every true top-k object o has UB(o) <=
//    (1 + eps) d(o), so the k-th smallest UB is <= (1 + eps) d_k); range = the objects with UB <= r: no
//    false positives, and every object with d <= r / (1 + eps) is in.
// Baselines without an index: one full VTP sweep per query vertex (exact distances to every object), and
// exact search with straight-line pruning (the stronger one): MMP from q with early stop at the k
// objects nearest in the chart gives tau = the largest of their distances; every object with
// lambda |q - o| <= tau is a candidate, and a second early-stopping MMP sweep to the candidates gives the
// exact k nearest (range: one sweep to the objects with lambda |q - o| <= r).
// usage: proximity mesh.off [--eps 0.2] [--gap-factor 1] [--kmax K] [--store code8] [--near-store sym]
//        [--threads 1] [--objects 10000 | --poi-file f] [--queries 100] [--k 1,10,100]
//        [--radius 0.01,0.03,0.1] (fractions of the chart diagonal) [--seed 7] [--truth-cache f] [--json out.json]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include "../src/certified.h"
#include "../src/mesh_io.h"
using namespace cci;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

static std::vector<double> parse_list(const std::string& s) {
    std::vector<double> v; std::stringstream ss(s); std::string x;
    while (std::getline(ss, x, ',')) v.push_back(std::stod(x));
    return v;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: proximity mesh.off [options]\n"); return 1; }
    std::string mesh = argv[1], poi_file, truth_cache, json, klist = "1,10,100", rlist = "0.01,0.03,0.1";
    Options opt; opt.eps = 0.2; opt.gap_factor = 1;
    int M = 10000, Q = 100, seed = 7;
    bool profile = false;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--eps") opt.eps = std::stod(nx());
        else if (a == "--gap-factor") opt.gap_factor = std::stod(nx());
        else if (a == "--kmax") opt.k_max = std::stoi(nx());
        else if (a == "--store") opt.store = nx();
        else if (a == "--near-store") opt.near_store = nx();
        else if (a == "--threads") opt.threads = std::stoi(nx());
        else if (a == "--solver-mem-gb") opt.solver_mem_mb = 1e3 * std::stod(nx());
        else if (a == "--objects") M = std::stoi(nx());
        else if (a == "--poi-file") poi_file = nx();
        else if (a == "--queries") Q = std::stoi(nx());
        else if (a == "--k") klist = nx();
        else if (a == "--radius") rlist = nx();
        else if (a == "--seed") seed = std::stoi(nx());
        else if (a == "--truth-cache") truth_cache = nx();
        else if (a == "--json") json = nx();
        else if (a == "--profile") profile = true; // per-call timing of the kNN candidate pass, by the answering cut's size
    }
    std::vector<double> V; std::vector<unsigned> F;
    read_off(mesh, V, F);
    Terrain T; T.init(V, F);
    std::mt19937_64 rng(seed);
    // objects and query vertices
    std::vector<int> obj;
    if (!poi_file.empty()) {
        std::ifstream in(poi_file); int v;
        while (in >> v) obj.push_back(v);
        if (obj.size() > 1 && obj[0] == (int)obj.size() - 1) obj.erase(obj.begin());
    } else {
        std::vector<int> all(T.n);
        for (int i = 0; i < T.n; ++i) all[i] = i;
        std::shuffle(all.begin(), all.end(), rng);
        obj.assign(all.begin(), all.begin() + std::min(M, T.n));
    }
    std::sort(obj.begin(), obj.end());
    obj.erase(std::unique(obj.begin(), obj.end()), obj.end());
    M = (int)obj.size();
    std::vector<int> qv(Q);
    for (auto& q : qv) q = (int)(rng() % T.n);
    // truth: exact distances from each query vertex to every object, one full VTP sweep each (the
    // on-the-fly baseline); its time is the baseline's query time
    std::vector<float> D((size_t)Q * M);
    std::vector<double> sweep_s(Q, 0.0);
    // the cache starts with (n, M, Q, seed): a cache made for other objects or queries is not used
    bool cached = false;
    const int64_t hdr[4] = {T.n, M, Q, seed};
    if (!truth_cache.empty()) {
        std::ifstream in(truth_cache, std::ios::binary);
        int64_t h[4] = {0, 0, 0, 0};
        if (in && in.read((char*)h, sizeof h) && std::equal(h, h + 4, hdr)) {
            in.read((char*)D.data(), D.size() * sizeof(float));
            in.read((char*)sweep_s.data(), Q * sizeof(double));
            cached = (bool)in;
        }
        if (!cached && in) std::fprintf(stderr, "truth cache %s: made for other objects / queries, recomputing\n", truth_cache.c_str());
    }
    std::string sweep_name;
    if (!cached) {
        std::vector<unsigned> F2(T.F.begin(), T.F.end());
        for (size_t f = 0; f < F2.size(); f += 3) { // counter-clockwise in the chart, as VTP needs
            const double *A = &T.V[3 * F2[f]], *B = &T.V[3 * F2[f + 1]], *C = &T.V[3 * F2[f + 2]];
            if ((B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]) < 0) std::swap(F2[f + 1], F2[f + 2]);
        }
        std::unique_ptr<Solver> sol(make_vtp(T.V, F2));
        if (!sol) sol.reset(make_mmp(T.V, F2, false));
        sweep_name = sol->name();
        std::vector<unsigned> tg(obj.begin(), obj.end());
        std::vector<double> out(M);
        for (int i = 0; i < Q; ++i) {
            auto t0 = clk::now();
            sol->dist(qv[i], tg, out.data());
            sweep_s[i] = since(t0);
            for (int j = 0; j < M; ++j) D[(size_t)i * M + j] = (float)out[j];
        }
        if (!truth_cache.empty()) {
            std::ofstream o(truth_cache, std::ios::binary);
            o.write((const char*)hdr, sizeof hdr);
            o.write((const char*)D.data(), D.size() * sizeof(float));
            o.write((const char*)sweep_s.data(), Q * sizeof(double));
        }
    }
    double sweep_mean = 0; for (double x : sweep_s) sweep_mean += x / Q;
    // the pruned on-the-fly baseline: one MMP solver (early stop) on the whole mesh, built once
    std::unique_ptr<Solver> mmp(make_mmp(T.V, T.F, true));
    // the index
    auto t0 = clk::now();
    CertifiedIndex ix(T, opt);
    double tb = since(t0);
    // object grid (chart coordinates)
    std::vector<double> oxy(2 * M);
    for (int j = 0; j < M; ++j) { oxy[2 * j] = T.xy[2 * obj[j]]; oxy[2 * j + 1] = T.xy[2 * obj[j] + 1]; }
    const double diag = std::hypot(T.hi[0] - T.lo[0], T.hi[1] - T.lo[1]);
    PointGrid og; og.init(oxy, std::max(diag / std::sqrt((double)M) , T.edge));
    const double lam = ix.lam;
    auto chart = [&](int q, int j) { return std::hypot(T.xy[2 * q] - oxy[2 * j], T.xy[2 * q + 1] - oxy[2 * j + 1]); };

    std::string argv_s;
    for (int i = 0; i < argc; ++i) { for (const char* c = argv[i]; *c; ++c) { if (*c == '"' || *c == '\\') argv_s += '\\'; argv_s += *c; } if (i + 1 < argc) argv_s += ' '; }
    std::ostringstream js;
    js << "{\"argv\":\"" << argv_s << "\",\"mesh\":\"" << mesh << "\",\"n\":" << T.n << ",\"objects\":" << M << ",\"queries\":" << Q << ",\"eps\":" << opt.eps
       << ",\"gap_factor\":" << opt.gap_factor << ",\"k_max\":" << opt.k_max << ",\"store\":\"" << opt.store << "\",\"near_store\":\"" << opt.near_store
       << "\",\"build_s\":" << tb << ",\"index_bytes\":" << ix.nbytes() << ",\"sweep_ms_mean\":" << sweep_mean * 1e3 << ",\"knn\":[";
    // rewritten after every k and radius, so a job that dies keeps what it measured
    auto flush = [&](int stage) { // 0: inside knn, 1: inside range, 2: done
        if (json.empty()) return;
        std::ofstream f(json);
        f << js.str() << (stage == 0 ? "],\"range\":[],\"complete\":false}" : stage == 1 ? "],\"complete\":false}" : "") << "\n";
    };
    std::printf("%s: n=%d, %d objects, %d queries | index build %.1f s, %.2f MB | exact on the fly (one full sweep): %.2f ms per query\n",
                mesh.c_str(), T.n, M, Q, tb, ix.nbytes() / 1e6, sweep_mean * 1e3);

    struct Cand { int j; double lb, ub, est; };
    auto pct = [](std::vector<double> v, double p) { std::sort(v.begin(), v.end()); return v.empty() ? 0.0 : v[(size_t)(p * (v.size() - 1))]; };
    // ---------------------------------------------------------------- kNN
    bool firstk = true;
    for (double kd : parse_list(klist)) {
        const int k = std::min((int)kd, M);
        std::vector<double> tc, ta, ncand, nref, tb2, tg_;
        double prof_grid = 0, prof_bounds = 0; size_t prof_calls = 0;
        double prof_bin_t[5] = {0, 0, 0, 0, 0}, prof_worst = 0; long prof_bin_n[5] = {0, 0, 0, 0, 0};
        int prof_worst_K = 0, prof_worst_cid = 0; size_t prof_worst_nm = 0;
        int correct = 0, refined_q = 0, bad2 = 0, ap_bad = 0;
        double ap_ratio = 0, ap_recall = 0;
        double recall = 0, ratio_max = 0;
        for (int qi = 0; qi < Q; ++qi) {
            const int q = qv[qi];
            const float* Dq = &D[(size_t)qi * M];
            std::vector<float> sorted(Dq, Dq + M);
            std::nth_element(sorted.begin(), sorted.begin() + (k - 1), sorted.end());
            const double dk = sorted[k - 1]; // true k-th distance
            // certified
            auto c0 = clk::now();
            double r = std::max(T.edge, diag / std::sqrt((double)M));
            std::vector<int> near;
            auto p0 = clk::now();
            for (;;) { near = og.radius(T.xy[2 * q], T.xy[2 * q + 1], r); if ((int)near.size() >= k || r > 2 * diag) break; r *= 2; }
            prof_grid += since(p0);
            std::vector<Cand> C;
            double o3[3];
            p0 = clk::now();
            for (int j : near) {
                auto b0 = clk::now();
                ix.bounds(q, obj[j], o3); C.push_back({j, o3[0], o3[1], o3[2]});
                double tt = since(b0);
                if (!profile) continue;
                int cid = ix.answering_cut(q, obj[j]);
                int Kc = cid >= 0 ? ix.cuts[cid].K() : 0;
                int bin = cid < 0 ? 0 : Kc < 64 ? 1 : Kc < 256 ? 2 : Kc < 1024 ? 3 : 4;
                prof_bin_t[bin] += tt; prof_bin_n[bin]++;
                if (tt > prof_worst) { prof_worst = tt; prof_worst_K = Kc; prof_worst_cid = cid; prof_worst_nm = cid >= 0 ? ix.cuts[cid].M.size() : 0; }
            }
            prof_bounds += since(p0); prof_calls += near.size();
            std::vector<double> ubs; for (auto& c : C) ubs.push_back(c.ub);
            std::nth_element(ubs.begin(), ubs.begin() + (k - 1), ubs.end());
            double tau = ubs[k - 1];
            // every object with lambda |q - o| <= tau may be among the k nearest
            std::vector<char> have(M, 0);
            for (auto& c : C) have[c.j] = 1;
            p0 = clk::now();
            std::vector<int> more = og.radius(T.xy[2 * q], T.xy[2 * q + 1], tau / lam);
            prof_grid += since(p0);
            p0 = clk::now();
            for (int j : more)
                if (!have[j]) { ix.bounds(q, obj[j], o3); C.push_back({j, o3[0], o3[1], o3[2]}); prof_calls++; }
            prof_bounds += since(p0);
            const double t_cand = since(c0);
            {   // certified approximate kNN: the k candidates of smallest UB
                std::vector<std::pair<double, int>> us;
                for (auto& c : C) us.push_back({c.ub, c.j});
                std::partial_sort(us.begin(), us.begin() + std::min<size_t>(k, us.size()), us.end());
                double wr = 0; int hit = 0;
                for (int a = 0; a < k && a < (int)us.size(); ++a) { wr = std::max(wr, Dq[us[a].second] / dk); hit += Dq[us[a].second] <= dk * (1 + 1e-6); }
                ap_ratio = std::max(ap_ratio, wr); ap_recall += (double)hit / k / Q;
                if (wr > 1 + ix.eps + 1e-9) ap_bad++;
                ta.push_back(since(c0));
            }
            (void)t_cand;
            // classify: sorted LB and UB lists of the candidates
            std::vector<double> L, U;
            for (auto& c : C) { L.push_back(c.lb); U.push_back(c.ub); }
            std::sort(L.begin(), L.end()); std::sort(U.begin(), U.end());
            std::vector<int> in, amb;
            for (size_t i = 0; i < C.size(); ++i) {
                // others with LB < my UB (minus myself if my LB < my UB), others with UB < my LB
                int lb_below = (int)(std::lower_bound(L.begin(), L.end(), C[i].ub) - L.begin()) - (C[i].lb < C[i].ub ? 1 : 0);
                int ub_below = (int)(std::lower_bound(U.begin(), U.end(), C[i].lb) - U.begin());
                if (lb_below < k) in.push_back((int)i);
                else if (ub_below >= k) continue;
                else amb.push_back((int)i);
            }
            std::vector<int> answer;
            for (int i : in) answer.push_back(C[i].j);
            int need = k - (int)answer.size();
            double t_glob = 0, t_sub = 0;
            if (need > 0 && !amb.empty()) {
                std::vector<int> tg; double ubm = 0;
                for (int i : amb) { tg.push_back(obj[C[i].j]); ubm = std::max(ubm, C[i].ub); }
                std::vector<double> ex(tg.size()), ex2(tg.size());
                auto r0 = clk::now();
                ix.exact_from(q, tg, ubm, ex.data());             // index only: ellipse sub-mesh, fresh solver
                t_sub = since(r0);
                std::vector<unsigned> tgu(tg.begin(), tg.end());
                r0 = clk::now();
                mmp->dist((unsigned)q, tgu, ex2.data());           // a whole-mesh solver kept in memory
                t_glob = since(r0);
                std::vector<std::pair<double, int>> e;
                for (size_t a = 0; a < amb.size(); ++a) e.push_back({ex[a], C[amb[a]].j});
                std::sort(e.begin(), e.end());
                for (int a = 0; a < need && a < (int)e.size(); ++a) answer.push_back(e[a].second);
                nref.push_back((double)amb.size()); refined_q++;
            } else nref.push_back(0);
            double t_all = since(c0) - t_glob; // the certified query with the sub-mesh refine
            tc.push_back(t_all); tg_.push_back(t_all - t_sub + t_glob); ncand.push_back((double)C.size());
            // baseline: exact search with straight-line pruning (no index)
            {
                auto b0 = clk::now();
                std::vector<std::pair<double, int>> byc;
                for (int j : near) byc.push_back({chart(q, j), j});
                std::partial_sort(byc.begin(), byc.begin() + k, byc.end());
                std::vector<unsigned> tg; for (int a = 0; a < k; ++a) tg.push_back((unsigned)obj[byc[a].second]);
                std::vector<double> d1(tg.size());
                mmp->dist((unsigned)q, tg, d1.data());
                double tau2 = *std::max_element(d1.begin(), d1.end());
                std::vector<int> cj = og.radius(T.xy[2 * q], T.xy[2 * q + 1], tau2 / lam);
                std::vector<unsigned> tg2; for (int j : cj) tg2.push_back((unsigned)obj[j]);
                std::vector<double> d2(tg2.size());
                mmp->dist((unsigned)q, tg2, d2.data());
                std::vector<std::pair<double, int>> ex;
                for (size_t a = 0; a < cj.size(); ++a) ex.push_back({d2[a], cj[a]});
                std::partial_sort(ex.begin(), ex.begin() + std::min<size_t>(k, ex.size()), ex.end());
                tb2.push_back(since(b0));
                double w2 = 0; for (int a = 0; a < k && a < (int)ex.size(); ++a) w2 = std::max(w2, (double)Dq[ex[a].second]);
                if (w2 > dk * (1 + 1e-6) + 1e-9) bad2++;
            }
            double worst = 0;
            for (int j : answer) worst = std::max(worst, (double)Dq[j]);
            if ((int)answer.size() == k && worst <= dk * (1 + 1e-6) + 1e-9) correct++;
            // estimate-only (same candidates, no refine)
            auto a0 = clk::now();
            std::vector<std::pair<double, int>> es;
            for (auto& c : C) es.push_back({c.est, c.j});
            std::partial_sort(es.begin(), es.begin() + std::min<size_t>(k, es.size()), es.end());
            double t_est = since(a0);
            (void)t_est;
            int hit = 0; double wr = 0;
            for (int a = 0; a < k && a < (int)es.size(); ++a) { if (Dq[es[a].second] <= dk * (1 + 1e-6)) hit++; wr = std::max(wr, Dq[es[a].second] / dk); }
            recall += (double)hit / k / Q; ratio_max = std::max(ratio_max, wr);
        }
        double tm = 0; for (double x : tc) tm += x / Q;
        double cm = 0; for (double x : ncand) cm += x / Q;
        double rm = 0; for (double x : nref) rm += x / Q;
        double bm = 0; for (double x : tb2) bm += x / Q;
        double gm = 0; for (double x : tg_) gm += x / Q;
        double am = 0; for (double x : ta) am += x / Q;
        if (profile) std::printf("  [profile k=%d] per query: grid %.1f us, %.1f bound calls %.1f us (%.2f us per call)\n", k, prof_grid / Q * 1e6,
                    (double)prof_calls / Q, prof_bounds / Q * 1e6, prof_calls ? prof_bounds / prof_calls * 1e6 : 0.0);
        if (profile) {
        const char* bn[5] = {"near", "K<64", "K<256", "K<1024", "K>=1024"};
        for (int b = 0; b < 5; ++b) if (prof_bin_n[b]) std::printf("    %s: %ld calls, %.2f us per call\n", bn[b], prof_bin_n[b], prof_bin_t[b] / prof_bin_n[b] * 1e6);
        std::printf("    slowest call %.1f us: cut %d, K %d, %zu members\n", prof_worst * 1e6, prof_worst_cid, prof_worst_K, prof_worst_nm);
        }
        std::printf("  kNN k=%d: approximate (UB, guaranteed <= 1+eps): %.1f us mean, p50 %.1f us, recall %.4f, worst d/d_k %.4f (%d over 1+eps) | certified exact %.1f us mean (%.1f us refining with a whole-mesh solver), p50 %.1f, p99 %.1f us | %.0f candidates, refined %d/%d queries (%.1f objects each) | exact %d/%d | estimate-only: recall %.4f, worst d/d_k %.4f | on the fly: pruned MMP %.2f ms (%.0fx, %d wrong), full sweep %.2f ms (%.0fx)\n",
                    k, am * 1e6, pct(ta, 0.5) * 1e6, ap_recall, ap_ratio, ap_bad, tm * 1e6, gm * 1e6, pct(tc, 0.5) * 1e6, pct(tc, 0.99) * 1e6, cm, refined_q, Q, refined_q ? rm * Q / refined_q : 0.0, correct, Q, recall, ratio_max, bm * 1e3, bm / tm, bad2, sweep_mean * 1e3, sweep_mean / tm);
        js << (firstk ? "" : ",") << "{\"k\":" << k << ",\"approx_us_mean\":" << am * 1e6 << ",\"approx_us_p50\":" << pct(ta, 0.5) * 1e6
           << ",\"approx_recall\":" << ap_recall << ",\"approx_worst_ratio\":" << ap_ratio << ",\"approx_over_bound\":" << ap_bad << ",\"cert_us_mean\":" << tm * 1e6 << ",\"cert_global_us_mean\":" << gm * 1e6 << ",\"cert_us_p50\":" << pct(tc, 0.5) * 1e6 << ",\"cert_us_p99\":" << pct(tc, 0.99) * 1e6
           << ",\"candidates_mean\":" << cm << ",\"refined_queries\":" << refined_q << ",\"refined_objects_mean\":" << rm << ",\"exact_answers\":" << correct
           << ",\"est_recall\":" << recall << ",\"est_worst_ratio\":" << ratio_max << ",\"pruned_ms_mean\":" << bm * 1e3
           << ",\"pruned_ms_p50\":" << pct(tb2, 0.5) * 1e3 << ",\"sweep_ms_mean\":" << sweep_mean * 1e3 << "}";
        firstk = false;
        flush(0);
    }
    js << "],\"range\":[";
    // ---------------------------------------------------------------- range
    bool firstr = true;
    for (double rf : parse_list(rlist)) {
        const double r = rf * diag;
        std::vector<double> tc, tb2, tg_, ta;
        int correct = 0, refined_q = 0, ap_fp = 0, ap_inner_miss = 0;
        double res_mean = 0, prec = 0, rec = 0, amb_mean = 0, ap_rec = 0;
        for (int qi = 0; qi < Q; ++qi) {
            const int q = qv[qi];
            const float* Dq = &D[(size_t)qi * M];
            auto c0 = clk::now();
            std::vector<int> answer, amb_o;
            std::vector<double> amb_ub;
            std::vector<std::pair<int, double>> est_list;
            double o3[3];
            std::vector<int> ub_in; // certified approximate: UB <= r
            for (int j : og.radius(T.xy[2 * q], T.xy[2 * q + 1], r / lam)) {
                ix.bounds(q, obj[j], o3);
                est_list.push_back({j, o3[2]});
                if (o3[1] <= r) ub_in.push_back(j);
                if (o3[1] <= r) answer.push_back(j);
                else if (o3[0] > r) continue;
                else { amb_o.push_back(j); amb_ub.push_back(o3[1]); }
            }
            ta.push_back(since(c0));
            {
                int fp = 0, tp = 0, inner = 0, inner_hit = 0;
                for (int j : ub_in) { if (Dq[j] > r * (1 + 1e-6)) fp++; else tp++; }
                std::set<int> got(ub_in.begin(), ub_in.end());
                for (int j = 0; j < M; ++j) if (Dq[j] <= r / (1 + ix.eps)) { inner++; inner_hit += got.count(j) ? 1 : 0; }
                ap_fp += fp; ap_inner_miss += inner - inner_hit;
                int tn = 0; for (int j = 0; j < M; ++j) tn += Dq[j] <= r;
                ap_rec += (tn ? (double)tp / tn : 1.0) / Q;
            }
            double t_glob = 0, t_sub = 0;
            if (!amb_o.empty()) {
                std::vector<int> tg; for (int j : amb_o) tg.push_back(obj[j]);
                std::vector<double> ex(tg.size()), ex2(tg.size());
                auto r0 = clk::now();
                ix.exact_from(q, tg, *std::max_element(amb_ub.begin(), amb_ub.end()), ex.data());
                t_sub = since(r0);
                std::vector<unsigned> tgu(tg.begin(), tg.end());
                r0 = clk::now();
                mmp->dist((unsigned)q, tgu, ex2.data());
                t_glob = since(r0);
                for (size_t a = 0; a < amb_o.size(); ++a) if (ex[a] <= r) answer.push_back(amb_o[a]);
                refined_q++; amb_mean += (double)amb_o.size() / Q;
            }
            double t_all = since(c0) - t_glob;
            tc.push_back(t_all); tg_.push_back(t_all - t_sub + t_glob);
            {   // baseline: one early-stopping MMP sweep to the objects with lambda |q - o| <= r
                auto b0 = clk::now();
                std::vector<int> cj = og.radius(T.xy[2 * q], T.xy[2 * q + 1], r / lam);
                std::vector<unsigned> tg2; for (int j : cj) tg2.push_back((unsigned)obj[j]);
                std::vector<double> d2(tg2.size());
                if (!tg2.empty()) mmp->dist((unsigned)q, tg2, d2.data());
                int cnt = 0; for (double d : d2) cnt += d <= r;
                tb2.push_back(since(b0));
                (void)cnt;
            }
            // check against the truth (ties within 1e-6 relative are not counted)
            std::set<int> got(answer.begin(), answer.end());
            int truth_n = 0; bool ok = true;
            for (int j = 0; j < M; ++j) {
                bool in_t = Dq[j] <= r;
                truth_n += in_t;
                if (in_t != (bool)got.count(j) && std::fabs(Dq[j] / r - 1) > 1e-6) ok = false;
            }
            correct += ok; res_mean += (double)truth_n / Q;
            int tp = 0, ne = 0;
            for (auto& e : est_list) if (e.second <= r) { ne++; tp += Dq[e.first] <= r; }
            prec += (ne ? (double)tp / ne : 1.0) / Q; rec += (truth_n ? (double)tp / truth_n : 1.0) / Q;
        }
        double tm = 0; for (double x : tc) tm += x / Q;
        double bm = 0; for (double x : tb2) bm += x / Q;
        double gm = 0; for (double x : tg_) gm += x / Q;
        double am = 0; for (double x : ta) am += x / Q;
        std::printf("  range r=%.3g diag (%.0f m, %.1f objects on average): approximate (UB <= r) %.1f us mean, recall %.4f, %d false positives, %d missed within r/(1+eps) | certified exact %.1f us mean (%.1f us with a whole-mesh solver), p99 %.1f us | refined %d/%d queries (%.1f objects) | exact %d/%d | estimate-only: precision %.4f recall %.4f | on the fly: pruned MMP %.2f ms (%.0fx), full sweep %.2f ms\n",
                    rf, r, res_mean, am * 1e6, ap_rec, ap_fp, ap_inner_miss, tm * 1e6, gm * 1e6, pct(tc, 0.99) * 1e6, refined_q, Q, refined_q ? amb_mean * Q / refined_q : 0.0, correct, Q, prec, rec, bm * 1e3, bm / tm, sweep_mean * 1e3);
        js << (firstr ? "" : ",") << "{\"radius_frac\":" << rf << ",\"radius\":" << r << ",\"results_mean\":" << res_mean
           << ",\"approx_us_mean\":" << am * 1e6 << ",\"approx_recall\":" << ap_rec << ",\"approx_false_pos\":" << ap_fp << ",\"approx_inner_missed\":" << ap_inner_miss
           << ",\"cert_us_mean\":" << tm * 1e6 << ",\"cert_global_us_mean\":" << gm * 1e6
           << ",\"cert_us_p99\":" << pct(tc, 0.99) * 1e6 << ",\"refined_queries\":" << refined_q << ",\"refined_objects_mean\":" << amb_mean
           << ",\"exact_answers\":" << correct << ",\"est_precision\":" << prec << ",\"est_recall\":" << rec
           << ",\"pruned_ms_mean\":" << bm * 1e3 << ",\"sweep_ms_mean\":" << sweep_mean * 1e3 << "}";
        firstr = false;
        flush(1);
    }
    js << "],\"complete\":true}";
    flush(2);
    return 0;
}
