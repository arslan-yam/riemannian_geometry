// Build the Certified Cut Index on a mesh and evaluate it on vertex pairs against exact MMP.
// usage: cci_bench mesh.off [--eps 0.1] [--gap-factor 0.5] [--solver mmp|vtp|hybrid] [--ratio 2]
//                  [--pairs file | --random N] [--seed 7] [--no-truth] [--json out.json] [--verbose]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include "../src/certified.h"
#include "../src/mesh_io.h"
using namespace cci;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: cci_bench mesh.off [options]\n"); return 1; }
    std::string mesh = argv[1], pairs_file, json;
    Options opt;
    int nrand = 1000, seed = 7;
    bool truth = true;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--eps") opt.eps = std::stod(nx());
        else if (a == "--gap-factor") opt.gap_factor = std::stod(nx());
        else if (a == "--solver") opt.solver = nx();
        else if (a == "--ratio") opt.hybrid_ratio = std::stod(nx());
        else if (a == "--pairs") pairs_file = nx();
        else if (a == "--random") nrand = std::stoi(nx());
        else if (a == "--seed") seed = std::stoi(nx());
        else if (a == "--no-truth") truth = false;
        else if (a == "--json") json = nx();
        else if (a == "--verbose") opt.verbose = true;
        else if (a == "--interp") opt.interp = std::stoi(nx());
    }
    std::vector<double> V; std::vector<unsigned> F;
    read_off(mesh, V, F);
    Terrain T; T.init(V, F);
    std::vector<std::pair<int, int>> qp;
    if (!pairs_file.empty()) {
        std::ifstream in(pairs_file); int s, t;
        while (in >> s >> t) qp.push_back({s, t});
    } else {
        std::mt19937_64 rng(seed);
        while ((int)qp.size() < nrand) { int s = rng() % T.n, t = rng() % T.n; if (s != t) qp.push_back({s, t}); }
    }
    auto t0 = clk::now();
    CertifiedIndex ix(T, opt);
    double tb = std::chrono::duration<double>(clk::now() - t0).count();
    size_t m = qp.size();
    std::vector<double> B(3 * m);
    for (size_t i = 0; i < m; ++i) ix.bounds(qp[i].first, qp[i].second, &B[3 * i]);
    // query latency: best of 5 passes over all pairs
    double best = 1e30, sink = 0;
    for (int rep = 0; rep < 5; ++rep) {
        auto q0 = clk::now();
        double o[3];
        for (auto& p : qp) { ix.bounds(p.first, p.second, o); sink += o[2]; }
        best = std::min(best, std::chrono::duration<double>(clk::now() - q0).count() / m);
    }
    std::vector<double> tr(m, NAN);
    double ttruth = 0;
    if (truth) { // exact MMP on the whole mesh, one early-stopped sweep per distinct source
        auto t1 = clk::now();
        std::unique_ptr<Solver> sol(make_mmp(T.V, T.F, true));
        std::map<int, std::vector<size_t>> by;
        for (size_t i = 0; i < m; ++i) by[qp[i].first].push_back(i);
        for (auto& kv : by) {
            std::vector<unsigned> tg; for (size_t i : kv.second) tg.push_back(qp[i].second);
            std::vector<double> d(tg.size());
            sol->dist(kv.first, tg, d.data());
            for (size_t k = 0; k < tg.size(); ++k) tr[kv.second[k]] = d[k];
        }
        ttruth = std::chrono::duration<double>(clk::now() - t1).count();
    }
    double emean = 0, emax = 0, ubmax = 0, ubmean = 0, lbmin = 1e30, wmed = 0;
    bool sound = true;
    std::vector<double> errs, widths;
    for (size_t i = 0; i < m; ++i) {
        double lb = B[3 * i], ub = B[3 * i + 1], est = B[3 * i + 2];
        widths.push_back(lb > 0 ? (ub - lb) / lb : 0);
        if (!truth) continue;
        double d = tr[i];
        if (!(lb <= d && d <= ub)) sound = false;
        double e = std::fabs(est / d - 1);
        errs.push_back(e); emean += e; emax = std::max(emax, e);
        ubmax = std::max(ubmax, ub / d - 1); ubmean += ub / d - 1; lbmin = std::min(lbmin, lb / d);
    }
    std::sort(widths.begin(), widths.end()); wmed = widths[widths.size() / 2];
    double p99 = 0;
    if (truth) { emean /= m; ubmean /= m; std::sort(errs.begin(), errs.end()); p99 = errs[(size_t)(0.99 * (m - 1))]; }
    std::printf("%s eps=%.2f g=%.2fe solver=%s: build %.2f s (labels %.2f, near %.2f; sweeps MMP %lld VTP %lld) size %.2f MB "
                "cuts %zu portals %lld | query %.3f us",
                mesh.c_str(), opt.eps, opt.gap_factor, opt.solver.c_str(), tb, ix.timing["labels"], ix.timing["near_tables"],
                (long long)ix.sweeps_mmp, (long long)ix.sweeps_vtp, ix.nbytes() / 1e6, ix.cuts.size(), (long long)ix.portals(), best * 1e6);
    if (truth)
        std::printf(" | %zu pairs: sound %s, max UB err %.3f%%, UB mean %.4f%%, est mean %.4f%% p99 %.3f%% max %.3f%%, min LB/d %.4f, width median %.2f%% (truth %.1f s)",
                    m, sound ? "yes" : "NO", 100 * ubmax, 100 * ubmean, 100 * emean, 100 * p99, 100 * emax, lbmin, 100 * wmed, ttruth);
    std::printf("\n");
    if (!json.empty()) {
        std::ofstream o(json);
        o << "{\"mesh\":\"" << mesh << "\",\"n\":" << T.n << ",\"faces\":" << T.m << ",\"eps\":" << opt.eps << ",\"gap_factor\":" << opt.gap_factor
          << ",\"solver\":\"" << opt.solver << "\",\"build_s\":" << tb << ",\"index_bytes\":" << ix.nbytes() << ",\"query_us\":" << best * 1e6
          << ",\"cuts\":" << ix.cuts.size() << ",\"portals\":" << ix.portals() << ",\"labels\":" << ix.labels() << ",\"near_entries\":" << ix.near_ids.size()
          << ",\"sweeps_mmp\":" << ix.sweeps_mmp << ",\"sweeps_vtp\":" << ix.sweeps_vtp << ",\"timing\":{";
        bool first = true;
        for (auto& kv : ix.timing) { o << (first ? "" : ",") << "\"" << kv.first << "\":" << kv.second; first = false; }
        o << "},\"pairs\":" << m;
        if (truth) o << ",\"sound\":" << (sound ? "true" : "false") << ",\"max_ub_err\":" << ubmax << ",\"ub_mean\":" << ubmean << ",\"est_mean\":" << emean
                     << ",\"est_p99\":" << p99 << ",\"est_max\":" << emax << ",\"min_lb_ratio\":" << lbmin << ",\"width_median\":" << wmed;
        o << ",\"values\":[";
        for (size_t i = 0; i < m; ++i) o << (i ? "," : "") << "[" << qp[i].first << "," << qp[i].second << "," << B[3 * i] << "," << B[3 * i + 1] << "," << B[3 * i + 2] << "]";
        o << "]}\n";
    }
    return sink < 0; // keep the timing loop alive
}
