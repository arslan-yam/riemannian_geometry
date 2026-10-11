// Experiments on the size of the certified index: one build (kd or shifted hierarchy, optional portal
// cap k_max) evaluated with several storages (dense float32, 16-bit coding, eps_c-covers), with or
// without certify-or-refine, against exact distances.
// usage: cci_exp mesh.off [--eps 0.1] [--gap-factor 0.5] [--solver hybrid] [--hierarchy kd|shifted]
//                [--kmax K] [--refine] [--variants dense,code16,cover0.01] [--pairs f | --random N]
//                [--seed 7] [--truth-cache f] [--json out.json]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include "../src/certified.h"
#include "../src/mesh_io.h"
using namespace cci;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: cci_exp mesh.off [options]\n"); return 1; }
    std::string mesh = argv[1], pairs_file, json, truth_cache, variants = "dense";
    Options opt;
    int nrand = 1000, seed = 7;
    bool refine = false;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--eps") opt.eps = std::stod(nx());
        else if (a == "--gap-factor") opt.gap_factor = std::stod(nx());
        else if (a == "--solver") opt.solver = nx();
        else if (a == "--ratio") opt.hybrid_ratio = std::stod(nx());
        else if (a == "--hierarchy") opt.hierarchy = nx();
        else if (a == "--copies") opt.copies = std::stoi(nx());
        else if (a == "--kmax") opt.k_max = std::stoi(nx());
        else if (a == "--refine") refine = true;
        else if (a == "--variants") variants = nx();
        else if (a == "--pairs") pairs_file = nx();
        else if (a == "--random") nrand = std::stoi(nx());
        else if (a == "--seed") seed = std::stoi(nx());
        else if (a == "--truth-cache") truth_cache = nx();
        else if (a == "--json") json = nx();
        else if (a == "--verbose") opt.verbose = true;
    }
    std::vector<double> V; std::vector<unsigned> F;
    read_off(mesh, V, F);
    Terrain T; T.init(V, F);
    std::vector<std::pair<int, int>> qp;
    std::vector<double> given;
    if (!pairs_file.empty()) {
        std::ifstream in(pairs_file); std::string line;
        while (std::getline(in, line)) {
            std::istringstream ls(line); int s, t; double d;
            if (!(ls >> s >> t)) continue;
            qp.push_back({s, t});
            if (ls >> d) given.push_back(d);
        }
        if (given.size() != qp.size()) given.clear();
    } else {
        std::mt19937_64 rng(seed);
        while ((int)qp.size() < nrand) { int s = rng() % T.n, t = rng() % T.n; if (s != t) qp.push_back({s, t}); }
    }
    size_t m = qp.size();
    // truth (exact MMP, cached)
    std::vector<double> tr(m, NAN);
    if (!given.empty()) tr = given;
    else if (!truth_cache.empty() && std::ifstream(truth_cache).good()) { std::ifstream in(truth_cache); for (auto& d : tr) in >> d; }
    else {
        std::unique_ptr<Solver> sol(make_mmp(T.V, T.F, true));
        std::map<int, std::vector<size_t>> by;
        for (size_t i = 0; i < m; ++i) by[qp[i].first].push_back(i);
        for (auto& kv : by) {
            std::vector<unsigned> tg; for (size_t i : kv.second) tg.push_back(qp[i].second);
            std::vector<double> d(tg.size());
            sol->dist(kv.first, tg, d.data());
            for (size_t k = 0; k < tg.size(); ++k) tr[kv.second[k]] = d[k];
        }
        if (!truth_cache.empty()) { std::ofstream o(truth_cache); o.precision(17); for (double d : tr) o << d << "\n"; }
    }
    auto t0 = clk::now();
    CertifiedIndex ix(T, opt);
    double tb = since(t0);
    std::ostringstream js;
    js << "{\"mesh\":\"" << mesh << "\",\"n\":" << T.n << ",\"eps\":" << opt.eps << ",\"gap_factor\":" << opt.gap_factor
       << ",\"hierarchy\":\"" << opt.hierarchy << "\",\"k_max\":" << opt.k_max << ",\"solver\":\"" << opt.solver
       << "\",\"build_s\":" << tb << ",\"cuts\":" << ix.cuts.size() << ",\"portals\":" << ix.portals() << ",\"labels\":" << ix.labels()
       << ",\"near_entries\":" << ix.near_ids.size() << ",\"timing\":{";
    bool first = true;
    for (auto& kv : ix.timing) { js << (first ? "" : ",") << "\"" << kv.first << "\":" << kv.second; first = false; }
    js << "},\"variants\":[";
    std::stringstream vs(variants);
    std::string var;
    bool firstv = true;
    while (std::getline(vs, var, ',')) {
        auto tv = clk::now();
        if (var == "dense") ix.storage = DENSE;
        else if (var == "code16") { ix.encode16(16); ix.storage = CODE16; }
        else if (var == "code8") { ix.encode16(8); ix.storage = CODE8; }
        else if (var.rfind("cover", 0) == 0) { ix.build_covers(std::stod(var.substr(5))); ix.storage = COVER; }
        else { std::fprintf(stderr, "unknown variant %s\n", var.c_str()); continue; }
        double tprep = since(tv);
        // latency of the plain certified query: best of 5 passes
        double best = 1e30, sink = 0, o[3];
        for (int rep = 0; rep < 5; ++rep) {
            auto q0 = clk::now();
            for (auto& p : qp) { ix.bounds(p.first, p.second, o); sink += o[2]; }
            best = std::min(best, since(q0) / m);
        }
        // answers (with certify-or-refine if asked), one pass, timed per query
        std::vector<double> B(3 * m), qt(m);
        int refined = 0;
        for (size_t i = 0; i < m; ++i) {
            auto q0 = clk::now();
            bool r = refine ? ix.bounds_refined(qp[i].first, qp[i].second, &B[3 * i]) : (ix.bounds(qp[i].first, qp[i].second, &B[3 * i]), false);
            qt[i] = since(q0);
            refined += r;
        }
        bool sound = true;
        double emean = 0, emax = 0, ubmax = 0, ubmean = 0, worst_ratio = 0;
        std::vector<double> errs, widths;
        for (size_t i = 0; i < m; ++i) {
            double lb = B[3 * i], ub = B[3 * i + 1], est = B[3 * i + 2], d = tr[i];
            if (!(lb <= d * (1 + 1e-9) && d <= ub * (1 + 1e-9))) sound = false;
            double e = std::fabs(est / d - 1);
            errs.push_back(e); emean += e; emax = std::max(emax, e);
            ubmax = std::max(ubmax, ub / d - 1); ubmean += ub / d - 1;
            worst_ratio = std::max(worst_ratio, ub / std::max(lb, 1e-300));
            widths.push_back(lb > 0 ? (ub - lb) / lb : 0);
        }
        emean /= m; ubmean /= m;
        std::sort(errs.begin(), errs.end()); std::sort(widths.begin(), widths.end());
        std::vector<double> qs(qt); std::sort(qs.begin(), qs.end());
        double qmean = 0; for (double x : qt) qmean += x; qmean /= m;
        double p99 = errs[(size_t)(0.99 * (m - 1))];
        std::printf("%s eps=%.2f g=%.2fe %s kmax=%d %s%s: build %.1f s size %.2f MB (prep %.1f s) | query %.3f us, with refine mean %.2f us p99 %.1f us max %.0f us, refined %d/%zu | sound %s, max UB err %.3f%%, est mean %.4f%% p99 %.3f%% max %.3f%%, width median %.2f%%, max UB/LB %.4f\n",
                    mesh.c_str(), opt.eps, opt.gap_factor, opt.hierarchy.c_str(), opt.k_max, var.c_str(), refine ? "+refine" : "",
                    tb, ix.nbytes() / 1e6, tprep, best * 1e6, qmean * 1e6, qs[(size_t)(0.99 * (m - 1))] * 1e6, qs.back() * 1e6,
                    refined, m, sound ? "yes" : "NO", 100 * ubmax, 100 * emean, 100 * p99, 100 * emax, 100 * widths[m / 2], worst_ratio);
        std::fflush(stdout);
        js << (firstv ? "" : ",") << "{\"variant\":\"" << var << "\",\"refine\":" << (refine ? "true" : "false")
           << ",\"index_bytes\":" << ix.nbytes() << ",\"prep_s\":" << tprep << ",\"query_us\":" << best * 1e6
           << ",\"query_us_refine_mean\":" << qmean * 1e6 << ",\"query_us_refine_p99\":" << qs[(size_t)(0.99 * (m - 1))] * 1e6
           << ",\"refined\":" << refined << ",\"pairs\":" << m << ",\"sound\":" << (sound ? "true" : "false")
           << ",\"max_ub_err\":" << ubmax << ",\"ub_mean\":" << ubmean << ",\"est_mean\":" << emean << ",\"est_p99\":" << p99
           << ",\"est_max\":" << emax << ",\"width_median\":" << widths[m / 2] << ",\"max_ub_over_lb\":" << worst_ratio
           << ",\"cover_entries\":" << ix.cover_entries << "}";
        firstv = false;
        (void)sink;
    }
    js << "]}\n";
    if (!json.empty()) { std::ofstream f(json); f << js.str(); }
    return 0;
}
