// exactness and speed of the exact solvers: MMP full sweep, MMP stopped at the targets, VTP full sweep
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <memory>
#include <random>
#include "../src/mesh_io.h"
#include "../src/solver.h"
using namespace cci;
using clk = std::chrono::steady_clock;

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: solver_check mesh.off [sources=5] [--json out.json]\n"); return 1; }
    std::string json;
    for (int i = 2; i + 1 < argc; ++i) if (std::string(argv[i]) == "--json") json = argv[i + 1];
    std::vector<double> V; std::vector<unsigned> F;
    read_off(argv[1], V, F);
    int ns = argc > 2 && std::string(argv[2]) != "--json" ? std::atoi(argv[2]) : 5;
    unsigned n = V.size() / 3;
    std::vector<unsigned> all(n); for (unsigned i = 0; i < n; ++i) all[i] = i;
    std::mt19937 rng(1);
    std::vector<unsigned> src(ns); for (auto& s : src) s = rng() % n;
    // nearby targets: the 2% of vertices closest (chart) to the source -- the early-stop pattern
    auto t0 = clk::now();
    std::unique_ptr<Solver> mmp(make_mmp(V, F, false)), mmpe(make_mmp(V, F, true)), vtp(make_vtp(V, F));
    double tinit = std::chrono::duration<double>(clk::now() - t0).count();
    std::vector<double> a(n), b(n), c(n);
    double ta = 0, tb = 0, tc = 0, tnear_e = 0, tnear_v = 0, maxdiff_vtp = 0, maxdiff_early = 0, maxrel = 0;
    for (unsigned s : src) {
        auto t1 = clk::now(); mmp->dist(s, all, a.data()); auto t2 = clk::now();
        mmpe->dist(s, all, b.data()); auto t3 = clk::now();
        vtp->dist(s, all, c.data()); auto t4 = clk::now();
        ta += std::chrono::duration<double>(t2 - t1).count();
        tb += std::chrono::duration<double>(t3 - t2).count();
        tc += std::chrono::duration<double>(t4 - t3).count();
        for (unsigned i = 0; i < n; ++i) {
            maxdiff_early = std::max(maxdiff_early, std::fabs(a[i] - b[i]));
            maxdiff_vtp = std::max(maxdiff_vtp, std::fabs(a[i] - c[i]));
            if (a[i] > 0) maxrel = std::max(maxrel, std::fabs(a[i] - c[i]) / a[i]);
        }
        std::vector<std::pair<double, unsigned>> d2(n);
        for (unsigned i = 0; i < n; ++i) {
            double dx = V[3 * i] - V[3 * s], dy = V[3 * i + 1] - V[3 * s + 1];
            d2[i] = {dx * dx + dy * dy, i};
        }
        unsigned k = std::max(10u, n / 50);
        std::partial_sort(d2.begin(), d2.begin() + k, d2.end());
        std::vector<unsigned> near(k); for (unsigned i = 0; i < k; ++i) near[i] = d2[i].second;
        std::vector<double> o(k);
        auto t5 = clk::now(); mmpe->dist(s, near, o.data()); auto t6 = clk::now(); vtp->dist(s, near, o.data()); auto t7 = clk::now();
        tnear_e += std::chrono::duration<double>(t6 - t5).count();
        tnear_v += std::chrono::duration<double>(t7 - t6).count();
    }
    std::printf("%s: n=%u faces=%zu | setup %.0f ms | full sweep: MMP %.1f ms, MMP-early(all targets) %.1f ms, VTP %.1f ms "
                "(VTP speedup %.2fx) | nearest 2%%: MMP-early %.2f ms, VTP %.1f ms | max |MMP-VTP| %.2e (rel %.2e), max |MMP-early| %.2e\n",
                argv[1], n, F.size() / 3, 1e3 * tinit, 1e3 * ta / ns, 1e3 * tb / ns, 1e3 * tc / ns, ta / tc,
                1e3 * tnear_e / ns, 1e3 * tnear_v / ns, maxdiff_vtp, maxrel, maxdiff_early);
    if (!json.empty()) {
        std::ofstream o(json);
        o << "{\"tool\":\"solver_check\",\"mesh\":\"" << argv[1] << "\",\"n\":" << n << ",\"faces\":" << F.size() / 3 << ",\"sources\":" << ns
          << ",\"setup_ms\":" << 1e3 * tinit << ",\"mmp_full_ms\":" << 1e3 * ta / ns << ",\"mmp_early_all_ms\":" << 1e3 * tb / ns
          << ",\"vtp_full_ms\":" << 1e3 * tc / ns << ",\"mmp_early_near2pct_ms\":" << 1e3 * tnear_e / ns << ",\"vtp_near2pct_ms\":" << 1e3 * tnear_v / ns
          << ",\"max_rel_mmp_vtp\":" << maxrel << ",\"max_abs_mmp_early\":" << maxdiff_early << "}\n";
    }
}
