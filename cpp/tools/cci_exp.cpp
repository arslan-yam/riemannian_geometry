// Experiments on the size of the certified index: one build (kd or shifted hierarchy, optional portal
// cap k_max) evaluated with several storages (dense float32, 16-bit coding, eps_c-covers), with or
// without certify-or-refine, against exact distances.
// usage: cci_exp mesh.off [--eps 0.1] [--gap-factor 0.5] [--solver hybrid] [--hierarchy kd|shifted]
//                [--kmax K] [--refine] [--variants dense,code16,cover0.01] [--pairs f | --random N]
//                [--seed 7] [--truth-cache f] [--json out.json]
//                [--threads 1] [--store dense|code16|code8]   (store: coded during the build; the only variant)
//                [--poi-random P [--poi-seed 11] | --poi-file f] [--poi-table] [--truth-only]
//                [--near-store sym] [--near-variants dense,sym,sym16,sym8] [--near-pairs N [--near-truth-cache f]]
// --random 0 with --near-pairs N evaluates the near pairs only.
//                [--sources S --per-source 10] [--buckets k --bucket-sources 100] [--a2a-sources S --a2a-per-source 5]
//                [--export prefix]   (prefix.{V2V,A2A}.{query,truth} for the competitors' harness, prefix.V2V.detour)
//                [--detour-sets]   (the random pairs again, by detour ratio d / |s - t|: sets detour0..3)
//                [--exact vtp]   (VTP for the truth and every exact solve of the index: meshes with holes)
//                [--weights f --solver steiner --steiner-m 5 --truth-steiner 24]   (weighted terrain)
// POI mode indexes only the POIs and draws the query pairs among them; --poi-table also builds the plain
// exact baseline (one full VTP sweep per POI, a P x P float32 table) and times it.
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
#include <sys/resource.h>
#include "../src/certified.h"
#include "../src/mesh_io.h"
using namespace cci;
using clk = std::chrono::steady_clock;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }
static double peak_rss_mb() {
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return ru.ru_maxrss / 1e6;  // bytes
#else
    return ru.ru_maxrss / 1e3;  // kilobytes
#endif
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: cci_exp mesh.off [options]\n"); return 1; }
    std::string mesh = argv[1], pairs_file, json, truth_cache, variants = "dense";
    Options opt;
    int nrand = 1000, seed = 7, poi_random = 0, poi_seed = 11, near_pairs = 0;
    int n_sources = 0, per_source = 10, buckets = 0, bucket_sources = 100, a2a_sources = 0, a2a_per_source = 5;
    std::string poi_file, near_variants, near_truth_cache, export_prefix, interp_variants, weights_file;
    int truth_steiner = 0;
    bool refine = false, poi_table = false, truth_only = false, detour_sets = false;
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
        else if (a == "--threads") opt.threads = std::stoi(nx());
        else if (a == "--solver-mem-gb") opt.solver_mem_mb = 1e3 * std::stod(nx());
        else if (a == "--store") opt.store = nx();
        else if (a == "--near-store") opt.near_store = nx();
        else if (a == "--near-variants") near_variants = nx();  // dense,sym,sym16: evaluated on one build
        else if (a == "--near-pairs") near_pairs = std::stoi(nx()); // + a set of pairs closer than the near radius
        else if (a == "--near-truth-cache") near_truth_cache = nx();
        else if (a == "--sources") n_sources = std::stoi(nx());       // random set = sources x --per-source targets
        else if (a == "--per-source") per_source = std::stoi(nx());
        else if (a == "--buckets") buckets = std::stoi(nx());         // targets per distance band and source
        else if (a == "--bucket-sources") bucket_sources = std::stoi(nx());
        else if (a == "--a2a-sources") a2a_sources = std::stoi(nx()); // A2A set = points x --a2a-per-source
        else if (a == "--a2a-per-source") a2a_per_source = std::stoi(nx());
        else if (a == "--export") export_prefix = nx();               // harness query files + exact truths
        else if (a == "--detour-sets") detour_sets = true;             // random pairs by detour ratio (obstacles)
        else if (a == "--exact") opt.exact = nx();                    // mmp | vtp: the exact solver (truth and index)
        else if (a == "--interp-variants") interp_variants = nx();     // 0,1,2: estimate interpolation, one build
        else if (a == "--weights") weights_file = nx();               // face weights, one per line (weighted terrain)
        else if (a == "--steiner-m") opt.steiner_m = std::stoi(nx()); // weighted: Steiner points per edge of the labels
        else if (a == "--truth-steiner") truth_steiner = std::stoi(nx()); // weighted: reference = Steiner graph, M points per edge
        else if (a == "--poi-random") poi_random = std::stoi(nx());
        else if (a == "--poi-seed") poi_seed = std::stoi(nx());
        else if (a == "--poi-file") poi_file = nx();
        else if (a == "--poi-table") poi_table = true;
        else if (a == "--truth-only") truth_only = true;   // compute (and cache) the exact distances, then exit
    }
    if (opt.store != "dense") variants = opt.store; // the dense labels are gone after the build
    std::vector<double> V; std::vector<unsigned> F;
    double origin[2] = {0, 0};
    read_off(mesh, V, F, origin);
    std::vector<double> weights;
    if (!weights_file.empty()) {
        std::ifstream in(weights_file); double x;
        while (in >> x) weights.push_back(x);
        if (weights.size() != F.size() / 3) { std::fprintf(stderr, "%zu weights for %zu faces\n", weights.size(), F.size() / 3); return 1; }
    }
    Terrain T; T.init(V, F, weights.empty() ? nullptr : &weights);
    const bool weighted = !weights.empty();
    if (weighted && truth_steiner <= 0) { std::fprintf(stderr, "weighted terrain: no exact solver; give --truth-steiner M (reference Steiner graph)\n"); return 1; }
    // the solver of the truth: exact MMP (early stop), or on a weighted terrain a fine Steiner graph (a
    // reference: graph distances are upper bounds that converge to the weighted distance)
    // (--exact vtp: full VTP sweeps; Kirsanov's MMP returns too-long distances on some meshes with holes)
    auto truth_solver = [&]() -> Solver* {
        if (truth_steiner > 0) return make_steiner(T.V, T.F, T.w, truth_steiner);
        if (opt.exact == "vtp") { if (Solver* s = make_vtp(T.V, T.F, true)) return s; std::fprintf(stderr, "VTP rejects the mesh: MMP truth\n"); }
        return make_mmp(T.V, T.F, true);
    };
    int truth_fallbacks = 0; // exact = vtp: A2A truths that fell back to MMP
    if (poi_random > 0) {
        std::mt19937_64 rng(poi_seed);
        std::vector<int> all(T.n);
        for (int i = 0; i < T.n; ++i) all[i] = i;
        std::shuffle(all.begin(), all.end(), rng);
        opt.poi.assign(all.begin(), all.begin() + std::min(poi_random, T.n));
    } else if (!poi_file.empty()) {
        std::ifstream in(poi_file); int v;
        while (in >> v) opt.poi.push_back(v);
        // UP-Oracle's POI files start with the count
        if (opt.poi.size() > 1 && opt.poi[0] == (int)opt.poi.size() - 1) opt.poi.erase(opt.poi.begin());
    }
    std::sort(opt.poi.begin(), opt.poi.end());
    // query sets (vertex pairs):
    //   "random": a pair file, --random N independent pairs, or --sources S x --per-source k (one exact sweep
    //             per source for the truth); among the POIs in POI mode
    //   "near":   --near-pairs N pairs closer than the near radius (answered by the near tables)
    //   "d<i>":   --buckets k: k targets per distance band and source (the first --bucket-sources sources),
    //             bands by chart distance [0, R], (R, 2R], (2R, 4R], ... up to the diameter, R the near radius
    // and one set of arbitrary surface points, --a2a-sources S x --a2a-per-source k ("a2a")
    struct PairSet { std::string name, cache; std::vector<std::pair<int, int>> qp; std::vector<double> tr; double lo = -1, hi = -1; };
    std::vector<PairSet> sets;
    const double g_ = opt.gap > 0 ? opt.gap : opt.gap_factor * T.edge, Rnear = g_ / (opt.eps * T.lam) * (1.0 + 1e-9);
    std::vector<int> srcs; // the random sources (sources mode), reused by the bands
    // the near pairs and the bands depend on the near radius: it is part of their cache names
    char rtag[32]; std::snprintf(rtag, sizeof rtag, ".R%.3g", Rnear / T.edge);
    {
        PairSet R; R.name = "random"; R.cache = truth_cache;
        std::vector<double> given;
        std::mt19937_64 rng(seed);
        const bool poi = !opt.poi.empty();
        auto pick = [&]() { return poi ? opt.poi[rng() % opt.poi.size()] : (int)(rng() % T.n); };
        if (!pairs_file.empty()) {
            std::ifstream in(pairs_file); std::string line;
            while (std::getline(in, line)) {
                std::istringstream ls(line); int s, t; double d;
                if (!(ls >> s >> t)) continue;
                R.qp.push_back({s, t});
                if (ls >> d) given.push_back(d);
            }
            if (given.size() == R.qp.size()) R.tr = given;
        } else if (n_sources > 0) {
            for (int i = 0; i < n_sources; ++i) {
                int s = pick();
                srcs.push_back(s);
                for (int k = 0; k < per_source; ++k) { int t = pick(); if (t != s) R.qp.push_back({s, t}); }
            }
        } else {
            while ((int)R.qp.size() < nrand) { int s = pick(), t = pick(); if (s != t) R.qp.push_back({s, t}); }
        }
        if (!R.qp.empty()) sets.push_back(std::move(R));
    }
    PointGrid pg; pg.init(T.xy, std::max(Rnear, T.edge));
    if (near_pairs > 0) {
        PairSet N; N.name = "near";
        N.cache = !near_truth_cache.empty() ? near_truth_cache : (truth_cache.empty() ? "" : truth_cache + ".near" + std::to_string(near_pairs) + rtag);
        std::mt19937_64 rng(seed + 1);
        const bool poi = !opt.poi.empty();
        auto is_poi = [&](int v) { return std::binary_search(opt.poi.begin(), opt.poi.end(), v); };
        int tries = 0;
        while ((int)N.qp.size() < near_pairs && ++tries < 100 * near_pairs + 1000) {
            int s = poi ? opt.poi[rng() % opt.poi.size()] : (int)(rng() % T.n);
            std::vector<int> nb;
            for (int t : pg.radius(T.xy[2 * s], T.xy[2 * s + 1], Rnear)) if (t != s && (!poi || is_poi(t))) nb.push_back(t);
            if (nb.empty()) continue;
            N.qp.push_back({s, nb[rng() % nb.size()]});
        }
        if (!N.qp.empty()) sets.push_back(std::move(N));
    }
    if (buckets > 0 && opt.poi.empty()) {
        // band sources: the first --bucket-sources random sources (or fresh ones), k targets per band
        std::mt19937_64 rng(seed + 2);
        std::vector<int> bs(srcs.begin(), srcs.begin() + std::min<size_t>(srcs.size(), bucket_sources));
        while ((int)bs.size() < bucket_sources) bs.push_back((int)(rng() % T.n));
        const double diam = std::hypot(T.hi[0] - T.lo[0], T.hi[1] - T.lo[1]);
        std::vector<std::pair<double, double>> bands{{0.0, Rnear}};
        for (double lo = Rnear; lo < diam; lo *= 2) bands.push_back({lo, std::min(2 * lo, diam)});
        for (size_t bi = 0; bi < bands.size(); ++bi) {
            PairSet B; B.name = "d" + std::to_string(bi); B.lo = bands[bi].first; B.hi = bands[bi].second;
            B.cache = truth_cache.empty() ? "" : truth_cache + ".b" + std::to_string(bucket_sources) + "x" + std::to_string(buckets) + rtag + "." + B.name;
            for (int s : bs) {
                const double sx = T.xy[2 * s], sy = T.xy[2 * s + 1];
                auto in_band = [&](int t) { double d = std::hypot(T.xy[2 * t] - sx, T.xy[2 * t + 1] - sy); return t != s && d > B.lo && d <= B.hi; };
                std::vector<int> cand;
                if (B.hi <= diam / 8) { for (int t : pg.radius(sx, sy, B.hi)) if (in_band(t)) cand.push_back(t); }
                int got = 0;
                for (int tries = 0; got < buckets && tries < 4000; ++tries) {
                    int t;
                    if (!cand.empty()) t = cand[rng() % cand.size()];
                    else { t = (int)(rng() % T.n); if (!in_band(t)) continue; }
                    B.qp.push_back({s, t}); ++got;
                }
            }
            if (!B.qp.empty()) sets.push_back(std::move(B));
        }
    }
    // truth for every vertex set without a cache: one exact MMP sweep per distinct source, all sets at once
    {
        std::map<int, std::vector<std::pair<size_t, size_t>>> by;
        for (size_t si = 0; si < sets.size(); ++si) {
            auto& S = sets[si];
            size_t m = S.qp.size();
            if (S.tr.size() == m) continue;
            S.tr.assign(m, NAN);
            if (!S.cache.empty() && std::ifstream(S.cache).good()) { std::ifstream in(S.cache); for (auto& d : S.tr) in >> d; continue; }
            for (size_t i = 0; i < m; ++i) by[S.qp[i].first].push_back({si, i});
        }
        if (!by.empty()) {
            std::unique_ptr<Solver> sol(truth_solver());
            for (auto& kv : by) {
                std::vector<unsigned> tg;
                for (auto& e : kv.second) tg.push_back(sets[e.first].qp[e.second].second);
                std::vector<double> d(tg.size());
                sol->dist(kv.first, tg, d.data());
                for (size_t k = 0; k < tg.size(); ++k) sets[kv.second[k].first].tr[kv.second[k].second] = d[k];
            }
            for (auto& S : sets)
                if (!S.cache.empty() && !std::ifstream(S.cache).good()) { std::ofstream o(S.cache); o.precision(17); for (double d : S.tr) o << d << "\n"; }
        }
    }
    // the detour ratio of a pair: d over the 3-d straight line (a lower bound of d, so >= 1); relief alone
    // keeps it near 1, obstacles push it up
    auto detour = [&](int s, int t, double d) {
        double e = 0;
        for (int j = 0; j < 3; ++j) { double u = T.V[3 * s + j] - T.V[3 * t + j]; e += u * u; }
        return d / std::max(std::sqrt(e), 1e-12);
    };
    // --detour-sets: the random pairs again, split by detour ratio [1, 1.01), [1.01, 1.05), [1.05, 1.2), >= 1.2
    if (detour_sets && !sets.empty() && sets[0].name == "random") {
        const PairSet R = sets[0];
        const double cut[5] = {1.0, 1.01, 1.05, 1.2, 1e9};
        for (int b = 0; b < 4; ++b) {
            PairSet D; D.name = "detour" + std::to_string(b); D.lo = cut[b]; D.hi = cut[b + 1];
            for (size_t i = 0; i < R.qp.size(); ++i) {
                double r = detour(R.qp[i].first, R.qp[i].second, R.tr[i]);
                if (r >= D.lo && r < D.hi) { D.qp.push_back(R.qp[i]); D.tr.push_back(R.tr[i]); }
            }
            if (!D.qp.empty()) sets.push_back(std::move(D));
        }
    }
    // arbitrary surface points: uniform random face, uniform barycentric coordinates (as the competitors'
    // generator); truth by MMP from the source point on the whole mesh, early stop at its targets
    std::vector<std::pair<APoint, APoint>> aq;
    std::vector<double> atr;
    if (a2a_sources > 0) {
        std::mt19937_64 rng(seed + 3);
        std::uniform_real_distribution<double> U01(0.0, 1.0);
        auto rand_point = [&]() {
            int f = (int)(rng() % T.m);
            double r1 = std::sqrt(U01(rng)), r2 = U01(rng), l[3] = {1 - r1, r1 * (1 - r2), r1 * r2};
            double x = 0, y = 0;
            for (int k = 0; k < 3; ++k) { x += l[k] * T.xy[2 * T.F[3 * f + k]]; y += l[k] * T.xy[2 * T.F[3 * f + k] + 1]; }
            return APoint{f, x, y};
        };
        for (int i = 0; i < a2a_sources; ++i) { APoint p = rand_point(); for (int k = 0; k < a2a_per_source; ++k) aq.push_back({p, rand_point()}); }
        std::string cache = truth_cache.empty() ? "" : truth_cache + ".a2a" + std::to_string(a2a_sources) + "x" + std::to_string(a2a_per_source);
        atr.assign(aq.size(), NAN);
        if (!cache.empty() && std::ifstream(cache).good()) { std::ifstream in(cache); for (auto& d : atr) in >> d; }
        else {
            std::unique_ptr<Solver> sol(truth_solver());
            for (size_t i = 0; i < aq.size(); i += a2a_per_source) {
                const APoint& p = aq[i].first;
                SurfPt sp{(unsigned)p.f, p.x, p.y, T.z_at(p.f, p.x, p.y)};
                std::vector<SurfPt> tg;
                for (int k = 0; k < a2a_per_source; ++k) { const APoint& q = aq[i + k].second; tg.push_back({(unsigned)q.f, q.x, q.y, T.z_at(q.f, q.x, q.y)}); }
                if (!sol->point_dists(sp, tg, &atr[i])) { // VTP rejected the mesh with the points inserted
                    ++truth_fallbacks;
                    std::unique_ptr<Solver> m(make_mmp(T.V, T.F, true));
                    m->point_dists(sp, tg, &atr[i]);
                }
            }
            if (!cache.empty()) { std::ofstream o(cache); o.precision(17); for (double d : atr) o << d << "\n"; }
        }
    }
    // the same queries for the competitors' harness ("x y z face x y z face", original coordinates; a vertex
    // query is moved 5 cm into its first face, as the harness cannot locate a point exactly at a vertex)
    if (!export_prefix.empty()) {
        std::vector<int> first_face(T.n, -1);
        for (int f = 0; f < T.m; ++f) for (int k = 0; k < 3; ++k) if (first_face[T.F[3 * f + k]] < 0) first_face[T.F[3 * f + k]] = f;
        auto vpoint = [&](int v, double* o) {
            int f = first_face[v];
            double c[3] = {0, 0, 0};
            for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j) c[j] += T.V[3 * T.F[3 * f + k] + j] / 3;
            double u[3], nu = 0;
            for (int j = 0; j < 3; ++j) { u[j] = c[j] - T.V[3 * v + j]; nu += u[j] * u[j]; }
            nu = std::max(std::sqrt(nu), 1e-12);
            for (int j = 0; j < 3; ++j) o[j] = T.V[3 * v + j] + u[j] * (0.05 / nu);
            return f;
        };
        auto put = [&](std::ofstream& o, const double* p, int f) {
            char buf[160];
            std::snprintf(buf, sizeof buf, "%.6f %.6f %.6f %d", p[0] + origin[0], p[1] + origin[1], p[2], f);
            o << buf;
        };
        if (!sets.empty() && sets[0].name == "random") {
            std::ofstream q(export_prefix + ".V2V.query"), t(export_prefix + ".V2V.truth"), pr(export_prefix + ".V2V.pairs"),
                dt(export_prefix + ".V2V.detour");
            t.precision(17);
            for (size_t i = 0; i < sets[0].qp.size(); ++i) {
                pr << sets[0].qp[i].first << " " << sets[0].qp[i].second << "\n"; // vertex pairs (DGG)
                dt << detour(sets[0].qp[i].first, sets[0].qp[i].second, sets[0].tr[i]) << "\n"; // per pair (bucketed reports)
                double a[3], b[3];
                int fa = vpoint(sets[0].qp[i].first, a), fb = vpoint(sets[0].qp[i].second, b);
                put(q, a, fa); q << " "; put(q, b, fb); q << "\n";
                t << sets[0].tr[i] << "\n";
            }
        }
        if (!aq.empty()) {
            std::ofstream q(export_prefix + ".A2A.query"), t(export_prefix + ".A2A.truth");
            t.precision(17);
            for (size_t i = 0; i < aq.size(); ++i) {
                const APoint &p = aq[i].first, &r = aq[i].second;
                double a[3] = {p.x, p.y, T.z_at(p.f, p.x, p.y)}, b[3] = {r.x, r.y, T.z_at(r.f, r.x, r.y)};
                put(q, a, p.f); q << " "; put(q, b, r.f); q << "\n";
                t << atr[i] << "\n";
            }
        }
    }
    if (truth_only) {
        for (auto& S : sets) std::printf("%s: %zu exact distances (%s)%s\n", mesh.c_str(), S.qp.size(), S.name.c_str(), S.cache.empty() ? "" : (" -> " + S.cache).c_str());
        if (!aq.empty()) std::printf("%s: %zu exact A2A distances\n", mesh.c_str(), aq.size());
        return 0;
    }
    if (sets.empty()) { std::fprintf(stderr, "no query pairs\n"); return 1; }
    const auto& qp = sets[0].qp;
    const auto& tr = sets[0].tr;
    const size_t m = qp.size();
    // plain exact POI baseline: one full sweep per POI (VTP, MMP if VTP rejects the mesh), P x P table
    std::string table_js;
    if (poi_table && !opt.poi.empty()) {
        auto tt = clk::now();
        std::vector<unsigned> F2(T.F.begin(), T.F.end());
        for (size_t f = 0; f < F2.size(); f += 3) { // counter-clockwise in the chart, as VTP needs
            const double *A = &T.V[3 * F2[f]], *B = &T.V[3 * F2[f + 1]], *C = &T.V[3 * F2[f + 2]];
            if ((B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]) < 0) std::swap(F2[f + 1], F2[f + 2]);
        }
        std::unique_ptr<Solver> sol(make_vtp(T.V, F2));
        if (!sol) sol.reset(make_mmp(T.V, F2, true));
        size_t P = opt.poi.size();
        std::vector<unsigned> tg(opt.poi.begin(), opt.poi.end());
        std::vector<float> tab(P * P);
        std::vector<double> out(P);
        for (size_t i = 0; i < P; ++i) {
            sol->dist(tg[i], tg, out.data());
            for (size_t j = 0; j < P; ++j) tab[i * P + j] = (float)out[j];
        }
        double tbt = since(tt);
        auto idx = [&](int v) { return (size_t)(std::lower_bound(opt.poi.begin(), opt.poi.end(), v) - opt.poi.begin()); };
        double best = 1e30;
        volatile double sink = 0;
        for (int rep = 0; rep < 5; ++rep) {
            auto q0 = clk::now();
            for (auto& pq : qp) sink = sink + tab[idx(pq.first) * P + idx(pq.second)];
            best = std::min(best, since(q0) / m);
        }
        double emax = 0, emean = 0;
        for (size_t i = 0; i < m; ++i) {
            double e = std::fabs(tab[idx(qp[i].first) * P + idx(qp[i].second)] / tr[i] - 1);
            emax = std::max(emax, e); emean += e / m;
        }
        std::printf("%s POI table P=%zu (%s): build %.1f s size %.2f MB | query %.3f us | mean err %.2e max %.2e\n",
                    mesh.c_str(), P, sol->name(), tbt, P * P * 4 / 1e6, best * 1e6, emean, emax);
        std::ostringstream o;
        o << ",\"poi_table\":{\"build_s\":" << tbt << ",\"index_bytes\":" << P * P * 4 << ",\"query_us\":" << best * 1e6
          << ",\"err_mean\":" << emean << ",\"err_max\":" << emax << ",\"solver\":\"" << sol->name() << "\"}";
        table_js = o.str();
        if (sink == -1.0) std::puts(""); // keeps the timed lookups
    }
    // near-table variants: several on one build keep the dense tables until every variant is encoded
    std::vector<std::string> vlist, nlist;
    std::vector<int> ilist;
    { std::stringstream is(interp_variants.empty() ? std::to_string(opt.interp) : interp_variants); std::string x; while (std::getline(is, x, ',')) ilist.push_back(std::stoi(x)); }
    { std::stringstream vs(variants); std::string x; while (std::getline(vs, x, ',')) vlist.push_back(x); }
    if (near_variants.empty()) near_variants = opt.near_store;
    { std::stringstream ns(near_variants); std::string x; while (std::getline(ns, x, ',')) nlist.push_back(x); }
    for (auto& nv : nlist) CertifiedIndex::near_mode_of(nv); // validates
    const std::string near_build = opt.near_store;
    if (nlist.size() > 1 || nlist[0] != opt.near_store) opt.near_store = "dense";
    auto t0 = clk::now();
    CertifiedIndex ix(T, opt);
    double tb = since(t0);
    double rss = peak_rss_mb();
    std::string argv_s;
    for (int i = 0; i < argc; ++i) { for (const char* c = argv[i]; *c; ++c) { if (*c == '"' || *c == '\\') argv_s += '\\'; argv_s += *c; } if (i + 1 < argc) argv_s += ' '; }
    std::ostringstream js;
    js << "{\"argv\":\"" << argv_s << "\",\"mesh\":\"" << mesh << "\",\"n\":" << T.n << ",\"eps\":" << opt.eps << ",\"gap_factor\":" << opt.gap_factor
       << ",\"hierarchy\":\"" << opt.hierarchy << "\",\"k_max\":" << opt.k_max << ",\"solver\":\"" << opt.solver << "\",\"exact\":\"" << opt.exact
       << "\",\"truth_mmp_fallbacks\":" << truth_fallbacks << ",\"build_mmp_fallbacks\":" << ix.mmp_fallbacks.load() << ",\"sweeps_vtp\":" << ix.sweeps_vtp << ",\"sweeps_mmp\":" << ix.sweeps_mmp
       << ",\"threads\":" << opt.threads << ",\"store\":\"" << opt.store << "\",\"near_store\":\"" << near_build
       << "\",\"pois\":" << opt.poi.size() << ",\"peak_rss_mb\":" << rss << table_js
       << ",\"build_s\":" << tb << ",\"cuts\":" << ix.cuts.size() << ",\"portals\":" << ix.portals() << ",\"labels\":" << ix.labels()
       << ",\"near_entries\":" << ix.near_total << ",\"timing\":{";
    bool first = true;
    for (auto& kv : ix.timing) { js << (first ? "" : ",") << "\"" << kv.first << "\":" << kv.second; first = false; }
    js << "},\"variants\":[";
    // the JSON is rewritten after every result (and once right after the build), so a job that dies
    // or times out keeps what it measured; "complete" says whether it finished
    auto flush = [&](bool done) {
        if (json.empty()) return;
        std::ofstream f(json);
        f << js.str() << "],\"complete\":" << (done ? "true" : "false") << "}\n";
    };
    flush(false);
    bool firstv = true;
    for (const auto& nv : nlist) {
      if (opt.near_store == "dense") { // built with dense tables: encode this variant (keeping them)
          auto mode = CertifiedIndex::near_mode_of(nv);
          if (mode == CertifiedIndex::NEAR_DENSE) ix.use_near(mode); else ix.encode_near(mode, true);
      }
      for (int iv : ilist) {
      ix.opt.interp = iv; // the estimate only: 0 none, 1 parabola, 2 model (bounds unchanged)
      for (const auto& var : vlist) {
        const std::string vname = (nlist.size() == 1 ? var : var + "/" + nv) + (ilist.size() > 1 ? "/i" + std::to_string(iv) : "");
        auto tv = clk::now();
        if (var == opt.store && opt.store != "dense") { /* coded during the build: storage already set */ }
        else if (var == "dense") ix.storage = DENSE;
        else if (var == "code16") { ix.encode16(16); ix.storage = CODE16; }
        else if (var == "code8") { ix.encode16(8); ix.storage = CODE8; }
        else if (var.rfind("cover", 0) == 0) { ix.build_covers(std::stod(var.substr(5))); ix.storage = COVER; }
        else { std::fprintf(stderr, "unknown variant %s\n", var.c_str()); continue; }
        double tprep = since(tv);
        // statistics of one query set, printed and appended to the JSON
        auto report = [&](const std::string& set, const std::vector<double>& B, const std::vector<double>& D,
                          const std::vector<double>& qt, int refined, double best, double lo, double hi) {
            const size_t mq = D.size();
            // weighted: the reference D is a fine Steiner graph (itself an upper bound), so soundness checks
            // the LB only, and an UB below D is a path the reference missed (counted, not an error)
            bool sound = true;
            int ub_below = 0;
            double emean = 0, emax = 0, ubmax = 0, ubmean = 0, worst_ratio = 0;
            std::vector<double> errs, widths;
            for (size_t i = 0; i < mq; ++i) {
                double l = B[3 * i], u = B[3 * i + 1], est = B[3 * i + 2], d = D[i];
                if (!(l <= d * (1 + 1e-9) && (weighted || d <= u * (1 + 1e-9)))) {
                    if (sound) std::fprintf(stderr, "UNSOUND %s %s pair %zu: LB %.9g d %.9g UB %.9g (LB/d - 1 = %.3g, d/UB - 1 = %.3g)\n",
                                            vname.c_str(), set.c_str(), i, l, d, u, l / d - 1, d / u - 1);
                    sound = false;
                }
                ub_below += weighted && u < d * (1 - 1e-9);
                double e = std::fabs(est / d - 1);
                errs.push_back(e); emean += e; emax = std::max(emax, e);
                ubmax = std::max(ubmax, u / d - 1); ubmean += u / d - 1;
                worst_ratio = std::max(worst_ratio, u / std::max(l, 1e-300));
                widths.push_back(l > 0 ? (u - l) / l : 0);
            }
            emean /= mq; ubmean /= mq;
            std::sort(errs.begin(), errs.end()); std::sort(widths.begin(), widths.end());
            std::vector<double> qs(qt); std::sort(qs.begin(), qs.end());
            double qmean = 0; for (double x : qt) qmean += x; qmean /= mq;
            double p99 = errs[(size_t)(0.99 * (mq - 1))];
            char band[64] = "";
            if (lo >= 0) std::snprintf(band, sizeof band, " %.0f-%.0f m", lo, hi);
            std::printf("%s eps=%.2f g=%.2fe %s kmax=%d threads=%d P=%zu %s%s [%s pairs%s, %zu]: build %.1f s (peak RSS %.0f MB) size %.2f MB (near %.2f MB) | query %.3f us, with refine mean %.2f us p99 %.1f us max %.0f us, refined %d/%zu | sound%s %s, max UB err %.3f%%, est mean %.4f%% p99 %.3f%% max %.3f%%, width median %.2f%%, max UB/LB %.4f\n",
                        mesh.c_str(), opt.eps, opt.gap_factor, opt.hierarchy.c_str(), opt.k_max, opt.threads, opt.poi.size(), vname.c_str(), refine ? "+refine" : "", set.c_str(), band, mq,
                        tb, rss, ix.nbytes() / 1e6, ix.near_bytes() / 1e6, best * 1e6, qmean * 1e6, qs[(size_t)(0.99 * (mq - 1))] * 1e6, qs.back() * 1e6,
                        refined, mq, weighted ? " (LB vs reference)" : "", sound ? "yes" : "NO", 100 * ubmax, 100 * emean, 100 * p99, 100 * emax, 100 * widths[mq / 2], worst_ratio);
            std::fflush(stdout);
            js << (firstv ? "" : ",") << "{\"variant\":\"" << vname << "\",\"set\":\"" << set << "\",\"near\":\"" << nv << "\",\"near_bytes\":" << ix.near_bytes();
            if (lo >= 0) js << ",\"band_lo\":" << lo << ",\"band_hi\":" << hi;
            js << ",\"refine\":" << (refine ? "true" : "false")
               << ",\"index_bytes\":" << ix.nbytes() << ",\"prep_s\":" << tprep << ",\"query_us\":" << best * 1e6
               << ",\"query_us_refine_mean\":" << qmean * 1e6 << ",\"query_us_refine_p99\":" << qs[(size_t)(0.99 * (mq - 1))] * 1e6
               << ",\"refined\":" << refined << ",\"mmp_fallbacks\":" << ix.mmp_fallbacks.load() << ",\"pairs\":" << mq << ",\"sound\":" << (sound ? "true" : "false") << ",\"ub_below_reference\":" << ub_below
               << ",\"max_ub_err\":" << ubmax << ",\"ub_mean\":" << ubmean << ",\"est_mean\":" << emean << ",\"est_p99\":" << p99
               << ",\"est_max\":" << emax << ",\"width_median\":" << widths[mq / 2] << ",\"max_ub_over_lb\":" << worst_ratio
               << ",\"cover_entries\":" << ix.cover_entries << "}";
            firstv = false;
            flush(false);
        };
        for (const auto& S : sets) {
            const auto& Q = S.qp;
            const size_t mq = Q.size();
            // latency of the plain certified query: best of 5 passes
            double best = 1e30, sink = 0, o[3];
            for (int rep = 0; rep < 5; ++rep) {
                auto q0 = clk::now();
                for (auto& p : Q) { ix.bounds(p.first, p.second, o); sink += o[2]; }
                best = std::min(best, since(q0) / mq);
            }
            // answers (with certify-or-refine if asked), one pass, timed per query
            std::vector<double> B(3 * mq), qt(mq);
            int refined = 0;
            for (size_t i = 0; i < mq; ++i) {
                auto q0 = clk::now();
                bool r = refine ? ix.bounds_refined(Q[i].first, Q[i].second, &B[3 * i]) : (ix.bounds(Q[i].first, Q[i].second, &B[3 * i]), false);
                qt[i] = since(q0);
                refined += r;
            }
            report(S.name, B, S.tr, qt, refined, best, S.lo, S.hi);
            (void)sink;
        }
        if (!aq.empty()) { // A2A: always certified (exact below the near radius, refine when wider than 1 + 2 eps)
            const size_t ma = aq.size();
            std::vector<double> B(3 * ma), qt(ma);
            int refined = 0;
            for (size_t i = 0; i < ma; ++i) {
                auto q0 = clk::now();
                refined += ix.a2a_bounds(aq[i].first, aq[i].second, &B[3 * i]);
                qt[i] = since(q0);
            }
            std::vector<double> qs(qt); std::sort(qs.begin(), qs.end());
            report("a2a", B, atr, qt, refined, qs[ma / 2], -1, -1); // "query" = the median, refines included
        }
      }
      }
    }
    flush(true);
    return 0;
}
