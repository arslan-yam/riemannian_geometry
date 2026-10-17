// Training data for NeuroGF (Zhang et al., NeurIPS 2023) on a terrain, in the layout of their per-model
// pipeline (code/neurogf_ovft): for each training source one exact MMP sweep (Kirsanov, the algorithm of
// the gdist package their code imports) gives
//   global   G random targets                         rows (source xyz, target xyz, distance)
//   local_1  the targets within r1 = 3 mean edges      ("several tens" in their data)
//   local_2  the targets within r2 = 25 mean edges     ("several thousands"), at most 8192
//   paths    the exact shortest paths to the first P global targets (128 points) and to up to 256
//            local_2 targets (32 points), resampled by arc length, source -> target
// and for each validation source the distances to 20,000 random targets (their checkpoint selection).
// Vertices of the test pairs (--exclude, "s t" per line) are never sources or targets.
// usage: ngf_data mesh.off out_dir [--sources 512] [--val 8] [--global 4096] [--paths 1024] [--seed 1]
//                 [--threads 1] [--exclude pairs_file ...]
// writes out_dir/{train_global, train_l1, train_l1_off, train_l2, train_l2_off, train_paths_global,
// train_paths_local, train_paths_local_off, val}.npy (float32 / int64) and meta.json
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "geodesic_algorithm_exact.h"  // third_party/kirsanov
#include "../src/mesh_io.h"
using clk = std::chrono::steady_clock;

static void npy(const std::string& path, const std::string& descr, const std::vector<size_t>& shape, const void* data, size_t bytes) {
    std::string sh = "(";
    for (size_t i = 0; i < shape.size(); ++i) sh += std::to_string(shape[i]) + (shape.size() == 1 || i + 1 < shape.size() ? "," : "");
    sh += ")";
    std::string h = "{'descr': '" + descr + "', 'fortran_order': False, 'shape': " + sh + ", }";
    while ((10 + h.size() + 1) % 64) h += ' ';
    h += '\n';
    std::ofstream o(path, std::ios::binary);
    o.write("\x93NUMPY\x01\x00", 8);
    unsigned short len = (unsigned short)h.size();
    o.write((const char*)&len, 2);
    o << h;
    o.write((const char*)data, (std::streamsize)bytes);
}

struct SourceData {
    std::vector<float> global, l1, l2, pg, pl, val;
    double secs = 0;
};

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: ngf_data mesh.off out_dir [options]\n"); return 1; }
    std::string mesh_file = argv[1], out = argv[2];
    int S = 512, NV = 8, G = 4096, P = 1024, PL = 256, seed = 1, threads = 1, NVT = 20000;
    const int PPG = 128, PPL = 32, L2MAX = 8192, L1MAX = 256;
    std::vector<std::string> excl;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--sources") S = std::stoi(nx());
        else if (a == "--val") NV = std::stoi(nx());
        else if (a == "--global") G = std::stoi(nx());
        else if (a == "--paths") P = std::stoi(nx());
        else if (a == "--seed") seed = std::stoi(nx());
        else if (a == "--threads") threads = std::stoi(nx());
        else if (a == "--exclude") excl.push_back(nx());
    }
    P = std::min(P, G);
    std::vector<double> V; std::vector<unsigned> F;
    double origin[2];
    cci::read_off(mesh_file, V, F, origin);
    for (size_t i = 0; i < V.size(); i += 3) { V[i] += origin[0]; V[i + 1] += origin[1]; } // the file's coordinates
    const unsigned n = (unsigned)(V.size() / 3);
    double el = 0; size_t ne = 0;
    for (size_t f = 0; f < F.size(); f += 3) for (int k = 0; k < 3; ++k) {
        unsigned a = F[f + k], b = F[f + (k + 1) % 3]; double d = 0;
        for (int j = 0; j < 3; ++j) d += (V[3 * a + j] - V[3 * b + j]) * (V[3 * a + j] - V[3 * b + j]);
        el += std::sqrt(d); ++ne;
    }
    el /= ne;
    const double r1 = 3 * el, r2 = 25 * el;
    std::vector<char> bad(n, 0);
    for (auto& f : excl) { std::ifstream in(f); long s, t; while (in >> s >> t) { if (s >= 0 && s < n) bad[s] = 1; if (t >= 0 && t < n) bad[t] = 1; } }
    std::vector<unsigned> pool;
    for (unsigned v = 0; v < n; ++v) if (!bad[v]) pool.push_back(v);
    std::mt19937_64 rng(seed);
    std::shuffle(pool.begin(), pool.end(), rng);
    if ((int)pool.size() < S + NV) { std::fprintf(stderr, "not enough vertices\n"); return 1; }
    std::vector<unsigned> src(pool.begin(), pool.begin() + S + NV); // S training sources, then NV validation sources
    std::vector<SourceData> res(S + NV);
    std::atomic<int> next{0};
    auto t0 = clk::now();
    auto worker = [&]() {
        geodesic::Mesh mesh;
        { std::vector<double> p(V); std::vector<unsigned> f(F); mesh.initialize_mesh_data(p, f); }
        geodesic::GeodesicAlgorithmExact alg(&mesh);
        std::vector<double> d(n);
        for (int i; (i = next++) < S + NV;) {
            auto ts = clk::now();
            SourceData& R = res[i];
            unsigned s = src[i];
            std::vector<geodesic::SurfacePoint> sp(1, geodesic::SurfacePoint(&mesh.vertices()[s]));
            alg.propagate(sp, geodesic::GEODESIC_INF);
            for (unsigned v = 0; v < n; ++v) { geodesic::SurfacePoint p(&mesh.vertices()[v]); alg.best_source(p, d[v]); }
            std::mt19937_64 r(seed * 1000003ull + i);
            auto row = [&](std::vector<float>& o, unsigned t) {
                for (int j = 0; j < 3; ++j) o.push_back((float)V[3 * s + j]);
                for (int j = 0; j < 3; ++j) o.push_back((float)V[3 * t + j]);
                o.push_back((float)d[t]);
            };
            auto ok = [&](unsigned t) { return t != s && !bad[t] && std::isfinite(d[t]) && d[t] < geodesic::GEODESIC_INF; };
            if (i >= S) { // validation source
                for (int k = 0; k < NVT;) { unsigned t = pool[r() % pool.size()]; if (ok(t)) { row(R.val, t); ++k; } }
                R.secs = std::chrono::duration<double>(clk::now() - ts).count();
                continue;
            }
            std::vector<unsigned> gl;
            while ((int)gl.size() < G) { unsigned t = pool[r() % pool.size()]; if (ok(t)) gl.push_back(t); }
            for (unsigned t : gl) row(R.global, t);
            std::vector<unsigned> n1, n2;
            for (unsigned v = 0; v < n; ++v) if (ok(v)) { if (d[v] <= r1) n1.push_back(v); if (d[v] <= r2) n2.push_back(v); }
            std::shuffle(n1.begin(), n1.end(), r); if (n1.size() > (size_t)L1MAX) n1.resize(L1MAX);
            std::shuffle(n2.begin(), n2.end(), r); if (n2.size() > (size_t)L2MAX) n2.resize(L2MAX);
            for (unsigned t : n1) row(R.l1, t);
            for (unsigned t : n2) row(R.l2, t);
            // shortest paths, resampled to k points by arc length, source -> target
            auto path = [&](unsigned t, int k, std::vector<float>& o) {
                geodesic::SurfacePoint tp(&mesh.vertices()[t]);
                std::vector<geodesic::SurfacePoint> pts;
                alg.trace_back(tp, pts);
                std::reverse(pts.begin(), pts.end());
                std::vector<double> acc(pts.size(), 0.0);
                for (size_t j = 1; j < pts.size(); ++j) acc[j] = acc[j - 1] + pts[j].distance(&pts[j - 1]);
                double L = acc.back();
                size_t seg = 0;
                for (int q = 0; q < k; ++q) {
                    double x = pts.size() < 2 ? 0 : L * q / (k - 1);
                    while (seg + 2 < pts.size() && acc[seg + 1] < x) ++seg;
                    double w = 0;
                    if (pts.size() >= 2 && acc[seg + 1] > acc[seg]) w = std::min(1.0, std::max(0.0, (x - acc[seg]) / (acc[seg + 1] - acc[seg])));
                    geodesic::SurfacePoint& A = pts[seg];
                    geodesic::SurfacePoint& B = pts[std::min(seg + 1, pts.size() - 1)];
                    o.push_back((float)(A.x() + w * (B.x() - A.x())));
                    o.push_back((float)(A.y() + w * (B.y() - A.y())));
                    o.push_back((float)(A.z() + w * (B.z() - A.z())));
                }
            };
            for (int j = 0; j < P; ++j) path(gl[j], PPG, R.pg);
            for (size_t j = 0; j < n2.size() && (int)j < PL; ++j) path(n2[j], PPL, R.pl);
            R.secs = std::chrono::duration<double>(clk::now() - ts).count();
        }
    };
    std::vector<std::thread> th;
    for (int k = 0; k < threads; ++k) th.emplace_back(worker);
    for (auto& t : th) t.join();
    double wall = std::chrono::duration<double>(clk::now() - t0).count(), cpu = 0;
    for (auto& R : res) cpu += R.secs;
    // gather and write
    auto cat = [&](std::vector<float> SourceData::*m, int lo, int hi, std::vector<int64_t>* off, size_t rowlen) {
        std::vector<float> all;
        if (off) off->assign(1, 0);
        for (int i = lo; i < hi; ++i) { all.insert(all.end(), (res[i].*m).begin(), (res[i].*m).end()); if (off) off->push_back((int64_t)(all.size() / rowlen)); }
        return all;
    };
    std::vector<int64_t> o1, o2, op;
    auto g_ = cat(&SourceData::global, 0, S, nullptr, 7);
    npy(out + "/train_global.npy", "<f4", {(size_t)S, (size_t)G, 7}, g_.data(), g_.size() * 4);
    auto l1 = cat(&SourceData::l1, 0, S, &o1, 7);
    npy(out + "/train_l1.npy", "<f4", {l1.size() / 7, 7}, l1.data(), l1.size() * 4);
    npy(out + "/train_l1_off.npy", "<i8", {o1.size()}, o1.data(), o1.size() * 8);
    auto l2 = cat(&SourceData::l2, 0, S, &o2, 7);
    npy(out + "/train_l2.npy", "<f4", {l2.size() / 7, 7}, l2.data(), l2.size() * 4);
    npy(out + "/train_l2_off.npy", "<i8", {o2.size()}, o2.data(), o2.size() * 8);
    auto pg = cat(&SourceData::pg, 0, S, nullptr, PPG * 3);
    npy(out + "/train_paths_global.npy", "<f4", {(size_t)S, (size_t)P, (size_t)PPG, 3}, pg.data(), pg.size() * 4);
    auto pl = cat(&SourceData::pl, 0, S, &op, PPL * 3);
    npy(out + "/train_paths_local.npy", "<f4", {pl.size() / (PPL * 3), (size_t)PPL, 3}, pl.data(), pl.size() * 4);
    npy(out + "/train_paths_local_off.npy", "<i8", {op.size()}, op.data(), op.size() * 8);
    auto va = cat(&SourceData::val, S, S + NV, nullptr, 7);
    npy(out + "/val.npy", "<f4", {va.size() / 7, 7}, va.data(), va.size() * 4);
    std::ofstream m(out + "/meta.json");
    m << "{\"tool\":\"ngf_data\",\"mesh\":\"" << mesh_file << "\",\"n\":" << n << ",\"faces\":" << F.size() / 3 << ",\"sources\":" << S
      << ",\"val_sources\":" << NV << ",\"global\":" << G << ",\"paths\":" << P << ",\"r1\":" << r1 << ",\"r2\":" << r2
      << ",\"mean_edge\":" << el << ",\"excluded_vertices\":" << (n - pool.size()) << ",\"threads\":" << threads
      << ",\"wall_s\":" << wall << ",\"cpu_s\":" << cpu << ",\"l1_rows\":" << l1.size() / 7 << ",\"l2_rows\":" << l2.size() / 7
      << ",\"local_paths\":" << pl.size() / (PPL * 3) << "}\n";
    std::printf("%s: %d sources (+%d validation), %.1f s wall, %.1f CPU-s (%d threads), r1 %.1f r2 %.1f, l1 %zu l2 %zu rows, %zu local paths\n",
                mesh_file.c_str(), S, NV, wall, cpu, threads, r1, r2, l1.size() / 7, l2.size() / 7, pl.size() / (PPL * 3));
}
