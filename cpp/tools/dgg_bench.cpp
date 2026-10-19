// DGG baseline (Adikusuma, Fang, He, TOG 2020; github.com/GeodesicGraph/DGG_and_SVG, macOS port in
// external/dgg_port): loads a precomputed discrete geodesic graph once and answers each query with
// one Dijkstra from the source (their query algorithm, double precision), scored against exact
// distances. usage: dgg_bench graph.binary pairs_file truth_file
// Also times a point-to-point Dijkstra on the same graph that stops once the target is settled (same answer,
// the standard early exit; their code always computes the whole distance field): query_ms_p2p in the JSON.
#include "dgg_dijkstra.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <string>
#include <map>
#include <queue>
#include <limits>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: dgg_bench graph.binary pairs truth [--json out.json]\n"); return 1; }
    std::string json;
    for (int i = 4; i + 1 < argc; ++i) if (std::string(argv[i]) == "--json") json = argv[i + 1];
    Dijkstra<double> g;
    g.readDGGFile(argv[1]);
    std::ifstream pf(argv[2]), tf(argv[3]);
    std::vector<std::pair<int, int>> qp; int s, t;
    while (pf >> s >> t) { qp.push_back({s, t}); std::string rest; std::getline(pf, rest); }
    std::vector<double> tr(qp.size());
    for (auto& d : tr) tf >> d;
    // point-to-point: the same Dijkstra (same graph, same double weights), stopped when the target is popped
    const int N = g.NodeNum();
    std::vector<double> pd(N, std::numeric_limits<double>::max()); std::vector<char> done(N, 0); std::vector<int> touched;
    auto p2p = [&](int src, int dst) {
        for (int v : touched) { pd[v] = std::numeric_limits<double>::max(); done[v] = 0; }
        touched.clear();
        typedef std::pair<double, int> QN;
        std::priority_queue<QN, std::vector<QN>, std::greater<QN>> pq;
        pd[src] = 0; touched.push_back(src); pq.push({0, src});
        while (!pq.empty()) {
            int u = pq.top().second; pq.pop();
            if (done[u]) continue;
            done[u] = 1;
            if (u == dst) return pd[u];
            const auto& nb = g.graphNeighbor(u); const auto& nw = g.graphNeighborDis(u);
            for (size_t i = 0; i < nb.size(); ++i) {
                int v = nb[i]; double w = pd[u] + nw[i];
                if (!done[v] && pd[v] > w) { if (pd[v] == std::numeric_limits<double>::max()) touched.push_back(v); pd[v] = w; pq.push({w, v}); }
            }
        }
        return pd[dst];
    };
    double emean = 0, emax = -1e30, emin = 1e30, tsum = 0, tp2p = 0, p2p_diff = 0;
    std::vector<double> aerr;
    for (size_t i = 0; i < qp.size(); ++i) {
        auto t0 = std::chrono::steady_clock::now();
        g.findShortestDistance(qp[i].first);
        double d = g.distanceToSource(qp[i].second);
        tsum += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        auto t1 = std::chrono::steady_clock::now();
        double d2 = p2p(qp[i].first, qp[i].second);
        tp2p += std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
        p2p_diff = std::max(p2p_diff, std::fabs(d2 - d) / std::max(d, 1e-300));
        double e = d / tr[i] - 1;
        emean += e; emax = std::max(emax, e); emin = std::min(emin, e); aerr.push_back(std::fabs(e));
    }
    size_t m = qp.size();
    std::printf("DGG %s: %zu queries, query %.3f ms (point-to-point with early exit %.3f ms, answers differ by %.1e), err mean %.4f%% max %.4f%% min %.4f%%\n",
                argv[1], m, 1e3 * tsum / m, 1e3 * tp2p / m, p2p_diff, 100 * emean / m, 100 * emax, 100 * emin);
    if (!json.empty()) {
        std::sort(aerr.begin(), aerr.end());
        double am = 0; for (double x : aerr) am += x / m;
        std::ofstream o(json);
        o << "{\"tool\":\"dgg_bench\",\"graph\":\"" << argv[1] << "\",\"queries\":" << m << ",\"query_ms\":" << 1e3 * tsum / m << ",\"query_ms_p2p\":" << 1e3 * tp2p / m << ",\"p2p_max_rel_diff\":" << p2p_diff
          << ",\"err_mean_signed\":" << emean / m << ",\"err_mean\":" << am << ",\"err_p99\":" << aerr[(size_t)(0.99 * (m - 1))]
          << ",\"err_max\":" << emax << ",\"err_min\":" << emin << "}\n";
    }
}
