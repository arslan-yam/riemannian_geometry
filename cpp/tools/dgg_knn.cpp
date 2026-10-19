// DGG (Adikusuma, Fang, He, TOG 2020; authors' graph, port in external/dgg_port) answering kNN and range
// queries the way a graph method naturally does: one Dijkstra from the query vertex on the precomputed
// geodesic graph, stopped once the k nearest other objects are settled (kNN) or the next vertex is beyond
// r (range). Scored on the pairs written by cpp/tools/proximity --export (the same queries, objects and
// exact distances as our certified proximity queries), with the metrics of scripts/competitors/knn_eval.py.
// usage: dgg_knn graph.binary mesh.off PREFIX.pairs [--k 1,10,25] [--radius 0.01,0.03,0.1] [--json out.json]
//   PREFIX.pairs: "query object distance" per line, grouped by query; --radius: fractions of the chart diagonal
#include "dgg_dijkstra.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

typedef std::chrono::steady_clock clk;
static double ms_since(clk::time_point t) { return std::chrono::duration<double, std::milli>(clk::now() - t).count(); }

static std::vector<double> parse_list(const std::string& s) {
    std::vector<double> v; std::stringstream ss(s); std::string x;
    while (std::getline(ss, x, ',')) if (!x.empty()) v.push_back(std::stod(x));
    return v;
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: dgg_knn graph.binary mesh.off pairs [--k 1,10,25] [--radius 0.01,0.03,0.1] [--json out.json]\n"); return 1; }
    std::string ks_s = "1,10,25", rs_s = "0.01,0.03,0.1", json;
    for (int i = 4; i + 1 < argc; ++i) {
        std::string a = argv[i];
        if (a == "--k") ks_s = argv[i + 1]; else if (a == "--radius") rs_s = argv[i + 1]; else if (a == "--json") json = argv[i + 1];
    }
    std::vector<int> ks; for (double k : parse_list(ks_s)) ks.push_back((int)k);
    std::vector<double> rfs = parse_list(rs_s);

    // chart diagonal of the mesh (as knn_eval.py and proximity)
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    {
        std::ifstream in(argv[2]); std::string head; in >> head; long nv, nf, ne; in >> nv >> nf >> ne;
        for (long i = 0; i < nv; ++i) { double x, y, z; in >> x >> y >> z; lo[0] = std::min(lo[0], x); lo[1] = std::min(lo[1], y); hi[0] = std::max(hi[0], x); hi[1] = std::max(hi[1], y); }
    }
    const double diag = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);

    // queries and their objects with exact distances, in file order
    std::vector<int> qv; std::vector<std::vector<std::pair<int, double>>> objs;
    {
        std::ifstream in(argv[3]); long q, o; double d; std::unordered_map<long, int> at;
        while (in >> q >> o >> d) {
            auto it = at.find(q);
            if (it == at.end()) { at[q] = (int)qv.size(); qv.push_back((int)q); objs.emplace_back(); it = at.find(q); }
            objs[it->second].push_back({(int)o, d});
        }
    }
    const int Q = (int)qv.size();
    if (Q == 0) { std::fprintf(stderr, "dgg_knn: no pairs in %s\n", argv[3]); return 1; }

    Dijkstra<double> g;
    g.readDGGFile(argv[1]);
    const int N = g.NodeNum();
    const double INF = std::numeric_limits<double>::max();
    std::vector<double> dist(N, INF); std::vector<char> done(N, 0); std::vector<int> touched;
    std::vector<int> objmark(N, -1);  // index of the vertex among the current query's objects, or -1
    typedef std::pair<double, int> QN;
    // one early-stopping Dijkstra from src; stop(settled object count, popped distance) ends it before
    // settling a vertex. Returns the settled objects with their graph distances, in settle order.
    auto sweep = [&](int src, auto stop) {
        for (int v : touched) { dist[v] = INF; done[v] = 0; }
        touched.clear();
        std::vector<std::pair<int, double>> got;
        std::priority_queue<QN, std::vector<QN>, std::greater<QN>> pq;
        dist[src] = 0; touched.push_back(src); pq.push({0.0, src});
        while (!pq.empty()) {
            QN top = pq.top();
            if (stop((int)got.size(), top.first)) break;
            pq.pop();
            int u = top.second;
            if (done[u]) continue;
            done[u] = 1;
            if (objmark[u] >= 0) got.push_back({objmark[u], dist[u]});
            const auto& nb = g.graphNeighbor(u); const auto& nw = g.graphNeighborDis(u);
            for (size_t i = 0; i < nb.size(); ++i) {
                int v = nb[i]; double w = dist[u] + nw[i];
                if (!done[v] && dist[v] > w) { if (dist[v] == INF) touched.push_back(v); dist[v] = w; pq.push({w, v}); }
            }
        }
        return got;
    };

    std::ostringstream knn_js, range_js;
    std::printf("DGG kNN/range %s: %d queries, %.1f objects per query, diag %.1f\n", argv[1], Q, Q ? (double)objs[0].size() : 0.0, diag);
    // kNN: the k other objects settled first (the query vertex itself is not among its objects)
    for (size_t ki = 0; ki < ks.size(); ++ki) {
        int k = ks[ki]; double rec_sum = 0, worst = 0; std::vector<double> tq;
        for (int i = 0; i < Q; ++i) {
            auto& ob = objs[i]; int kk = std::min(k, (int)ob.size());
            for (size_t j = 0; j < ob.size(); ++j) objmark[ob[j].first] = (int)j;
            auto t0 = clk::now();
            auto got = sweep(qv[i], [&](int n, double) { return n >= kk; });
            tq.push_back(ms_since(t0));
            for (auto& o : ob) objmark[o.first] = -1;
            std::vector<double> d; for (auto& o : ob) d.push_back(o.second);
            std::vector<int> ord(ob.size()); for (size_t j = 0; j < ord.size(); ++j) ord[j] = (int)j;
            std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return d[a] < d[b]; });
            std::set<int> truth(ord.begin(), ord.begin() + kk);
            double dk = d[ord[kk - 1]], wq = 0; int hit = 0;
            for (int j = 0; j < (int)got.size() && j < kk; ++j) { hit += truth.count(got[j].first); wq = std::max(wq, d[got[j].first]); }
            rec_sum += (double)hit / kk;
            worst = std::max(worst, dk > 0 ? wq / dk : 1.0);
        }
        double mean = 0; for (double t : tq) mean += t / Q;
        std::vector<double> s = tq; std::sort(s.begin(), s.end());
        if (ki) knn_js << ",";
        knn_js << "{\"k\":" << k << ",\"recall\":" << rec_sum / Q << ",\"worst_ratio\":" << worst << ",\"query_ms_mean\":" << mean << ",\"query_ms_p50\":" << s[Q / 2] << "}";
        std::printf("  k=%d: recall %.4f, worst ratio %.5f, %.3f ms per query\n", k, rec_sum / Q, worst, mean);
    }
    // range r: every object settled before the first vertex beyond r
    for (size_t ri = 0; ri < rfs.size(); ++ri) {
        double r = rfs[ri] * diag, prec_sum = 0, rec_sum = 0, mean = 0; long fp = 0, fn = 0;
        for (int i = 0; i < Q; ++i) {
            auto& ob = objs[i];
            for (size_t j = 0; j < ob.size(); ++j) objmark[ob[j].first] = (int)j;
            auto t0 = clk::now();
            auto got = sweep(qv[i], [&](int, double dd) { return dd > r; });
            mean += ms_since(t0) / Q;
            for (auto& o : ob) objmark[o.first] = -1;
            std::vector<char> in(ob.size(), 0); for (auto& x : got) in[x.first] = 1;
            int inter = 0, ns = 0, nt = 0;
            for (size_t j = 0; j < ob.size(); ++j) { bool t = ob[j].second <= r, sel = in[j]; inter += t && sel; ns += sel; nt += t; fp += sel && !t; fn += t && !sel; }
            prec_sum += ns ? (double)inter / ns : 1.0; rec_sum += nt ? (double)inter / nt : 1.0;
        }
        if (ri) range_js << ",";
        range_js << "{\"radius_frac\":" << rfs[ri] << ",\"radius\":" << r << ",\"precision\":" << prec_sum / Q << ",\"recall\":" << rec_sum / Q
                 << ",\"false_positives\":" << fp << ",\"missed\":" << fn << ",\"query_ms_mean\":" << mean << "}";
        std::printf("  r=%g: precision %.4f, recall %.4f (%ld false positives, %ld missed), %.3f ms per query\n", rfs[ri], prec_sum / Q, rec_sum / Q, fp, fn, mean);
    }
    if (!json.empty()) {
        std::ofstream o(json);
        o << "{\"tool\":\"dgg_knn\",\"graph\":\"" << argv[1] << "\",\"pairs\":\"" << argv[3] << "\",\"queries\":" << Q << ",\"diag\":" << diag
          << ",\"methods\":[{\"method\":\"DGG\",\"knn\":[" << knn_js.str() << "],\"range\":[" << range_js.str() << "]}]}\n";
    }
}
