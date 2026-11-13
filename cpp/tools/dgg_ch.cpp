// DGG + a graph distance index: a contraction hierarchy (Geisberger et al. 2008) on the discrete geodesic graph of
// DGG (Adikusuma, Fang, He, TOG 2020; authors' graph, port in external/dgg_port), the review's "good surface graph
// plus a good graph distance index". The graph is directed (an edge u->v is stored with u's neighbours), so the
// hierarchy is too. Dense graphs make contraction expensive (every contraction runs witness searches from all
// in-neighbours), so the contraction stops when the remaining graph gets too dense (--core-degree) or after
// --time-limit hours, and the remaining nodes form a core that both query directions search freely (core-based
// CH). Queries: upward searches to the core, then a bidirectional Dijkstra on the core. The same pairs are also answered by a plain Dijkstra stopped at the target (the baseline of
// review item 1) and checked against it, and both are scored against the exact distances.
// usage: dgg_ch graph.binary pairs truth [--json out.json] [--settle 500] [--core-degree 400 | --core-factor F] [--time-limit 24] [--ch 0]
//               [--pll degree|kd --mesh mesh.off] [--pll-time-limit 24]
//   pairs: "s t" per line (comp/<mesh>.V2V.pairs); truth: one exact distance per line
#include "dgg_dijkstra.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

typedef std::chrono::steady_clock clk;
static double since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }
static const double INF = std::numeric_limits<double>::infinity();
struct E { int v; double w; };
typedef std::pair<double, int> QN;
typedef std::priority_queue<QN, std::vector<QN>, std::greater<QN>> MinQ;

static long peak_rss_kb() {
    std::ifstream f("/proc/self/status"); std::string line;
    while (std::getline(f, line)) if (line.rfind("VmHWM:", 0) == 0) return std::stol(line.substr(6));
    return -1;
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: dgg_ch graph.binary pairs truth [--json out.json] [--settle 500] [--core-degree 400] [--time-limit 24]\n"); return 1; }
    std::string json, pll_order, mesh_path; int settle_limit = 500; double core_degree = 400, core_factor = 0, time_limit_h = 24, pll_time_limit_h = 24;
    bool do_ch = true;
    for (int i = 4; i + 1 < argc; ++i) {
        std::string a = argv[i];
        if (a == "--json") json = argv[i + 1];
        else if (a == "--settle") settle_limit = std::stoi(argv[i + 1]);
        else if (a == "--core-degree") core_degree = std::stod(argv[i + 1]);
        else if (a == "--core-factor") core_factor = std::stod(argv[i + 1]);  // stop at this multiple of the initial mean out-degree
        else if (a == "--time-limit") time_limit_h = std::stod(argv[i + 1]);
        else if (a == "--pll") pll_order = argv[i + 1];          // hub labels too, in this vertex order: degree | kd
        else if (a == "--mesh") mesh_path = argv[i + 1];         // the terrain (.off): chart coordinates for --pll kd
        else if (a == "--pll-time-limit") pll_time_limit_h = std::stod(argv[i + 1]);
        else if (a == "--ch") do_ch = std::string(argv[i + 1]) != "0";
    }
    std::vector<std::pair<int, int>> qp; std::vector<double> tr;
    { std::ifstream pf(argv[2]); int s, t; std::string rest; while (pf >> s >> t) { qp.push_back({s, t}); std::getline(pf, rest); } }
    { std::ifstream tf(argv[3]); double d; while (tf >> d) tr.push_back(d); }
    if (tr.size() < qp.size()) qp.resize(tr.size());
    const size_t m = qp.size();

    // ---- the graph
    auto t_load = clk::now();
    std::vector<std::vector<E>> out, in;
    long long m0 = 0;
    int n;
    {
        Dijkstra<double> g;
        g.readDGGFile(argv[1]);
        n = g.NodeNum();
        out.resize(n); in.resize(n);
        for (int u = 0; u < n; ++u) {
            const auto& nb = g.graphNeighbor(u); const auto& nw = g.graphNeighborDis(u);
            for (size_t i = 0; i < nb.size(); ++i) if (nb[i] != u) { out[u].push_back({nb[i], nw[i]}); in[nb[i]].push_back({u, nw[i]}); ++m0; }
        }
    }
    double load_s = since(t_load);
    // original graph in CSR, for the plain point-to-point Dijkstra
    std::vector<long long> o_off(n + 1, 0); std::vector<E> o_adj; o_adj.reserve(m0);
    for (int u = 0; u < n; ++u) { for (auto& e : out[u]) o_adj.push_back(e); o_off[u + 1] = (long long)o_adj.size(); }
    std::printf("DGG graph %s: %d nodes, %lld edges (mean out-degree %.1f), loaded in %.1f s\n", argv[1], n, m0, (double)m0 / n, load_s);
    std::fflush(stdout);

    // ---- contraction
    auto t_build = clk::now();
    std::vector<char> contracted(n, 0);
    std::vector<int> rank(n, n), deleted_nb(n, 0);
    std::vector<double> wd(n, INF); std::vector<int> wtouch; std::vector<char> wdone(n, 0);
    // witness search from u in the remaining graph without v, up to distance lim and settle_limit nodes
    auto witness = [&](int u, int v, double lim) {
        for (int x : wtouch) { wd[x] = INF; wdone[x] = 0; }
        wtouch.clear();
        MinQ q; wd[u] = 0; wtouch.push_back(u); q.push({0, u});
        int settled = 0;
        while (!q.empty()) {
            QN tp = q.top(); q.pop();
            int x = tp.second;
            if (wdone[x]) continue;
            if (tp.first > lim || ++settled > settle_limit) break;
            wdone[x] = 1;
            for (auto& e : out[x]) {
                if (contracted[e.v] || e.v == v) continue;
                double nd = tp.first + e.w;
                if (nd < wd[e.v]) { if (wd[e.v] == INF) wtouch.push_back(e.v); wd[e.v] = nd; q.push({nd, e.v}); }
            }
        }
    };
    // the shortcuts that contracting v needs (or only their number, when simulating)
    struct SC { int a, b; double w; };
    // direct edges of u (stamped), the cheapest witness: the DGG is close to transitive, so most pairs u -> x
    // around v already have an edge at most as long as the path through v, and need no search
    std::vector<double> dw(n, INF); std::vector<int> dstamp(n, -1); int stamp = 0;
    std::vector<int> need;
    auto shortcuts = [&](int v, std::vector<SC>* res) {
        int cnt = 0;
        for (auto& ei : in[v]) {
            int u = ei.v;
            if (contracted[u]) continue;
            ++stamp;
            for (auto& e : out[u]) if (!contracted[e.v]) { if (dstamp[e.v] != stamp || e.w < dw[e.v]) dw[e.v] = e.w; dstamp[e.v] = stamp; }
            need.clear();
            double lim = 0;
            for (size_t j = 0; j < out[v].size(); ++j) {
                int x = out[v][j].v;
                if (contracted[x] || x == u) continue;
                double c = ei.w + out[v][j].w;
                if (dstamp[x] == stamp && dw[x] <= c * (1 + 1e-12)) continue;
                need.push_back((int)j); lim = std::max(lim, c);
            }
            if (need.empty()) continue;
            witness(u, v, lim);
            for (int j : need) {
                int x = out[v][j].v;
                double c = ei.w + out[v][j].w;
                if (wd[x] <= c * (1 + 1e-12)) continue;  // a witness at most as long (rounding: keep the graph's answers)
                ++cnt;
                if (res) res->push_back({u, x, c});
            }
        }
        return cnt;
    };
    auto live_deg = [&](int v) {
        int d = 0;
        for (auto& e : out[v]) d += !contracted[e.v];
        for (auto& e : in[v]) d += !contracted[e.v];
        return d;
    };
    auto prio = [&](int v) { return (double)shortcuts(v, nullptr) - live_deg(v) + deleted_nb[v]; };
    // initial priorities with a cheap witness limit (the lazy updates recompute them with the full one)
    MinQ pq;
    if (do_ch) {
        int keep = settle_limit; settle_limit = std::min(settle_limit, 50);
        for (int v = 0; v < n; ++v) pq.push({prio(v), v});
        settle_limit = keep;
    }
    std::printf("initial priorities in %.1f s\n", since(t_build)); std::fflush(stdout);
    long long live_edges = m0, shortcuts_added = 0;
    if (core_factor > 0) core_degree = core_factor * (double)m0 / std::max(n, 1);
    int next_rank = 0, remaining = n;
    double contract_s = 0; bool stopped_time = false;
    std::vector<SC> sc;
    auto last_print = clk::now();
    while (do_ch && !pq.empty()) {
        if ((double)live_edges / std::max(remaining, 1) > core_degree) break;
        if (since(t_build) > time_limit_h * 3600) { stopped_time = true; break; }
        QN tp = pq.top(); pq.pop();
        int v = tp.second;
        if (contracted[v]) continue;
        double p = prio(v);
        if (!pq.empty() && p > pq.top().first) { pq.push({p, v}); continue; }  // lazy update
        sc.clear(); shortcuts(v, &sc);
        contracted[v] = 1; rank[v] = next_rank++; --remaining;
        for (auto& e : out[v]) if (!contracted[e.v]) { --live_edges; ++deleted_nb[e.v]; }
        for (auto& e : in[v]) if (!contracted[e.v]) { --live_edges; ++deleted_nb[e.v]; }
        for (auto& s : sc) { // add or improve a -> b
            bool found = false;
            for (auto& e : out[s.a]) if (e.v == s.b) { found = true; if (s.w < e.w) { e.w = s.w; for (auto& f : in[s.b]) if (f.v == s.a) f.w = s.w; } break; }
            if (!found) { out[s.a].push_back({s.b, s.w}); in[s.b].push_back({s.a, s.w}); ++live_edges; ++shortcuts_added; }
        }
        if (since(last_print) > 600) {
            last_print = clk::now();
            std::printf("  %.1f h: contracted %d of %d, live mean degree %.1f, shortcuts %lld\n", since(t_build) / 3600, next_rank, n,
                        (double)live_edges / std::max(remaining, 1), shortcuts_added);
            std::fflush(stdout);
        }
    }
    contract_s = since(t_build);
    const int core = remaining;
    // search graphs: up[v] = edges v -> x with rank[x] > rank[v] (core nodes all have rank n: core edges in both);
    // down[v] = edges x -> v with rank[x] > rank[v], reversed (the backward search goes up from t along them)
    std::vector<long long> u_off(n + 1, 0), d_off(n + 1, 0); std::vector<E> u_adj, d_adj;
    for (int v = 0; v < n; ++v) {
        for (auto& e : out[v]) if (rank[e.v] > rank[v] || (rank[v] == n && rank[e.v] == n)) u_adj.push_back(e);
        u_off[v + 1] = (long long)u_adj.size();
        for (auto& e : in[v]) if (rank[e.v] > rank[v] || (rank[v] == n && rank[e.v] == n)) d_adj.push_back(e);
        d_off[v + 1] = (long long)d_adj.size();
    }
    std::vector<std::vector<E>>().swap(out); std::vector<std::vector<E>>().swap(in);
    double build_s = since(t_build);
    const long long index_bytes = (long long)(u_adj.size() + d_adj.size()) * (long long)sizeof(E) + 2LL * (n + 1) * (long long)sizeof(long long);
    std::printf("CH: contracted %d of %d nodes in %.1f s (%s), core %d nodes, %lld shortcuts, %zu up + %zu down edges, %.1f MB\n",
                n - core, n, contract_s, stopped_time ? "time limit" : (core ? "core degree" : "complete"), core, shortcuts_added,
                u_adj.size(), d_adj.size(), index_bytes / 1e6);
    std::fflush(stdout);

    // ---- queries
    std::vector<double> df(n, INF), db(n, INF); std::vector<char> sf(n, 0), sb(n, 0); std::vector<int> tf_, tb_;
    long long settled_ch = 0;
    // query: (1) upward searches from s (up graph) and t (reversed down graph) over the contracted nodes, run to
    // completion, which give exact distances to the core nodes where they enter it (core nodes are not expanded);
    // (2) a bidirectional Dijkstra on the core from those entry distances, stopped when the two queue minima add up
    // to the best meeting distance. A shortest path goes up to the core, through it, and down, or meets below it.
    auto ch_query = [&](int s, int t) {
        for (int x : tf_) { df[x] = INF; sf[x] = 0; } for (int x : tb_) { db[x] = INF; sb[x] = 0; }
        tf_.clear(); tb_.clear();
        if (s == t) return 0.0;
        double best = INF;
        auto upward = [&](int src, std::vector<double>& d, std::vector<char>& st, std::vector<int>& tch, const std::vector<long long>& off, const std::vector<E>& adj) {
            MinQ q; d[src] = 0; tch.push_back(src); q.push({0, src});
            while (!q.empty()) {
                QN tp = q.top(); q.pop();
                int x = tp.second;
                if (st[x] || rank[x] == n) continue;  // core nodes wait for phase 2
                st[x] = 1; ++settled_ch;
                for (long long i = off[x]; i < off[x + 1]; ++i) {
                    double nd = tp.first + adj[i].w; int y = adj[i].v;
                    if (nd < d[y]) { if (d[y] == INF) tch.push_back(y); d[y] = nd; q.push({nd, y}); }
                }
            }
        };
        upward(s, df, sf, tf_, u_off, u_adj);
        upward(t, db, sb, tb_, d_off, d_adj);
        for (int x : tf_) if (db[x] < INF) best = std::min(best, df[x] + db[x]);  // meeting below the core (or at a core entry)
        MinQ qf, qb;
        for (int x : tf_) if (rank[x] == n) qf.push({df[x], x});
        for (int x : tb_) if (rank[x] == n) qb.push({db[x], x});
        while (!qf.empty() && !qb.empty() && qf.top().first + qb.top().first < best) {
            bool fwd = qf.top().first <= qb.top().first;
            MinQ& q = fwd ? qf : qb; auto& d = fwd ? df : db; auto& st = fwd ? sf : sb; auto& tch = fwd ? tf_ : tb_;
            const auto& other = fwd ? db : df; const auto& off = fwd ? u_off : d_off; const auto& adj = fwd ? u_adj : d_adj;
            QN tp = q.top(); q.pop();
            int x = tp.second;
            if (st[x] || tp.first > d[x]) continue;
            st[x] = 1; ++settled_ch;
            for (long long i = off[x]; i < off[x + 1]; ++i) {  // a core node's edges in the search graphs lead to core nodes
                double nd = tp.first + adj[i].w; int y = adj[i].v;
                if (nd < d[y]) { if (d[y] == INF) tch.push_back(y); d[y] = nd; q.push({nd, y}); if (other[y] < INF) best = std::min(best, nd + other[y]); }
            }
        }
        return best;
    };
    std::vector<double> pd(n, INF); std::vector<char> pdone(n, 0); std::vector<int> ptouch;
    long long settled_p2p = 0;
    auto p2p = [&](int s, int t) {
        for (int x : ptouch) { pd[x] = INF; pdone[x] = 0; }
        ptouch.clear();
        MinQ q; pd[s] = 0; ptouch.push_back(s); q.push({0, s});
        while (!q.empty()) {
            QN tp = q.top(); q.pop();
            int x = tp.second;
            if (pdone[x]) continue;
            pdone[x] = 1; ++settled_p2p;
            if (x == t) return tp.first;
            for (long long i = o_off[x]; i < o_off[x + 1]; ++i) {
                double nd = tp.first + o_adj[i].w; int y = o_adj[i].v;
                if (nd < pd[y]) { if (pd[y] == INF) ptouch.push_back(y); pd[y] = nd; q.push({nd, y}); }
            }
        }
        return INF;
    };
    // answers and per-query times (one pass each), then the best of three timed passes for the CH (warm)
    std::vector<double> a_ch(m), a_p2p(m), t_ch(m), t_p2p(m);
    for (size_t i = 0; i < m; ++i) { auto t0 = clk::now(); a_p2p[i] = p2p(qp[i].first, qp[i].second); t_p2p[i] = since(t0); }
    for (size_t i = 0; i < m; ++i) { auto t0 = clk::now(); a_ch[i] = ch_query(qp[i].first, qp[i].second); t_ch[i] = since(t0); }
    double ch_best = INF; volatile double sink = 0;
    for (int rep = 0; rep < 3; ++rep) { auto t0 = clk::now(); for (auto& q : qp) sink = sink + ch_query(q.first, q.second); ch_best = std::min(ch_best, since(t0) / m); }
    double diff = 0, emean = 0, emax = 0, emean_p = 0, emax_p = 0;
    std::vector<double> errs;
    for (size_t i = 0; i < m; ++i) {
        diff = std::max(diff, std::fabs(a_ch[i] - a_p2p[i]) / std::max(a_p2p[i], 1e-300));
        double e = std::fabs(a_ch[i] / tr[i] - 1), ep = std::fabs(a_p2p[i] / tr[i] - 1);
        emean += e / m; emax = std::max(emax, e); emean_p += ep / m; emax_p = std::max(emax_p, ep); errs.push_back(e);
    }
    auto stats = [&](std::vector<double> v) { std::sort(v.begin(), v.end()); double s = 0; for (double x : v) s += x; return std::make_pair(1e3 * s / v.size(), 1e3 * v[(size_t)(0.99 * (v.size() - 1))]); };
    auto sc_ = stats(t_ch), sp_ = stats(t_p2p);
    std::sort(errs.begin(), errs.end());
    std::printf("queries %zu: CH %.3f ms (best of 3; one pass mean %.3f ms, p99 %.3f ms, %.0f settled), point-to-point Dijkstra %.3f ms (p99 %.3f ms, %.0f settled), "
                "answers differ by %.1e; error mean %.5f%% max %.4f%% (Dijkstra %.5f%% / %.4f%%)\n",
                m, 1e3 * ch_best, sc_.first, sc_.second, (double)settled_ch / m, sp_.first, sp_.second, (double)settled_p2p / m, diff,
                100 * emean, 100 * emax, 100 * emean_p, 100 * emax_p);
    // ---- hub labels: pruned landmark labelling (Akiba, Iwata, Yoshida 2013), directed: L_out(v) holds (hub, d(v, hub)),
    // L_in(v) holds (hub, d(hub, v)); one forward and one backward pruned Dijkstra per hub, hubs in the given order.
    // Exact on the graph, like the CH. Orders: degree (the usual one) or kd (separator strips of a kd split of the
    // chart first, level by level: the vertex-separator order that keeps labels small on planar-like graphs).
    std::ostringstream pjs;
    if (!pll_order.empty()) {
        auto t_pll = clk::now();
        std::vector<long long> r_off(n + 1, 0); std::vector<E> r_adj(o_adj.size());
        for (int u = 0; u < n; ++u) for (long long i = o_off[u]; i < o_off[u + 1]; ++i) ++r_off[o_adj[i].v + 1];
        for (int v = 0; v < n; ++v) r_off[v + 1] += r_off[v];
        { std::vector<long long> pos(r_off.begin(), r_off.end() - 1);
          for (int u = 0; u < n; ++u) for (long long i = o_off[u]; i < o_off[u + 1]; ++i) r_adj[pos[o_adj[i].v]++] = {u, o_adj[i].w}; }
        std::vector<int> ord;
        auto deg = [&](int v) { return (o_off[v + 1] - o_off[v]) + (r_off[v + 1] - r_off[v]); };
        if (pll_order == "kd") {
            std::vector<double> X(n), Y(n);
            { std::ifstream in(mesh_path); std::string head; long nv, nf, ne; in >> head >> nv >> nf >> ne;
              if (!in || nv != n) { std::fprintf(stderr, "--pll kd needs --mesh with the graph's %d vertices\n", n); return 1; }
              for (int i = 0; i < n; ++i) { double z; in >> X[i] >> Y[i] >> z; } }
            // strip half-width: half the 99th percentile of the chart length of the graph's edges, so that almost no
            // edge jumps a strip (pruning keeps the labels exact whatever the order)
            std::vector<double> el; el.reserve(o_adj.size() / 16 + 1);
            for (int u = 0; u < n; ++u) for (long long i = o_off[u]; i < o_off[u + 1]; i += 16) el.push_back(std::hypot(X[u] - X[o_adj[i].v], Y[u] - Y[o_adj[i].v]));
            std::sort(el.begin(), el.end());
            const double band = 0.5 * el[(size_t)(0.99 * (el.size() - 1))];
            std::vector<std::vector<int>> level{std::vector<int>(n)};
            for (int i = 0; i < n; ++i) level[0][i] = i;
            while (!level.empty()) {
                std::vector<std::vector<int>> next;
                for (auto& S : level) {
                    if (S.size() <= 2) { for (int v : S) ord.push_back(v); continue; }
                    double lo[2] = {INF, INF}, hi[2] = {-INF, -INF};
                    for (int v : S) { lo[0] = std::min(lo[0], X[v]); hi[0] = std::max(hi[0], X[v]); lo[1] = std::min(lo[1], Y[v]); hi[1] = std::max(hi[1], Y[v]); }
                    const auto& C = hi[0] - lo[0] >= hi[1] - lo[1] ? X : Y;
                    std::vector<double> c; for (int v : S) c.push_back(C[v]);
                    std::nth_element(c.begin(), c.begin() + c.size() / 2, c.end());
                    double med = c[c.size() / 2];
                    std::vector<int> sep, L, R;
                    for (int v : S) (std::fabs(C[v] - med) <= band ? sep : (C[v] < med ? L : R)).push_back(v);
                    if (L.empty() && R.empty()) { std::sort(sep.begin(), sep.end(), [&](int a, int b) { return deg(a) > deg(b); }); for (int v : sep) ord.push_back(v); continue; }
                    std::sort(sep.begin(), sep.end(), [&](int a, int b) { return deg(a) > deg(b); });
                    for (int v : sep) ord.push_back(v);
                    if (!L.empty()) next.push_back(std::move(L));
                    if (!R.empty()) next.push_back(std::move(R));
                }
                level.swap(next);
            }
        } else {
            ord.resize(n); for (int i = 0; i < n; ++i) ord[i] = i;
            std::sort(ord.begin(), ord.end(), [&](int a, int b) { return deg(a) > deg(b); });
        }
        struct LE { int h; double d; };
        std::vector<std::vector<LE>> Lout(n), Lin(n);
        std::vector<double> tmp(n, INF), pdist(n, INF); std::vector<int> ptch2; std::vector<char> pset(n, 0);
        long long entries = 0; int done_hubs = 0; bool pll_stopped = false;
        auto prune_search = [&](int h, int hr, const std::vector<long long>& off, const std::vector<E>& adj,
                                std::vector<std::vector<LE>>& own, std::vector<std::vector<LE>>& target) {
            for (auto& e : own[h]) tmp[e.h] = e.d;  // own = the labels of h on the other side
            MinQ q; pdist[h] = 0; ptch2.push_back(h); q.push({0, h});
            while (!q.empty()) {
                QN tp = q.top(); q.pop();
                int x = tp.second;
                if (pset[x] || tp.first > pdist[x]) continue;
                pset[x] = 1;
                bool prune = false;
                for (auto& e : target[x]) if (tmp[e.h] + e.d <= tp.first) { prune = true; break; }
                if (prune) continue;
                target[x].push_back({hr, tp.first}); ++entries;
                for (long long i = off[x]; i < off[x + 1]; ++i) {
                    double nd = tp.first + adj[i].w; int y = adj[i].v;
                    if (!pset[y] && nd < pdist[y]) { if (pdist[y] == INF) ptch2.push_back(y); pdist[y] = nd; q.push({nd, y}); }
                }
            }
            for (int x : ptch2) { pdist[x] = INF; pset[x] = 0; }
            ptch2.clear();
            for (auto& e : own[h]) tmp[e.h] = INF;
        };
        auto last = clk::now();
        for (int r = 0; r < n; ++r) {
            int h = ord[r];
            prune_search(h, r, o_off, o_adj, Lout, Lin);  // forward from h: d(h, x), into L_in(x)
            prune_search(h, r, r_off, r_adj, Lin, Lout);  // backward to h: d(x, h), into L_out(x)
            ++done_hubs;
            if (since(last) > 600) {
                last = clk::now();
                std::printf("  hub labels: %d of %d hubs, %.1f entries per vertex, %.2f h\n", r + 1, n, (double)entries / n, since(t_pll) / 3600);
                std::fflush(stdout);
            }
            if (since(t_pll) > pll_time_limit_h * 3600) { pll_stopped = true; break; }
        }
        double pll_s = since(t_pll);
        const long long pll_bytes = entries * 8;  // a hub id and a float32 distance per entry
        auto hl = [&](int s, int t) {
            if (s == t) return 0.0;
            const auto &A = Lout[s], &B = Lin[t]; double best = INF; size_t i = 0, j = 0;
            while (i < A.size() && j < B.size()) {
                if (A[i].h == B[j].h) { best = std::min(best, A[i].d + B[j].d); ++i; ++j; }
                else if (A[i].h < B[j].h) ++i; else ++j;
            }
            return best;
        };
        double hl_best = INF, hdiff = 0, hemean = 0, hemax = 0; volatile double hs = 0;
        for (size_t i = 0; i < m; ++i) {
            double a = hl(qp[i].first, qp[i].second);
            hdiff = std::max(hdiff, std::fabs(a - a_p2p[i]) / std::max(a_p2p[i], 1e-300));
            double e = std::fabs(a / tr[i] - 1); hemean += e / m; hemax = std::max(hemax, e);
        }
        for (int rep = 0; rep < 4; ++rep) { auto t0 = clk::now(); for (auto& q : qp) hs = hs + hl(q.first, q.second); if (rep) hl_best = std::min(hl_best, since(t0) / m); }
        std::printf("hub labels (%s order): %s %d of %d hubs in %.1f s, %.1f entries per vertex (%.1f MB), query %.3f us, answers differ from Dijkstra by %.1e, error mean %.5f%% max %.4f%%\n",
                    pll_order.c_str(), pll_stopped ? "STOPPED at" : "all", done_hubs, n, pll_s, (double)entries / n, pll_bytes / 1e6, 1e6 * hl_best, hdiff, 100 * hemean, 100 * hemax);
        std::fflush(stdout);
        pjs << ",\"pll\":{\"order\":\"" << pll_order << "\",\"complete\":" << (pll_stopped ? "false" : "true") << ",\"hubs_done\":" << done_hubs
            << ",\"build_s\":" << pll_s << ",\"entries\":" << entries << ",\"entries_per_vertex\":" << (double)entries / n << ",\"index_bytes\":" << pll_bytes
            << ",\"query_us\":" << 1e6 * hl_best << ",\"vs_p2p_max_rel_diff\":" << hdiff << ",\"err_mean\":" << hemean << ",\"err_max\":" << hemax << "}";
    }
    if (!json.empty()) {
        std::ofstream o(json);
        o << "{\"tool\":\"dgg_ch\",\"graph\":\"" << argv[1] << "\",\"nodes\":" << n << ",\"edges\":" << m0 << ",\"load_s\":" << load_s
          << ",\"settle_limit\":" << settle_limit << ",\"core_degree\":" << core_degree << ",\"core_factor\":" << core_factor << ",\"time_limit_h\":" << time_limit_h
          << ",\"contract_s\":" << contract_s << ",\"build_s\":" << build_s << ",\"stopped\":\"" << (stopped_time ? "time" : (core ? "core_degree" : "complete"))
          << "\",\"contracted\":" << n - core << ",\"core\":" << core << ",\"shortcuts\":" << shortcuts_added
          << ",\"up_edges\":" << u_adj.size() << ",\"down_edges\":" << d_adj.size() << ",\"index_bytes\":" << index_bytes
          << ",\"queries\":" << m << ",\"ch_query_ms\":" << 1e3 * ch_best << ",\"ch_query_ms_mean\":" << sc_.first << ",\"ch_query_ms_p99\":" << sc_.second
          << ",\"ch_settled_mean\":" << (double)settled_ch / m << ",\"p2p_query_ms\":" << sp_.first << ",\"p2p_query_ms_p99\":" << sp_.second
          << ",\"p2p_settled_mean\":" << (double)settled_p2p / m << ",\"ch_vs_p2p_max_rel_diff\":" << diff
          << ",\"err_mean\":" << emean << ",\"err_p99\":" << errs[(size_t)(0.99 * (m - 1))] << ",\"err_max\":" << emax
          << ",\"p2p_err_mean\":" << emean_p << ",\"p2p_err_max\":" << emax_p << pjs.str() << ",\"peak_rss_kb\":" << peak_rss_kb() << "}\n";
    }
    return 0;
}
