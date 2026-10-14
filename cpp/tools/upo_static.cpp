// Static driver for the UP-Oracle authors' code (github.com/yanyinzhao/UpdatedStructureTerrainCode, TKDE 2025):
// builds their oracles ONCE on one terrain (no update) for a POI set, answers many POI-pair queries
// and scores every answer against exact MMP (their exact_distance). Their own main() measures a single
// query (POI 0 -> 1) on a pre/post terrain pair.
// Compile with -I<their src> (see cpp/Makefile, target upo_static).
// usage: upo_static mesh.off poi_file|all|randomK eps [queries=200] [methods=UP,WSPD,EAR,RC,WAV,ESP]
#include "algorithms.h"
#include <random>
#include <sstream>

static double ms_since(std::chrono::high_resolution_clock::time_point t) {
    return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t).count();
}

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: upo_static mesh.off poi_file|all|randomK eps [queries] [methods]\n"; return 1; }
    std::string mesh_file = argv[1], poi_arg = argv[2];
    double epsilon = std::stod(argv[3]);
    int nq = argc > 4 ? std::stoi(argv[4]) : 200;
    std::string methods = argc > 5 ? argv[5] : "UP,WSPD,EAR,RC,WAV,ESP";
    auto want = [&](const std::string& m) { return ("," + methods + ",").find("," + m + ",") != std::string::npos; };

    std::vector<double> points; std::vector<unsigned> faces;
    geodesic::read_mesh_from_file(&mesh_file[0], points, faces);
    geodesic::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);
    int n = mesh.vertices().size();
    std::vector<int> poi;
    if (poi_arg == "all") { for (int i = 0; i < n; ++i) poi.push_back(i); }
    else if (poi_arg.rfind("random", 0) == 0) {
        int k = std::stoi(poi_arg.substr(6)); std::mt19937 rng(3); std::vector<int> all(n);
        for (int i = 0; i < n; ++i) all[i] = i;
        std::shuffle(all.begin(), all.end(), rng); poi.assign(all.begin(), all.begin() + k);
    } else { std::ifstream in(poi_arg); int k; in >> k; poi.resize(k); for (auto& p : poi) in >> p; }
    int poi_num = poi.size();
    // query pairs among POIs; exact truth by MMP
    std::mt19937 rng(11);
    std::vector<std::pair<int, int>> qp;
    while ((int)qp.size() < nq) { int s = rng() % poi_num, t = rng() % poi_num; if (s != t) qp.push_back({s, t}); }
    std::vector<double> truth(nq);
    for (int i = 0; i < nq; ++i) exact_distance(&mesh, poi, qp[i].first, qp[i].second, truth[i]);
    std::cout << "mesh " << mesh_file << " vertices " << n << " faces " << mesh.faces().size() << " POIs " << poi_num << " eps " << epsilon << " queries " << nq << std::endl;

    // JSON (--json out.json after the positional arguments), rewritten after every method
    std::string json;
    for (int i = 1; i + 1 < argc; ++i) if (std::string(argv[i]) == "--json") json = argv[i + 1];
    std::vector<std::string> rows;
    auto report = [&](const std::string& name, double build_ms, double size_bytes, const std::vector<double>& ans, const std::vector<double>& qms) {
        double em = 0, ex = 0, emin = 1e30, qm = 0, ea = 0;
        for (int i = 0; i < nq; ++i) { double e = ans[i] / truth[i] - 1; em += e; ea += std::fabs(e); ex = std::max(ex, e); emin = std::min(emin, e); qm += qms[i]; }
        if (!json.empty()) {
            std::ostringstream r;
            r << "{\"method\":\"" << name << "\",\"build_ms\":" << build_ms << ",\"index_bytes\":" << size_bytes << ",\"query_ms\":" << qm / nq
              << ",\"err_mean_signed\":" << em / nq << ",\"err_mean\":" << ea / nq << ",\"err_max\":" << ex << ",\"err_min\":" << emin << "}";
            rows.push_back(r.str());
            std::ofstream o(json);
            o << "{\"tool\":\"upo_static\",\"mesh\":\"" << mesh_file << "\",\"vertices\":" << n << ",\"pois\":" << poi_num << ",\"eps\":" << epsilon
              << ",\"queries\":" << nq << ",\"methods\":[";
            for (size_t k = 0; k < rows.size(); ++k) o << (k ? "," : "") << rows[k];
            o << "]}\n";
        }
        std::printf("RESULT %-10s build %10.1f ms  size %9.3f MB  query %10.4f ms  err mean %.4f%% max %.4f%% min %.4f%%\n",
                    name.c_str(), build_ms, size_bytes / 1e6, qm / nq, 100 * em / nq, 100 * ex, 100 * emin);
        std::fflush(stdout);
    };
    int MST_dummy = 1;
    (void)MST_dummy;

    if (want("UP")) { // complete graph of exact POI distances + hierarchy greedy spanner, query = spanner Dijkstra
        std::vector<std::vector<std::vector<int>>> fsi; std::vector<std::vector<double>> pd; std::vector<std::vector<bool>> pdc;
        std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> paths;
        std::vector<std::vector<double>> pdv(poi_num, std::vector<double>(n, 0));
        double t_cg = 0, mem = 0;
        pre_complete_graph_construction(poi_num, &mesh, poi, fsi, pd, pdc, paths, pdv, t_cg, mem, 1);
        std::unordered_map<int, double> pdm; std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pathm, spm;
        double cg_size = 0, cg_w = 0, t_hash = 0, hs_size = 0, hs_w = 0, t_hgs = 0, cg_size2 = 0; int cg_e = 0, hs_e = 0;
        get_pairwise_distance_and_path_poi_to_poi_map(pd, pdm, paths, pathm, cg_size, cg_e, cg_w, t_hash);
        Graph graph(poi_num);
        hierarchy_greedy_spanner(epsilon, graph, pd, paths, spm, hs_size, hs_e, hs_w, t_hgs, cg_size2);
        std::vector<double> ans(nq), qms(nq);
        for (int i = 0; i < nq; ++i) {
            std::vector<geodesic::SurfacePoint> path; double qt = 0;
            spanner_query(poi_num, graph, spm, qp[i].first, qp[i].second, ans[i], path, qt);
            qms[i] = qt;
        }
        report("UP-Oracle", t_cg + t_hash + t_hgs, hs_size, ans, qms);
    }
    if (want("WSPD")) {
        int id = 1; std::vector<GeoNode*> all; std::unordered_map<int, GeoNode*> gm; std::unordered_map<int, GeoPair*> gp; std::unordered_map<int, int> pm;
        std::vector<std::vector<double>> pdv(poi_num, std::vector<double>(n, 0));
        double t = 0, mem = 0, size = 0, w = 0; int e = 0;
        pre_or_post_WSPD_Oracle_and_pre_WSPD_UP_Oracle_construction(poi_num, &mesh, poi, epsilon, id, all, gm, gp, pm, pdv, false, t, mem, size, e, w, 1);
        std::vector<double> ans(nq), qms(nq);
        for (int i = 0; i < nq; ++i) {
            std::vector<geodesic::SurfacePoint> path; double qt = 0;
            pre_or_post_WSPD_Oracle_and_pre_WSPD_UP_Oracle_query(&mesh, id, all, gm, gp, qp[i].first, qp[i].second, qt, ans[i], path);
            qms[i] = qt;
        }
        report("WSPD-Oracle", t, size, ans, qms);
    }
    if (want("RC")) {
        int id = 1; std::vector<GeoNode*> all; std::unordered_map<int, GeoNode*> gm; std::unordered_map<int, GeoPair*> gp; std::unordered_map<int, int> pm;
        std::vector<std::vector<double>> pdv(poi_num, std::vector<double>(n, 0));
        double t = 0, mem = 0, size = 0, w = 0; int e = 0;
        pre_or_post_RC_TIN_Oracle_and_pre_RC_TIN_UP_Oracle_construction(poi_num, &mesh, poi, epsilon, id, all, gm, gp, pm, pdv, false, t, mem, size, e, w, 1);
        std::vector<double> ans(nq), qms(nq);
        for (int i = 0; i < nq; ++i) {
            std::vector<geodesic::SurfacePoint> path; double qt = 0;
            pre_or_post_RC_TIN_Oracle_and_pre_RC_TIN_UP_Oracle_query(&mesh, id, all, gm, gp, qp[i].first, qp[i].second, qt, ans[i], path);
            qms[i] = 2 * qt; // their RC query divides nanoseconds by 2e6; report milliseconds
        }
        report("RC-Oracle", t, size, ans, qms);
    }
    if (want("EAR")) {
        int sq = 2, id = 1;
        std::vector<int> hw; std::unordered_map<int, std::unordered_map<int, int>> hwb; std::vector<GeoNode*> all;
        std::unordered_map<int, GeoNode*> gm; std::unordered_map<int, GeoPair*> gp; std::unordered_map<int, int> hm;
        std::vector<std::vector<std::vector<int>>> fsi; std::vector<std::vector<double>> pd; std::vector<std::vector<bool>> pdc;
        std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> paths;
        std::vector<std::vector<double>> pdv(poi_num, std::vector<double>(n, 0));
        std::unordered_map<int, double> dph; std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pph;
        double t = 0, mem = 0, size = 0, w = 0; int e = 0;
        pre_or_post_EAR_Oracle_and_pre_EAR_UP_Oracle_construction(sq, poi_num, &mesh, poi, epsilon, hw, hwb, id, all, gm, gp, hm, fsi, pd, pdc, paths, pdv,
                                                                  dph, pph, false, t, mem, size, e, w);
        std::vector<double> ans(nq), qms(nq);
        for (int i = 0; i < nq; ++i) {
            std::vector<geodesic::SurfacePoint> path; double qt = 0;
            pre_or_post_EAR_Oracle_and_pre_EAR_UP_Oracle_query(&mesh, poi, sq, id, hwb, all, gm, gp, dph, pph, qp[i].first, qp[i].second, qt, ans[i], path);
            qms[i] = qt;
        }
        report("EAR-Oracle", t, size, ans, qms);
    }
    if (want("WAV")) { // exact on-the-fly (MMP with path)
        std::vector<double> ans(nq), qms(nq);
        for (int i = 0; i < nq; ++i) {
            std::vector<geodesic::SurfacePoint> path; double qt = 0, mem = 0;
            pre_or_post_WAV_Fly_Algo_query(&mesh, poi, qp[i].first, qp[i].second, qt, mem, ans[i], path);
            qms[i] = qt;
        }
        report("WAV-Fly", 0, 0, ans, qms);
    }
    if (want("ESP")) { // eps-approximate on-the-fly (Steiner subdivision)
        std::vector<double> ans(nq), qms(nq);
        for (int i = 0; i < nq; ++i) {
            std::vector<geodesic::SurfacePoint> path; double qt = 0, mem = 0;
            pre_or_post_ESP_Fly_Algo_query(&mesh, poi, epsilon, qp[i].first, qp[i].second, qt, mem, ans[i], path);
            qms[i] = qt;
        }
        report("ESP-Fly", 0, 0, ans, qms);
    }
    // the pairs, for our index to answer the same queries (vertex ids)
    std::ofstream pf("upo_pairs.txt");
    for (int i = 0; i < nq; ++i) pf << poi[qp[i].first] << " " << poi[qp[i].second] << " " << std::setprecision(17) << truth[i] << "\n";
    return 0;
}
