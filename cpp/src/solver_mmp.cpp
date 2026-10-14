#include "solver.h"
#include "geodesic_algorithm_exact.h"  // third_party/kirsanov

namespace cci {

struct MMPSolver : Solver {
    geodesic::Mesh mesh;
    geodesic::GeodesicAlgorithmExact* alg = nullptr;
    bool early;
    MMPSolver(const std::vector<double>& V, const std::vector<unsigned>& F, bool early_stop) : early(early_stop) {
        std::vector<double> p(V);
        std::vector<unsigned> f(F);
        mesh.initialize_mesh_data(p, f);
        alg = new geodesic::GeodesicAlgorithmExact(&mesh);
    }
    ~MMPSolver() override { delete alg; }
    void dist(unsigned src, const std::vector<unsigned>& targets, double* out) override {
        std::vector<geodesic::SurfacePoint> sources(1, geodesic::SurfacePoint(&mesh.vertices()[src]));
        std::vector<geodesic::SurfacePoint> stops;
        stops.reserve(targets.size());
        for (unsigned t : targets) stops.push_back(geodesic::SurfacePoint(&mesh.vertices()[t]));
        alg->propagate(sources, early ? 0.0 : geodesic::GEODESIC_INF, early ? &stops : nullptr);
        for (size_t i = 0; i < stops.size(); ++i) alg->best_source(stops[i], out[i]);
    }
    // surface points: Kirsanov's stop rule takes each target's closest vertex v and stops once the queue
    // is beyond d(v) + |target - v| >= d(target), so the early stop stays exact
    bool point_dists(const SurfPt& src, const std::vector<SurfPt>& targets, double* out) override {
        auto sp = [&](const SurfPt& p) { return geodesic::SurfacePoint(&mesh.faces()[p.face], p.x, p.y, p.z); };
        std::vector<geodesic::SurfacePoint> sources(1, sp(src));
        std::vector<geodesic::SurfacePoint> stops;
        stops.reserve(targets.size());
        for (const auto& t : targets) stops.push_back(sp(t));
        alg->propagate(sources, early ? 0.0 : geodesic::GEODESIC_INF, early ? &stops : nullptr);
        for (size_t i = 0; i < stops.size(); ++i) alg->best_source(stops[i], out[i]);
        return true;
    }
    const char* name() const override { return early ? "MMP (early stop)" : "MMP"; }
};

Solver* make_mmp(const std::vector<double>& V, const std::vector<unsigned>& F, bool early_stop) {
    return new MMPSolver(V, F, early_stop);
}

}  // namespace cci
