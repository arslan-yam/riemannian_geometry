// VTP lives in its own namespace: its headers reuse Kirsanov's class names in namespace geodesic
#define geodesic vtpgeo
#define sscanf_s sscanf  // MSVC-only name in the OBJ reader, which is not used here
#include "stdafx.h"
#include "geodesic_mesh.h"
#include "geodesic_algorithm_exact.h"  // third_party/vtp
#undef geodesic
#include "solver.h"
#include <stdexcept>

namespace cci {

struct VTPSolver : Solver {
    vtpgeo::Mesh mesh;
    VTPSolver(const std::vector<double>& V, const std::vector<unsigned>& F) {
        std::vector<double> p(V);
        std::vector<unsigned> f(F);
        mesh.initialize_mesh_data(p, f);
        vtpgeo::GeodesicAlgorithmExact probe(&mesh); // runs VTP's mesh consistency checks (throws)
    }
    void dist(unsigned src, const std::vector<unsigned>& targets, double* out) override {
        // VTP keeps its per-sweep state on the mesh vertices (distance, wavefront state): reset it,
        // the window lists live in the algorithm object, which is fresh per sweep
        for (auto& v : mesh.vertices()) {
            v.geodesic_distance() = vtpgeo::GEODESIC_INF;
            v.state() = vtpgeo::Vertex::OUTSIDE;
            v.incident_face() = nullptr;
        }
        vtpgeo::GeodesicAlgorithmExact alg(&mesh);
        alg.propagate(src);
        for (size_t i = 0; i < targets.size(); ++i) out[i] = mesh.vertices()[targets[i]].geodesic_distance();
    }
    const char* name() const override { return "VTP"; }
};

// nullptr when VTP rejects the mesh (its checks threw); the caller falls back to MMP
Solver* make_vtp(const std::vector<double>& V, const std::vector<unsigned>& F) {
    try {
        return new VTPSolver(V, F);
    } catch (const std::exception&) {
        return nullptr;
    }
}

}  // namespace cci
