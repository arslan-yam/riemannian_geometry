// Exact geodesic solvers behind one interface: Kirsanov's MMP (the code pygeodesic wraps) and
// VTP (Qin et al., SIGGRAPH 2016). Both are exact polyhedral geodesic algorithms.
#pragma once
#include <vector>

namespace cci {

// a point on a face of the solver's mesh (x, y, z on the face's plane)
struct SurfPt { unsigned face; double x, y, z; };

struct Solver {
    virtual ~Solver() {}
    // exact geodesic distance from vertex src to every vertex in targets (out[i] for targets[i])
    virtual void dist(unsigned src, const std::vector<unsigned>& targets, double* out) = 0;
    // exact geodesic distance from a surface point to surface points; false if the solver cannot (VTP)
    virtual bool point_dists(const SurfPt&, const std::vector<SurfPt>&, double*) { return false; }
    virtual const char* name() const = 0;
};

// early_stop: MMP stops once every target's distance is final (Kirsanov's stop rule with a zero
// distance bound); the distances are identical to a full sweep
Solver* make_mmp(const std::vector<double>& V, const std::vector<unsigned>& F, bool early_stop);
// VTP sweeps the whole (sub-)mesh: its vertex queue does not make a popped vertex final, so there
// is no exact early stop. Returns nullptr if VTP's own mesh checks reject the mesh
Solver* make_vtp(const std::vector<double>& V, const std::vector<unsigned>& F);
// weighted terrains: Dijkstra on the Steiner graph with m points per edge (face weights W, one per face);
// upper bounds of the weighted geodesic distance (no exact solver exists for weighted regions)
Solver* make_steiner(const std::vector<double>& V, const std::vector<unsigned>& F, const std::vector<double>& W, int m);
double steiner_bytes_per_face(int m);

}  // namespace cci
