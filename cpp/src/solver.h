// Exact geodesic solvers behind one interface: Kirsanov's MMP (the code pygeodesic wraps) and
// VTP (Qin et al., SIGGRAPH 2016). Both are exact polyhedral geodesic algorithms.
#pragma once
#include <vector>

namespace cci {

struct Solver {
    virtual ~Solver() {}
    // exact geodesic distance from vertex src to every vertex in targets (out[i] for targets[i])
    virtual void dist(unsigned src, const std::vector<unsigned>& targets, double* out) = 0;
    virtual const char* name() const = 0;
};

// early_stop: MMP stops once every target's distance is final (Kirsanov's stop rule with a zero
// distance bound); the distances are identical to a full sweep
Solver* make_mmp(const std::vector<double>& V, const std::vector<unsigned>& F, bool early_stop);
// VTP sweeps the whole (sub-)mesh: its vertex queue does not make a popped vertex final, so there
// is no exact early stop
Solver* make_vtp(const std::vector<double>& V, const std::vector<unsigned>& F);

}  // namespace cci
