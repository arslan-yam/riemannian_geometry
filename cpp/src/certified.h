// Certified Cut Index in C++: the same construction and query as rgeo/oracles/certified.py
// (theorem and proofs in notes/guarantee.md), with the label sweeps done by Kirsanov's MMP (early
// stop at the members), by VTP (full sweeps of the label sub-mesh), or by a per-cut hybrid.
#pragma once
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "graph.h"
#include "solver.h"
#include "terrain.h"

namespace cci {

inline float f32_up(double d) {
    float x = (float)d;
    if ((double)x < d) x = std::nextafter(x, std::numeric_limits<float>::infinity());
    return x;
}
constexpr double F32_REL = 1.0 / (1 << 22); // two float32 ulps, relative

struct Run {
    std::vector<double> ub, ab; // breakpoints along the line, cumulative arc length
    double smax = 0;
};

struct Cut {
    int axis = 0, depth = 0, hub = -1;
    double v = 0, lo[2], hi[2], mu = 0;
    std::vector<int> M;              // members
    std::vector<Run> runs;
    std::vector<double> pos, arc, gap, X; // portal positions on the line, arc positions, gaps (inf across runs), chart coords
    std::vector<int64_t> run;
    std::vector<int64_t> end_idx;    // run-end stretches: portal index, 2 x stretch arc
    std::vector<double> end_gap;
    bool lb_ok = true;
    std::vector<float> L;            // |M| x K labels, rounded up
    std::vector<float> foot_r, foot_arc;
    std::vector<int64_t> foot_run;
    int submesh_faces = 0;
    bool used_vtp = false;
    int K() const { return (int)pos.size(); }
};

struct Options {
    double eps = 0.1, gap_factor = 0.5, gap = 0, hub_slack = 4.0;
    int interp = 2;                 // 0 none, 1 parabola, 2 model
    std::string solver = "hybrid";  // mmp | vtp | hybrid
    double hybrid_ratio = 2.0;      // hybrid: VTP when sub-mesh faces <= ratio x faces of the members' box
    bool verbose = false;
};

class CertifiedIndex {
  public:
    const Terrain& T;
    Options opt;
    double g = 0, Rn = 0, Rq = 0, lam = 1, eps = 0.1;
    bool exact = true; // labels are exact geodesic distances (MMP / VTP)
    FaceGrid grid;
    PointGrid pgrid;
    Graph G;
    std::vector<double> ucoord[2];
    int D = 0; // tree depth (columns of anc / row)
    std::vector<int32_t> anc, row, depth_of;
    std::vector<Cut> cuts;
    std::vector<std::vector<int>> leaves;
    std::vector<int64_t> near_ptr;
    std::vector<int32_t> near_ids;
    std::vector<float> near_d;
    std::map<std::string, double> timing;
    int64_t sweeps_mmp = 0, sweeps_vtp = 0;

    CertifiedIndex(const Terrain& t, const Options& o);
    void bounds(int s, int t, double* out) const; // LB, UB, estimate
    size_t nbytes() const;
    int64_t portals() const { int64_t k = 0; for (auto& c : cuts) k += c.K(); return k; }
    int64_t labels() const { int64_t k = 0; for (auto& c : cuts) k += (int64_t)c.L.size(); return k; }

  private:
    void build_tree();
    double split_value(std::vector<double> c, int axis) const;
    std::vector<Run> line_pieces(int a, double v, double lo, double hi) const;
    std::vector<double> place(const Run& r, double r_lo, double r_hi) const;
    void cut_geometry();
    void labels_phase();
    void finalize_cut(Cut& C, std::vector<float>& Lc, std::vector<double>& Xu, int K);
    void near_phase();
    void feet_phase();
    void submesh(const std::vector<int>& faces, const std::vector<double>& X, const double* along,
                 std::vector<double>& V2, std::vector<unsigned>& F2, std::vector<int>& remap_used,
                 std::vector<int>& pidx, std::vector<double>& Xu) const;
    double interp(const Cut& C, const double* f, int K, int i, int rs, int rt) const;
};

}  // namespace cci
