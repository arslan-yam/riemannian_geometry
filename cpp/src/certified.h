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
    std::vector<int32_t> run_beg;    // first portal of every run, then K
    std::vector<double> lbsub;       // LB term j >= 1: 0.5 (f[j-1] + f[j]) - lbsub[j]; -inf when j-1, j are in different runs
    std::vector<int64_t> end_idx;    // run-end stretches: portal index, 2 x stretch arc
    std::vector<double> end_gap;
    bool lb_ok = true;
    std::vector<float> L;            // |M| x K labels, rounded up
    std::vector<float> foot_r, foot_arc;
    std::vector<int64_t> foot_run;
    int submesh_faces = 0;
    bool used_vtp = false;
    double gap_box = 0;              // portal gap inside the box (fixed g, capped, or scale-dependent)
    double cell_side = 0;            // shifted hierarchy: side of the cell this line splits
    // 16-bit coding (encode16): closed-loop deltas along each row, an anchor per run, decoded >= L
    std::vector<int16_t> Lq;         // |M| x K
    std::vector<int8_t> Lq8;         // |M| x K (8-bit coding)
    std::vector<float> La;           // |M| x (runs with portals): the anchors, rounded up
    std::vector<double> qstep, qerr; // per portal: quantum, and bound of (decoded - true) per label
    std::vector<double> lbsub_q;     // lbsub + the decoding errors of the two labels (coded storages)
    int n_anchor = 0;
    // Thorup-style covers (build_covers): per row a sorted list of (portal, label)
    std::vector<int32_t> cov_ptr;    // |M| + 1
    std::vector<uint16_t> cov_idx;
    std::vector<float> cov_val;
    int K() const { return (int)pos.size(); }
};

// one regular, shifted quadtree (the shifted hierarchy); node ids per level and vertex
struct ShiftNode {
    int cut[2] = {-1, -1};           // cut ids of the vertical (axis 0) and horizontal (axis 1) line
    double side = 0, mid[2] = {0, 0};
};
struct ShiftCopy {
    int D = 0;
    std::vector<int32_t> anc, row;   // n x D: node at each level, row of the vertex in the node's members
    std::vector<int32_t> depth_of;
    std::vector<ShiftNode> nodes;
};

struct Options {
    double eps = 0.1, gap_factor = 0.5, gap = 0, hub_slack = 4.0;
    int interp = 2;                 // 0 none, 1 parabola, 2 model
    std::string solver = "hybrid";  // mmp | vtp | hybrid
    double hybrid_ratio = 2.0;      // hybrid: VTP when sub-mesh faces <= ratio x faces of the members' box
    bool verbose = false;
    int k_max = 0;                  // version A: at most k_max portals per cut inside the box (0 = no cap)
    std::string hierarchy = "kd";   // kd (median cuts) | shifted (3 shifted quadtrees, gap max(g, eps lam S / 6))
    int copies = 3;                 // shifted: number of shifted copies (3 for the 2-d shifting lemma)
    int threads = 1;                // label sweeps, near tables and cut geometry on this many threads
                                    // (results identical to one thread; benchmarks use one)
    std::string store = "dense";    // dense | code16 | code8: code each cut's labels as soon as it is built
                                    // and drop its float32 labels (peak memory ~ the coded index)
    std::vector<int> poi;           // POI mode: index only these vertices (queries between them); empty = all
    std::string near_store = "sym";  // dense | sym | sym16 | sym8: near tables as built (int32 id + float32),
                                    // or each pair once with varint id deltas and the float32 distance (sym,
                                    // lossless) or a 16 / 8-bit code of its ratio to the chord, rounded up
                                    // (sym16, sym8); encoded after the build, the dense tables dropped
};

enum Storage { DENSE = 0, CODE16 = 1, COVER = 2, CODE8 = 3 };

// an arbitrary surface point: chart coordinates inside face f
struct APoint { int f; double x, y; };

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
    int64_t near_total = 0;          // near entries as built (both directions)
    // compact near tables (encode_near): row of slot(u) holds the neighbours t > u only, as a byte stream
    // of (varint id delta, distance) entries; distance = float32 (NEAR_SYM) or a uint16 / uint8 code c
    // with decoded value chord(u,t) (1 + c nz_q[row]) >= the float32 distance (NEAR_SYM16 / NEAR_SYM8)
    enum NearMode { NEAR_DENSE = 0, NEAR_SYM = 1, NEAR_SYM16 = 2, NEAR_SYM8 = 3 };
    static NearMode near_mode_of(const std::string& name); // "dense" | "sym" | "sym16" | "sym8"
    NearMode near_mode = NEAR_DENSE;
    std::vector<uint32_t> nz_ptr;    // NP + 1 byte offsets
    std::vector<uint8_t> nz_bytes;
    std::vector<float> nz_q;         // NP relative quanta (sym16, sym8)
    std::map<std::string, double> timing;
    int64_t sweeps_mmp = 0, sweeps_vtp = 0;

    std::vector<int32_t> slot_of;    // vertex -> row in anc / row / depth_of / near tables (-1: not indexed)
    std::vector<int> poi_list;       // POI mode: the indexed vertices, ascending
    int NP = 0;                      // indexed vertices
    std::vector<ShiftCopy> copies;   // shifted hierarchy
    Storage storage = DENSE;
    double cover_eps = 0;
    int64_t cover_entries = 0;

    CertifiedIndex(const Terrain& t, const Options& o);
    void bounds(int s, int t, double* out) const; // LB, UB, estimate, with the current storage
    // certify-or-refine: if UB > (1 + eps) LB the pair is computed exactly (local MMP sweep on the
    // ellipse sub-mesh); returns true when it refined
    bool bounds_refined(int s, int t, double* out) const;
    double exact_pair(int s, int t, double ub) const;
    // exact distances from vertex s to vertices T (each <= ub): one early-stopping MMP sweep on the faces
    // meeting the bounding box of s and T dilated by ub / (2 lambda) (ellipse lemma)
    void exact_from(int s, const std::vector<int>& T, double ub, double* out) const;
    // arbitrary points (A2A), certified: bounds from the vertex certificates of the two faces' corners;
    // pairs within the near radius, and pairs whose interval is wider than UB > (1 + 2 eps) LB, are
    // computed exactly (MMP from the surface point on the ellipse sub-mesh). Returns true if refined
    // because the interval was wide (near pairs are always exact and not counted)
    bool a2a_bounds(const APoint& p, const APoint& q, double* out) const;
    double exact_points(const APoint& p, const APoint& q, double ub) const;
    void encode16(int bits = 16);    // builds the 16-bit (or 8-bit) coded labels (keeps the dense ones);
                                     // cuts coded during the build (Options::store) are left as they are
    void build_covers(double eps_c); // builds eps_c-covers per row (keeps the dense labels; needs them)
    // builds the compact near tables (NEAR_SYM, NEAR_SYM16 or NEAR_SYM8) and switches the queries to them; the dense
    // tables are kept only if keep_dense (to evaluate several modes on one build)
    void encode_near(NearMode mode, bool keep_dense);
    void use_near(NearMode mode) { near_mode = mode; } // switch between modes that are present
    size_t near_bytes() const;       // bytes of the near tables in the current mode
    size_t nbytes() const;           // bytes of what the current storage mode reads
    int64_t portals() const { int64_t k = 0; for (auto& c : cuts) k += c.K(); return k; }
    int64_t labels() const { int64_t k = 0; for (auto& c : cuts) k += (int64_t)c.M.size() * c.K(); return k; }

  private:
    void build_tree();
    void build_shifted();
    double split_value(std::vector<double> c, int axis) const;
    std::vector<Run> line_pieces(int a, double v, double lo, double hi) const;
    std::vector<double> place(const Run& r, double r_lo, double r_hi, double gap) const;
    bool select_cut(int s, int t, int& cid, int& rs, int& rt) const;
    void eval_cut(const Cut& C, int s, int t, int rs, int rt, double* out) const;
    void cut_geometry();
    void labels_phase();
    // one cut's label sub-mesh (portals inserted), its members' local ids and the solver choice
    struct CutWork {
        std::vector<double> V2, Xu;
        std::vector<unsigned> F2, tgt;
        std::vector<int> used, pidx;
        size_t nfaces = 0;
        bool vtp = false;
    };
    void prepare_cut(const Cut& C, Dijkstra& dij, CutWork& W) const;
    Solver* cut_solver(const CutWork& W, bool& vtp) const;
    void finish_cut(Cut& C, std::vector<float>& Lc, CutWork& W, int K);
    bool encode_cut(Cut& C, int bits);
    void finalize_cut(Cut& C, std::vector<float>& Lc, std::vector<double>& Xu, int K);
    void near_phase();
    void feet_phase();
    void submesh(const std::vector<int>& faces, const std::vector<double>& X, const double* along,
                 std::vector<double>& V2, std::vector<unsigned>& F2, std::vector<int>& remap_used,
                 std::vector<int>& pidx, std::vector<double>& Xu) const;
    double interp(const Cut& C, const double* f, int K, int i, int rs, int rt) const;
    bool near_lookup(int s, int t, double* out) const; // near pair from the tables of the current mode
    double near_chord(int u, int v) const;
};

}  // namespace cci
