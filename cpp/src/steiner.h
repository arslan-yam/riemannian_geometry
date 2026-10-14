// Steiner graph of a weighted triangle mesh (the skeleton of the weighted-region literature, Aleksandrov,
// Maheshwari, Sack, JACM 2005; rgeo/terrain/steiner.py): m points on every edge; along an edge the points
// are chained at the cheaper adjacent face cost, and inside every face each pair of its boundary points
// that do not lie on a common edge is joined by a straight segment costing w_f x 3-d length. A graph path
// is a surface curve, so graph distances are upper bounds of the weighted geodesic distance and converge
// to it as m grows. One graph is built per terrain and shared (read-only) by every search; a search keeps
// its state in a State (one per thread) and stops once every target is settled, so it only explores the
// ball of the farthest target.
#pragma once
#include <cstdint>
#include <utility>
#include <vector>

namespace cci {

struct SteinerGraph {
    int n = 0, m = 5, nf = 0;
    int64_t N = 0;
    std::vector<double> P;           // node coordinates: the mesh vertices (ids 0..n-1), then m points per edge
    std::vector<int64_t> ptr;        // CSR
    std::vector<int32_t> adj;
    std::vector<double> wt;
    std::vector<int32_t> fpts;       // nf x (3 + 3m): the boundary nodes of every face
    std::vector<double> fw;          // face weights

    struct State {
        std::vector<double> D;
        std::vector<int32_t> touched;
        std::vector<char> done, want;
    };
    using Seeds = std::vector<std::pair<int32_t, double>>;

    void build(const std::vector<double>& V, const std::vector<unsigned>& F, const std::vector<double>& W, int m);
    int B() const { return 3 + 3 * m; }
    void init(State& S) const;
    // Dijkstra from weighted seeds; stops once every node of `stop` is settled (S.D valid until the next search)
    void search(State& S, const Seeds& seeds, const std::vector<int32_t>& stop) const;
    // a point (x, y, z) of face f: seeds at the face's boundary nodes (w_f x straight distance)
    void point_seeds(int f, double x, double y, double z, Seeds& out) const;
    // after a search that settled the boundary nodes of face f: the cost to reach (x, y, z) in f through them
    double point_dist(const State& S, int f, double x, double y, double z) const;
    size_t bytes() const { return P.size() * 8 + ptr.size() * 8 + adj.size() * 4 + wt.size() * 8 + fpts.size() * 4 + fw.size() * 8; }
};

}  // namespace cci
