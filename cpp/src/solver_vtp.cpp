// VTP lives in its own namespace: its headers reuse Kirsanov's class names in namespace geodesic.
// third_party/vtp carries two fixes to the upstream code (vertex queue ties, the pairwise window check): notes/vtp_exactness.md
#define geodesic vtpgeo
#define sscanf_s sscanf  // MSVC-only name in the OBJ reader, which is not used here
#include "stdafx.h"
#include "geodesic_mesh.h"
#include "geodesic_algorithm_exact.h"  // third_party/vtp
#undef geodesic
#include "solver.h"
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>

namespace cci {

// a copy of the mesh with surface points inserted as vertices: each point splits the triangle containing
// it into three (a point on an edge is first moved 1e-6 of the way to the triangle's centroid; a point on
// a vertex is that vertex). ids[i] is the vertex of pts[i]
static void insert_points(std::vector<double>& V, std::vector<unsigned>& F, const std::vector<SurfPt>& pts, std::vector<unsigned>& ids) {
    std::map<unsigned, std::vector<unsigned>> parts; // original face -> its current triangles
    ids.clear();
    for (const SurfPt& p0 : pts) {
        auto& L = parts[p0.face];
        if (L.empty()) L.push_back(p0.face);
        double P[3] = {p0.x, p0.y, p0.z};
        auto bary = [&](unsigned f, const double* Q, double* l) {
            const double *A = &V[3 * F[3 * f]], *B = &V[3 * F[3 * f + 1]], *C = &V[3 * F[3 * f + 2]];
            double e1[3], e2[3], n[3];
            for (int k = 0; k < 3; ++k) { e1[k] = B[k] - A[k]; e2[k] = C[k] - A[k]; }
            n[0] = e1[1] * e2[2] - e1[2] * e2[1]; n[1] = e1[2] * e2[0] - e1[0] * e2[2]; n[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
            auto area = [&](const double* X, const double* Y) { // n . ((X - Q) x (Y - Q))
                double u[3], v[3];
                for (int k = 0; k < 3; ++k) { u[k] = X[k] - Q[k]; v[k] = Y[k] - Q[k]; }
                return n[0] * (u[1] * v[2] - u[2] * v[1]) + n[1] * (u[2] * v[0] - u[0] * v[2]) + n[2] * (u[0] * v[1] - u[1] * v[0]);
            };
            l[0] = area(B, C) / nn; l[1] = area(C, A) / nn; l[2] = area(A, B) / nn;
        };
        unsigned best = L[0]; double bl[3], bmin = -1e300;
        for (unsigned f : L) { double l[3]; bary(f, P, l); double m = std::min({l[0], l[1], l[2]}); if (m > bmin) { bmin = m; best = f; std::copy(l, l + 3, bl); } }
        int at_vertex = -1;
        for (int k = 0; k < 3; ++k) if (bl[k] > 1 - 1e-9) at_vertex = k;
        if (at_vertex >= 0) { ids.push_back(F[3 * best + at_vertex]); continue; }
        if (bmin < 1e-7) { // on an edge: move inside
            for (int k = 0; k < 3; ++k) {
                double g = (V[3 * F[3 * best] + k] + V[3 * F[3 * best + 1] + k] + V[3 * F[3 * best + 2] + k]) / 3;
                P[k] += 1e-6 * (g - P[k]);
            }
        }
        unsigned a = F[3 * best], b = F[3 * best + 1], c = F[3 * best + 2], v = (unsigned)(V.size() / 3);
        V.insert(V.end(), P, P + 3);
        F[3 * best + 2] = v;                                     // (a, b, v)
        unsigned f2 = (unsigned)(F.size() / 3);
        F.insert(F.end(), {b, c, v, c, a, v});                   // (b, c, v), (c, a, v): orientation kept
        L.push_back(f2); L.push_back(f2 + 1);
        ids.push_back(v);
    }
}

// VTP links its window lists by a search from the first edge, so it rejects a mesh in several pieces. With
// split_parts each edge-connected piece gets its own VTP mesh, and a target in another piece than the
// source is unreachable (INF). That is exact only if no vertex is shared by two pieces (a path could pass
// through it): such a mesh is rejected
struct VTPSolver : Solver {
    struct Part { vtpgeo::Mesh mesh; };
    std::vector<std::unique_ptr<Part>> parts;
    std::vector<int> part_of, local_of; // vertex -> its piece and its id there
    std::vector<double> V0; std::vector<unsigned> F0; // for point queries (points become vertices of a copy)
    bool split;
    VTPSolver(const std::vector<double>& V, const std::vector<unsigned>& F, bool split_parts) : V0(V), F0(F), split(split_parts) {
        const size_t n = V.size() / 3, m = F.size() / 3;
        if (!split) { // one VTP mesh, exactly as given
            part_of.assign(n, 0); local_of.resize(n);
            for (size_t i = 0; i < n; ++i) local_of[i] = (int)i;
            std::vector<double> p(V); std::vector<unsigned> f(F);
            parts.emplace_back(new Part);
            parts.back()->mesh.initialize_mesh_data(p, f);
            vtpgeo::GeodesicAlgorithmExact probe(&parts.back()->mesh); // VTP's mesh consistency checks (throw)
            return;
        }
        std::vector<int> comp(m, 0);
        int nc = 1;
        { // pieces: union of faces across shared edges
            std::vector<int> par(m);
            for (size_t f = 0; f < m; ++f) par[f] = (int)f;
            auto find = [&](int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
            std::map<std::pair<unsigned, unsigned>, int> first;
            for (size_t f = 0; f < m; ++f) for (int k = 0; k < 3; ++k) {
                unsigned u = F[3 * f + k], v = F[3 * f + (k + 1) % 3];
                auto key = std::make_pair(std::min(u, v), std::max(u, v));
                auto it = first.find(key);
                if (it == first.end()) first[key] = (int)f; else par[find((int)f)] = find(it->second);
            }
            std::map<int, int> id;
            for (size_t f = 0; f < m; ++f) { int r = find((int)f); auto it = id.find(r); comp[f] = it == id.end() ? (id[r] = (int)id.size()) : it->second; }
            nc = (int)id.size();
        }
        part_of.assign(n, -1); local_of.assign(n, -1);
        std::vector<std::vector<double>> pv(nc); std::vector<std::vector<unsigned>> pf(nc);
        for (size_t f = 0; f < m; ++f) {
            int c = comp[f];
            for (int k = 0; k < 3; ++k) {
                unsigned v = F[3 * f + k];
                if (part_of[v] < 0) { part_of[v] = c; local_of[v] = (int)(pv[c].size() / 3); pv[c].insert(pv[c].end(), &V[3 * v], &V[3 * v] + 3); }
                else if (part_of[v] != c) throw std::runtime_error("VTP: a vertex shared by two pieces");
                pf[c].push_back((unsigned)local_of[v]);
            }
        }
        for (int c = 0; c < nc; ++c) {
            parts.emplace_back(new Part);
            parts.back()->mesh.initialize_mesh_data(pv[c], pf[c]);
            vtpgeo::GeodesicAlgorithmExact probe(&parts.back()->mesh); // VTP's mesh consistency checks (throw)
        }
    }
    void dist(unsigned src, const std::vector<unsigned>& targets, double* out) override {
        int c = part_of[src];
        if (c < 0) { for (size_t i = 0; i < targets.size(); ++i) out[i] = targets[i] == src ? 0.0 : vtpgeo::GEODESIC_INF; return; }
        vtpgeo::Mesh& mesh = parts[c]->mesh;
        // VTP keeps its per-sweep state on the mesh vertices (distance, wavefront state): reset it,
        // the window lists live in the algorithm object, which is fresh per sweep
        for (auto& v : mesh.vertices()) {
            v.geodesic_distance() = vtpgeo::GEODESIC_INF;
            v.state() = vtpgeo::Vertex::OUTSIDE;
            v.incident_face() = nullptr;
        }
        vtpgeo::GeodesicAlgorithmExact alg(&mesh);
        alg.propagate((unsigned)local_of[src]);
        for (size_t i = 0; i < targets.size(); ++i)
            out[i] = part_of[targets[i]] == c ? mesh.vertices()[local_of[targets[i]]].geodesic_distance() : vtpgeo::GEODESIC_INF;
    }
    // surface points: inserted as vertices of a copy of the mesh, then one VTP sweep from the source's vertex
    // (false if VTP rejects the copy)
    bool point_dists(const SurfPt& src, const std::vector<SurfPt>& targets, double* out) override {
        std::vector<double> V(V0); std::vector<unsigned> F(F0), ids;
        std::vector<SurfPt> all(1, src);
        all.insert(all.end(), targets.begin(), targets.end());
        insert_points(V, F, all, ids);
        std::unique_ptr<VTPSolver> s;
        try { s.reset(new VTPSolver(V, F, split)); } catch (const std::exception&) { return false; }
        std::vector<unsigned> tg(ids.begin() + 1, ids.end());
        s->dist(ids[0], tg, out);
        return true;
    }
    const char* name() const override { return "VTP"; }
};

// nullptr when VTP rejects the mesh (its checks threw); the caller falls back to MMP
Solver* make_vtp(const std::vector<double>& V, const std::vector<unsigned>& F, bool split_parts) {
    try {
        return new VTPSolver(V, F, split_parts);
    } catch (const std::exception&) {
        if (const char* dir = std::getenv("CCI_DUMP_REJECTED")) { // debugging: the meshes VTP rejects
            static std::atomic<int> k{0};
            int i = k++;
            if (i < 20) { std::ofstream o(std::string(dir) + "/rej" + std::to_string(i) + ".off"); o.precision(17);
                o << "OFF\n" << V.size() / 3 << " " << F.size() / 3 << " 0\n";
                for (size_t j = 0; j < V.size(); j += 3) o << V[j] << " " << V[j + 1] << " " << V[j + 2] << "\n";
                for (size_t j = 0; j < F.size(); j += 3) o << "3 " << F[j] << " " << F[j + 1] << " " << F[j + 2] << "\n"; }
        }
        return nullptr;
    }
}

}  // namespace cci
