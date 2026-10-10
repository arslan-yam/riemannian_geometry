#pragma once
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cci {

// OFF reader; xy is translated so the chart minimum is at the origin (float-safe coordinates)
inline void read_off(const std::string& path, std::vector<double>& V, std::vector<unsigned>& F, double* origin = nullptr) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string tok;
    in >> tok;
    if (tok != "OFF") throw std::runtime_error("not an OFF file: " + path);
    size_t nv, nf, ne;
    in >> nv >> nf >> ne;
    V.resize(3 * nv);
    for (size_t i = 0; i < 3 * nv; ++i) in >> V[i];
    F.resize(3 * nf);
    for (size_t i = 0; i < nf; ++i) {
        int k;
        in >> k >> F[3 * i] >> F[3 * i + 1] >> F[3 * i + 2];
        if (k != 3) throw std::runtime_error("non-triangle face");
    }
    double mx = 1e300, my = 1e300;
    for (size_t i = 0; i < nv; ++i) { mx = std::min(mx, V[3 * i]); my = std::min(my, V[3 * i + 1]); }
    for (size_t i = 0; i < nv; ++i) { V[3 * i] -= mx; V[3 * i + 1] -= my; }
    if (origin) { origin[0] = mx; origin[1] = my; }
}

}  // namespace cci
