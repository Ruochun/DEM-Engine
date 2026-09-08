// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_FENG_TEST_UTILS_HPP
#define DEME_FENG_TEST_UTILS_HPP
#include "DEM/utils/FengGeometry.cuh"
#include <array>
#include <vector>

namespace feng_test {
using Triangle = std::array<double3, 3>;

// Outward-wound unit cube, with a selectable diagonal and 2*n*n triangles on each physical face.
inline std::vector<Triangle> cube(int n, bool flip, double3 offset = make_double3(0.0, 0.0, 0.0)) {
    const double3 origins[] = {{0, 0, 0}, {1, 0, 0}, {0, 0, 0}, {0, 1, 0}, {0, 0, 0}, {0, 0, 1}};
    const double3 u[] = {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}};
    const double3 v[] = {{0, 1, 0}, {0, 0, 1}, {0, 0, 1}, {1, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    std::vector<Triangle> result;
    for (int face = 0; face < 6; ++face) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                const double3 p = origins[face] + offset + u[face] * (double(i) / n) + v[face] * (double(j) / n);
                const double3 q = p + u[face] / double(n), r = q + v[face] / double(n), s = p + v[face] / double(n);
                if (flip) {
                    result.push_back({p, q, s});
                    result.push_back({q, r, s});
                } else {
                    result.push_back({p, q, r});
                    result.push_back({p, r, s});
                }
            }
        }
    }
    return result;
}
}  // namespace feng_test
#endif
