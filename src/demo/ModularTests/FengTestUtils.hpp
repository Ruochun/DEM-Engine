// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_FENG_TEST_UTILS_HPP
#define DEME_FENG_TEST_UTILS_HPP
#include "DEM/utils/FengGeometry.cuh"
#include "DEM/BdrsAndObjs.h"
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

// Weld vertices for host validation and the solver's mesh topology tables.
inline deme::DEMMesh cubeMesh(int refinement = 1, bool flip = false) {
    deme::DEMMesh mesh;
    for (const auto& triangle : cube(refinement, flip)) {
        int ids[3];
        for (int k = 0; k < 3; ++k) {
            const float3 vertex = to_float3(triangle[k]);
            size_t index = 0;
            for (; index < mesh.m_vertices.size(); ++index)
                if (length(mesh.m_vertices[index] - vertex) < 1e-7f)
                    break;
            if (index == mesh.m_vertices.size())
                mesh.m_vertices.push_back(vertex);
            ids[k] = static_cast<int>(index);
        }
        mesh.m_face_v_indices.push_back(make_int3(ids[0], ids[1], ids[2]));
    }
    mesh.nTri = mesh.m_face_v_indices.size();
    mesh.SetMass(1.f);
    mesh.SetMOI(make_float3(1.f));
    return mesh;
}
}  // namespace feng_test
#endif
