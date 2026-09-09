// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_FENG_MESH_VALIDATION_HPP
#define DEME_FENG_MESH_VALIDATION_HPP

#include "DEM/BdrsAndObjs.h"
#include <map>

namespace deme {
namespace feng {

// The first force implementation supports small, closed, outward convex solids with one material/mesh patch.
// Verify the geometry rather than trusting SetConvex(). Welding exact coincident vertices supports OBJ seams;
// near-coincident seams conservatively fall back. Plane tests also reject inward winding and concave surfaces.
inline bool validatedSolid(const DEMMesh& mesh) {
    if (mesh.IsShell() || mesh.GetNumPatches() != 1 || mesh.nTri < 4 || mesh.nTri > 256 ||
        mesh.m_face_v_indices.size() != mesh.nTri || mesh.m_vertices.empty() || mesh.m_vertices.size() > 3 * mesh.nTri)
        return false;
    std::vector<double3> vertices;
    std::vector<size_t> canonical;
    for (const auto& v : mesh.m_vertices) {
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z))
            return false;
        size_t j = 0;
        for (; j < vertices.size(); ++j)
            if (v.x == vertices[j].x && v.y == vertices[j].y && v.z == vertices[j].z)
                break;
        if (j == vertices.size())
            vertices.push_back(to_double3(v));
        canonical.push_back(j);
    }
    double scale = 0;
    double3 center = make_double3(0.0, 0.0, 0.0);
    for (const auto& v : vertices) {
        center += v;
        scale = std::max(scale, length(v - vertices[0]));
    }
    if (!(scale > 0))
        return false;
    center = center / double(vertices.size());
    std::map<std::pair<size_t, size_t>, unsigned int> edges;
    double volume6 = 0;
    for (const auto& f : mesh.m_face_v_indices) {
        const int ids[] = {f.x, f.y, f.z};
        size_t c[3];
        for (int k = 0; k < 3; ++k) {
            if (ids[k] < 0 || size_t(ids[k]) >= canonical.size())
                return false;
            c[k] = canonical[ids[k]];
        }
        const auto a = vertices[c[0]], b = vertices[c[1]], d = vertices[c[2]];
        const auto normal = cross(b - a, d - a);
        const double norm = length(normal);
        if (!(norm > 1e-12 * scale * scale))
            return false;
        for (const auto& v : vertices)
            if (dot(normal, v - a) > 1e-7 * scale * norm)
                return false;
        volume6 += dot(a - center, cross(b - center, d - center));
        for (int k = 0; k < 3; ++k)
            if (++edges[{c[k], c[(k + 1) % 3]}] != 1)
                return false;
    }
    for (const auto& edge : edges)
        if (!edges.count({edge.first.second, edge.first.first}))
            return false;
    // Euler characteristic excludes disconnected closed components and unsupported topology.
    return double(vertices.size()) - double(edges.size()) / 2 + mesh.nTri == 2 &&
           volume6 > 1e-12 * scale * scale * scale;
}
}  // namespace feng
}  // namespace deme
#endif
