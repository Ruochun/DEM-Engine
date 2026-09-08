// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_FENG_GEOMETRY_CUH
#define DEME_FENG_GEOMETRY_CUH

#include "kernel/CUDAMathHelpers.cuh"
#include <cmath>

namespace deme {
namespace feng {

// Ambiguous features are deliberately excluded: their multiplicity cannot be resolved from one triangle pair.
enum class Intersection { NONE, SEGMENT, AMBIGUOUS };

// Intersect a triangle with a plane, retaining on-plane vertices once. Edge-on-plane cases are handled by the caller.
__host__ __device__ inline int planeCut(const double3* v, const double* d, double tolerance, double3* points) {
    int count = 0;
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        if (fabs(d[i]) <= tolerance) {
            points[count++] = v[i];
        } else if ((d[i] < -tolerance && d[j] > tolerance) || (d[i] > tolerance && d[j] < -tolerance)) {
            points[count++] = v[i] + (v[j] - v[i]) * (d[i] / (d[i] - d[j]));
        }
    }
    return count;
}

// Plane cuts followed by interval overlap give the true zero-thickness surface segment. Coordinates should be
// owner-relative doubles. Tolerances scale with triangle size, rather than the world origin or solver margins.
__host__ __device__ inline Intersection intersectionSegment(const double3* a,
                                                            const double3* b,
                                                            double3& start,
                                                            double3& end) {
    double scale = 0.0;
    for (int i = 0; i < 3; ++i) {
        scale = fmax(scale, length(a[(i + 1) % 3] - a[i]));
        scale = fmax(scale, length(b[(i + 1) % 3] - b[i]));
    }
    double3 na = cross(a[1] - a[0], a[2] - a[0]);
    double3 nb = cross(b[1] - b[0], b[2] - b[0]);
    const double la = length(na), lb = length(nb);
    if (!(scale > 0.0) || !(la > 1e-12 * scale * scale) || !(lb > 1e-12 * scale * scale)) {
        return Intersection::AMBIGUOUS;
    }
    na = na / la;
    nb = nb / lb;
    const double tolerance = 1e-10 * scale;
    double da[3], db[3];
    int za = 0, zb = 0;
    for (int i = 0; i < 3; ++i) {
        da[i] = dot(a[i] - b[0], nb);
        db[i] = dot(b[i] - a[0], na);
        za += fabs(da[i]) <= tolerance;
        zb += fabs(db[i]) <= tolerance;
    }
    if ((fmin(da[0], fmin(da[1], da[2])) > tolerance) || (fmax(da[0], fmax(da[1], da[2])) < -tolerance) ||
        (fmin(db[0], fmin(db[1], db[2])) > tolerance) || (fmax(db[0], fmax(db[1], db[2])) < -tolerance)) {
        return Intersection::NONE;
    }
    double3 tangent = cross(na, nb);
    const double sine = length(tangent);
    if (sine < 1e-8 || za >= 2 || zb >= 2) {
        return Intersection::AMBIGUOUS;
    }
    tangent = tangent / sine;
    double3 pa[3], pb[3];
    if (planeCut(a, da, tolerance, pa) != 2 || planeCut(b, db, tolerance, pb) != 2) {
        return Intersection::AMBIGUOUS;
    }
    // Relative interval coordinates avoid cancellation when the triangle pair is far from its owner center.
    const double3 reference = pa[0];
    double a0 = 0.0, a1 = dot(pa[1] - reference, tangent);
    double b0 = dot(pb[0] - reference, tangent), b1 = dot(pb[1] - reference, tangent);
    const double lower = fmax(fmin(a0, a1), fmin(b0, b1));
    const double upper = fmin(fmax(a0, a1), fmax(b0, b1));
    if (upper < lower - tolerance) {
        return Intersection::NONE;
    }
    if (!(upper - lower > tolerance)) {
        return Intersection::AMBIGUOUS;
    }
    start = reference + tangent * lower;
    end = reference + tangent * upper;
    return Intersection::SEGMENT;
}

// Stokes' theorem for G = integral r cross dS gives G = -1/2 integral |r|^2 dr on a CLOSED boundary.
// Integrating the quadratic along a straight segment yields the expression below. The proposal's -1/3 prefactor
// would give only 2/3 of this moment (and the wrong planar centroid). S and G must use the SAME local origin.
__host__ __device__ inline void contribution(const double3& start, const double3& end, double3& s, double3& g) {
    const double3 dx = end - start;
    s = cross(start, end) * 0.5;
    g = dx * (-0.5 * (dot(start, end) + dot(dx, dx) / 3.0));
}

// Choose the point on the boundary-derived line nearest a supplied local anchor. Unlike world lambda=0, this
// choice translates with the bodies and keeps the tangential-force lever arm near the legacy contact location.
// The orientation nA cross nB produces S toward B for an outward-wound penetrating A; DEME uses B2A = -S/|S|.
__host__ __device__ inline bool contactLine(const double3& s,
                                            const double3& g,
                                            const double3& anchor,
                                            double& area,
                                            double3& b2a,
                                            double3& point) {
    area = length(s);
    if (!(area > 0.0) || !std::isfinite(area)) {
        return false;
    }
    const double3 n = s / area;
    const double3 base = cross(n, g) / area;
    b2a = n * -1.0;
    point = base + n * dot(anchor - base, n);
    return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z);
}

}  // namespace feng
}  // namespace deme
#endif
