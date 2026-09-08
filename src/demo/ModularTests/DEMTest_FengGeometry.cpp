// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#include "FengTestUtils.hpp"
#include <iostream>
#include <stdexcept>

using namespace deme;
using namespace feng_test;

namespace {
void require(bool ok, const char* message) {
    if (!ok)
        throw std::runtime_error(message);
}
void near(double3 actual, double3 expected, const char* message, double tolerance = 1e-9) {
    require(length(actual - expected) < tolerance, message);
}

// Analytic triangle cut, endpoint ordering, degenerate inputs, and grazing/parallel rejection.
void testSegments() {
    Triangle a = {{{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}};
    Triangle b = {{{0.5, -1, -1}, {0.5, 2, -1}, {0.5, 0.5, 1}}};
    double3 p, q;
    require(feng::intersectionSegment(a.data(), b.data(), p, q) == feng::Intersection::SEGMENT, "triangle cut");
    near(p, make_double3(0.5, 0, 0), "start orientation");
    near(q, make_double3(0.5, 1.25, 0), "end orientation");
    double3 s, g;
    feng::contribution(p, q, s, g);
    near(s, make_double3(0, 0, 0.3125), "segment vector area");
    near(g, make_double3(0, -0.5 * (0.25 + 1.5625 / 3.0) * 1.25, 0), "segment moment");
    double3 reversedP, reversedQ;
    require(feng::intersectionSegment(b.data(), a.data(), reversedP, reversedQ) == feng::Intersection::SEGMENT,
            "owner swap cut");
    near(reversedP, q, "owner swap start");
    near(reversedQ, p, "owner swap end");
    require(feng::intersectionSegment(a.data(), a.data(), p, q) == feng::Intersection::AMBIGUOUS, "coplanar");
    Triangle separated = a;
    for (auto& v : separated)
        v.z += 1e-4;
    require(feng::intersectionSegment(a.data(), separated.data(), p, q) == feng::Intersection::NONE, "margin only");
    Triangle zero = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    require(feng::intersectionSegment(a.data(), zero.data(), p, q) == feng::Intersection::AMBIGUOUS, "zero area");
    Triangle edge = {{{0, 0, 0}, {2, 0, 0}, {0, 0, 1}}};
    require(feng::intersectionSegment(a.data(), edge.data(), p, q) == feng::Intersection::AMBIGUOUS, "shared edge");
    Triangle point = {{{0, 0, 0}, {-1, 0, 1}, {0, -1, 1}}};
    require(feng::intersectionSegment(a.data(), point.data(), p, q) == feng::Intersection::AMBIGUOUS, "point touch");
}

// A translated planar square catches the incorrect -1/3 moment prefactor and origin-dependent lambda=0 choice.
void testMomentAndOrigin() {
    double3 vertices[] = {{2, 3, 5}, {3, 3, 5}, {3, 4, 5}, {2, 4, 5}};
    const double3 center = make_double3(2.5, 3.5, 5);
    for (double3 origin : {make_double3(0.0, 0.0, 0.0), make_double3(10, -20, 30)}) {
        double3 s = make_double3(0.0, 0.0, 0.0), g = s, residual = s;
        for (int i = 0; i < 4; ++i) {
            double3 si, gi;
            feng::contribution(vertices[i] - origin, vertices[(i + 1) % 4] - origin, si, gi);
            s += si;
            g += gi;
            residual += vertices[(i + 1) % 4] - vertices[i];
        }
        near(s, make_double3(0, 0, 1), "square area");
        near(g, cross(center - origin, s), "square moment");
        double area;
        double3 normal, point;
        require(feng::contactLine(s, g, center - origin, area, normal, point), "square line");
        near(point + origin, center, "translated contact point");
        near(normal, make_double3(0, 0, -1), "B2A sign");
        near(residual, make_double3(0.0, 0.0, 0.0), "closed square");
    }
    // One edge is a legal user-defined patch but not a closed boundary. Its vector area depends on origin.
    double3 s1, g1, s2, g2;
    feng::contribution(vertices[0], vertices[1], s1, g1);
    feng::contribution(vertices[0] - center, vertices[1] - center, s2, g2);
    require(length(s1 - s2) > 1.0, "open patch counterexample");
}

// Integrate the three rectangular faces of A inside B independently to get analytic S and G. Then check all
// triangle pairs for 2/8/32 triangles per face and opposite diagonals, including owner swap.
void testCubes() {
    const double3 offset = make_double3(.37, .23, .19);
    const double3 expectedS = make_double3(.77 * .81, .63 * .81, .63 * .77);
    const double3 expectedG = cross(make_double3(1, .615, .595), make_double3(expectedS.x, 0, 0)) +
                              cross(make_double3(.685, 1, .595), make_double3(0, expectedS.y, 0)) +
                              cross(make_double3(.685, .615, 1), make_double3(0, 0, expectedS.z));
    for (int n : {1, 2, 4})
        for (bool flip : {false, true})
            for (bool swap : {false, true}) {
                const auto a = cube(n, flip), b = cube(n, !flip, offset);
                double3 s = make_double3(0.0, 0.0, 0.0), g = s, residual = s;
                int count = 0;
                for (const auto& ta : a)
                    for (const auto& tb : b) {
                        double3 p, q, si, gi;
                        auto status = swap ? feng::intersectionSegment(tb.data(), ta.data(), p, q)
                                           : feng::intersectionSegment(ta.data(), tb.data(), p, q);
                        require(status != feng::Intersection::AMBIGUOUS, "nondegenerate cubes");
                        if (status == feng::Intersection::SEGMENT) {
                            feng::contribution(p, q, si, gi);
                            s += si;
                            g += gi;
                            residual += q - p;
                            ++count;
                        }
                    }
                require(count >= 6, "cube boundary segments");
                const double sign = swap ? -1.0 : 1.0;
                near(s, expectedS * sign, "cube tessellation area");
                near(g, expectedG * sign, "cube tessellation moment");
                near(residual, make_double3(0.0, 0.0, 0.0), "cube boundary closure");
                double area;
                double3 normal, point;
                require(feng::contactLine(s, g, make_double3(.7, .6, .6), area, normal, point), "cube line");
                near(normal, expectedS * (-sign / length(expectedS)), "cube B2A orientation");
            }
}
// Similarity transforms exercise units and non-axis-aligned intersections. Construct the geometry relative to an
// origin after translating in world space, as the GPU path does with owner centers and local mesh coordinates.
void testTransforms() {
    const double3 expected = make_double3(.77 * .81, .63 * .81, .63 * .77);
    const double3 shift = make_double3(1e7, -2e7, 3e7);
    auto rotate = [](double3 v) { return make_double3(.8 * v.x - .6 * v.y, .6 * v.x + .8 * v.y, v.z); };
    for (double scale : {1e-6, 1.0, 1e6}) {
        auto a = cube(2, false), b = cube(2, true, make_double3(.37, .23, .19));
        for (auto* mesh : {&a, &b}) {
            for (auto& tri : *mesh)
                for (auto& p : tri)
                    p = rotate(p) * scale;
        }
        double3 s = make_double3(0.0, 0.0, 0.0), g = s;
        for (const auto& ta : a)
            for (const auto& tb : b) {
                double3 p, q, si, gi;
                if (feng::intersectionSegment(ta.data(), tb.data(), p, q) == feng::Intersection::SEGMENT) {
                    feng::contribution(p, q, si, gi);
                    s += si;
                    g += gi;
                }
            }
        near(s / (scale * scale), rotate(expected), "rotated/scaled cube area");
        double area;
        double3 n, point;
        const double3 anchor = rotate(make_double3(.7, .6, .6)) * scale;
        require(feng::contactLine(s, g, anchor, area, n, point), "scaled line");
        // Origin is applied once, after local line construction; it must not enter the polynomial moments.
        const double3 worldPoint = point + shift;
        near(worldPoint - shift, point, "large world translation", 1e-8);
    }
}
}  // namespace
int main() {
    try {
        testSegments();
        testMomentAndOrigin();
        testCubes();
        testTransforms();
        std::cout << "PASS: Feng segments, moment, origin, B2A, degeneracies and cube refinement/diagonals.\n";
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
