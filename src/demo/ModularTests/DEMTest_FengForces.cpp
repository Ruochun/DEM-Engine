// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#include "DEM/API.h"
#include "FengTestUtils.hpp"
#include <iostream>
#include <stdexcept>

using namespace deme;
namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

// Evaluate one initially stationary contact: CoR=1 and mu=Crr=0 isolate the elastic Hertz term and its lever arm.
// Keep a sphere pair active too, checking that the geometry selector only changes mesh-mesh records.
std::vector<float3> run(bool fengMode,
                        bool aligned = false,
                        bool split = false,
                        float overlap = .005f,
                        bool sliding = false,
                        bool edge = false) {
    DEMSolver sim;
    sim.SetVerbosity("ERROR");
    sim.InstructBoxDomainDimension(10.f, 10.f, 10.f);
    sim.SetGravitationalAcceleration(make_float3(0.f));
    sim.SetMeshUniversalContact(true);
    sim.SetMeshParticlesLowPoly(true);
    sim.SetSimplePatchCombination(true);
    sim.SetMeshMeshContactGeometry(fengMode ? "feng" : "default");
    // Force mode must produce snapshots even without the separate diagnostic switch.
    sim.SetMeshMeshFengDiagnostics(!fengMode);
    // Include opposing faces in the FIRST candidate list of this initially overlapping configuration.
    sim.SetExpandFactor(.02f);
    sim.SetCDUpdateFreq(1);
    sim.DisableAdaptiveUpdateFreq();
    sim.SetTimeStepSize(1e-7f);
    auto material =
        sim.LoadMaterial({{"E", 1e5f}, {"nu", .3f}, {"CoR", 1.f}, {"mu", sliding ? .3f : 0.f}, {"Crr", 0.f}});
    auto sphere = sim.LoadSphereType(1.f, .1f, material);
    sim.AddClumps(sphere, std::vector<float3>{make_float3(-2.f, 0.f, 0.f), make_float3(-1.85f, 0.f, 0.f)});
    auto a = feng_test::cubeMesh(), b = feng_test::cubeMesh();
    // A small tilted cube penetrating one face of a larger slab avoids the opposed-face voting of two deep cubes.
    // This is a local particle/wall contact, matching the hopper's intended force-validation regime.
    for (auto& v : a.m_vertices)
        v = (v - make_float3(.5f)) * .1f;
    for (auto& v : b.m_vertices)
        v = make_float3(2 * v.x - 1, 2 * v.y - 1, .1f * v.z - .1f);
    const float angle = aligned ? 0.f : .31f;
    const auto orientation = make_float4(0.f, std::sin(angle / 2), 0.f, std::cos(angle / 2));
    a.SetInitQuat(orientation);
    float bottom = 0;
    for (auto v : a.m_vertices) {
        applyOriQToVector3(v, orientation);
        bottom = std::min(bottom, v.z);
    }
    a.SetInitPos(make_float3(edge ? .975f : (aligned ? 0.f : .013f), aligned ? 0.f : .007f, -bottom - overlap));
    // Meshes default to the fixed family; sliding requires a dynamic owner across both evaluations.
    a.SetFamily(0);
    a.SetConvex(true);
    b.SetConvex(true);
    a.SetMaterial(material);
    b.SetMaterial(material);
    if (split) {
        a.SetEachTriangleAsPatch();
        b.SetEachTriangleAsPatch();
    }
    auto ah = sim.AddMesh(a), bh = sim.AddMesh(b);
    sim.Initialize();
    if (sliding)
        sim.SetOwnerVelocity(2, std::vector<float3>{make_float3(0.f, 1.f, 0.f)});
    sim.DoDynamicsThenSync(1e-7);
    const auto records = sim.GetMeshMeshFengDiagnostics();
    require(!records.empty(), "mesh force snapshot missing");
    std::vector<float3> forces, torques;
    sim.GetOwnerContactWrench(forces, torques, 0, 4);
    require(length(forces[2]) > 0 && length(forces[3]) > 0, "positive mesh Hertz forces required");
    require(length(forces[2] + forces[3]) < 1e-4 * length(forces[2]), "equal/opposite mesh forces");
    if (!split && !aligned) {
        require(records.size() == 1, "one owner-pair contact");
        const auto& r = records[0];
        std::cout << "mode=" << fengMode << " depth=" << r.legacyPenetration << " area=" << r.legacyArea
                  << " Feng=" << r.fengArea << " eligible=" << r.fengEligible
                  << " fallback=" << static_cast<unsigned int>(r.fallbackReason)
                  << " pointDifference=" << length(r.fengContactPoint - r.legacyContactPoint) << '\n';
        require(r.fengEligible && r.boundaryValidated, "active closed mesh contact must be Feng eligible");
        require(r.usedFeng == fengMode, "force selector did not use requested geometry");
        if (!edge)
            require(length(cross(r.fengContactPoint - r.legacyContactPoint, r.fengNormal)) > 1e-4,
                    "fixture must distinguish the selected torque lever arm");
        const double area = fengMode ? r.fengArea : r.legacyArea;
        const double3 normal = fengMode ? r.fengNormal : to_double3(r.legacyNormal);
        const double3 point = fengMode ? r.fengContactPoint : r.legacyContactPoint;
        const double magnitude = (4.0 / 3.0) * (1e5 / (2 * (1 - .3 * .3))) * std::sqrt(area / PI) * r.legacyPenetration;
        auto expected = normal * magnitude;
        double tangentMagnitude = 0;
        if (sliding) {
            const double3 relative = r.ownerA == 2 ? make_double3(0, 1, 0) : make_double3(0, -1, 0);
            const auto tangent = relative - normal * dot(relative, normal);
            const double effectiveG = 1e5 / (4 * (2 - .3) * (1 + .3));
            const auto tangentForce = tangent * (-8 * effectiveG * std::sqrt(area / PI) * 1e-7);
            tangentMagnitude = length(tangentForce);
            const auto observed = to_double3(forces[r.ownerA]) - normal * dot(to_double3(forces[r.ownerA]), normal);
            require(length(observed - tangentForce) < .05 * tangentMagnitude, "first-step Hertz tangential force");
            expected += tangentForce;
        }
        require(length(to_double3(forces[r.ownerA]) - expected) < 2e-4 * magnitude, "Hertz force formula mismatch");
        require(length(to_double3(torques[r.ownerA]) - cross(point - r.referenceOrigin, expected)) < 1e-5 * magnitude,
                "Feng contact point not used for torque");
        if (sliding) {
            // A second force evaluation must accumulate the same patch's elastic tangential history.
            sim.DoDynamicsThenSync(1e-7);
            const auto next = sim.GetMeshMeshFengDiagnostics().at(0);
            sim.GetOwnerContactWrench(forces, torques, 0, 4);
            const auto force = to_double3(forces[next.ownerA]);
            const double tangential = length(force - next.fengNormal * dot(force, next.fengNormal));
            require(
                next.usedFeng == fengMode && tangential > 1.5 * tangentMagnitude && tangential < 2.5 * tangentMagnitude,
                "tangential history did not accumulate across Feng steps");
        }
        if (fengMode) {
            sim.DoDynamicsThenSync(0.0);
            // An explicit mesh edit invalidates its certificate, even when the coordinates happen to be unchanged.
            sim.SetTriNodeRelPos(ah->owner, 0, ah->m_vertices);
            sim.DoDynamicsThenSync(1e-7);
            for (const auto& next : sim.GetMeshMeshFengDiagnostics())
                require(!next.usedFeng && !next.ownersValidated, "deformation must trigger whole-patch fallback");
        }
    } else {
        for (const auto& r : records)
            require(!r.usedFeng, "ambiguous/partial boundary must use legacy geometry");
    }
    return forces;
}
}  // namespace

int main() {
    try {
        const auto baseline = run(false), feng = run(true);
        require(length(baseline[0] - feng[0]) == 0.f && length(baseline[1] - feng[1]) == 0.f, "sphere forces changed");
        // A planar wall can give the same area in both flavors; formula and selection checks apply separately.
        run(true, false, false, .008f);
        run(false, false, false, .005f, true);
        run(true, false, false, .005f, true);
        // Crossing a slab edge also checks force selection on a nonplanar boundary.
        run(false, false, false, .005f, false, true);
        run(true, false, false, .005f, false, true);

        for (bool split : {false, true}) {
            const auto legacy = run(false, true, split), fallback = run(true, true, split);
            for (size_t i = 0; i < legacy.size(); ++i)
                require(length(legacy[i] - fallback[i]) == 0.f, "fallback changed Hertz force");
        }
        std::cout
            << "PASS: Hertz forces, selected torque, tangential history, unchanged spheres and whole-patch fallback.\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
