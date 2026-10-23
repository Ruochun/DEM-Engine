// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause

#include "DEM/API.h"
#include "kernel/DEMHelperKernels.cuh"

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Execute the actual built-in model snippets on the host. A shared translational
// boost must not erase the smaller contact velocity due to rocking. This catches
// float arithmetic before promotion as well as mixed-vector narrowing overloads.
int main(int argc, char* argv[]) {
    const bool gpu = argc == 2 && std::string(argv[1]) == "--gpu";
    if (argc > 1 && !gpu) {
        std::cout << "Usage: DEMTest_HertzianPrecision [--gpu]\n";
        return std::string(argv[1]) == "--help" ? 0 : 1;
    }
    auto evaluate = [](bool frictional, float boost, float mass, bool tangential, bool contact, float restitution) {
        const unsigned int bodyAMatType = 0, bodyBMatType = 0;
        const float E[] = {1e6f}, nu[] = {0.3f};
        const float CoR[1][1] = {{restitution}}, mu[1][1] = {{0.3f}}, Crr[1][1] = {{0.2f}};
        const float AOwnerMass = mass, BOwnerMass = mass, ts = 1e-5f;
        const float3 ALinVel = make_float3(boost, 0, 0), BLinVel = ALinVel;
        const float3 ARotVel = make_float3(0, 0, 1), BRotVel = make_float3(0, 0, 0);
        const float3 locCPA = make_float3(0, -1, 0), locCPB = make_float3(0, 0, 0);
        const float4 AOriQ = make_float4(0, 0, 0, 1), BOriQ = AOriQ;
        const float3 B2A = tangential ? make_float3(0, 0, 1) : make_float3(1, 0, 0);
        const double overlapDepth = contact ? 1e-4 : -1e-4;
        const double overlapArea = deme::PI * 1e-6;
        float3 force = make_float3(0, 0, 0), torque_only_force = make_float3(0, 0, 0);
        float delta_tan_x = 0, delta_tan_y = 0, delta_tan_z = 0, delta_time = 1;
        if (frictional) {
#include "kernel/DEMCustomizablePolicies/FullHertzianForceModel.cu"
        } else {
#include "kernel/DEMCustomizablePolicies/FrictionlessHertzianForceModel.cu"
        }
        return std::array<float, 10>{
            force.x,     force.y,     force.z,     torque_only_force.x, torque_only_force.y, torque_only_force.z,
            delta_tan_x, delta_tan_y, delta_tan_z, delta_time};
    };
    bool passed = true;
    // The unboosted and boosted cases describe identical contact motion; compare
    // force, rolling resistance, and stored friction history after the same step.
    for (bool frictional : {false, true}) {
        for (bool tangential : {false, true}) {
            for (float restitution : {0.0f, 0.2f, 1.0f}) {
                const auto reference = evaluate(frictional, 0, 1, tangential, true, restitution);
                const auto boosted = evaluate(frictional, 1e8f, 1, tangential, true, restitution);
                for (size_t i = 0; i < reference.size(); ++i) {
                    passed = passed && std::isfinite(boosted[i]) &&
                             std::abs(boosted[i] - reference[i]) <= 1e-6f * std::max(1e-6f, std::abs(reference[i]));
                }
            }
        }
        // Both masses are representable in float, but their product is not.
        // Effective mass must be formed in double before evaluating the dashpot.
        const auto massive = evaluate(frictional, 0, 1e25f, false, true, 0.2f);
        passed = passed && std::isfinite(massive[0]) && massive[0] < -1e10f;
        const auto separated = evaluate(frictional, 0, 1, false, false, 0.2f);
        for (size_t i = 0; i < (frictional ? separated.size() : 6); ++i)
            passed = passed && separated[i] == 0;
    }
    // Ensure the tangential/rolling branches were exercised, not merely equal zeros.
    const auto sliding = evaluate(true, 0, 1, true, true, 0.2f);
    passed = passed && sliding[0] < 0 && sliding[2] > 0 && sliding[3] < 0 && sliding[6] != 0;
    std::cout << (passed ? "PASS" : "FAIL")
              << ": Hertzian contact-velocity precision, mass overflow, friction/rolling, and separation\n";
    if (!passed || !gpu)
        return passed ? 0 : 1;

    // Exercise actual runtime compilation and force writeback for primitive and
    // mesh-patch contacts in both models. These short overlaps are deliberate;
    // they test nonzero finite action/reaction forces, not settling behavior.
    using namespace deme;
    for (bool frictional : {false, true}) {
        DEMSolver sim;
        sim.SetVerbosity("INFO");
        if (frictional)
            sim.UseFrictionalHertzianModel();
        else
            sim.UseFrictionlessHertzianModel();
        sim.SetMeshUniversalContact(true);
        sim.SetMeshParticlesLowPoly(true);
        sim.InstructBoxDomainDimension(0.2, 0.2, 0.2);
        sim.SetInitBinSize(0.01);
        sim.SetTimeStepSize(1e-6);
        sim.SetGravitationalAcceleration(make_float3(0, 0, 0));
        sim.SetFamilyFixed(10);
        auto material = sim.LoadMaterial({{"E", 1e6}, {"nu", 0.3}, {"CoR", 0.2}, {"mu", 0.3}, {"Crr", 0.2}});
        auto mesh_type = sim.LoadMeshType(GetDEMEDataFile("mesh/cube.obj"), material, false);
        mesh_type->Scale(0.01);
        mesh_type->SetMass(1);
        mesh_type->SetMOI(make_float3(1e-4, 1e-4, 1e-4));
        mesh_type->SetConvex(true);
        auto base = sim.AddMeshFromTemplate(mesh_type, make_float3(0, 0, 0));
        base->SetFamily(10);
        auto top = sim.AddMeshFromTemplate(mesh_type, make_float3(0, 0, 0.0099));
        auto base_tracker = sim.Track(base), top_tracker = sim.Track(top);
        auto sphere = sim.LoadSphereType(1, 0.001, material);
        auto spheres = sim.AddClumps(sphere, {make_float3(0.04, 0, 0), make_float3(0.0419, 0, 0)});
        auto sphere_tracker = sim.Track(spheres);
        sim.Initialize();
        sim.DoDynamicsThenSync(2e-5);
        auto wrench = [&](bodyID_t owner) {
            std::vector<float3> forces, torques;
            sim.GetOwnerContactWrench(forces, torques, owner);
            return forces[0];
        };
        const auto mesh_a = wrench(base_tracker->GetOwnerID()), mesh_b = wrench(top_tracker->GetOwnerID());
        const auto sphere_a = wrench(sphere_tracker->GetOwnerID()), sphere_b = wrench(sphere_tracker->GetOwnerID(1));
        bool finite = true;
        for (const auto force : {mesh_a, mesh_b, sphere_a, sphere_b})
            finite = finite && std::isfinite(force.x) && std::isfinite(force.y) && std::isfinite(force.z);
        const bool valid = finite && mesh_b.z > 0 && sphere_a.x < 0 &&
                           length(mesh_a + mesh_b) < 1e-5f * length(mesh_b) &&
                           length(sphere_a + sphere_b) < 1e-5f * length(sphere_a);
        passed = passed && valid;
        std::cout << (valid ? "PASS" : "FAIL") << ": " << (frictional ? "frictional" : "frictionless")
                  << " GPU primitive and mesh-patch force paths; mesh Fz=" << mesh_b.z << ", sphere Fx=" << sphere_a.x
                  << '\n';
    }
    return passed ? 0 : 1;
}
