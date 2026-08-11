// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "DEM/API.h"

using namespace deme;

namespace {

constexpr float kBodyMass = 2.0f;
constexpr float kCohesion = 50.0f;

// Run one fixed, overlapping sphere pair and return the force reported for the first sphere. Fixing both bodies keeps
// the Hertzian contribution identical between the zero-cohesion and cohesive runs, isolating the added force term.
float3 contactForce(bool frictionless, float cohesion) {
    DEMSolver solver(1);
    solver.SetVerbosity("ERROR");
    solver.InstructBoxDomainDimension(4.0f, 4.0f, 4.0f);
    solver.SetGravitationalAcceleration(make_float3(0.0f));
    solver.SetInitTimeStep(1.0e-5f);
    if (frictionless) {
        solver.UseFrictionlessHertzianModel();
    }

    std::unordered_map<std::string, float> properties{
        {"E", 1.0e5f}, {"nu", 0.3f}, {"CoR", 0.5f}, {"mu", 0.0f}, {"Crr", 0.0f}, {"Cohesion", cohesion},
    };
    auto material = solver.LoadMaterial(properties);
    auto sphere = solver.LoadSphereType(kBodyMass, 0.5f, material);

    auto body_a = solver.AddClumps(sphere, std::vector<float3>{make_float3(-0.45f, 0.0f, 0.0f)});
    auto body_b = solver.AddClumps(sphere, std::vector<float3>{make_float3(0.45f, 0.0f, 0.0f)});
    body_a->SetFamily(1);
    body_b->SetFamily(1);
    solver.SetFamilyFixed(1);
    auto tracker_a = solver.Track(body_a);

    solver.Initialize();
    solver.DoStepDynamics();

    std::vector<float3> points;
    std::vector<float3> forces;
    tracker_a->GetContactForcesForAll(points, forces);
    if (forces.size() != 1) {
        std::cerr << "FAIL: expected one sphere-sphere contact, found " << forces.size() << std::endl;
        return make_float3(NAN);
    }
    return forces.front();
}

bool verifyModel(bool frictionless) {
    const float3 baseline = contactForce(frictionless, 0.0f);
    const float3 cohesive = contactForce(frictionless, kCohesion);
    const float3 added = cohesive - baseline;

    // Equal masses give m_eff = m / 2, so this setup's expected cohesion magnitude is exactly kCohesion.
    const float expected = kCohesion * (kBodyMass * kBodyMass) / (kBodyMass + kBodyMass);
    const float actual = length(added);
    const float tolerance = 1.0e-2f * expected;
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::cerr << "FAIL: " << (frictionless ? "frictionless" : "frictional") << " Hertzian cohesion magnitude was "
                  << actual << ", expected " << expected << std::endl;
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!verifyModel(false) || !verifyModel(true)) {
        return 1;
    }
    std::cout << "PASS: default frictional and frictionless Hertzian models apply pairwise cohesion." << std::endl;
    return 0;
}
