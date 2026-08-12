// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include <iostream>

#include "DEM/API.h"

using namespace deme;

namespace {

bool near(float actual, float expected, float tolerance = 1.0e-5f) {
    return std::abs(actual - expected) <= tolerance;
}

}  // namespace

// Verify that visualization snapshots expose current transformed geometry and honor category filtering.
int main() {
    DEMSolver solver(1);
    solver.SetVerbosity("ERROR");
    solver.InstructBoxDomainDimension(8.0f, 8.0f, 8.0f);

    auto material = solver.LoadMaterial({{"E", 1.0e7}, {"nu", 0.3}, {"CoR", 0.5}, {"mu", 0.4}, {"Crr", 0.0}});
    auto sphere = solver.LoadSphereType(1.0f, 0.25f, material);
    solver.AddClumps(sphere, std::vector<float3>{make_float3(1.0f, 2.0f, 3.0f)});

    auto mesh = solver.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/plane_20by20.obj").string(), material);
    mesh->Scale(0.1f);
    mesh->SetFamily(1);
    solver.SetFamilyFixed(1);
    solver.Initialize();

    const auto sphere_snapshot = solver.GetVisualizationSnapshot(true, false);
    if (sphere_snapshot.spheres.size() != 1 || !sphere_snapshot.triangles.empty()) {
        std::cerr << "FAIL: sphere-only snapshot returned the wrong geometry categories." << std::endl;
        return 1;
    }
    const auto& rendered_sphere = sphere_snapshot.spheres.front();
    if (!near(rendered_sphere.position.x, 1.0f) || !near(rendered_sphere.position.y, 2.0f) ||
        !near(rendered_sphere.position.z, 3.0f) || !near(rendered_sphere.radius, 0.25f)) {
        std::cerr << "FAIL: visualization sphere does not match its initialized world transform." << std::endl;
        return 1;
    }

    const auto triangle_snapshot = solver.GetVisualizationSnapshot(false, true);
    if (!triangle_snapshot.spheres.empty() || triangle_snapshot.triangles.size() != mesh->GetNumTriangles()) {
        std::cerr << "FAIL: triangle-only snapshot returned the wrong geometry categories." << std::endl;
        return 1;
    }

    std::cout << "PASS: visualization snapshots expose and filter current solver geometry." << std::endl;
    return 0;
}
