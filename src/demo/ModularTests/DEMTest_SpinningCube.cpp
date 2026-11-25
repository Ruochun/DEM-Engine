// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "DEM/API.h"

using namespace deme;
namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

// Compare two free cubes in the same scene: only Cspin differs. Zero mu/Crr isolates the spin couple.
// Zero intrinsic curvature deliberately exercises spin-only contacts with an unresolved RR bound.
void run(bool mesh_plane, bool rr_error = false) {
    DEMSolver solver(1);
    solver.SetVerbosity("ERROR");
    solver.InstructBoxDomainDimension(8, 8, 8);
    solver.SetMeshUniversalContact(true);
    solver.SetTimeStepSize(2e-5);
    solver.SetCDUpdateFreq(1);
    solver.SetInitBinNumTarget(16);
    solver.SetGravitationalAcceleration(make_float3(0, 0, -9.81f));
    auto wall =
        solver.LoadMaterial({{"E", 1e5f}, {"nu", .25f}, {"CoR", .3f}, {"mu", 0.f}, {"Crr", 0.f}, {"Cspin", 0.f}});
    if (mesh_plane) {
        DEMMesh plane;
        // Both footprints stay inside one triangle, avoiding edges and triangulation seams.
        plane.m_vertices = {make_float3(-3, -2, 0), make_float3(3, -2, 0), make_float3(0, 3, 0)};
        plane.m_face_v_indices = {make_int3(0, 1, 2)};
        plane.nTri = 1;
        plane.SetPatchIDs({0});
        plane.SetMass(1000.f);
        plane.SetMOI(make_float3(1000.f));
        plane.SetMaterial(wall);
        plane.SetFamily(10);
        plane.SetContactRadius(0.f);
        solver.AddMesh(plane);
        solver.SetFamilyFixed(10);
    } else {
        solver.AddBCPlane(make_float3(0), make_float3(0, 0, 1), wall);
    }
    std::vector<std::shared_ptr<DEMTracker>> cubes;
    for (int i = 0; i < 2; ++i) {
        auto mat =
            solver.LoadMaterial({{"E", 1e5f}, {"nu", .25f}, {"CoR", .3f}, {"mu", 0.f}, {"Crr", 0.f}, {"Cspin", 0.f}});
        solver.SetMaterialPropertyPair("Cspin", mat, wall, i == 0 ? .1f : 0.f);
        solver.SetMaterialPropertyPair("Crr", mat, wall, rr_error ? .1f : 0.f);
        auto cube = solver.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat);
        cube->Scale(.2f);
        cube->SetEachTriangleAsPatch();
        cube->SetMass(1.f);
        cube->SetMOI(make_float3(.2f * .2f / 6.f));
        cube->SetContactRadius(0.f);
        // Loaded meshes default to the fixed wall family; these two bodies must evolve freely.
        cube->SetFamily(i + 1);
        cube->SetInitPos(make_float3(i == 0 ? -.5f : .5f, 0, .1f));
        cubes.push_back(solver.Track(cube));
    }
    solver.Initialize();
    solver.DoDynamicsThenSync(.15);
    for (auto& cube : cubes)
        cube->SetAngVel(make_float3(0, 0, 1));
    // Observe a sustained low-speed interval, not a single zero crossing. The law is regularized near rest.
    int quiet_samples = 0;
    for (int sample = 0; sample < 60; ++sample) {
        solver.DoDynamicsThenSync(.01);
        for (int i = 0; i < 2; ++i) {
            const auto w = cubes[i]->AngVelGlobal();
            const auto p = cubes[i]->Pos();
            const auto v = cubes[i]->Vel();
            require(std::isfinite(length(w)) && std::isfinite(length(p)) && std::isfinite(length(v)),
                    "non-finite cube state");
            require(std::abs(p.x - (i == 0 ? -.5f : .5f)) < .02f && std::abs(p.y) < .02f &&
                        std::abs(p.z - .1f) < .02f && length(v) < .2f,
                    "spin test developed excessive translational motion");
        }
        const float spin = length(cubes[0]->AngVelGlobal());
        quiet_samples = spin < .005f ? quiet_samples + 1 : 0;
        if (sample % 10 == 9)
            std::cout << "mesh_plane=" << mesh_plane << " t=" << (sample + 1) * .01 << " spin=" << spin
                      << " control=" << length(cubes[1]->AngVelGlobal()) << std::endl;
    }
    require(quiet_samples >= 10, "spinning friction did not stop and hold the cube below tolerance");
    require(cubes[1]->AngVelGlobal().z > .8f, "zero-Cspin control lost substantial spin");
    std::cout << "PASS: spinning cube, mesh_plane=" << mesh_plane << std::endl;
}
}  // namespace
int main(int argc, char** argv) {
    try {
        const std::string option = argc > 1 ? argv[1] : "";
        if (option == "--rr-error") {
            run(true, true);  // Deliberate GPU trap: run after ordinary tests on WSL.
            return 2;
        }
        run(false);
        run(true);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << std::endl;
        return 1;
    }
}
