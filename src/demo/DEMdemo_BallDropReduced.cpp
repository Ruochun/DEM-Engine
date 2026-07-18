//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//  SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Reduced meshed-ball drop test for debugging sphere angular velocity growth.
// A single deterministic layer of spherical particles rests on a ground plane,
// then the same meshed ball used by DEMdemo_BallDrop is released just above it.
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <vector>

using namespace deme;
using namespace std::filesystem;

int main() {
    std::cout << "==== DEME demo/test: DEMdemo_BallDropReduced ====" << std::endl;
    std::cout << "========================================" << std::endl;
    DEMSolver DEMSim;
    DEMSim.SetVerbosity("METRIC");
    DEMSim.SetOutputFormat("CSV");
    DEMSim.SetOutputContent({"ABSV", "ANG_VEL", "FAMILY"});
    DEMSim.SetMeshOutputFormat("VTK");
    DEMSim.SetMeshParticlesLowPoly(true);
    // Current staged branch uses the simple patch-combination route by default.
    // DEMSim.SetSimplePatchCombination(true);
    // DEMSim.UseFrictionlessHertzianModel();

    path out_dir = current_path();
    out_dir /= "DemoOutput_BallDropReduced";
    create_directory(out_dir);

    constexpr double ball_rad = 0.0254 / 2.;
    constexpr float ball_density = 6.2e3f;
    constexpr double particle_rad = 0.0025 / 2.;
    constexpr double layer_halfwidth = 0.045;
    constexpr double drop_clearance = 0.006;
    constexpr float step_size = 2e-6f;
    constexpr float settle_time = 0.05f;
    constexpr float sim_time = 0.35f;
    constexpr unsigned int fps = 200;
    constexpr unsigned int projectile_family = 2;

    auto mat_type_ball = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});
    auto mat_type_particles = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});
    auto mat_type_ground = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});

    DEMSim.InstructBoxDomainDimension({-0.075, 0.075}, {-0.075, 0.075}, {0, 0.18});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_ground);
    DEMSim.AddBCPlane(make_float3(0.f, 0.f, 0.f), make_float3(0.f, 0.f, 1.f), mat_type_ground);

    auto projectile = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/sphere.obj").string(), mat_type_ball);
    projectile->Scale(ball_rad);
    projectile->SetInitPos(make_float3(0.f, 0.f, 0.14f));
    const float ball_mass = ball_density * 4.f / 3.f * PI * ball_rad * ball_rad * ball_rad;
    projectile->SetMass(ball_mass);
    projectile->SetMOI(make_float3(ball_mass * 2.f / 5.f * ball_rad * ball_rad,
                                   ball_mass * 2.f / 5.f * ball_rad * ball_rad,
                                   ball_mass * 2.f / 5.f * ball_rad * ball_rad));
    projectile->SetFamily(projectile_family);
    DEMSim.SetFamilyFixed(projectile_family);
    auto projectile_tracker = DEMSim.Track(projectile);

    const float particle_mass = 2.5e3f * 4.f / 3.f * PI * particle_rad * particle_rad * particle_rad;
    auto particle_template = DEMSim.LoadSphereType(particle_mass, particle_rad, mat_type_particles);

    std::vector<float3> particle_xyz;
    const double spacing_x = 2.02 * particle_rad;
    const double spacing_y = std::sqrt(3.0) * particle_rad * 1.01;
    unsigned int row = 0;
    for (double y = -layer_halfwidth; y <= layer_halfwidth; y += spacing_y, row++) {
        const double row_offset = (row % 2 == 0) ? 0.0 : particle_rad * 1.01;
        for (double x = -layer_halfwidth + row_offset; x <= layer_halfwidth; x += spacing_x) {
            particle_xyz.push_back(
                make_float3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(particle_rad)));
        }
    }
    auto particles = DEMSim.AddClumps(particle_template, particle_xyz);
    particles->SetVel(make_float3(0.f, 0.f, 0.f));
    particles->SetAngVel(make_float3(0.f, 0.f, 0.f));

    std::cout << "Meshed ball triangles: " << projectile->GetNumTriangles() << std::endl;
    std::cout << "Particle layer count: " << particle_xyz.size() << std::endl;
    std::cout << "Particle radius: " << particle_rad << std::endl;
    std::cout << "Ball radius: " << ball_rad << std::endl;
    std::cout << "Drop clearance above particle tops: " << drop_clearance << std::endl;

    auto max_angvel_finder = DEMSim.CreateInspector("max_absangvel");
    auto max_absv_finder = DEMSim.CreateInspector("max_absv");

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetMaxVelocity(30.);
    DEMSim.SetGravitationalAcceleration(make_float3(0.f, 0.f, -9.81f));
    DEMSim.Initialize();

    const float frame_time = 1.f / static_cast<float>(fps);
    unsigned int currframe = 0;
    auto write_frame = [&](const char* stage) {
        char sphere_filename[100], mesh_filename[100];
        std::snprintf(sphere_filename, sizeof(sphere_filename), "DEMdemo_output_%04u.csv", currframe);
        std::snprintf(mesh_filename, sizeof(mesh_filename), "DEMdemo_mesh_%04u.vtk", currframe);
        DEMSim.WriteSphereFile(out_dir / sphere_filename);
        DEMSim.WriteMeshFile(out_dir / mesh_filename);
        std::cout << stage << " frame " << currframe << ", max |omega| = " << max_angvel_finder->GetValue()
                  << ", max |v| = " << max_absv_finder->GetValue() << std::endl;
        currframe++;
    };

    for (float t = 0.f; t < settle_time; t += frame_time) {
        write_frame("settle");
        DEMSim.DoDynamicsThenSync(frame_time);
    }

    const float release_z = static_cast<float>(2.0 * particle_rad + ball_rad + drop_clearance);
    DEMSim.ChangeFamily(projectile_family, 0);
    projectile_tracker->SetPos(make_float3(0.f, 0.f, release_z));
    projectile_tracker->SetVel(make_float3(0.f, 0.f, 0.f));
    projectile_tracker->SetAngVel(make_float3(0.f, 0.f, 0.f));

    std::cout << "Released mesh ball at z = " << release_z << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (float t = 0.f; t < sim_time; t += frame_time) {
        write_frame("drop");
        DEMSim.DoDynamics(frame_time);
        DEMSim.ShowThreadCollaborationStats();
    }
    DEMSim.DoDynamicsThenSync(0.f);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << elapsed.count() << " seconds (wall time) to finish reduced ball drop" << std::endl;
    DEMSim.ShowTimingStats();
    std::cout << "DEMdemo_BallDropReduced exiting..." << std::endl;
    return 0;
}
