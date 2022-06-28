//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//  All rights reserved.

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <core/utils/chpf/particle_writer.hpp>
#include <DEM/ApiSystem.h>
#include <DEM/HostSideHelpers.hpp>

#include <cstdio>
#include <time.h>
#include <filesystem>

using namespace sgps;
using namespace std::filesystem;
const double PI = 3.141592653589793;

void EllpsiodFallingOver(DEMSolver& DEM_sim) {
    // An ellipsoid a,b,c = 0.2,0.2,0.5, represented several sphere components
    std::vector<float> radii = {0.095, 0.136, 0.179, 0.204, 0.204, 0.179, 0.136, 0.095};
    std::vector<float3> relPos = {make_float3(0, 0, 0.4),    make_float3(0, 0, 0.342),  make_float3(0, 0, 0.228),
                                  make_float3(0, 0, 0.071),  make_float3(0, 0, -0.071), make_float3(0, 0, -0.228),
                                  make_float3(0, 0, -0.342), make_float3(0, 0, -0.4)};
    // Then calculate mass and MOI
    float mass = 5.0;
    // E, nu, CoR, mu, Crr
    auto mat_type_1 = DEM_sim.LoadMaterialType(1e8, 0.3, 0.5, 0.25, 0.2);
    float3 MOI = make_float3(1. / 5. * mass * (0.2 * 0.2 + 0.5 * 0.5), 1. / 5. * mass * (0.2 * 0.2 + 0.5 * 0.5),
                             1. / 5. * mass * (0.2 * 0.2 + 0.2 * 0.2));
    auto ellipsoid_template = DEM_sim.LoadClumpType(mass, MOI, radii, relPos, mat_type_1);

    // Add the ground
    float3 normal_dir = make_float3(0, 0, 1);
    float3 tang_dir = make_float3(0, 1, 0);
    DEM_sim.AddBCPlane(make_float3(0, 0, 0), normal_dir, mat_type_1);

    // Add an ellipsoid with init vel
    auto ellipsoid = DEM_sim.AddClumps(ellipsoid_template, normal_dir * 0.5);
    ellipsoid->SetVel(tang_dir * 0.3);
    auto ellipsoid_tracker = DEM_sim.Track(ellipsoid);

    DEM_sim.SetTimeStepSize(1e-3);
    DEM_sim.Initialize();

    float frame_time = 1e-1;
    path out_dir = current_path();
    out_dir += "/DEMdemo_TestPack";
    create_directory(out_dir);
    for (int i = 0; i < 6.0 / frame_time; i++) {
        char filename[100];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), i);
        DEM_sim.WriteClumpFile(std::string(filename));
        std::cout << "Frame: " << i << std::endl;
        float4 oriQ = ellipsoid_tracker->OriQ();
        float3 angVel = ellipsoid_tracker->AngVel();
        std::cout << "Time: " << frame_time * i << std::endl;
        std::cout << "Quaternion of the ellipsoid: " << oriQ.x << ", " << oriQ.y << ", " << oriQ.z << ", " << oriQ.w
                  << std::endl;
        std::cout << "Angular velocity of the ellipsoid: " << angVel.x << ", " << angVel.y << ", " << angVel.z
                  << std::endl;

        DEM_sim.DoDynamics(frame_time);
    }
}

void SphereRollUpIncline(DEMSolver& DEM_sim) {
    // Clump initial profile
    std::vector<float3> input_xyz;
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_template_type;
    std::vector<unsigned int> family_code;
    std::vector<float3> input_vel;

    auto mat_type_1 = DEM_sim.LoadMaterialType(1e7, 0.3, 0.5, 0.5, 0.3);
    // A ball
    float sphere_rad = 0.2;
    float mass = 5.0;
    auto sphere_template = DEM_sim.LoadClumpSimpleSphere(mass, sphere_rad, mat_type_1);

    // Incline angle
    float alpha = 35.;
    // Add the incline
    float3 normal_dir = make_float3(-std::sin(2. * PI * (alpha / 360.)), 0., std::cos(2. * PI * (alpha / 360.)));
    float3 tang_dir = make_float3(std::cos(2. * PI * (alpha / 360.)), 0., std::sin(2. * PI * (alpha / 360.)));
    DEM_sim.AddBCPlane(make_float3(0, 0, 0), normal_dir, mat_type_1);

    // Add a ball rolling
    auto sphere = DEM_sim.AddClumps(sphere_template, normal_dir * sphere_rad);
    sphere->SetVel(tang_dir * 0.5);
    auto sphere_tracker = DEM_sim.Track(sphere);

    float step_time = 1e-5;
    DEM_sim.SetTimeStepSize(step_time);
    DEM_sim.Initialize();

    path out_dir = current_path();
    out_dir += "/DEMdemo_TestPack";
    create_directory(out_dir);
    for (int i = 0; i < 0.15 / step_time; i++) {
        char filename[100];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), i);
        // DEM_sim.WriteClumpFile(std::string(filename));
        std::cout << "Frame: " << i << std::endl;
        float3 vel = sphere_tracker->Vel();
        float3 angVel = sphere_tracker->AngVel();
        std::cout << "Time: " << step_time * i << std::endl;
        std::cout << "Velocity of the sphere: " << vel.x << ", " << vel.y << ", " << vel.z << std::endl;
        std::cout << "Angular velocity of the sphere: " << angVel.x << ", " << angVel.y << ", " << angVel.z
                  << std::endl;

        DEM_sim.DoStepDynamics();
    }
}

int main() {
    DEMSolver DEM_sim;
    DEM_sim.SetVerbosity(DEBUG);
    DEM_sim.SetOutputFormat(DEM_OUTPUT_FORMAT::CSV);

    DEM_sim.InstructBoxDomainNumVoxel(22, 22, 20, 7.5e-11);
    DEM_sim.CenterCoordSys();
    DEM_sim.SetGravitationalAcceleration(make_float3(0, 0, -9.8));
    DEM_sim.SetCDUpdateFreq(0);

    // Validation tests
    // SphereRollUpIncline(DEM_sim);
    EllpsiodFallingOver(DEM_sim);

    std::cout << "DEMdemo_TestPack exiting..." << std::endl;
    // TODO: add end-game report APIs

    return 0;
}