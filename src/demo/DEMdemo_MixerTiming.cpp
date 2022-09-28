//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <cmath>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

int main() {
    float granular_rad = 0.005;
    unsigned int num_particles = 0;

    while (num_particles < 3e8) {
        DEMSolver DEMSim;
        DEMSim.SetVerbosity(ERROR);
        DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
        DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
        DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);

        // E, nu, CoR, mu, Crr...
        auto mat_type_mixer = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.2}, {"mu", 0.5}, {"Crr", 0.0}});
        auto mat_type_granular =
            DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.2}, {"mu", 0.5}, {"Crr", 0.0}});

        float step_size = 1e-5;
        const double world_size = 1;
        const float chamber_height = world_size / 3.;
        const float fill_height = chamber_height;
        const float chamber_bottom = -world_size / 2.;
        const float fill_bottom = chamber_bottom + chamber_height;

        DEMSim.InstructBoxDomainDimension(world_size, world_size, world_size);
        DEMSim.InstructBoxDomainBoundingBC("all", mat_type_granular);
        DEMSim.SetCoordSysOrigin("center");

        // Now add a cylinderical boundary
        auto walls = DEMSim.AddExternalObject();
        walls->AddCylinder(make_float3(0), make_float3(0, 0, 1), world_size / 2., mat_type_mixer, 0);

        auto mixer =
            DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/internal_mixer.obj").string(), mat_type_mixer);
        std::cout << "Total num of triangles: " << mixer->GetNumTriangles() << std::endl;
        mixer->Scale(make_float3(world_size / 2, world_size / 2, chamber_height));
        mixer->SetFamily(10);
        // Define the prescribed motion of mixer
        DEMSim.SetFamilyPrescribedAngVel(10, "0", "0", "2 * 3.14159");

        auto template_granular = DEMSim.LoadSphereType(
            granular_rad * granular_rad * granular_rad * 2.8e3 * 4 / 3 * 3.14, granular_rad, mat_type_granular);

        // Track the mixer
        auto mixer_tracker = DEMSim.Track(mixer);

        // Sampler to use
        HCPSampler sampler(2.1f * granular_rad);
        float3 fill_center = make_float3(0, 0, fill_bottom + fill_height / 2);
        const float fill_radius = world_size / 2. - 2. * granular_rad;
        auto input_xyz = sampler.SampleCylinderZ(fill_center, fill_radius, fill_height / 2);
        DEMSim.AddClumps(template_granular, input_xyz);
        num_particles = input_xyz.size();
        std::cout << "Particle size: " << granular_rad << std::endl;
        std::cout << "Total num of particles: " << num_particles << std::endl;

        DEMSim.SetInitTimeStep(step_size);
        DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
        // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
        DEMSim.SetCDUpdateFreq(10);
        // DEMSim.SetExpandFactor(1e-3);
        DEMSim.SetMaxVelocity(10.);
        DEMSim.SetExpandSafetyParam(1.0);
        DEMSim.SetInitBinSize(4 * granular_rad);
        DEMSim.Initialize();

        path out_dir = current_path();
        out_dir += "/DemoOutput_Mixer";
        create_directory(out_dir);

        float sim_end = 3.0;
        unsigned int fps = 20;

        mixer_tracker->SetPos(make_float3(0, 0, chamber_bottom + chamber_height / 2.0));
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

        DEMSim.DoDynamicsThenSync(sim_end);
        DEMSim.ShowThreadCollaborationStats();
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        std::cout << (time_sec.count()) / sim_end << " seconds (wall time) to finish 1 second's simulation"
                  << std::endl;

        granular_rad *= std::pow(0.5, 1. / 3.);
    }

    return 0;
}
