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
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

inline void buildContactMap(std::vector<std::vector<bodyID_t>>& map,
                            std::vector<std::vector<float3>>& relative_pos,
                            const DEMSolver& DEMSim,
                            std::shared_ptr<DEMTracker>& particle_tracker,
                            unsigned int num_particles) {
    auto cnt_pairs = DEMSim.GetClumpContacts();
    map.clear();
    relative_pos.clear();
    map.resize(num_particles);
    relative_pos.resize(num_particles);
    std::vector<float3> particle_xyz(num_particles);
    // At system-level, the clump's ID may not start from 0; but a batch of clumps loaded together have consecutive IDs.
    size_t clump_ID_offset = particle_tracker->GetOwnerID();
    for (unsigned int i = 0; i < num_particles; i++) {
        particle_xyz[i] = particle_tracker->Pos(i);
    }
    for (unsigned int i = 0; i < cnt_pairs.size(); i++) {
        const auto& pair = cnt_pairs.at(i);
        map[pair.first - clump_ID_offset].push_back(pair.second);
        map[pair.second - clump_ID_offset].push_back(pair.first);
    }
    for (unsigned int i = 0; i < num_particles; i++) {
        // Main particle location
        float3 main_loc = particle_xyz[i];
        std::vector<float3> init_rel_pos;
        // Compute all this guy's partners' relative positions wrt to itself
        for (const auto& ID : map[i]) {
            init_rel_pos.push_back(DEMSim.GetOwnerPosition(ID) - main_loc);
        }
        relative_pos[i] = init_rel_pos;
    }
}

int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    // DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV | OUTPUT_CONTENT::VEL);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    // You can enforce owner wildcard output by the following call, or directly include OUTPUT_CONTENT::OWNER_WILDCARD
    // in SetOutputContent
    DEMSim.EnableOwnerWildcardOutput();

    path out_dir = current_path();
    out_dir += "/DemoOutput_Indentation";
    create_directory(out_dir);

    // E, nu, CoR, mu, Crr...
    auto mat_type_cube = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}});
    auto mat_type_granular = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}});

    float granular_rad = 0.001;  // 0.002;
    auto template_granular = DEMSim.LoadSphereType(granular_rad * granular_rad * granular_rad * 2.6e3 * 4 / 3 * 3.14,
                                                   granular_rad, mat_type_granular);

    float step_size = 1e-6;
    const double world_size = 0.5;
    const float fill_height = 0.2;
    const float chamber_bottom = -world_size / 2.;
    const float fill_bottom = chamber_bottom + granular_rad;

    DEMSim.InstructBoxDomainDimension(world_size, world_size, world_size);
    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_granular);
    DEMSim.SetCoordSysOrigin("center");

    // Now add a cylinderical boundary
    auto walls = DEMSim.AddExternalObject();
    walls->AddCylinder(make_float3(0), make_float3(0, 0, 1), world_size / 2., mat_type_cube, 0);

    auto cube = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/cube.obj"), mat_type_cube);
    std::cout << "Total num of triangles: " << cube->GetNumTriangles() << std::endl;
    // Make the cube about 10cm by 2cm
    float cube_width = 0.1;
    float cube_height = 0.04;
    double cube_speed = 0.05;  // 0.1 and 0.02, try them too... very similar though
    cube->Scale(make_float3(cube_width, cube_width, cube_height));
    cube->SetFamily(10);
    DEMSim.SetFamilyFixed(10);
    DEMSim.SetFamilyPrescribedLinVel(11, "0", "0", to_string_with_precision(-cube_speed));
    // Track the cube
    auto cube_tracker = DEMSim.Track(cube);

    // Sampler to use
    const float spacing = 2.0005f * granular_rad;
    const float fill_radius = world_size / 2. - 2. * granular_rad;

    PDSampler sampler(spacing);
    std::vector<float3> input_xyz;
    float layer_z = 0;
    while (layer_z < fill_height) {
        float3 sample_center = make_float3(0, 0, fill_bottom + layer_z + spacing / 2);
        auto layer_xyz = sampler.SampleCylinderZ(sample_center, fill_radius, 0);
        input_xyz.insert(input_xyz.end(), layer_xyz.begin(), layer_xyz.end());
        layer_z += spacing;
    }

    // HCPSampler sampler(spacing);
    // float3 fill_center = make_float3(0, 0, fill_bottom + fill_height / 2);
    // auto input_xyz = sampler.SampleCylinderZ(fill_center, fill_radius, fill_height / 2);

    // Calling AddClumps a second time will just add more clumps to the system
    auto particles = DEMSim.AddClumps(template_granular, input_xyz);
    particles->SetFamily(1);
    // Initially, no contact between the brick and the granular material
    DEMSim.DisableContactBetweenFamilies(1, 10);

    // Use a owner wildcard to record tangential displacement compared to initial pos
    DEMSim.ReadContactForceModel("SampleCustomForceModel.cu");
    auto force_model = DEMSim.GetContactForceModel();
    force_model->SetPerOwnerWildcards({"gran_strain", "mu_custom"});
    force_model->SetPerContactWildcards({"delta_tan_x", "delta_tan_y", "delta_tan_z"});
    particles->AddOwnerWildcard("gran_strain", 0.0);
    // Or simply DEMSim.SetOwnerWildcards({"gran_strain"}); it does the job too

    // Low mu at start. This will make the terrain settle into a more densely-packed configuration.
    particles->AddOwnerWildcard("mu_custom", 0.0);

    unsigned int num_particles = input_xyz.size();
    std::cout << "Total num of particles: " << num_particles << std::endl;
    auto particle_tracker = DEMSim.Track(particles);

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.SetCDUpdateFreq(15);
    DEMSim.SetInitBinSize(4 * granular_rad);
    DEMSim.Initialize();

    float sim_end = cube_height * 1.0 / cube_speed;  // 3.0;
    unsigned int fps = 60;                           // 20;
    float frame_time = 1.0 / fps;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    // Keep tab of some sim quantities
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    // Settle
    for (double t = 0; t < 0.9; t += frame_time) {
        // char filename[200], meshname[200];
        // std::cout << "Outputting frame: " << currframe << std::endl;
        // sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        // sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
        // DEMSim.WriteSphereFile(std::string(filename));
        // DEMSim.WriteMeshFile(std::string(meshname));
        DEMSim.ShowThreadCollaborationStats();

        DEMSim.DoDynamicsThenSync(frame_time);
    }
    double init_max_z = max_z_finder->GetValue();
    std::cout << "After settling, max particle Z coord is " << init_max_z << std::endl;

    // Record init positions of the particles
    std::vector<std::vector<bodyID_t>> particle_cnt_map;
    std::vector<std::vector<float3>> particle_init_relative_pos;
    // Build contact map (contact partner owner IDs) for all particles
    buildContactMap(particle_cnt_map, particle_init_relative_pos, DEMSim, particle_tracker, num_particles);

    // Ready to start indentation
    std::cout << "Simulation starts..." << std::endl;
    // Let the brick sink with a downward velocity.
    DEMSim.ChangeFamily(10, 11);
    // Add some friction which is physical... perhaps 0 friction works too,
    DEMSim.SetFamilyOwnerWildcardValue(1, "mu_custom", 0.4);
    double cube_zpos = max_z_finder->GetValue() + cube_height / 2;
    cube_tracker->SetPos(make_float3(0, 0, cube_zpos));
    std::cout << "Initially the cube is at Z = " << cube_zpos << std::endl;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (double t = 0; t < sim_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            // Compute relative displacement
            std::vector<float> gran_strain(num_particles);
            for (unsigned int i = 0; i < num_particles; i++) {
                float3 main_loc = particle_tracker->Pos(i);
                // Compute contact partners' new locations
                std::vector<float3> rel_pos;
                for (auto& ID : particle_cnt_map.at(i)) {
                    rel_pos.push_back(DEMSim.GetOwnerPosition(ID) - main_loc);
                }
                // How large is the strain?
                // float3 strains = make_float3(0);
                float strains = 0.;
                int num_neighbors = particle_init_relative_pos.at(i).size();
                for (int j = 0; j < num_neighbors; j++) {
                    // strains += particle_init_relative_pos.at(i).at(j) - rel_pos.at(j);
                    strains += length(particle_init_relative_pos.at(i).at(j) - rel_pos.at(j));
                }
                gran_strain[i] = (num_neighbors > 0) ? (strains / num_neighbors) : 0.0;
            }
            // Re-build contact map, for the next output step
            buildContactMap(particle_cnt_map, particle_init_relative_pos, DEMSim, particle_tracker, num_particles);
            std::cout << "A new contact map constructed..." << std::endl;

            // Feed displacement info to wildcard, then leverage the output method to output it to the file
            DEMSim.SetFamilyOwnerWildcardValue(1, "gran_strain", gran_strain);
            char filename[200], meshname[200];
            std::cout << "Outputting frame: " << currframe << std::endl;
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshname));
            DEMSim.ShowThreadCollaborationStats();
        }

        DEMSim.DoDynamics(step_size);
        // cube_zpos -= cube_speed * step_size;
        // cube_tracker->SetPos(make_float3(0, 0, cube_zpos));
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count()) / sim_end / (1e-5 / step_size)
              << " seconds (wall time) to finish 1e5 steps' simulation" << std::endl;

    std::cout << "DEMdemo_Indentation exiting..." << std::endl;
    return 0;
}