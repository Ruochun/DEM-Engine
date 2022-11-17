//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

using namespace deme;

const double math_PI = 3.14159;

int main() {
    std::filesystem::path out_dir = std::filesystem::current_path();
    out_dir += "/DEMdemo_WheelDP_Slip_mu0.3";
    std::filesystem::create_directory(out_dir);

    // `World'
    float G_mag = 9.81;
    float step_size = 5e-7;
    double world_size_y = 0.52;
    double world_size_x = 1.53;
    double world_size_z = 4.0;

    // Define the wheel geometry
    float wheel_rad = 0.22;  // 1.5 cm grousers
    float wheel_width = 0.25;
    float wheel_mass = 8.7;
    float total_pressure = 480.0;
    float added_pressure = (total_pressure - wheel_mass * G_mag);
    float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);

    float DPs[] = {200, 150, 100, 50, 30, 20, 10};
    unsigned int run_mode = 0;
    unsigned int currframe = 0;

    for (float DP : DPs) {
        DEMSolver DEMSim;
        DEMSim.SetVerbosity(INFO);
        DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
        DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
        DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
        DEMSim.SetContactOutputContent(OWNER | FORCE | POINT);

        // E, nu, CoR, mu, Crr...
        // 0.8 is the default
        auto mat_type_wheel = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.3}, {"Crr", 0.00}});
        auto mat_type_terrain =
            DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.3}, {"Crr", 0.00}});

        DEMSim.InstructBoxDomainDimension(world_size_x, world_size_y, world_size_z);
        DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);
        DEMSim.SetCoordSysOrigin("center");
        float bottom = -0.5;
        auto bot_wall = DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);
        auto bot_wall_tracker = DEMSim.Track(bot_wall);

        auto wheel = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/rover_wheels/curiosity_wheel_surface.obj"),
                                                   mat_type_wheel);
        wheel->SetMass(wheel_mass);
        wheel->SetMOI(make_float3(wheel_IXX, wheel_IYY, wheel_IXX));
        // Give the wheel a family number so we can potentially add prescription
        wheel->SetFamily(1);
        // Track it
        auto wheel_tracker = DEMSim.Track(wheel);

        // Define the GRC terrain particle templates
        DEMClumpTemplate shape_template;
        shape_template.ReadComponentFromFile(GetDEMEDataFile("clumps/triangular_flat.csv"));
        // Calculate its mass and MOI
        float terrain_density = 2.6e3;
        double clump_vol = 5.5886717;
        float mass = terrain_density * clump_vol;  // in kg or g
        float3 MOI = make_float3(1.8327927, 2.1580013, 0.77010059) * (double)2.6e3;
        // Scale the template we just created
        std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates;
        std::vector<double> scales = {0.00063, 0.00033, 0.00022, 0.00015, 0.00009};
        std::for_each(scales.begin(), scales.end(), [](double& r) { r *= 20.; });
        for (double scaling : scales) {
            auto this_template = shape_template;
            this_template.mass = (double)mass * scaling * scaling * scaling;
            this_template.MOI.x = (double)MOI.x * (double)(scaling * scaling * scaling * scaling * scaling);
            this_template.MOI.y = (double)MOI.y * (double)(scaling * scaling * scaling * scaling * scaling);
            this_template.MOI.z = (double)MOI.z * (double)(scaling * scaling * scaling * scaling * scaling);
            std::cout << "Mass: " << this_template.mass << std::endl;
            std::cout << "MOIX: " << this_template.MOI.x << std::endl;
            std::cout << "MOIY: " << this_template.MOI.y << std::endl;
            std::cout << "MOIZ: " << this_template.MOI.z << std::endl;
            std::cout << "=====================" << std::endl;
            std::for_each(this_template.radii.begin(), this_template.radii.end(),
                          [scaling](float& r) { r *= scaling; });
            std::for_each(this_template.relPos.begin(), this_template.relPos.end(),
                          [scaling](float3& r) { r *= scaling; });
            this_template.materials = std::vector<std::shared_ptr<DEMMaterial>>(this_template.nComp, mat_type_terrain);
            this_template.SetVolume(clump_vol * scaling * scaling * scaling);
            ground_particle_templates.push_back(DEMSim.LoadClumpType(this_template));
        }

        // Now we load clump locations from a checkpointed file
        {
            std::cout << "Making terrain..." << std::endl;
            auto clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC_20e6.csv");
            auto clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC_20e6.csv");
            std::vector<float3> in_xyz;
            std::vector<float4> in_quat;
            std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
            unsigned int t_num = 0;
            for (int i = 0; i < scales.size(); i++) {
                char t_name[20];
                sprintf(t_name, "%04d", t_num);

                auto this_type_xyz = clump_xyz[std::string(t_name)];
                auto this_type_quat = clump_quaternion[std::string(t_name)];

                size_t n_clump_this_type = this_type_xyz.size();
                std::cout << "Loading clump " << std::string(t_name) << " which has particle num: " << n_clump_this_type
                          << std::endl;
                // Prepare clump type identification vector for loading into the system (don't forget type 0 in
                // ground_particle_templates is the template for rover wheel)
                std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                         ground_particle_templates.at(t_num));

                // Add them to the big long vector
                in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
                in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
                in_types.insert(in_types.end(), this_type.begin(), this_type.end());
                std::cout << "Added clump type " << t_num << std::endl;
                // Our template names are 0000, 0001 etc.
                t_num++;
            }

            // Now, we don't need all particles loaded...
            std::vector<notStupidBool_t> elem_to_remove(in_xyz.size(), 0);
            for (size_t i = 0; i < in_xyz.size(); i++) {
                if (std::abs(in_xyz.at(i).y) > (world_size_y - 0.03) / 2 ||
                    std::abs(in_xyz.at(i).x) > (world_size_x - 0.06) / 2)
                    elem_to_remove.at(i) = 1;
            }
            in_xyz.erase(std::remove_if(in_xyz.begin(), in_xyz.end(),
                                        [&elem_to_remove, &in_xyz](const float3& i) {
                                            return elem_to_remove.at(&i - in_xyz.data());
                                        }),
                         in_xyz.end());
            in_quat.erase(std::remove_if(in_quat.begin(), in_quat.end(),
                                         [&elem_to_remove, &in_quat](const float4& i) {
                                             return elem_to_remove.at(&i - in_quat.data());
                                         }),
                          in_quat.end());
            in_types.erase(std::remove_if(in_types.begin(), in_types.end(),
                                          [&elem_to_remove, &in_types](const auto& i) {
                                              return elem_to_remove.at(&i - in_types.data());
                                          }),
                           in_types.end());
            DEMClumpBatch base_batch(in_xyz.size());
            base_batch.SetTypes(in_types);
            base_batch.SetPos(in_xyz);
            base_batch.SetOriQ(in_quat);

            /*
            std::vector<float> x_shift_dist = {-0.5, 0.5};
            std::vector<float> y_shift_dist = {0};
            // Add some patches of such graular bed
            for (float x_shift : x_shift_dist) {
                for (float y_shift : y_shift_dist) {
                    DEMClumpBatch another_batch = base_batch;
                    std::vector<float3> my_xyz = in_xyz;
                    std::for_each(my_xyz.begin(), my_xyz.end(), [x_shift, y_shift](float3& xyz) {
                        xyz.x += x_shift;
                        xyz.y += y_shift;
                    });
                    another_batch.SetPos(my_xyz);
                    DEMSim.AddClumps(another_batch);
                }
            }
            */

            DEMSim.AddClumps(base_batch);
        }

        // Now add a plane to compress the sample
        // auto compressor = DEMSim.AddExternalObject();
        // compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, -1), mat_type_terrain);
        // compressor->SetFamily(2);
        // auto compressor_tracker = DEMSim.Track(compressor);

        // Families' prescribed motions
        float w_r = math_PI / 12.;
        float v_ref = w_r * wheel_rad;

        double sim_end = 3.;
        // Note: this wheel is not `dictated' by our prescrption of motion because it can still fall onto the ground
        // (move freely linearly)
        DEMSim.SetFamilyPrescribedAngVel(1, "0", to_string_with_precision(w_r), "0", false);
        DEMSim.AddFamilyPrescribedAcc(1, "none", "none", to_string_with_precision(-added_pressure / wheel_mass));

        // `Real sim' family number
        DEMSim.SetFamilyPrescribedAngVel(2, "0", to_string_with_precision(w_r), "0", false);
        DEMSim.AddFamilyPrescribedAcc(2, to_string_with_precision(-DP / wheel_mass), "none",
                                      to_string_with_precision(-added_pressure / wheel_mass));

        // Some inspectors
        auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
        auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
        auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
        auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

        DEMSim.SetInitTimeStep(step_size);
        DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -G_mag));
        DEMSim.SetCDUpdateFreq(20);
        DEMSim.SetMaxVelocity(30.);
        DEMSim.SetExpandSafetyParam(1.1);
        DEMSim.SetInitBinSize(2 * scales.at(2));
        DEMSim.Initialize();

        // Compress until dense enough
        unsigned int curr_step = 0;
        unsigned int fps = 10;
        unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
        double frame_time = 1.0 / fps;
        unsigned int report_ps = 1000;
        unsigned int report_steps = (unsigned int)(1.0 / (report_ps * step_size));
        std::cout << "Output at " << fps << " FPS" << std::endl;

        // Put the wheel in place, then let the wheel sink in initially
        float max_z = max_z_finder->GetValue();
        float init_x = -0.4;
        // if (DP > 140) {
        //     init_x = 0.;
        // }
        wheel_tracker->SetPos(make_float3(init_x, 0, max_z + 0.04 + wheel_rad));
        for (double t = 0; t < 0.6; t += frame_time) {
            char filename[200], meshname[200];
            std::cout << "Outputting frame: " << currframe << std::endl;
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshname));

            // float3 pos = bot_wall_tracker->Pos();
            // float3 force = bot_wall_tracker->ContactAcc();
            // std::cout << "Bottom wall pos: " << pos.x << ", " << pos.y << ", " << pos.z << std::endl;
            // std::cout << "Bottom wall force: " << force.x << ", " << force.y << ", " << force.z << std::endl;
            DEMSim.DoDynamicsThenSync(frame_time);
        }

        // Switch wheel from free fall into DP test
        DEMSim.ChangeFamily(1, 2);

        // DEMSim.SetInitTimeStep(step_size * 2);
        // DEMSim.UpdateSimParams();

        // float x_at_start;
        bool start_measure = false;
        float time_at_start = 0.5;
        for (double t = 0; t < sim_end; t += step_size, curr_step++) {
            if (curr_step % out_steps == 0) {
                char filename[200], meshname[200];
                std::cout << "Outputting frame: " << currframe << std::endl;
                sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
                sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
                DEMSim.WriteSphereFile(std::string(filename));
                DEMSim.WriteMeshFile(std::string(meshname));
                DEMSim.ShowThreadCollaborationStats();
            }

            if (t >= time_at_start && !start_measure) {
                start_measure = true;
                // Get X location now
                // x_at_start = wheel_tracker->Pos().x;
            }

            if (curr_step % report_steps == 0 && start_measure) {
                // float dist_moved = wheel_tracker->Pos().x - x_at_start;
                // float slip = 1.0 - dist_moved / (wheel_rad * w_r * (t - time_at_start));
                float Vx = wheel_tracker->Vel().x;
                float slip = 1.0 - Vx / (w_r * wheel_rad);
                std::cout << "Current run mode: " << run_mode << std::endl;
                std::cout << "Time: " << t << std::endl;
                // std::cout << "Distance: " << dist_moved << std::endl;
                std::cout << "Vx: " << Vx << std::endl;
                std::cout << "Slip: " << slip << std::endl;
                std::cout << "Max system velocity: " << max_v_finder->GetValue() << std::endl;
            }

            DEMSim.DoDynamics(step_size);
        }

        run_mode++;
        DEMSim.ShowTimingStats();
    }

    std::cout << "DEMdemo_WheelDP_Slip exiting..." << std::endl;
    return 0;
}