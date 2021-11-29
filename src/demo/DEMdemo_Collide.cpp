//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//  All rights reserved.

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <core/utils/chpf/particle_writer.hpp>
#include <granular/ApiSystem.h>

#include <cstdio>
#include <time.h>

using namespace sgps;

int main() {
    DEMSolver DEM_sim(1.f);

    srand(time(NULL));

    // total number of random clump templates to generate
    int num_template = 2;

    int min_sphere = 1;
    int max_sphere = 1;

    float min_rad = 0.08;
    float max_rad = 0.2;

    float min_relpos = -0.1;
    float max_relpos = 0.1;

    auto mat_type_1 = DEM_sim.LoadMaterialType(1, 10);

    for (int i = 0; i < num_template; i++) {
        // first decide the number of spheres that live in this clump
        int num_sphere = rand() % (max_sphere - min_sphere + 1) + 1;

        // then allocate the clump template definition arrays
        float mass = (float)rand() / RAND_MAX;
        float3 MOI = make_float3((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX);
        std::vector<float> radii;
        std::vector<float3> relPos;
        std::vector<unsigned int> mat;

        // randomly generate clump template configurations

        // the relPos of a sphere is always seeded from one of the already-generated sphere
        float3 seed_pos = make_float3(0);
        for (int j = 0; j < num_sphere; j++) {
            radii.push_back(((float)rand() / RAND_MAX) * (max_rad - min_rad) + min_rad);
            /*
            if (i == 0)
                radii.push_back(0.16138418018817902);
            else
                radii.push_back(0.16321393847465515);
            */
            float3 tmp;
            if (j == 0) {
                tmp.x = 0;
                tmp.y = 0;
                tmp.z = 0;
            } else {
                tmp.x = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.y = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.z = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
            }
            tmp += seed_pos;
            relPos.push_back(tmp);
            mat.push_back(mat_type_1);

            // seed relPos from one of the previously generated spheres
            int choose_from = rand() % (j + 1);
            seed_pos = relPos.at(choose_from);
        }

        // it returns the numbering of this clump template (although here we don't care)
        auto template_num = DEM_sim.LoadClumpType(mass, MOI, radii, relPos, mat);
    }

    std::vector<unsigned int> input_template_num;
    std::vector<float3> input_xyz;
    std::vector<float3> input_vel;

    // show one for each template configuration
    for (int i = 0; i < num_template; i++) {
        input_template_num.push_back(i);
        input_xyz.push_back(make_float3(i * 0.3, 0, 0));
        float sgn_vel = (i % 2 == 0) ? 1.0 : -1.0;
        input_vel.push_back(make_float3(sgn_vel * 50.0, 0, 0));
    }
    DEM_sim.SetClumps(input_template_num, input_xyz);
    DEM_sim.SetClumpVels(input_vel);

    DEM_sim.InstructBoxDomainNumVoxel(22, 21, 21, 1e-10);

    DEM_sim.CenterCoordSys();
    DEM_sim.SetTimeStepSize(1e-4);
    DEM_sim.SetGravitationalAcceleration(make_float3(0, 0, 0));

    DEM_sim.Initialize();

    for (int i = 0; i < 10; i++) {
        std::cout << "Iteration: " << i + 1 << std::endl;
        DEM_sim.LaunchThreads();

        char filename[100];
        sprintf(filename, "./DEMdemo_collide_output_%04d.csv", i);
        DEM_sim.WriteFileAsSpheres(std::string(filename));
    }

    std::cout << "DEMdemo_Collide exiting..." << std::endl;
    return 0;
}