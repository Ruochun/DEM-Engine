// Copyright (c) 2021, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause

// C++ regression counterpart of deme3_tests/05_geometry4_100p_3s.py. Two inward-facing wall meshes deliberately remain
// unsplit so the test can detect and diagnose abnormal particle acceleration caused by merged concave-patch contacts.

#include <core/ApiVersion.h>
#include <DEM/API.h>
#include <DEM/utils/Samplers.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace deme;
using namespace std::filesystem;

namespace {

constexpr float ABNORMAL_VELOCITY_MM_S = 1000.0f;
constexpr float SPHERE_RADIUS_MM = 1.0f;

int fail(const std::string& message) {
    std::cerr << "DEMTest_InwardBoxCornerContact failed: " << message << std::endl;
    return 1;
}

bool finite(const float3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

// Match the Python case's evenly distributed down-selection from the full HCP sample.
std::vector<float3> distributePositions(const std::vector<float3>& candidates, size_t limit) {
    if (candidates.size() <= limit) {
        return candidates;
    }
    std::vector<float3> selected;
    selected.reserve(limit);
    for (size_t index = 0; index < limit; index++) {
        const double location =
            static_cast<double>(index) * static_cast<double>(candidates.size() - 1) / static_cast<double>(limit - 1);
        selected.push_back(candidates[static_cast<size_t>(std::llround(location))]);
    }
    return selected;
}

size_t writeAbnormalContacts(std::ofstream& output,
                             double time,
                             bodyID_t owner,
                             float speed,
                             const float3& position,
                             const float3& velocity,
                             const std::shared_ptr<ContactInfoContainer>& contacts) {
    size_t related_contacts = 0;
    for (size_t contact = 0; contact < contacts->Size(); contact++) {
        if (contacts->GetAOwner()[contact] != owner && contacts->GetBOwner()[contact] != owner) {
            continue;
        }
        related_contacts++;
        const float3 point = contacts->GetPoint()[contact];
        const float3 normal = contacts->GetNormal()[contact];
        const float3 force = contacts->GetForce()[contact];
        const float3 particle_to_contact = point - position;
        const float contact_distance = length(particle_to_contact);
        const float penetration = contacts->Get<float>("diagnostic_penetration")[contact];
        const float area = contacts->Get<float>("diagnostic_area")[contact];
        output << time << ',' << owner << ',' << speed << ',' << position.x << ',' << position.y << ',' << position.z
               << ',' << velocity.x << ',' << velocity.y << ',' << velocity.z << ','
               << contacts->GetContactType()[contact] << ',' << contacts->GetAOwner()[contact] << ','
               << contacts->GetBOwner()[contact] << ',' << contacts->GetAGeo()[contact] << ','
               << contacts->GetBGeo()[contact] << ',' << point.x << ',' << point.y << ',' << point.z << ','
               << particle_to_contact.x << ',' << particle_to_contact.y << ',' << particle_to_contact.z << ','
               << contact_distance << ',' << contact_distance - SPHERE_RADIUS_MM << ',' << penetration << ',' << area
               << ',' << normal.x << ',' << normal.y << ',' << normal.z << ',' << force.x << ',' << force.y << ','
               << force.z << '\n';
    }
    return related_contacts;
}

}  // namespace

int main() {
    constexpr float density_kg_mm3 = 2.6e-9f;
    constexpr float youngs_modulus_kg_mm_s2 = 100.0f;
    constexpr float gravity_mm_s2 = 9810.0f;
    constexpr float spacing_mm = 2.2f;
    constexpr size_t particle_limit = 100;
    constexpr double time_step = 1.0e-6;
    constexpr double diagnostic_start_time = 2.0;
    constexpr double prediagnostic_interval = 0.01;
    // Sample every few integration steps so a short-lived contact that produces an abnormal impulse is not skipped.
    constexpr unsigned int diagnostic_steps = 5;
    constexpr double diagnostic_interval = diagnostic_steps * time_step;
    constexpr double output_interval = 0.20;
    constexpr double simulation_end = 3.0;

    const auto candidates = DEMBoxHCPSampler(make_float3(0.0f, 31.1f, 180.0f), make_float3(40.0f, 1.1f, 20.0f),
                                             spacing_mm * SPHERE_RADIUS_MM);
    const std::vector<float3> positions = distributePositions(candidates, particle_limit);
    if (positions.empty()) {
        return fail("HCP sampling produced no particle positions");
    }

    DEMSolver solver;
    solver.SetVerbosity("INFO");
    solver.SetOutputFormat(OUTPUT_FORMAT::CSV);
    solver.SetOutputContent({"ABSV", "XYZ", "VEL", "FAMILY"});
    solver.SetMeshOutputFormat(MESH_FORMAT::VTK);
    solver.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    solver.SetContactOutputContent({"OWNER", "GEO_ID", "FORCE", "POINT", "NORMAL", "CNT_WILDCARD"});
    solver.InstructBoxDomainDimension({-93.0f, 73.0f}, {-126.5f, 96.5f}, {146.0f, 249.25f});

    auto particle_material = solver.LoadMaterial(
        {{"E", youngs_modulus_kg_mm_s2}, {"nu", 0.30f}, {"CoR", 0.80f}, {"mu", 0.20f}, {"Crr", 0.0f}});
    auto wall_material = solver.LoadMaterial(
        {{"E", youngs_modulus_kg_mm_s2}, {"nu", 0.30f}, {"CoR", 0.50f}, {"mu", 0.50f}, {"Crr", 0.0f}});
    solver.SetMaterialPropertyPair("mu", particle_material, particle_material, 0.20f);
    solver.SetMaterialPropertyPair("mu", particle_material, wall_material, 0.35f);
    solver.SetMaterialPropertyPair("CoR", particle_material, particle_material, 0.80f);
    solver.SetMaterialPropertyPair("CoR", particle_material, wall_material, 0.50f);
    solver.SetMaterialPropertyPair("Crr", particle_material, particle_material, 0.0f);
    solver.SetMaterialPropertyPair("Crr", particle_material, wall_material, 0.0f);

    // Retain the stock Hertzian model's patch-level penetration and area for abnormal-contact diagnosis. Per-contact
    // wildcards make the values available through both GetContactDetailedInfo and WriteContactFile.
    auto diagnostic_force_model = solver.DefineContactForceModel(HERTZIAN_FORCE_MODEL() + R"(
        if (overlapDepth > 0) {
            diagnostic_penetration = (float)overlapDepth;
            diagnostic_area = (float)overlapArea;
        } else {
            diagnostic_penetration = 0.f;
            diagnostic_area = 0.f;
        }
    )");
    diagnostic_force_model->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr"});
    diagnostic_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr"});
    diagnostic_force_model->SetPerContactWildcards(
        {"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "diagnostic_penetration", "diagnostic_area"});

    struct MeshSpec {
        const char* filename;
        float3 translation;
    };
    const MeshSpec mesh_specs[] = {
        {"WALL_OTR_WALL_OTR-1_global.obj", make_float3(-10.0f, -15.0649415f, 193.5f)},
        {"WALL_INPUT_WALL_INPUT-1_global.obj", make_float3(-2.2601318f, 16.5789261f, 200.625f)},
    };
    for (size_t mesh_index = 0; mesh_index < 2; mesh_index++) {
        const path mesh_path = GET_DATA_PATH() / "mesh/geometry4" / mesh_specs[mesh_index].filename;
        auto mesh = solver.AddWavefrontMeshObject(mesh_path.string(), wall_material);
        mesh->SetInitPos(mesh_specs[mesh_index].translation);
        mesh->SetFamily(static_cast<unsigned int>(10 + mesh_index));
        // Do not call SplitIntoConvexPatches: this regression targets the default all-triangles-in-patch-0 behavior.
        if (mesh->GetNumPatches() != 1 || mesh->ArePatchesExplicitlySet()) {
            return fail("a geometry4 wall did not retain its default single implicit patch");
        }
        solver.SetFamilyFixed(static_cast<unsigned int>(10 + mesh_index));
    }

    const float sphere_mass = (4.0f / 3.0f) * static_cast<float>(PI) * SPHERE_RADIUS_MM * SPHERE_RADIUS_MM *
                              SPHERE_RADIUS_MM * density_kg_mm3;
    auto sphere_type = solver.LoadSphereType(sphere_mass, SPHERE_RADIUS_MM, particle_material);
    auto sphere_batch = solver.AddClumps(sphere_type, positions);
    auto sphere_tracker = solver.Track(sphere_batch);

    solver.SetTimeStepSize(time_step);
    solver.SetErrorOutVelocity(5.0e4f);
    solver.SetErrorOutAvgContacts(200.0f);
    // Match the active reference configuration, which uses ten percent of standard gravity.
    solver.SetGravitationalAcceleration(make_float3(0.0f, -0.1f * gravity_mm_s2, 0.0f));
    solver.Initialize();

    const path output_dir = current_path() / "modular_test_output" / "DEMTest_InwardBoxCornerContact";
    std::error_code output_error;
    create_directories(output_dir, output_error);
    if (output_error) {
        return fail("could not create output directory: " + output_error.message());
    }
    // Older versions of this test emitted full-system abnormal_contacts_*.csv snapshots. Remove only those obsolete
    // generated files so a rerun cannot make stale snapshots look like current diagnostic output.
    for (const auto& entry : directory_iterator(output_dir)) {
        const std::string filename = entry.path().filename().string();
        if (entry.is_regular_file() && filename.rfind("abnormal_contacts_", 0) == 0 &&
            entry.path().extension() == ".csv") {
            remove(entry.path(), output_error);
            if (output_error) {
                return fail("could not remove obsolete contact snapshot: " + output_error.message());
            }
        }
    }
    solver.WriteMeshFile(output_dir / "walls_static.vtk");

    std::ofstream events(output_dir / "abnormal_velocity_contacts.csv");
    events << std::setprecision(9)
           << "time,owner,speed_mm_s,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,contact_type,a_owner,b_owner,a_geo,b_geo,"
              "point_x,point_y,point_z,particle_to_point_x,particle_to_point_y,particle_to_point_z,"
              "particle_to_point_distance,particle_to_point_surface_offset,penetration,area,normal_x,normal_y,normal_z,"
              "force_x,force_y,force_z\n";

    bool abnormal_velocity_observed = false;
    std::set<bodyID_t> reported_owners;
    unsigned int output_frame = 0;
    double next_output_time = 0.0;
    while (solver.GetSimTime() < simulation_end - 0.5 * time_step) {
        if (solver.GetSimTime() + 0.5 * diagnostic_interval >= next_output_time) {
            char sphere_filename[64];
            std::snprintf(sphere_filename, sizeof(sphere_filename), "sphere_%04u.csv", output_frame++);
            solver.WriteSphereFile(output_dir / sphere_filename);
            std::cout << "t=" << solver.GetSimTime() << " s" << std::endl;
            next_output_time += output_interval;
        }

        const double current_time = solver.GetSimTime();
        const double remaining = simulation_end - current_time;
        const double advance_interval = current_time < diagnostic_start_time
                                            ? std::min(prediagnostic_interval, diagnostic_start_time - current_time)
                                            : diagnostic_interval;
        // DoDynamics lets kT and dT retain their normal asynchronous collaboration. The dT state queries below are
        // safe after this call returns and do not impose the full worker synchronization of DoDynamicsThenSync.
        solver.DoDynamics(std::min(advance_interval, remaining));

        // The reference failure occurs after 2.1 s. Avoid high-frequency device-to-host diagnostic queries during the
        // known-good prefix, then inspect every few integration steps once the failure window begins.
        if (solver.GetSimTime() <= diagnostic_start_time) {
            continue;
        }

        const auto particle_positions = sphere_tracker->Positions();
        const auto particle_velocities = sphere_tracker->Velocities();
        const auto owner_ids = sphere_tracker->GetOwnerIDs();
        float max_speed = 0.0f;
        std::vector<size_t> abnormal_particles;
        for (size_t particle = 0; particle < particle_velocities.size(); particle++) {
            const float3 velocity = particle_velocities[particle];
            if (!finite(particle_positions[particle]) || !finite(velocity)) {
                return fail("particle state became non-finite");
            }
            const float speed = length(velocity);
            max_speed = std::max(max_speed, speed);
            if (speed > ABNORMAL_VELOCITY_MM_S) {
                abnormal_particles.push_back(particle);
            }
        }

        if (abnormal_particles.empty()) {
            continue;
        }

        abnormal_velocity_observed = true;
        auto contacts = solver.GetContactDetailedInfo(DEME_TINY_FLOAT);
        for (size_t particle : abnormal_particles) {
            const bodyID_t owner = owner_ids[particle];
            const float speed = length(particle_velocities[particle]);
            const size_t related_contacts =
                writeAbnormalContacts(events, solver.GetSimTime(), owner, speed, particle_positions[particle],
                                      particle_velocities[particle], contacts);
            if (reported_owners.insert(owner).second) {
                std::cerr << "Abnormal particle owner " << owner << " reached " << speed
                          << " mm/s at t=" << solver.GetSimTime() << " s with " << related_contacts
                          << " related active contacts (" << contacts->Size() << " system-wide)." << std::endl;
            }
        }
    }

    // Always write the terminal particle state, even when it does not align exactly with the regular output interval.
    solver.WriteSphereFile(output_dir / "sphere_final.csv");
    if (abnormal_velocity_observed) {
        return fail("one or more particles exceeded the 1000 mm/s abnormal-velocity threshold; see " +
                    (output_dir / "abnormal_velocity_contacts.csv").string());
    }

    std::cout << "DEMTest_InwardBoxCornerContact passed: " << positions.size()
              << " particles remained below 1000 mm/s for " << simulation_end << " s." << std::endl;
    return 0;
}
