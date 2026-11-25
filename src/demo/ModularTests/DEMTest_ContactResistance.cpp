// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "DEM/API.h"
#include "DEM/utils/ForceModelHelpers.hpp"
#include "DEM/utils/ContactOutput.hpp"
#include "kernel/DEMHelperKernels.cuh"

using namespace deme;

namespace {
void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

void near(float3 actual, float3 expected, const std::string& message, float tolerance = 2e-3f) {
    if (!(length(actual - expected) <= tolerance * (1.f + length(expected)))) {
        std::cerr << message << ": actual=(" << actual.x << "," << actual.y << "," << actual.z << ") expected=("
                  << expected.x << "," << expected.y << "," << expected.z << ")\n";
        require(false, message);
    }
}

// Distinct A/B values expose swapped sides and verify that legacy accessors alias A storage.
void testContactOutputFields() {
    const unsigned int flags = CNT_OUTPUT_CONTENT::FORCE | CNT_OUTPUT_CONTENT::TORQUE;
    ContactInfoContainer info(flags, {});
    info.ResizeAll(1);
    info.GetContactType()[0] = "SS";
    info.GetForceA()[0] = make_float3(1, 2, 3);
    info.GetForceB()[0] = make_float3(4, 5, 6);
    info.GetTorqueA()[0] = make_float3(7, 8, 9);
    info.GetTorqueB()[0] = make_float3(10, 11, 12);
    require(&info.GetForce() == &info.GetForceA(), "legacy force accessor must alias A");
    require(&info.GetTorque() == &info.GetTorqueA(), "legacy torque accessor must alias A");
    require(&info.Get<float3>("Force") == &info.GetForceA(), "legacy force key must alias A");
    require(&info.Get<float3>("Torque") == &info.GetTorqueA(), "legacy torque key must alias A");
    std::ostringstream output;
    writeContactSnapshotAsCsv(output, info, flags, {});
    require(output.str().find("contact_type,f_a_x,f_a_y,f_a_z,f_b_x,f_b_y,f_b_z,"
                              "torque_a_x,torque_a_y,torque_a_z,torque_b_x,torque_b_y,torque_b_z\n") == 0,
            "CSV headers must identify both sides explicitly");
    for (const auto& column : {"f_a_x", "f_b_z", "torque_a_y", "torque_b_z"})
        require(CNT_FILE_KNOWN_COL_NAMES.count(column) != 0, "built-in A/B columns must not become wildcards");
    require(output.str().find("SS,1,2,3,4,5,6,7,8,9,10,11,12\n") != std::string::npos,
            "CSV must preserve both independent forces and torques");
}

// CPU checks exercise the exact helper compiled into the force kernels, including its rest/regularization branches.
void testConstitutiveLaw() {
    const float3 n = make_float3(1, 0, 0);
    auto torque = [&](float3 omega, float load = 10.f, float overlap = .02f) {
        return contactResistanceTorque<float>(omega, n, load, .1f, overlap, .2f, .3f, .5f);
    };
    // Double instantiation must preserve precision rather than silently taking a float-only math path.
    const double3 double_torque =
        contactResistanceTorque<double>(make_double3(2, 3, 0), make_double3(1, 0, 0), 10., .1, .02, .2, .3, .5);
    require(length(double_torque - make_double3(-(3. * PI / 16.) * .3, -1., 0.)) < 1e-12,
            "double-precision resistance magnitude");
    const double3 double_zero =
        contactResistanceTorque<double>(make_double3(2, 3, 0), make_double3(1, 0, 0), 0., .1, .02, .2, .3, .5);
    require(length(double_zero) == 0., "double-precision zero load");
    const float spin = (3.f * PI / 16.f) * .3f;
    near(torque(make_float3(2, 3, 0)), make_float3(-spin, -1, 0), "mixed rolling/spinning magnitude", 1e-6f);
    near(torque(make_float3(0)), make_float3(0), "zero motion", 1e-6f);
    near(torque(make_float3(1), 0), make_float3(0), "zero normal load", 1e-6f);
    near(torque(make_float3(1), 10, 0), make_float3(0), "zero overlap", 1e-6f);
    near(torque(make_float3(5e-4f, 0, 0)), make_float3(-.5f * spin, 0, 0), "regularized spin", 1e-6f);
    for (float sign : {-1.f, 1.f}) {
        const float3 w = sign * make_float3(2, 3, 4);
        require(dot(torque(w), w) < 0.f, "resistance must dissipate power");
        near(torque(-1.f * w), -1.f * torque(w), "reversal symmetry", 1e-6f);
    }
}

// Deprecation is rejected when defining the model, before any GPU setup. Comments and diagnostics are legal.
void testDeprecation() {
    DEMForceModel model;
    model.DefineCustomModel(R"code(
        // torque_only_force is removed.
        /* torque_only_force */
        const char* message = "torque_only_force";
        const char* raw = R"note(torque_only_force)note";
        float torque_only_force_suffix = 0;
        A_only_torque = make_float3(0);
    )code");
    for (bool prerequisites : {false, true}) {
        bool rejected = false;
        try {
            if (prerequisites)
                model.DefineCustomModelPrerequisites("torque_only_force = make_float3(0);");
            else
                model.DefineCustomModel("torque_only_force/*comment*/=make_float3(0);");
        } catch (const SolverException& error) {
            rejected = std::string(error.what()).find("deprecated") != std::string::npos;
        }
        require(rejected, "removed output must produce deprecation error");
    }
}

// One isolated sphere pair makes recorded forces, tracker resultants, and applied acceleration independently
// comparable. Rotating the local inertia frames catches accidentally applying global torques as local vectors.
void testSolver(bool in_place, bool custom, bool no_record = false) {
    std::cout << "Contact outputs: in_place=" << in_place << " custom=" << custom << " no_record=" << no_record
              << std::endl;
    DEMSolver solver(1);
    solver.SetVerbosity("ERROR");
    solver.InstructBoxDomainDimension(4, 4, 4);
    solver.SetGravitationalAcceleration(make_float3(0));
    solver.SetTimeStepSize(1e-7);
    solver.SetCDUpdateFreq(1);
    solver.SetCollectAccRightAfterForceCalc(in_place);
    solver.SetNoForceRecord(no_record);
    solver.SetContactOutputContent({"OWNER", "FORCE", "TORQUE", "POINT", "NORMAL", "CNT_WILDCARD"});
    if (custom) {
        auto model = solver.DefineContactForceModel(R"(
            // The B-only torque case must survive a zero A force/torque and a zero contact force.
            if (overlapDepth > 0) {
                if (mode[AOwner] < 0.5f) {
                    force = make_float3(2, 3, 5);
                    A_only_force = make_float3(7, 11, 13);
                    B_only_force = make_float3(17, 19, 23);
                    A_only_torque = make_float3(29, 31, 37);
                    B_only_torque = make_float3(41, 43, 47);
                } else {
                    B_only_torque = make_float3(41, 43, 47);
                }
            }
        )");
        model->SetPerOwnerWildcards({"mode"});
    }
    auto material = solver.LoadMaterial({{"E", 1e5f}, {"nu", .25f}, {"CoR", 1.f}, {"mu", 0.f}, {"Crr", .2f}});
    // Exercise both the optional zero default and the explicit material-pair override.
    if (!no_record)
        solver.SetMaterialPropertyPair("Cspin", material, material, .3f);
    auto sphere = solver.LoadSphereType(2.f, .5f, material);
    sphere->SetMOI(make_float3(.2f));
    auto batch = solver.AddClumps(sphere, {make_float3(-.49f, 0, 0), make_float3(.49f, 0, 0)});
    const float q = std::sqrt(.5f);
    batch->SetOriQ({make_float4(0, 0, q, q), make_float4(q, 0, 0, q)});
    auto tracker = solver.Track(batch);
    solver.Initialize();
    const std::vector<float3> motions = {make_float3(2, 0, 0), make_float3(0, 3, 0), make_float3(2, 3, 0),
                                         make_float3(0), make_float3(5e-4f, 0, 0)};
    const int cases = custom ? 2 : static_cast<int>(motions.size());
    for (int c = 0; c < cases; ++c) {
        tracker->SetPos({make_float3(-.49f, 0, 0), make_float3(.49f, 0, 0)});
        tracker->SetVel(std::vector<float3>(2, make_float3(0)));
        // Both bodies share a common angular velocity in case 3: resistance must be exactly zero.
        const float3 wa = custom ? make_float3(0) : motions[c];
        const float3 wb = (!custom && c == 3) ? make_float3(1, 2, 3) : make_float3(0);
        const float3 global_a = (!custom && c == 3) ? wb : wa;
        for (size_t owner = 0; owner < 2; ++owner) {
            // Use an exactly representable common frame for the zero-relative-motion case. Different float
            // quaternion round trips otherwise leave a tiny relative velocity, amplified by the regularized slope.
            if (!custom && c == 3)
                tracker->SetOriQ(make_float4(0, 0, 0, 1), owner);
            float3 local = owner == 0 ? global_a : wb;
            const float4 orientation = tracker->OriQ(owner);
            applyOriQToVector3(local, make_float4(-orientation.x, -orientation.y, -orientation.z, orientation.w));
            tracker->SetAngVel(local, owner);
        }
        if (custom)
            tracker->SetOwnerWildcardValues("mode", {static_cast<float>(c), static_cast<float>(c)});
        solver.RequestContactUpdate();
        solver.DoDynamicsThenSync(1e-7);
        if (no_record) {
            // Omitted Cspin disables pure spin; rolling must still work without stored contact results.
            if (c == 0)
                near(tracker->ContactAngAccGlobal(0), make_float3(0), "default Cspin is zero");
            if (c == 1)
                require(tracker->ContactAngAccGlobal(0).y < 0.f, "no-record rolling decelerates A");
            if (c == 3)
                near(tracker->ContactAngAccGlobal(0), make_float3(0), "no-record common rotation");
            continue;
        }
        auto info = solver.GetContactDetailedInfo();
        require(info->Size() == 1, "expected one contact record");
        const size_t a = info->GetAOwner()[0] - tracker->GetOwnerID();
        const size_t b = info->GetBOwner()[0] - tracker->GetOwnerID();
        const float3 fa = info->GetForce()[0], fb = info->GetForceB()[0];
        const float3 ta = info->GetTorque()[0], tb = info->GetTorqueB()[0];
        if (custom) {
            near(fa, c == 0 ? make_float3(9, 14, 18) : make_float3(0), "independent A force");
            near(fb, c == 0 ? make_float3(15, 16, 18) : make_float3(0), "independent B force");
            near(ta, c == 0 ? make_float3(29, 31, 37) : make_float3(0), "independent A torque");
            near(tb, make_float3(41, 43, 47), "independent B torque");
        } else {
            near(fa + fb, make_float3(0), "force balance");
            near(ta + tb, make_float3(0), "couple balance");
            const float3 relative = a == 0 ? global_a - wb : wb - global_a;
            const float3 n = info->GetNormal()[0];
            // Existing sphere geometry supplies area pi*(R*delta-delta^2/4), not pi*R_eff*delta.
            const float radius = std::sqrt(.5f * .02f - .02f * .02f / 4.f);
            const float load = dot(fa, n);
            near(ta, contactResistanceTorque<float>(relative, n, load, radius, .02f, .2f, .3f, .25f),
                 "Hertzian resistance magnitude and direction");
            require(dot(ta, relative) <= 1e-5f, "contact couple dissipates energy");
        }
        for (size_t owner : {a, b}) {
            const float3 f = owner == a ? fa : fb, t = owner == a ? ta : tb;
            near(tracker->ContactAcc(owner), f / 2.f, "applied linear acceleration");
            const float3 arm = info->GetPoint()[0] - tracker->Pos(owner);
            near(tracker->ContactAngAccGlobal(owner), (cross(arm, f) + t) / .2f, "applied angular acceleration", 4e-3f);
            // Tracker force-pair queries are explicitly disabled by the existing in-place reduction API.
            // Detailed records and applied accelerations above still verify that path's independent outputs.
            if (in_place)
                continue;
            std::vector<float3> points, forces, torques;
            tracker->GetContactForcesAndGlobalTorque(points, forces, torques, owner);
            if (length(f) + length(t) > 0.f) {
                require(forces.size() == 1, "owner query must retain torque-only contact");
                near(forces[0], f, "owner force reporting");
                near(torques[0], t, "owner direct torque reporting");
                tracker->GetContactForcesAndLocalTorque(points, forces, torques, owner);
                float3 local_torque = t;
                const float4 orientation = tracker->OriQ(owner);
                applyOriQToVector3(local_torque,
                                   make_float4(-orientation.x, -orientation.y, -orientation.z, orientation.w));
                near(torques[0], local_torque, "owner-local direct torque reporting");
            }
        }
    }
    // Separation must not retain either previous torque or force.
    tracker->SetPos(make_float3(1.5f, 0, 0), 1);
    solver.RequestContactUpdate();
    solver.DoDynamicsThenSync(1e-7);
    require(solver.GetNumContacts() == 0, "separation must clear patch contact count");
    for (size_t owner = 0; owner < 2; ++owner) {
        near(tracker->ContactAcc(owner), make_float3(0), "separated force");
        near(tracker->ContactAngAccGlobal(owner), make_float3(0), "separated torque");
    }
    if (!custom && !no_record) {
        // A new contact after an empty step must not recover the old sliding displacement or direct torque.
        tracker->SetPos({make_float3(-.49f, 0, 0), make_float3(.49f, 0, 0)});
        tracker->SetVel(std::vector<float3>(2, make_float3(0)));
        tracker->SetAngVel(std::vector<float3>(2, make_float3(0)));
        solver.RequestContactUpdate();
        solver.DoDynamicsThenSync(1e-7);
        const auto fresh = solver.GetContactDetailedInfo();
        require(fresh->Size() == 1, "expected re-contact");
        near(make_float3(fresh->Get<float>("delta_tan_x")[0], fresh->Get<float>("delta_tan_y")[0],
                         fresh->Get<float>("delta_tan_z")[0]),
             make_float3(0), "re-contact history reset", 1e-7f);
        near(fresh->GetTorque()[0], make_float3(0), "re-contact torque reset");
    }
}
// Mesh contacts follow a separate patch kernel. A B-only normal force also checks that triangle pressure
// is derived from the mesh's own force, even when the sphere's force is zero.
void testPatchOutputs(bool in_place) {
    std::cout << "Patch outputs: in_place=" << in_place << std::endl;
    DEMSolver solver(1);
    solver.SetVerbosity("ERROR");
    solver.InstructBoxDomainDimension(2, 2, 2);
    solver.SetGravitationalAcceleration(make_float3(0));
    solver.SetTimeStepSize(1e-7);
    solver.SetCDUpdateFreq(1);
    solver.SetInitBinNumTarget(8);
    solver.SetCollectAccRightAfterForceCalc(in_place);
    solver.SetContactOutputContent({"OWNER", "FORCE", "TORQUE", "POINT", "NORMAL", "CNT_WILDCARD"});
    solver.DefineContactForceModel(R"(
        if (overlapDepth > 0) {
            B_only_force = -7.f * B2A;
            A_only_torque = make_float3(.1f, .2f, .3f);
            B_only_torque = make_float3(.4f, .5f, .6f);
        }
    )");
    auto mat = solver.LoadMaterial({{"E", 1e5f}, {"nu", .25f}, {"CoR", 1.f}});
    auto mesh = solver.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat);
    mesh->Scale(.2f);
    mesh->SetMass(2.f);
    mesh->SetMOI(make_float3(.2f));
    auto mesh_tracker = solver.Track(mesh);
    auto sphere = solver.LoadSphereType(2.f, .1f, mat);
    sphere->SetMOI(make_float3(.2f));
    auto batch = solver.AddClumps(sphere, make_float3(0, 0, .19f));
    auto sphere_tracker = solver.Track(batch);
    solver.Initialize();
    solver.SetTrianglePVTrackingOwners({mesh_tracker->GetOwnerID()});
    solver.DoDynamicsThenSync(1e-7);
    const auto info = solver.GetContactDetailedInfo();
    require(info->Size() > 0, "expected sphere-mesh patch contact");
    for (const auto& tracker : {mesh_tracker, sphere_tracker}) {
        float3 total_force = make_float3(0), total_torque = make_float3(0);
        for (size_t i = 0; i < info->Size(); ++i) {
            const bool a = info->GetAOwner()[i] == tracker->GetOwnerID();
            const float3 force = a ? info->GetForce()[i] : info->GetForceB()[i];
            const float3 torque = a ? info->GetTorque()[i] : info->GetTorqueB()[i];
            total_force += force;
            total_torque += cross(info->GetPoint()[i] - tracker->Pos(), force) + torque;
        }
        near(tracker->ContactAcc(), total_force / 2.f, "patch linear acceleration");
        near(tracker->ContactAngAccGlobal(), total_torque / .2f, "patch angular acceleration");
    }
    near(sphere_tracker->ContactAcc(), make_float3(0), "patch A force remains zero");
    require(length(mesh_tracker->ContactAcc()) > 0.f, "patch B force applied");
    std::vector<float> pressure, speed, pv;
    require(solver.GetTrackedOwnerTrianglePV(mesh_tracker->GetOwnerID(), pressure, speed, pv, false),
            "mesh pressure query");
    float total_pressure = 0.f;
    for (float value : pressure)
        total_pressure += value;
    require(total_pressure > 0.f, "B-only force must contribute mesh pressure");
}
}  // namespace

int main(int argc, char** argv) {
    try {
        testContactOutputFields();
        testConstitutiveLaw();
        testDeprecation();
        if (argc > 1 && std::string(argv[1]) == "--host-only") {
            std::cout << "PASS: contact resistance constitutive law and deprecation\n";
            return 0;
        }
        for (bool in_place : {false, true}) {
            testSolver(in_place, true);
            testSolver(in_place, false);
            testPatchOutputs(in_place);
        }
        testSolver(true, false, true);
        std::cout << "PASS: direct contact outputs, rolling/spinning, both reduction paths and no-record mode\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
