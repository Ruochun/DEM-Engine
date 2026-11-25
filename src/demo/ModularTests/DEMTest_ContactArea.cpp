// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "DEM/API.h"
#include "kernel/DEMHelperKernels.cuh"

using namespace deme;
namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}
void near(double actual,
          double expected,
          const char* message,
          double tolerance = 2e-4,
          double absolute_tolerance = 0.) {
    if (std::abs(actual - expected) > tolerance * std::max(std::abs(expected), 1e-10) + absolute_tolerance) {
        std::cerr << message << ": " << actual << " != " << expected << '\n';
        throw std::runtime_error(message);
    }
}
void setup(DEMSolver& solver) {
    solver.SetVerbosity("ERROR");
    solver.InstructBoxDomainDimension(8, 8, 8);
    solver.SetGravitationalAcceleration(make_float3(0));
    solver.SetTimeStepSize(1e-8);
    solver.SetCDUpdateFreq(1);
    solver.SetInitBinNumTarget(16);
    solver.SetContactOutputContent({"OWNER", "FORCE"});
    // Encode independently observable geometry in recorded force, keeping motion negligible.
    solver.DefineContactForceModel(R"(
        if (overlapDepth > 0.)
            force = make_float3(overlapArea, overlapDepth, rollingResistanceLength);
        else
            force = make_float3(overlapArea, overlapDepth, 1.f);
    )");
}

// Exercise radius policies without constructing a CUDA solver, so inference can also be checked on build hosts.
void testMeshRadiusHost() {
    DEMMesh mesh((GET_DATA_PATH() / "mesh/cube.obj").string());
    double volume;
    float3 center, inertia;
    mesh.ComputeMassProperties(volume, center, inertia);
    const double radius = std::cbrt(3. * volume / (4. * PI));
    near(mesh.GetContactRadius(), 0., "mass query leaves automatic radius unresolved");
    mesh.ResolveContactRadius(false);
    near(mesh.GetContactRadius(), 0., "legacy mesh defaults to wall");
    mesh.ResolveContactRadius(true);
    near(mesh.GetContactRadius(), radius, "closed solid infers radius");
    mesh.SetContactRadius(.4f);
    mesh.Scale(2.f);
    mesh.ResolveContactRadius(true);
    near(mesh.GetContactRadius(), .4f, "explicit radius is independent of geometry scaling");
    mesh.ClearContactRadius();
    mesh.ResolveContactRadius(true);
    near(mesh.GetContactRadius(), 2. * radius, "automatic radius follows geometry scaling");
    mesh.SetShellThickness(.01f);
    mesh.ResolveContactRadius(true);
    near(mesh.GetContactRadius(), 0., "shell does not infer solid curvature");
    mesh.DisableShell();
    mesh.m_face_v_indices.pop_back();
    mesh.ResolveContactRadius(true);
    near(mesh.GetContactRadius(), 0., "open solid falls back to wall");
    mesh.m_face_v_indices.clear();
    mesh.ResolveContactRadius(true);
    near(mesh.GetContactRadius(), 0., "failed volume falls back to wall");
}

// Radius inference is geometry-based, independent of mass queries, and respects copied template overrides.
void testRadiusLifecycle() {
    std::cout << "Testing mesh radius lifecycle" << std::endl;
    DEMSolver solver(1);
    setup(solver);
    solver.SetMeshUniversalContact(true);
    auto mat = solver.LoadMaterial({{"E", 1e5f}, {"nu", .25f}, {"CoR", 1.f}});
    auto mesh = solver.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat);
    mesh->Scale(.2f);
    double volume;
    float3 center, inertia;
    mesh->ComputeMassProperties(volume, center, inertia);
    near(mesh->GetContactRadius(), 0., "mass query is pure");
    const double radius = std::cbrt(3. * volume / (4. * PI));
    auto template_mesh = solver.LoadMeshType((GET_DATA_PATH() / "mesh/cube.obj").string(), mat);
    template_mesh->SetContactRadius(.4f);
    auto first_instance = solver.AddMeshFromTemplate(template_mesh, make_float3(2, 0, 0));
    template_mesh->SetContactRadius(.6f);
    auto later_instance = solver.AddMeshFromTemplate(template_mesh, make_float3(-2, 0, 0));
    auto tracker = solver.Track(mesh);
    solver.Initialize();
    near(first_instance->GetContactRadius(), .4f, "existing instance retains template override");
    near(later_instance->GetContactRadius(), .6f, "later instance inherits changed template override");
    near(mesh->GetContactRadius(), radius, "automatic volume radius");
    auto vertices = tracker->GetMesh()->GetCoordsVertices();
    for (auto& v : vertices)
        v *= 2.f;
    // Leave the CPU geometry cache dirty: full Update must synchronize it before radius inference.
    float3* device_vertices = nullptr;
    require(cudaMalloc(reinterpret_cast<void**>(&device_vertices), vertices.size() * sizeof(float3)) == cudaSuccess,
            "allocate deformation input");
    require(cudaMemcpy(device_vertices, vertices.data(), vertices.size() * sizeof(float3), cudaMemcpyHostToDevice) ==
                cudaSuccess,
            "upload deformation input");
    tracker->UpdateMeshFromDevice(device_vertices, solver.GetGPUDeviceIDs().front(), true, false);
    require(cudaFree(device_vertices) == cudaSuccess, "release deformation input");
    near(mesh->GetContactRadius(), radius, "lightweight deformation preserves radius");
    solver.Update();
    near(mesh->GetContactRadius(), 2. * radius, "full update refreshes radius");
    mesh->SetContactRadius(.7f);
    solver.Update();
    near(mesh->GetContactRadius(), .7f, "explicit override survives update");
    mesh->SetContactRadius(0.f);
    solver.Update();
    near(mesh->GetContactRadius(), 0.f, "explicit wall survives universal contact");
    mesh->ClearContactRadius();
    solver.Update();
    near(mesh->GetContactRadius(), 2. * radius, "clear override restores automatic radius");
    for (float invalid : {-1.f, std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()}) {
        bool failed = false;
        try {
            mesh->SetContactRadius(invalid);
        } catch (const SolverException&) {
            failed = true;
        }
        require(failed, "invalid radius must be rejected");
    }
    DEMMesh copy = *mesh;
    copy.SetContactRadius(.4f);
    DEMMesh instance = copy;
    copy.SetContactRadius(.6f);
    instance.ResolveContactRadius(false);
    near(instance.GetContactRadius(), .4f, "template instance owns override copy");
    instance.ClearContactRadius();
    instance.ResolveContactRadius(false);
    near(instance.GetContactRadius(), 0., "universal disabled defaults to wall");
    instance.SetShellThickness(.01f);
    instance.ResolveContactRadius(true);
    near(instance.GetContactRadius(), 0., "shell does not infer solid radius");
    instance.DisableShell();
    instance.m_face_v_indices.pop_back();
    instance.ResolveContactRadius(true);
    near(instance.GetContactRadius(), 0., "open mesh defaults to wall");
}

// Build one open planar patch; refinement changes only its triangulation, not its footprint or radius.
DEMMesh makeFlatPatch(int divisions, bool reverse, const std::shared_ptr<DEMMaterial>& mat) {
    DEMMesh mesh;
    for (int y = 0; y <= divisions; ++y)
        for (int x = 0; x <= divisions; ++x)
            mesh.m_vertices.push_back(make_float3(.2f * x / divisions - .1f, .2f * y / divisions - .1f, 0));
    for (int y = 0; y < divisions; ++y)
        for (int x = 0; x < divisions; ++x) {
            const int a = y * (divisions + 1) + x, b = a + 1, c = a + divisions + 1, d = c + 1;
            mesh.m_face_v_indices.push_back(make_int3(a, b, d));
            mesh.m_face_v_indices.push_back(make_int3(a, d, c));
        }
    if (reverse)
        std::reverse(mesh.m_face_v_indices.begin(), mesh.m_face_v_indices.end());
    mesh.nTri = mesh.m_face_v_indices.size();
    mesh.nPatches = 1;
    mesh.m_patch_ids.assign(mesh.nTri, 0);
    mesh.SetFamily(1);
    mesh.SetMass(1.f);
    mesh.SetMOI(make_float3(.1f));
    mesh.SetContactRadius(.2f);
    mesh.SetMaterial(mat);
    return mesh;
}

// A refined sphere-mesh patch precedes a sphere-plane contact in type order. The plane's compacted
// patch index lies inside the mesh's primitive range: wholesale primitive dispatch corrupts that range.
// Compare enabled/disabled plane contacts to preserve the exact same mesh geometry reference.
void testTypeSweep(bool in_place, bool no_record = false) {
    std::cout << "Testing type sweep, in_place=" << in_place << " no_record=" << no_record << std::endl;
    DEMSolver solver(1);
    setup(solver);
    solver.SetMeshUniversalContact(true);  // This path sorts mesh contacts before sphere-analytical contacts.
    solver.SetCollectAccRightAfterForceCalc(in_place);
    solver.SetNoForceRecord(no_record);
    solver.SetContactOutputContent({"OWNER", "FORCE", "TORQUE", "POINT", "NORMAL"});
    solver.DefineContactForceModel(R"(
        if (overlapDepth > 0.) {
            force = make_float3(overlapArea, overlapDepth, B2A.z);
            B_only_force = make_float3(.3f, .2f, .1f);
            A_only_torque = make_float3(.01f, .02f, .03f);
            B_only_torque = make_float3(.04f, .05f, .06f);
        }
    )");
    auto mat = solver.LoadMaterial({{"E", 1e5f}, {"nu", .25f}, {"CoR", 1.f}});
    auto mesh = makeFlatPatch(8, false, mat);
    mesh.SetInitPos(make_float3(0, 0, 1));
    auto mesh_tracker = solver.Track(solver.AddMesh(mesh));
    auto sphere = solver.LoadSphereType(1.f, .2f, mat);
    auto batch = solver.AddClumps(sphere, make_float3(0, 0, 1.19f));
    batch->SetFamily(1);
    auto sphere_tracker = solver.Track(batch);
    auto other = solver.AddClumps(sphere, make_float3(2, 0, .19f));
    other->SetFamily(2);
    auto plane = solver.AddBCPlane(make_float3(0), make_float3(0, 0, 1), mat);
    plane->SetFamily(3);
    solver.DisableContactBetweenFamilies(2, 3);
    solver.Initialize();
    auto check_vector = [](float3 actual, float3 expected) {
        require(std::isfinite(actual.x) && std::isfinite(actual.y) && std::isfinite(actual.z),
                "type sweep produced non-finite output");
        require(length(actual - expected) < 1e-4f * (1.f + length(expected)),
                "later contact type corrupted earlier mesh geometry or force");
    };
    {
        solver.DoDynamicsThenSync(1e-8);
        const auto accA = sphere_tracker->ContactAcc(), accB = mesh_tracker->ContactAcc();
        require(length(accA) > .5f && length(accB) > .5f, "regression requires active mesh force");
        std::shared_ptr<ContactInfoContainer> reference;
        if (!no_record) {
            reference = solver.GetContactDetailedInfo();
            require(reference->Size() == 1, "refined primitives must collapse to one mesh patch");
        }
        solver.EnableContactBetweenFamilies(2, 3);
        solver.RequestContactUpdate();
        solver.DoDynamicsThenSync(1e-8);
        require(solver.GetNumContacts() == 2, "mixed regression needs mesh and analytical patches");
        check_vector(sphere_tracker->ContactAcc(), accA);
        check_vector(mesh_tracker->ContactAcc(), accB);
        if (!no_record) {
            const auto mixed = solver.GetContactDetailedInfo();
            require(mixed->Size() == 2 && mixed->GetAOwner()[0] == sphere_tracker->GetOwnerID(),
                    "mesh patch must precede analytical contact");
            check_vector(mixed->GetForceA()[0], reference->GetForceA()[0]);
            check_vector(mixed->GetForceB()[0], reference->GetForceB()[0]);
            check_vector(mixed->GetTorqueA()[0], reference->GetTorqueA()[0]);
            check_vector(mixed->GetTorqueB()[0], reference->GetTorqueB()[0]);
            check_vector(mixed->GetPoint()[0], reference->GetPoint()[0]);
        }
        solver.DisableContactBetweenFamilies(2, 3);
        solver.RequestContactUpdate();
        solver.DoDynamicsThenSync(1e-8);
    }
    std::cout << "PASS: type sweep preserves mesh geometry and independent outputs" << std::endl;
}

// Check all contact types and the runtime RR policy independently of area-mode selection.
void testProjectedGeometry(bool in_place) {
    DEMSolver solver(1);
    setup(solver);
    solver.InstructBoxDomainDimension(12, 12, 8);
    solver.SetMeshUniversalContact(true);
    solver.SetCollectAccRightAfterForceCalc(in_place);
    solver.SetFamilyFixed(1);  // Keep the diagnostic geometry fixed, independent of diagnostic forces.
    solver.SetFamilyFixed(2);
    solver.SetFamilyExtraMargin(2, .02f);
    auto mat = solver.LoadMaterial({{"E", 1e5f}, {"nu", .25f}, {"CoR", 1.f}});
    auto sphere = solver.LoadSphereType(1.f, .2f, mat);
    // This separated pair is deliberately retained for margin-aware cohesive/custom models.
    auto margin_a = solver.AddClumps(sphere, make_float3(0, -2, 2));
    margin_a->SetFamily(2);
    solver.AddClumps(sphere, make_float3(0, -2, 2.41f))->SetFamily(2);
    auto margin_tracker = solver.Track(margin_a);
    std::vector<std::shared_ptr<DEMTracker>> primary;
    for (int kind = 0; kind < 5; ++kind) {
        const float x = 2.f * (kind - 2);
        const float z = kind == 1 || kind == 3 ? .19f : 2.39f;
        if (kind >= 3) {
            auto a = solver.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat);
            a->Scale(.4f);
            a->SetFamily(1);
            a->SetInitPos(make_float3(x, 0, z));
            a->SetContactRadius(.2f);
            primary.push_back(solver.Track(a));
        } else {
            auto a = solver.AddClumps(sphere, make_float3(x, 0, z));
            a->SetFamily(1);
            primary.push_back(solver.Track(a));
        }
        if (kind == 0)
            solver.AddClumps(sphere, make_float3(x, 0, 2))->SetFamily(1);
        if (kind == 2 || kind == 4) {
            auto b = solver.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat);
            b->Scale(.4f);
            b->SetInitPos(make_float3(x, 0, 2));
            b->SetContactRadius(0.f);  // Closed wall override must win over automatic volume inference.
        }
    }
    std::vector<std::shared_ptr<DEMTracker>> flat;
    for (int i = 0; i < 4; ++i) {
        auto mesh = makeFlatPatch(i < 2 ? 1 : (i < 4 ? 4 : 16), i % 2, mat);
        mesh.SetInitPos(make_float3(1.2f * (i - 2), 2.2f, -.001f));
        flat.push_back(solver.Track(solver.AddMesh(mesh)));
    }
    solver.AddBCPlane(make_float3(0), make_float3(0, 0, 1), mat);
    solver.Initialize();
    std::cout << "Initialized shared contact scene" << std::endl;

    for (int policy = 0; policy < 5; ++policy) {
        if (policy == 3)
            primary[4]->GetMesh()->SetContactRadius(0.f);
        const bool bounded = policy != 1 && policy != 4;
        solver.SetRollingResistanceLengthBoundsEnabled(bounded);
        solver.SetRollingResistanceLengthLimit(1.f, policy == 2 ? .03f : 0.f);
        solver.Update();
        solver.DoDynamicsThenSync(1e-8);
        const auto info = solver.GetContactDetailedInfo();
        std::vector<bool> found(5, false);
        std::vector<unsigned int> flat_count(flat.size(), 0);
        for (size_t j = 0; j < info->Size(); ++j) {
            const auto value = info->GetForce()[j];
            if (value.y <= 0.)
                continue;
            for (size_t kind = 0; kind < primary.size(); ++kind) {
                const auto owner = primary[kind]->GetOwnerID();
                if (info->GetAOwner()[j] != owner && info->GetBOwner()[j] != owner)
                    continue;
                found[kind] = true;
                near(value.y, .01, "known projected penetration", 2e-3);
                const double raw = value.x / (PI * value.y);
                double expected = !bounded ? raw : std::min(raw, kind == 0 ? .1 : .2);
                if (policy == 2)
                    expected = std::min(expected, .03);
                if (policy == 3 && kind == 4)
                    require(value.z < 0.f, "unresolved bound must expose a negative ingredient");
                else
                    near(value.z, expected, "runtime rolling-length policy");
            }
        }
        // A flat footprint must retain its area under refinement and triangle-order reversal.
        for (size_t j = 0; j < info->Size(); ++j)
            for (size_t i = 0; i < flat.size(); ++i)
                if (info->GetAOwner()[j] == flat[i]->GetOwnerID()) {
                    ++flat_count[i];
                    const auto value = info->GetForce()[j];
                    near(value.x, .04, "refinement-independent projected footprint");
                    near(value.y, .001, "flat penetration", 2e-3);
                    const double raw = value.x / (PI * value.y);
                    double expected = bounded ? std::min(raw, .2) : raw;
                    if (policy == 2)
                        expected = std::min(expected, .03);
                    near(value.z, expected, "flat rolling length");
                }
        for (auto count : flat_count)
            require(count == 1, "flat primitives must collapse to one patch");
        for (bool hit : found)
            require(hit, "missing projected contact type");
    }
    std::cout << "PASS: all projected contact types and runtime rolling bounds" << std::endl;
}

// Exercise the producer contract independently of GPU state, including bounded overflow avoidance.
void testRollingLimits() {
    const auto spin_only =
        contactResistanceTorque<float>(make_float3(0, 0, 1), make_float3(0, 0, 1), 1.f, 1.f, 1e-40f, 0.f, .2f, -1.f);
    require(std::isfinite(spin_only.z) && spin_only.z < 0.f, "spin-only torque needs no rolling scale");
    for (double depth : {1e-3, 1e-6, 1e-12, 1e-30}) {
        near(resolveRollingResistanceLength<double>(.1, 1., depth, .5, 0.), .05, "intrinsic cap");
        near(resolveRollingResistanceLength<double>(0., 1., depth, 1., .03), .03, "absolute wall cap");
        near(resolveRollingResistanceLength<double>(0., 1., depth, 1., 0., false), 1. / (PI * depth),
             "unbounded wall length");
    }
    near(resolveRollingResistanceLength<float>(.1f, 1.f, 1e-40f, 1.f, 0.f), .1, "cap before float narrowing");
    require(resolveRollingResistanceLength<double>(0., 1., .1, 1., 0.) < 0., "unresolved bound marker");
    for (double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        bool rejected = false;
        try {
            resolveRollingResistanceLength<double>(.1, invalid, .01, 1., 0.);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        require(rejected, "producer must reject non-finite geometry");
    }
    bool overflow_rejected = false;
    try {
        resolveRollingResistanceLength<float>(0.f, 1.f, 1e-40f, 1.f, 0.f, false);
    } catch (const std::runtime_error&) {
        overflow_rejected = true;
    }
    require(overflow_rejected, "producer must reject unbounded float overflow");
}
}  // namespace
int main(int argc, char** argv) {
    try {
        testRollingLimits();
        testMeshRadiusHost();
        if (argc > 1 && std::string(argv[1]) == "--host-only") {
            std::cout << "PASS: mesh radius policies and rolling length bounds" << std::endl;
            return 0;
        }
        if (argc > 1 && std::string(argv[1]) == "--type-sweep") {
            testTypeSweep(false);
            testTypeSweep(true);
            testTypeSweep(true, true);
            return 0;
        }
        testTypeSweep(false);
        testTypeSweep(true);
        testTypeSweep(true, true);
        // Immediate setup must not race background allocation when rapidly constructing successive solvers.
        for (int i = 0; i < 8; ++i) {
            DEMSolver startup(1);
            startup.SetFamilyExtraMargin(2, .02f);
        }
        testRadiusLifecycle();
        testProjectedGeometry(false);
        testProjectedGeometry(true);
        std::cout << "PASS: projected contact sweep, radius lifecycle and rolling bounds\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
