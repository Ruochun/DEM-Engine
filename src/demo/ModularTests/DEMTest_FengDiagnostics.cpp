// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#include "DEM/API.h"
#include "FengTestUtils.hpp"
#include <iostream>
#include <stdexcept>

using namespace deme;
namespace {
void require(bool ok, const char* message) {
    if (!ok)
        throw std::runtime_error(message);
}

// Weld shared vertices so DEME's watertight check sees a closed cube, even after face refinement.
DEMMesh cubeMesh(int refinement, bool flip) {
    DEMMesh mesh;
    for (const auto& triangle : feng_test::cube(refinement, flip)) {
        int ids[3];
        for (int k = 0; k < 3; ++k) {
            const float3 vertex = to_float3(triangle[k]);
            size_t index = 0;
            for (; index < mesh.m_vertices.size(); ++index) {
                if (length(mesh.m_vertices[index] - vertex) < 1e-7f)
                    break;
            }
            if (index == mesh.m_vertices.size())
                mesh.m_vertices.push_back(vertex);
            ids[k] = static_cast<int>(index);
        }
        mesh.m_face_v_indices.push_back(make_int3(ids[0], ids[1], ids[2]));
    }
    mesh.nTri = mesh.m_face_v_indices.size();
    mesh.SetMass(1.f);
    mesh.SetMOI(make_float3(1.f));
    return mesh;
}

// Compare enabled/disabled solver motion and validate the GPU reduction against an analytic cube contact. Add a
// sphere-sphere contact too, so mesh keys have a nonzero global offset (a common patch-index integration pitfall).
std::vector<float3> run(bool enable, int refinement, bool flip, bool splitPatches = false) {
    DEMSolver sim;
    sim.SetVerbosity("ERROR");
    sim.InstructBoxDomainDimension(10.f, 10.f, 10.f);
    sim.SetGravitationalAcceleration(make_float3(0.f));
    sim.SetMeshUniversalContact(true);
    sim.SetSimplePatchCombination(true);
    sim.SetMeshMeshFengDiagnostics(enable);
    sim.SetTimeStepSize(1e-8f);
    auto material = sim.LoadMaterial({{"E", 1e5f}, {"nu", .3f}, {"CoR", 1.f}, {"mu", 0.f}});
    auto sphere = sim.LoadSphereType(1.f, .1f, material);
    sim.AddClumps(sphere, std::vector<float3>{make_float3(-2.f, 0.f, 0.f), make_float3(-1.85f, 0.f, 0.f)});
    auto a = cubeMesh(refinement, flip), b = cubeMesh(refinement, !flip);
    a.SetMaterial(material);
    b.SetMaterial(material);
    if (splitPatches) {
        a.SetEachTriangleAsPatch();
        b.SetEachTriangleAsPatch();
    }
    auto ah = sim.AddMesh(a), bh = sim.AddMesh(b);
    ah->SetInitPos(make_float3(0.f));
    bh->SetInitPos(make_float3(.37f, .23f, .19f));
    sim.Initialize();
    require(sim.GetMeshMeshFengDiagnostics().empty(), "snapshot before force evaluation");
    sim.DoDynamics(1e-8);
    const auto records = sim.GetMeshMeshFengDiagnostics();
    const auto velocities = sim.GetOwnerVelocity(0, 4);
    if (enable) {
        require(!records.empty(), "GPU diagnostic records");
        if (!splitPatches) {
            require(records.size() == 1, "one cube patch pair");
            const auto& r = records[0];
            require(r.patchContact > 0, "mixed contact type global patch offset");
            require(r.ownersWatertight, "watertight cube owners");
            require(r.closurePassed, "GPU cube boundary closure");
            require(r.hasContactLine, "GPU cube contact line");
            const double3 expected = make_double3(.77 * .81, .63 * .81, .63 * .77);
            require(length(r.vectorArea - expected) < 2e-5, "GPU cube vector area");
            const double3 surfaceMoment = cross(make_double3(1, .615, .595), make_double3(expected.x, 0, 0)) +
                                          cross(make_double3(.685, 1, .595), make_double3(0, expected.y, 0)) +
                                          cross(make_double3(.685, .615, 1), make_double3(0, 0, expected.z));
            const double3 expectedMoment = surfaceMoment - cross(r.referenceOrigin, expected);
            require(length(r.geometricMoment - expectedMoment) < 2e-5, "GPU cube surface moment");
            const double3 n = expected / length(expected);
            const double3 localPoint = r.fengContactPoint - r.referenceOrigin;
            require(length(cross(localPoint, expected) - (expectedMoment - n * dot(n, expectedMoment))) < 2e-5,
                    "GPU contact line reproduces transverse moment");
            require(fabs(dot(r.fengContactPoint - r.legacyContactPoint, n)) < 2e-5,
                    "GPU contact point nearest legacy anchor");
            require(length(r.fengNormal - expected * (-1.0 / length(expected))) < 2e-5, "GPU cube B2A");
            require(std::isfinite(r.legacyArea) && std::isfinite(r.legacyPenetration), "legacy geometry snapshot");
            std::cout << "cube n=" << refinement << " flip=" << flip << " segments=" << r.segmentCount
                      << " legacy area=" << r.legacyArea << " boundary area=" << r.fengArea << '\n';
        } else {
            for (const auto& r : records)
                require(!r.hasContactLine, "open per-triangle patches must fail line gate");
        }
        sim.SetMeshMeshFengDiagnostics(false);
        require(sim.GetMeshMeshFengDiagnostics().empty(), "disable clears snapshot");
        sim.SetMeshMeshFengDiagnostics(true);
        require(sim.GetMeshMeshFengDiagnostics().empty(), "reenable does not expose stale snapshot");
        sim.SetOwnerPosition(3, std::vector<float3>{make_float3(3.f, 2.f, 2.f)});
        sim.RequestContactUpdate();
        sim.DoDynamicsThenSync(1e-8);
        // The first force evaluation may still use the old candidate list; it must nevertheless contain fresh
        // geometry at the moved pose. The next call consumes the contact detection result produced by the sync.
        for (const auto& row : sim.GetMeshMeshFengDiagnostics()) {
            require(row.segmentCount == 0.0 && !row.hasContactLine, "stale geometry after moving mesh");
        }
        sim.DoDynamicsThenSync(1e-8);
        require(sim.GetMeshMeshFengDiagnostics().empty(), "contact-free step clears mesh snapshot");
    } else
        require(records.empty(), "disabled diagnostics");
    return velocities;
}
}  // namespace
int main() {
    try {
        const auto baseline = run(false, 1, false);
        require(length(baseline[0]) > 0.f && length(baseline[1]) > 0.f, "active sphere forces in baseline");
        const auto diagnostic = run(true, 1, false);
        for (size_t i = 0; i < baseline.size(); ++i)
            require(length(baseline[i] - diagnostic[i]) == 0.f, "diagnostics changed motion");
        run(true, 2, true);
        run(true, 4, false);
        run(true, 1, false, true);
        std::cout << "PASS: GPU Feng diagnostics, refinement, open patches, mixed keys and unchanged motion.\n";
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
