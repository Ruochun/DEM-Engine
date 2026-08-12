// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#include <cuda_runtime_api.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "DEM/API.h"
#include "core/utils/DataMigrationHelper.hpp"

using namespace deme;

namespace {

constexpr float kTolerance = 2e-5f;

bool close(float a, float b) {
    return std::abs(a - b) <= kTolerance;
}

bool close(const float3& a, const float3& b) {
    return close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z);
}

bool close(const float4& a, const float4& b) {
    return close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z) && close(a.w, b.w);
}

template <typename T>
class TestDeviceBuffer {
  public:
    TestDeviceBuffer(size_t count, int device) : m_count(count), m_device(device) {
        ScopedCudaDevice scope(device);
        DEME_GPU_CALL(cudaMalloc(reinterpret_cast<void**>(&m_pointer), count * sizeof(T)));
    }

    ~TestDeviceBuffer() {
        ScopedCudaDevice scope(m_device);
        (void)cudaFree(m_pointer);
    }

    T* data() { return m_pointer; }
    size_t size() const { return m_count; }
    int device() const { return m_device; }

    std::vector<T> ToHost(size_t count = 0) const {
        const size_t output_count = count == 0 ? m_count : count;
        std::vector<T> output(output_count);
        ScopedCudaDevice scope(m_device);
        DEME_GPU_CALL(cudaMemcpy(output.data(), m_pointer, output_count * sizeof(T), cudaMemcpyDeviceToHost));
        return output;
    }

  private:
    T* m_pointer = nullptr;
    size_t m_count;
    int m_device;
};

template <typename T>
bool compareVectors(const std::vector<T>& actual, const std::vector<T>& expected, const std::string& label) {
    if (actual.size() != expected.size()) {
        std::cerr << "FAIL: " << label << " size mismatch." << std::endl;
        return false;
    }
    for (size_t i = 0; i < actual.size(); i++) {
        if (!close(actual[i], expected[i])) {
            std::cerr << "FAIL: " << label << " differs at index " << i << "." << std::endl;
            return false;
        }
    }
    return true;
}

template <>
bool compareVectors<unsigned int>(const std::vector<unsigned int>& actual,
                                  const std::vector<unsigned int>& expected,
                                  const std::string& label) {
    if (actual != expected) {
        std::cerr << "FAIL: " << label << " differs." << std::endl;
        return false;
    }
    return true;
}

bool testFixedOwnerData(int visible_devices) {
    DEMSolver solver(1);
    solver.SetVerbosity("ERROR");
    solver.InstructBoxDomainDimension(8, 8, 8);
    solver.SetGravitationalAcceleration(make_float3(0, 0, 0));
    solver.SetCDUpdateFreq(1);

    auto force_model = solver.GetContactForceModel();
    force_model->SetPerOwnerWildcards({"retrieval_tag"});
    auto material = solver.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.2}, {"mu", 0.4}, {"Crr", 0.0}});
    auto sphere = solver.LoadSphereType(2.5f, 0.5f, material);

    const std::vector<float3> positions = {make_float3(1.0f, 1.0f, 1.0f), make_float3(1.8f, 1.0f, 1.0f),
                                           make_float3(4.0f, 2.0f, 1.5f)};
    const std::vector<float3> velocities = {make_float3(1, 2, 3), make_float3(-1, 0.5f, 2), make_float3(0, -2, 1)};
    const std::vector<float3> angular_velocities = {make_float3(1, 0, 0), make_float3(0, 2, 0), make_float3(0, 0, 3)};
    const float half_sqrt_two = std::sqrt(0.5f);
    const std::vector<float4> orientations = {make_float4(0, 0, half_sqrt_two, half_sqrt_two),
                                              make_float4(half_sqrt_two, 0, 0, half_sqrt_two), make_float4(0, 0, 0, 1)};
    const std::vector<unsigned int> families = {2, 3, 4};
    const std::vector<float> wildcards = {10.5f, -2.0f, 7.25f};

    DEMClumpBatch input_batch(3);
    input_batch.SetTypes(sphere);
    input_batch.SetPos(positions);
    input_batch.SetVel(velocities);
    input_batch.SetAngVel(angular_velocities);
    input_batch.SetOriQ(orientations);
    input_batch.SetFamilies(families);
    input_batch.AddOwnerWildcard("retrieval_tag", wildcards);
    auto batch = solver.AddClumps(input_batch);
    auto tracker = solver.Track(batch);

    solver.SetInitTimeStep(1e-5);
    solver.Initialize();
    solver.DoDynamicsThenSync(1e-5);

    const size_t count = positions.size();
    TestDeviceBuffer<float3> float3_output(count, 0);
    TestDeviceBuffer<float4> float4_output(count, 0);
    TestDeviceBuffer<float> float_output(count, 0);
    TestDeviceBuffer<unsigned int> uint_output(count, 0);

    tracker->PositionsToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->Positions(), "positions"))
        return false;
    tracker->VelocitiesToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->Velocities(), "velocities"))
        return false;
    tracker->AngularVelocitiesLocalToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->AngularVelocitiesLocal(), "local angular velocities"))
        return false;
    tracker->AngularVelocitiesGlobalToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->AngularVelocitiesGlobal(), "global angular velocities"))
        return false;
    tracker->OrientationQuaternionsToDevice(float4_output.data(), count, 0);
    if (!compareVectors(float4_output.ToHost(), tracker->OrientationQuaternions(), "orientations"))
        return false;
    tracker->ContactAccelerationsToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->ContactAccelerations(), "contact accelerations"))
        return false;
    tracker->ContactAngularAccelerationsLocalToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->ContactAngularAccelerationsLocal(),
                        "local contact angular accelerations"))
        return false;
    tracker->ContactAngularAccelerationsGlobalToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->ContactAngularAccelerationsGlobal(),
                        "global contact angular accelerations"))
        return false;
    tracker->FamiliesToDevice(uint_output.data(), count, 0);
    if (!compareVectors(uint_output.ToHost(), tracker->GetFamilies(), "families"))
        return false;
    tracker->MassesToDevice(float_output.data(), count, 0);
    if (!compareVectors(float_output.ToHost(), tracker->Masses(), "masses"))
        return false;
    tracker->MOIsToDevice(float3_output.data(), count, 0);
    if (!compareVectors(float3_output.ToHost(), tracker->MOIs(), "moments of inertia"))
        return false;
    tracker->OwnerWildcardValuesToDevice("retrieval_tag", float_output.data(), count, 0);
    if (!compareVectors(float_output.ToHost(), tracker->GetOwnerWildcardValues("retrieval_tag"), "owner wildcards"))
        return false;

    TestDeviceBuffer<float3> range_output(2, 0);
    solver.GetOwnerVelocityToDevice(range_output.data(), 2, 0, tracker->GetOwnerID(1), 2);
    const auto host_velocities = tracker->Velocities();
    if (!compareVectors(range_output.ToHost(), std::vector<float3>(host_velocities.begin() + 1, host_velocities.end()),
                        "nonzero owner range"))
        return false;

    // Repeated output exercises reusable scratch/staging ownership and ensures no query leaves a claimed temp buffer.
    for (int iteration = 0; iteration < 3; iteration++) {
        tracker->PositionsToDevice(float3_output.data(), count, 0);
    }

    bool capacity_rejected = false;
    try {
        tracker->PositionsToDevice(float3_output.data(), count - 1, 0);
    } catch (const SolverException&) {
        capacity_rejected = true;
    }
    if (!capacity_rejected) {
        std::cerr << "FAIL: insufficient fixed-output capacity was not rejected." << std::endl;
        return false;
    }

    bool host_pointer_rejected = false;
    std::vector<float3> host_only_output(count);
    try {
        tracker->PositionsToDevice(host_only_output.data(), count, 0);
    } catch (const SolverException&) {
        host_pointer_rejected = true;
    }
    if (!host_pointer_rejected) {
        std::cerr << "FAIL: ordinary host memory was accepted as CUDA device output." << std::endl;
        return false;
    }

    if (visible_devices >= 2) {
        TestDeviceBuffer<float3> peer_output(count, 1);
        tracker->PositionsToDevice(peer_output.data(), count, 1);
        if (!compareVectors(peer_output.ToHost(), tracker->Positions(), "cross-device positions"))
            return false;

        bool wrong_device_rejected = false;
        try {
            tracker->PositionsToDevice(peer_output.data(), count, 0);
        } catch (const SolverException&) {
            wrong_device_rejected = true;
        }
        if (!wrong_device_rejected) {
            std::cerr << "FAIL: mismatched pointer device declaration was not rejected." << std::endl;
            return false;
        }
    }

    const size_t contact_capacity = solver.GetNumContacts();
    if (contact_capacity == 0) {
        std::cerr << "FAIL: contact-force retrieval scenario produced no contacts." << std::endl;
        return false;
    }
    TestDeviceBuffer<float3> contact_points(contact_capacity, 0);
    TestDeviceBuffer<float3> contact_forces(contact_capacity, 0);
    TestDeviceBuffer<float3> contact_torques(contact_capacity, 0);

    std::vector<float3> host_points;
    std::vector<float3> host_forces;
    std::vector<float3> host_torques;
    const size_t host_count = tracker->GetContactForcesAndGlobalTorqueForAll(host_points, host_forces, host_torques);
    const size_t device_count = tracker->GetContactForcesAndGlobalTorqueForAllToDevice(
        contact_points.data(), contact_forces.data(), contact_torques.data(), contact_capacity, 0);
    if (device_count != host_count ||
        !compareVectors(contact_points.ToHost(device_count), host_points, "contact points") ||
        !compareVectors(contact_forces.ToHost(device_count), host_forces, "contact forces") ||
        !compareVectors(contact_torques.ToHost(device_count), host_torques, "global contact torques")) {
        return false;
    }

    const size_t force_only_count =
        tracker->GetContactForcesForAllToDevice(contact_points.data(), contact_forces.data(), contact_capacity, 0);
    if (force_only_count != host_count ||
        !compareVectors(contact_forces.ToHost(force_only_count), host_forces, "force-only contact retrieval")) {
        return false;
    }

    const size_t local_count = tracker->GetContactForcesAndLocalTorqueForAll(host_points, host_forces, host_torques);
    const size_t local_device_count = tracker->GetContactForcesAndLocalTorqueForAllToDevice(
        contact_points.data(), contact_forces.data(), contact_torques.data(), contact_capacity, 0);
    if (local_device_count != local_count ||
        !compareVectors(contact_torques.ToHost(local_device_count), host_torques, "local contact torques")) {
        return false;
    }

    return true;
}

bool testNonJitifiedMassData() {
    DEMSolver solver(1);
    solver.SetVerbosity("ERROR");
    solver.DisableJitifyMassProperties();
    solver.InstructBoxDomainDimension(4, 4, 4);
    auto material = solver.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.2}, {"mu", 0.4}, {"Crr", 0.0}});
    auto sphere = solver.LoadSphereType(3.25f, 0.4f, material);
    auto batch = solver.AddClumps(sphere, std::vector<float3>{make_float3(1, 1, 1), make_float3(3, 1, 1)});
    auto tracker = solver.Track(batch);
    solver.SetInitTimeStep(1e-5);
    solver.Initialize();

    TestDeviceBuffer<float> masses(2, 0);
    TestDeviceBuffer<float3> moments(2, 0);
    tracker->MassesToDevice(masses.data(), 2, 0);
    tracker->MOIsToDevice(moments.data(), 2, 0);
    return compareVectors(masses.ToHost(), tracker->Masses(), "non-jitified masses") &&
           compareVectors(moments.ToHost(), tracker->MOIs(), "non-jitified moments of inertia");
}

}  // namespace

int main() {
    int visible_devices = 0;
    const cudaError_t status = cudaGetDeviceCount(&visible_devices);
    if (status != cudaSuccess) {
        std::cerr << "SKIP: CUDA runtime unavailable: " << cudaGetErrorString(status) << std::endl;
        return 0;
    }
    if (visible_devices == 0) {
        std::cout << "SKIP: DEMTest_DeviceDataRetrieval requires a visible CUDA device." << std::endl;
        return 0;
    }

    if (!testFixedOwnerData(visible_devices)) {
        return 1;
    }
    if (!testNonJitifiedMassData()) {
        return 1;
    }
    std::cout << "PASS: tracker device-data retrieval matches host APIs." << std::endl;
    return 0;
}
