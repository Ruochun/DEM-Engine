// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DeviceDataTransfer.hpp"

#include <algorithm>

namespace deme {
namespace device_data {

namespace {

class DeviceRestorer {
  public:
    cudaError_t Capture() { return cudaGetDevice(&m_device); }
    ~DeviceRestorer() {
        if (m_device >= 0) {
            (void)cudaSetDevice(m_device);
        }
    }

  private:
    int m_device = -1;
};

}  // namespace

cudaError_t ValidateOutputPointer(const void* pointer, size_t bytes, int device) {
    if (bytes == 0) {
        return cudaSuccess;
    }
    if (!pointer || device < 0) {
        return cudaErrorInvalidValue;
    }

    int device_count = 0;
    cudaError_t status = cudaGetDeviceCount(&device_count);
    if (status != cudaSuccess) {
        return status;
    }
    if (device >= device_count) {
        return cudaErrorInvalidDevice;
    }

    cudaPointerAttributes attributes{};
    status = cudaPointerGetAttributes(&attributes, pointer);
    if (status != cudaSuccess) {
        return status;
    }
    if (attributes.type == cudaMemoryTypeDevice && attributes.device != device) {
        return cudaErrorInvalidDevicePointer;
    }
    if (attributes.type != cudaMemoryTypeDevice && attributes.type != cudaMemoryTypeManaged) {
        return cudaErrorInvalidDevicePointer;
    }
    return cudaSuccess;
}

TransferBuffer::~TransferBuffer() {
    if (m_staging) {
        (void)cudaFreeHost(m_staging);
    }
}

cudaError_t TransferBuffer::EnsureStagingCapacity(size_t bytes) {
    if (bytes <= m_staging_capacity) {
        return cudaSuccess;
    }

    void* replacement = nullptr;
    cudaError_t status = cudaMallocHost(&replacement, bytes);
    if (status != cudaSuccess) {
        return status;
    }
    if (m_staging) {
        status = cudaFreeHost(m_staging);
        if (status != cudaSuccess) {
            (void)cudaFreeHost(replacement);
            return status;
        }
    }
    m_staging = replacement;
    m_staging_capacity = bytes;
    return cudaSuccess;
}

cudaError_t TransferBuffer::Copy(void* destination,
                                 int destination_device,
                                 const void* source,
                                 int source_device,
                                 size_t bytes) {
    if (bytes == 0) {
        return cudaSuccess;
    }
    if (!source) {
        return cudaErrorInvalidValue;
    }
    cudaError_t status = ValidateOutputPointer(destination, bytes, destination_device);
    if (status != cudaSuccess) {
        return status;
    }

    std::lock_guard<std::mutex> lock(m_mutex);
    DeviceRestorer restore;
    status = restore.Capture();
    if (status != cudaSuccess) {
        return status;
    }

    if (source_device == destination_device) {
        status = cudaSetDevice(source_device);
        if (status != cudaSuccess) {
            return status;
        }
        return cudaMemcpy(destination, source, bytes, cudaMemcpyDeviceToDevice);
    }

    int peer_access = 0;
    status = cudaDeviceCanAccessPeer(&peer_access, destination_device, source_device);
    if (status != cudaSuccess) {
        return status;
    }
    if (peer_access) {
        return cudaMemcpyPeer(destination, destination_device, source, source_device, bytes);
    }

    status = EnsureStagingCapacity(bytes);
    if (status != cudaSuccess) {
        return status;
    }
    status = cudaSetDevice(source_device);
    if (status != cudaSuccess) {
        return status;
    }
    status = cudaMemcpy(m_staging, source, bytes, cudaMemcpyDeviceToHost);
    if (status != cudaSuccess) {
        return status;
    }
    status = cudaSetDevice(destination_device);
    if (status != cudaSuccess) {
        return status;
    }
    return cudaMemcpy(destination, m_staging, bytes, cudaMemcpyHostToDevice);
}

}  // namespace device_data
}  // namespace deme
