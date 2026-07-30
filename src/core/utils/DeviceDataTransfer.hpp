// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DEME_DEVICE_DATA_TRANSFER_HPP
#define DEME_DEVICE_DATA_TRANSFER_HPP

#include <cuda_runtime_api.h>

#include <cstddef>
#include <mutex>

namespace deme {
namespace device_data {

/// Validate that a caller-provided output range is backed by CUDA-accessible memory on the declared device.
cudaError_t ValidateOutputPointer(const void* pointer, size_t bytes, int device);

/// Reusable synchronous device-transfer helper. Same-device and peer-capable transfers stay on GPU; systems without
/// peer access use a reusable pinned-host staging allocation.
///
/// This utility deliberately depends only on the CUDA runtime and the C++ standard library so it can be moved to a
/// standalone interoperability library without carrying DEM-Engine types or ownership rules with it.
class TransferBuffer {
  public:
    TransferBuffer() = default;
    ~TransferBuffer();

    TransferBuffer(const TransferBuffer&) = delete;
    TransferBuffer& operator=(const TransferBuffer&) = delete;

    cudaError_t Copy(void* destination, int destination_device, const void* source, int source_device, size_t bytes);

  private:
    cudaError_t EnsureStagingCapacity(size_t bytes);

    void* m_staging = nullptr;
    size_t m_staging_capacity = 0;
    std::mutex m_mutex;
};

}  // namespace device_data
}  // namespace deme

#endif
