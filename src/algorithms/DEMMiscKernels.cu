//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <algorithms/DEMStaticDeviceSubroutines.h>
#include <algorithms/DEMStaticDeviceUtilities.cuh>

#include <kernel/DEMHelperKernels.cuh>

namespace deme {

////////////////////////////////////////////////////////////////////////////////
// Misc kernels implementations
////////////////////////////////////////////////////////////////////////////////

__global__ void fillMarginValues_impl(DEMSimParams* simParams,
                                      DEMDataKT* granData,
                                      float* marginSizeArr,
                                      bodyID_t* ownerIDArr,
                                      size_t n) {
    size_t ID = blockIdx.x * blockDim.x + threadIdx.x;
    if (ID < n) {
        bodyID_t ownerID = ownerIDArr[ID];
        unsigned int my_family = granData->familyID[ownerID];
        marginSizeArr[ID] = simParams->dyn.beta + granData->familyExtraMarginSize[my_family];
    }
}

void fillMarginValues(DEMSimParams* simParams,
                      DEMDataKT* granData,
                      size_t nSphere,
                      size_t nTri,
                      size_t nAnal,
                      cudaStream_t& this_stream) {
    size_t blocks_needed = (nSphere + DEME_MAX_THREADS_PER_BLOCK - 1) / DEME_MAX_THREADS_PER_BLOCK;
    if (blocks_needed > 0) {
        fillMarginValues_impl<<<blocks_needed, DEME_MAX_THREADS_PER_BLOCK, 0, this_stream>>>(
            simParams, granData, granData->marginSizeSphere, granData->ownerClumpBody, nSphere);
        DEME_GPU_DEBUG_SYNC(this_stream);
    }
    blocks_needed = (nTri + DEME_MAX_THREADS_PER_BLOCK - 1) / DEME_MAX_THREADS_PER_BLOCK;
    if (blocks_needed > 0) {
        fillMarginValues_impl<<<blocks_needed, DEME_MAX_THREADS_PER_BLOCK, 0, this_stream>>>(
            simParams, granData, granData->marginSizeTriangle, granData->ownerTriMesh, nTri);
        DEME_GPU_DEBUG_SYNC(this_stream);
    }
    blocks_needed = (nAnal + DEME_MAX_THREADS_PER_BLOCK - 1) / DEME_MAX_THREADS_PER_BLOCK;
    if (blocks_needed > 0) {
        fillMarginValues_impl<<<blocks_needed, DEME_MAX_THREADS_PER_BLOCK, 0, this_stream>>>(
            simParams, granData, granData->marginSizeAnalytical, granData->ownerAnalBody, nAnal);
        DEME_GPU_DEBUG_SYNC(this_stream);
    }
}

}  // namespace deme
