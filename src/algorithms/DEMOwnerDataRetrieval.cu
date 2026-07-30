// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DEMOwnerDataRetrieval.h"

#include "DEM/Defines.h"
#include "kernel/DEMHelperKernels.cuh"

namespace deme {

namespace {

constexpr unsigned int OWNER_DATA_RETRIEVAL_BLOCK_SIZE = 256;

__global__ void PackOwnerDataKernel(void* output,
                                    OwnerDataField field,
                                    bodyID_t owner_begin,
                                    size_t count,
                                    const DEMSimParams* sim_params,
                                    const DEMDataDT* data,
                                    bool mass_properties_are_jitified,
                                    unsigned int wildcard_index) {
    const size_t output_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (output_index >= count) {
        return;
    }
    const bodyID_t owner = owner_begin + output_index;

    const float4 orientation =
        make_float4(data->oriQx[owner], data->oriQy[owner], data->oriQz[owner], data->oriQw[owner]);
    switch (field) {
        case OwnerDataField::POSITION: {
            double x, y, z;
            voxelIDToPosition<double, voxelID_t, subVoxelPos_t>(
                x, y, z, data->voxelID[owner], data->locX[owner], data->locY[owner], data->locZ[owner],
                sim_params->nvXp2, sim_params->nvYp2, sim_params->voxelSize, sim_params->l);
            static_cast<float3*>(output)[output_index] =
                make_float3(x + sim_params->LBFX, y + sim_params->LBFY, z + sim_params->LBFZ);
            break;
        }
        case OwnerDataField::VELOCITY:
            static_cast<float3*>(output)[output_index] = make_float3(data->vX[owner], data->vY[owner], data->vZ[owner]);
            break;
        case OwnerDataField::ANGULAR_VELOCITY_LOCAL:
            static_cast<float3*>(output)[output_index] =
                make_float3(data->omgBarX[owner], data->omgBarY[owner], data->omgBarZ[owner]);
            break;
        case OwnerDataField::ANGULAR_VELOCITY_GLOBAL: {
            float3 value = make_float3(data->omgBarX[owner], data->omgBarY[owner], data->omgBarZ[owner]);
            applyOriQToVector3(value, orientation);
            static_cast<float3*>(output)[output_index] = value;
            break;
        }
        case OwnerDataField::ORIENTATION:
            static_cast<float4*>(output)[output_index] = orientation;
            break;
        case OwnerDataField::CONTACT_ACCELERATION:
            static_cast<float3*>(output)[output_index] = make_float3(data->aX[owner], data->aY[owner], data->aZ[owner]);
            break;
        case OwnerDataField::CONTACT_ANGULAR_ACCELERATION_LOCAL:
            static_cast<float3*>(output)[output_index] =
                make_float3(data->alphaX[owner], data->alphaY[owner], data->alphaZ[owner]);
            break;
        case OwnerDataField::CONTACT_ANGULAR_ACCELERATION_GLOBAL: {
            float3 value = make_float3(data->alphaX[owner], data->alphaY[owner], data->alphaZ[owner]);
            applyOriQToVector3(value, orientation);
            static_cast<float3*>(output)[output_index] = value;
            break;
        }
        case OwnerDataField::FAMILY:
            static_cast<unsigned int*>(output)[output_index] = static_cast<unsigned int>(+(data->familyID[owner]));
            break;
        case OwnerDataField::MASS: {
            const inertiaOffset_t offset = mass_properties_are_jitified ? data->inertiaPropOffsets[owner] : owner;
            static_cast<float*>(output)[output_index] = data->massOwnerBody[offset];
            break;
        }
        case OwnerDataField::MOI: {
            const inertiaOffset_t offset = mass_properties_are_jitified ? data->inertiaPropOffsets[owner] : owner;
            static_cast<float3*>(output)[output_index] =
                make_float3(data->mmiXX[offset], data->mmiYY[offset], data->mmiZZ[offset]);
            break;
        }
        case OwnerDataField::WILDCARD:
            static_cast<float*>(output)[output_index] = data->ownerWildcards[wildcard_index][owner];
            break;
    }
}

}  // namespace

void PackOwnerData(void* output,
                   OwnerDataField field,
                   bodyID_t owner_begin,
                   size_t count,
                   const DEMSimParams* sim_params,
                   const DEMDataDT* data,
                   bool mass_properties_are_jitified,
                   unsigned int wildcard_index,
                   cudaStream_t stream) {
    if (count == 0) {
        return;
    }
    const size_t blocks = (count + OWNER_DATA_RETRIEVAL_BLOCK_SIZE - 1) / OWNER_DATA_RETRIEVAL_BLOCK_SIZE;
    PackOwnerDataKernel<<<blocks, OWNER_DATA_RETRIEVAL_BLOCK_SIZE, 0, stream>>>(
        output, field, owner_begin, count, sim_params, data, mass_properties_are_jitified, wildcard_index);
    DEME_GPU_DEBUG_SYNC(stream);
}

size_t OwnerDataElementSize(OwnerDataField field) {
    switch (field) {
        case OwnerDataField::ORIENTATION:
            return sizeof(float4);
        case OwnerDataField::FAMILY:
            return sizeof(unsigned int);
        case OwnerDataField::MASS:
        case OwnerDataField::WILDCARD:
            return sizeof(float);
        default:
            return sizeof(float3);
    }
}

}  // namespace deme
