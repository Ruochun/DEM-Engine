// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#include "algorithms/DEMFengDiagnostics.h"
#include "algorithms/DEMStaticDeviceSubroutines.h"
#include "DEM/utils/FengGeometry.cuh"
#include "kernel/DEMHelperKernels.cuh"

namespace deme {
namespace {

// Reconstruct owner pose once per triangle; rotate local vertices before adding the relative owner displacement.
__device__ double3 ownerPosition(DEMSimParams* params, DEMDataDT* data, bodyID_t owner) {
    double3 p;
    voxelIDToPosition<double, voxelID_t, subVoxelPos_t>(p.x, p.y, p.z, data->voxelID[owner], data->locX[owner],
                                                        data->locY[owner], data->locZ[owner], params->nvXp2,
                                                        params->nvYp2, params->voxelSize, params->l);
    return p + make_double3(params->LBFX, params->LBFY, params->LBFZ);
}

__device__ void triangleVertices(DEMDataDT* data, bodyID_t tri, const double3& displacement, double3* v) {
    const bodyID_t owner = data->ownerTriMesh[tri];
    const float4 q = make_float4(data->oriQx[owner], data->oriQy[owner], data->oriQz[owner], data->oriQw[owner]);
    v[0] = to_double3(data->relPosNode1[tri]);
    v[1] = to_double3(data->relPosNode2[tri]);
    v[2] = to_double3(data->relPosNode3[tri]);
    for (int i = 0; i < 3; ++i) {
        applyOriQToVector3(v[i], q);
        v[i] += displacement;
    }
}

// Include every broad-phase candidate, even if legacy projection rejected it. Filtering by legacy area/depth
// would discard transverse segments and can open an otherwise closed boundary.
__global__ void contributions(DEMSimParams* params,
                              DEMDataDT* data,
                              contactPairs_t start,
                              contactPairs_t count,
                              double3* s,
                              double3* g,
                              double3* residual,
                              double3* statistics,
                              double3* starts,
                              double3* ends) {
    const contactPairs_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count)
        return;
    s[i] = g[i] = residual[i] = statistics[i] = make_double3(0.0, 0.0, 0.0);
    starts[i] = ends[i] = make_double3(0.0, 0.0, 0.0);
    const bodyID_t ta = data->idPrimitiveA[start + i], tb = data->idPrimitiveB[start + i];
    const double3 origin = ownerPosition(params, data, data->ownerTriMesh[ta]);
    const double3 offset = ownerPosition(params, data, data->ownerTriMesh[tb]) - origin;
    double3 a[3], b[3], p, q;
    triangleVertices(data, ta, make_double3(0.0, 0.0, 0.0), a);
    triangleVertices(data, tb, offset, b);
    const feng::Intersection result = feng::intersectionSegment(a, b, p, q);
    if (result == feng::Intersection::SEGMENT) {
        starts[i] = p;
        ends[i] = q;
        feng::contribution(p, q, s[i], g[i]);
        residual[i] = q - p;
        statistics[i] = make_double3(length(q - p), 1.0, 0.0);
    } else if (result == feng::Intersection::AMBIGUOUS) {
        statistics[i].z = 1.0;
    }
}

// A group's first primitive supplies its owner pair and common origin. Reduction ordering matches legacy patch
// voting: geomToPatchMap keys are contiguous and global patch indices; reduced values are relative to this contact-type
// batch.
__global__ void finalize(DEMSimParams* params,
                         DEMDataDT* data,
                         contactPairs_t start,
                         contactPairs_t count,
                         contactPairs_t patchStart,
                         const double3* s,
                         const double3* g,
                         const double3* residual,
                         const double3* statistics,
                         double* areas,
                         float3* normals,
                         const double* penetrations,
                         double3* points,
                         MeshMeshFengDiagnostic* output,
                         const double3* starts,
                         const double3* ends,
                         bool useFeng,
                         bool simpleGrouping) {
    const contactPairs_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count)
        return;
    const contactPairs_t key = data->geomToPatchMap[start + i];
    if (i > 0 && data->geomToPatchMap[start + i - 1] == key)
        return;
    const contactPairs_t patch = key - patchStart;
    MeshMeshFengDiagnostic r{};
    r.patchContact = key;
    r.ownerA = data->ownerTriMesh[data->idPrimitiveA[start + i]];
    r.ownerB = data->ownerTriMesh[data->idPrimitiveB[start + i]];
    r.referenceOrigin = ownerPosition(params, data, r.ownerA);
    r.vectorArea = s[patch];
    r.geometricMoment = g[patch];
    r.boundaryResidual = residual[patch];
    r.boundaryLength = statistics[patch].x;
    r.segmentCount = statistics[patch].y;
    r.ambiguousPairCount = statistics[patch].z;
    r.legacyArea = areas[patch];
    r.legacyNormal = normals[patch];
    r.legacyPenetration = penetrations[patch];
    r.legacyContactPoint = points[patch];
    r.ownersWatertight =
        data->ownerMeshWatertight && data->ownerMeshWatertight[r.ownerA] && data->ownerMeshWatertight[r.ownerB];
    r.closurePassed =
        r.segmentCount >= 3.0 && r.boundaryLength > 0.0 && length(r.boundaryResidual) <= 1e-8 * r.boundaryLength;
    // This is a diagnostic gate, not a proof of topology. In particular, cancelling open chains can pass closure.
    if (r.closurePassed && r.ambiguousPairCount == 0.0 && r.ownersWatertight &&
        length(r.vectorArea) > 1e-12 * r.boundaryLength * r.boundaryLength) {
        double3 localPoint;
        r.hasContactLine = feng::contactLine(r.vectorArea, r.geometricMoment, points[patch] - r.referenceOrigin,
                                             r.fengArea, r.fengNormal, localPoint);
        if (r.hasContactLine)
            r.fengContactPoint = localPoint + r.referenceOrigin;
    }
    // A single mesh patch on each validated solid plus simple grouping keeps all owner-pair candidates together.
    // The endpoint audit rejects incomplete chains and multiple loops; candidate completeness still relies on CD.
    r.ownersValidated = data->ownerMeshFengValidated[r.ownerA] && data->ownerMeshFengValidated[r.ownerB];
    r.fallbackReason = FengFallback::NONE;
    if (!(r.legacyArea > 0 && r.legacyPenetration > 0))
        r.fallbackReason = FengFallback::INACTIVE_LEGACY;
    else if (!r.ownersValidated)
        r.fallbackReason = FengFallback::UNSUPPORTED_SOLID;
    else if (!simpleGrouping)
        r.fallbackReason = FengFallback::UNSUPPORTED_GROUPING;
    else if (!r.hasContactLine)
        r.fallbackReason = FengFallback::DIAGNOSTIC_GATE;
    else {
        contactPairs_t end = i + 1;
        while (end < count && data->geomToPatchMap[start + end] == key)
            ++end;
        r.boundaryValidated = feng::singleBoundary(starts + i, ends + i, end - i, r.boundaryLength);
        if (!r.boundaryValidated)
            r.fallbackReason = FengFallback::BOUNDARY;
        else if (!(dot(r.fengNormal, to_double3(r.legacyNormal)) > 0.0))
            r.fallbackReason = FengFallback::NORMAL;
    }
    r.fengEligible = r.fallbackReason == FengFallback::NONE;
    r.usedFeng = useFeng && r.fengEligible;
    if (r.usedFeng) {
        // Replace all three geometry outputs together; never mix segment and projection contributions in a patch.
        // The force kernel still consumes legacy penetration and transports the same friction-history variables.
        areas[patch] = r.fengArea;
        normals[patch] = to_float3(r.fengNormal);
        points[patch] = r.fengContactPoint;
    }
    output[patch] = r;
}
}  // namespace

// Reuse existing FP64 CUB sum-by-key instantiations, with no atomics or host readback in the diagnostic pipeline.
void computeMeshMeshFengDiagnostics(DEMSimParams* params,
                                    DEMDataDT* data,
                                    contactPairs_t* patchKeys,
                                    contactPairs_t primitiveStart,
                                    contactPairs_t primitiveCount,
                                    contactPairs_t patchStart,
                                    contactPairs_t patchCount,
                                    double* legacyAreas,
                                    float3* legacyNormals,
                                    const double* legacyPenetrations,
                                    double3* legacyPoints,
                                    MeshMeshFengDiagnostic* output,
                                    bool useFeng,
                                    bool simpleGrouping,
                                    cudaStream_t& stream,
                                    DEMSolverScratchData& scratch) {
    if (!primitiveCount || !patchCount)
        return;
    constexpr unsigned int FENG_DIAGNOSTIC_BLOCK = 128;
    const size_t blocks = (static_cast<size_t>(primitiveCount) + FENG_DIAGNOSTIC_BLOCK - 1) / FENG_DIAGNOSTIC_BLOCK;
    const char* names[4] = {"fengS", "fengG", "fengResidual", "fengStatistics"};
    const char* sums[4] = {"fengSumS", "fengSumG", "fengSumResidual", "fengSumStatistics"};
    double3* values[4];
    double3* totals[4];
    for (int i = 0; i < 4; ++i) {
        values[i] = reinterpret_cast<double3*>(scratch.allocateTempVector(names[i], primitiveCount * sizeof(double3)));
        totals[i] = reinterpret_cast<double3*>(scratch.allocateTempVector(sums[i], patchCount * sizeof(double3)));
    }
    auto* starts =
        reinterpret_cast<double3*>(scratch.allocateTempVector("fengStarts", primitiveCount * sizeof(double3)));
    auto* ends = reinterpret_cast<double3*>(scratch.allocateTempVector("fengEnds", primitiveCount * sizeof(double3)));
    contributions<<<blocks, FENG_DIAGNOSTIC_BLOCK, 0, stream>>>(params, data, primitiveStart, primitiveCount, values[0],
                                                                values[1], values[2], values[3], starts, ends);
    DEME_GPU_CALL(cudaGetLastError());
    auto* keys = reinterpret_cast<contactPairs_t*>(
        scratch.allocateTempVector("fengKeys", primitiveCount * sizeof(contactPairs_t)));
    scratch.allocateDualStruct("fengNumKeys");
    for (int i = 0; i < 4; ++i) {
        cubSumReduceByKey<contactPairs_t, double3>(patchKeys, keys, values[i], totals[i],
                                                   scratch.getDualStructDevice("fengNumKeys"), primitiveCount, stream,
                                                   scratch);
    }
    finalize<<<blocks, FENG_DIAGNOSTIC_BLOCK, 0, stream>>>(
        params, data, primitiveStart, primitiveCount, patchStart, totals[0], totals[1], totals[2], totals[3],
        legacyAreas, legacyNormals, legacyPenetrations, legacyPoints, output, starts, ends, useFeng, simpleGrouping);
    DEME_GPU_CALL(cudaGetLastError());
    for (int i = 0; i < 4; ++i) {
        scratch.finishUsingTempVector(names[i]);
        scratch.finishUsingTempVector(sums[i]);
    }
    scratch.finishUsingTempVector("fengStarts");
    scratch.finishUsingTempVector("fengEnds");
    scratch.finishUsingTempVector("fengKeys");
    scratch.finishUsingDualStruct("fengNumKeys");
    DEME_GPU_DEBUG_SYNC(stream);
}
}  // namespace deme
