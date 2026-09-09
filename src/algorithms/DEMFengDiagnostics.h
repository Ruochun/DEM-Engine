// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_ALGORITHMS_FENG_DIAGNOSTICS_H
#define DEME_ALGORITHMS_FENG_DIAGNOSTICS_H
#include "DEM/Structs.h"
#include "DEM/utils/FengDiagnostics.h"

namespace deme {
// Capture boundary/legacy diagnostics and optionally select validated Feng geometry. Penetration/history are unchanged.
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
                                    DEMSolverScratchData& scratch);
}  // namespace deme
#endif
