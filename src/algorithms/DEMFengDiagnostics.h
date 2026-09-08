// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_ALGORITHMS_FENG_DIAGNOSTICS_H
#define DEME_ALGORITHMS_FENG_DIAGNOSTICS_H
#include "DEM/Structs.h"
#include "DEM/utils/FengDiagnostics.h"

namespace deme {
// Read-only side path over the existing primitive candidates and patch keys. No force or history arrays are modified.
void computeMeshMeshFengDiagnostics(DEMSimParams* params,
                                    DEMDataDT* data,
                                    contactPairs_t* patchKeys,
                                    contactPairs_t primitiveStart,
                                    contactPairs_t primitiveCount,
                                    contactPairs_t patchStart,
                                    contactPairs_t patchCount,
                                    const double* legacyAreas,
                                    const float3* legacyNormals,
                                    const double* legacyPenetrations,
                                    const double3* legacyPoints,
                                    MeshMeshFengDiagnostic* output,
                                    cudaStream_t& stream,
                                    DEMSolverScratchData& scratch);
}  // namespace deme
#endif
