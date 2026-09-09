// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_FENG_DIAGNOSTICS_H
#define DEME_FENG_DIAGNOSTICS_H

#include "DEM/VariableTypes.h"
#include <cuda_runtime.h>

namespace deme {

// One snapshot per current mesh-mesh patch contact, evaluated at the last force calculation (before integration).
// The original line gates remain diagnostic. Force selection additionally validates solid shape, grouping and a
// single boundary cycle. Patch indices are transient, not persistent history IDs.
enum class FengFallback : unsigned int {
    NONE = 0,
    INACTIVE_LEGACY = 1,
    UNSUPPORTED_SOLID = 2,
    UNSUPPORTED_GROUPING = 3,
    DIAGNOSTIC_GATE = 4,
    BOUNDARY = 5,
    NORMAL = 6
};
struct MeshMeshFengDiagnostic {
    contactPairs_t patchContact;
    bodyID_t ownerA, ownerB;
    double3 referenceOrigin;
    double3 vectorArea, geometricMoment, boundaryResidual;
    double boundaryLength, segmentCount, ambiguousPairCount;
    double legacyArea, legacyPenetration;
    float3 legacyNormal;
    double3 legacyContactPoint;
    double fengArea;
    double3 fengNormal, fengContactPoint;
    bool ownersWatertight;
    bool closurePassed;
    bool hasContactLine;
    bool ownersValidated;
    bool boundaryValidated;
    bool fengEligible;
    bool usedFeng;
    FengFallback fallbackReason;
};

}  // namespace deme
#endif
