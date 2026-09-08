// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_FENG_DIAGNOSTICS_H
#define DEME_FENG_DIAGNOSTICS_H

#include "DEM/VariableTypes.h"
#include <cuda_runtime.h>

namespace deme {

// One snapshot per current mesh-mesh patch contact, evaluated at the last force calculation (before integration).
// These are diagnostic candidates, never force inputs. Closure is necessary but does NOT prove loop completeness,
// consistent winding, or absence of duplicate segments. Patch indices are transient, not persistent history IDs.
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
};

}  // namespace deme
#endif
