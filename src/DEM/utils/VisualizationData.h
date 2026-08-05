// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DEME_VISUALIZATION_DATA_H
#define DEME_VISUALIZATION_DATA_H

#include <cuda_runtime.h>

#include <vector>

#include "../VariableTypes.h"

namespace deme {

/// One sphere component in a host-side visualization snapshot.
struct DEMVisualizationSphere {
    float3 position;
    float radius;
    family_t family;
    bodyID_t owner;
};

/// One mesh facet in a host-side visualization snapshot.
struct DEMVisualizationTriangle {
    float3 a;
    float3 b;
    float3 c;
    family_t family;
    bodyID_t owner;
};

/// Host-side geometry captured at one solver time for rendering or external visualization.
struct DEMVisualizationSnapshot {
    double simulation_time = 0.0;
    std::vector<DEMVisualizationSphere> spheres;
    std::vector<DEMVisualizationTriangle> triangles;
};

}  // namespace deme

#endif
