//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#ifndef DEME_DYNAMIC_THREAD_HELPERS_HPP
#define DEME_DYNAMIC_THREAD_HELPERS_HPP

#include <cstddef>
#include <vector>

#include <DEM/Defines.h>
#include <DEM/VariableTypes.h>

namespace deme {

// Host-side snapshot of dT's old contact arrays, used only for DEBUG-verbosity lost-contact diagnostics.
struct LostContactDebugSnapshot {
    size_t nPatchContacts = 0;
    size_t nPrimitiveContacts = 0;
    std::vector<bodyID_t> idPatchA;
    std::vector<bodyID_t> idPatchB;
    std::vector<contact_t> contactTypePatch;
    std::vector<bodyID_t> contactPatchIsland;
    std::vector<bodyID_t> idPrimitiveA;
    std::vector<bodyID_t> idPrimitiveB;
    std::vector<contact_t> contactTypePrimitive;
    std::vector<contactPairs_t> geomToPatchMap;
    std::vector<float3> primitivePenetrationStorage;
    std::vector<float3> primitiveAreaStorage;
};

}  // namespace deme

#endif
