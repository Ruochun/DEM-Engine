//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//  All rights reserved.

#include <cstring>
#include <iostream>
#include <thread>

#include <core/ApiVersion.h>
#include <core/utils/Macros.h>
#include <core/utils/chpf/particle_writer.hpp>
#include <granular/GranularDefines.h>
#include <granular/dT.h>
#include <granular/kT.h>
#include <core/utils/JitHelper.h>
#include <granular/HostSideHelpers.cpp>
#include <helper_math.cuh>

#include <algorithms/DEMCubHelperFunctions.h>

namespace sgps {

// Put sim data array pointers in place
void DEMDynamicThread::packDataPointers() {
    granData->inertiaPropOffsets = inertiaPropOffsets.data();
    granData->familyID = familyID.data();
    granData->voxelID = voxelID.data();
    granData->locX = locX.data();
    granData->locY = locY.data();
    granData->locZ = locZ.data();
    granData->h2aX = h2aX.data();
    granData->h2aY = h2aY.data();
    granData->h2aZ = h2aZ.data();
    granData->hvX = hvX.data();
    granData->hvY = hvY.data();
    granData->hvZ = hvZ.data();
    granData->oriQ0 = oriQ0.data();
    granData->oriQ1 = oriQ1.data();
    granData->oriQ2 = oriQ2.data();
    granData->oriQ3 = oriQ3.data();
    granData->hOmgBarX = hOmgBarX.data();
    granData->hOmgBarY = hOmgBarY.data();
    granData->hOmgBarZ = hOmgBarZ.data();
    granData->h2AlphaX = h2AlphaX.data();
    granData->h2AlphaY = h2AlphaY.data();
    granData->h2AlphaZ = h2AlphaZ.data();
    granData->idGeometryA = idGeometryA.data();
    granData->idGeometryB = idGeometryB.data();
    granData->idGeometryA_buffer = idGeometryA_buffer.data();
    granData->idGeometryB_buffer = idGeometryB_buffer.data();
    granData->contactForces = contactForces.data();
    granData->contactPointGeometryA = contactPointGeometryA.data();
    granData->contactPointGeometryB = contactPointGeometryB.data();

    // The offset info that indexes into the template arrays
    granData->ownerClumpBody = ownerClumpBody.data();
    granData->clumpComponentOffset = clumpComponentOffset.data();
    granData->materialTupleOffset = materialTupleOffset.data();

    // Template array pointers, which will be removed after JIT is fully functional
    granTemplates->radiiSphere = radiiSphere.data();
    granTemplates->relPosSphereX = relPosSphereX.data();
    granTemplates->relPosSphereY = relPosSphereY.data();
    granTemplates->relPosSphereZ = relPosSphereZ.data();
    granTemplates->massClumpBody = massClumpBody.data();
    granTemplates->mmiXX = mmiXX.data();
    granTemplates->mmiYY = mmiYY.data();
    granTemplates->mmiZZ = mmiZZ.data();
    granTemplates->EProxy = EProxy.data();
    granTemplates->GProxy = GProxy.data();
    granTemplates->CoRProxy = CoRProxy.data();
}
void DEMDynamicThread::packTransferPointers(DEMKinematicThread* kT) {
    // These are the pointers for sending data to dT
    granData->pKTOwnedBuffer_voxelID = kT->granData->voxelID_buffer;
    granData->pKTOwnedBuffer_locX = kT->granData->locX_buffer;
    granData->pKTOwnedBuffer_locY = kT->granData->locY_buffer;
    granData->pKTOwnedBuffer_locZ = kT->granData->locZ_buffer;
    granData->pKTOwnedBuffer_oriQ0 = kT->granData->oriQ0_buffer;
    granData->pKTOwnedBuffer_oriQ1 = kT->granData->oriQ1_buffer;
    granData->pKTOwnedBuffer_oriQ2 = kT->granData->oriQ2_buffer;
    granData->pKTOwnedBuffer_oriQ3 = kT->granData->oriQ3_buffer;
}

void DEMDynamicThread::setSimParams(unsigned char nvXp2,
                                    unsigned char nvYp2,
                                    unsigned char nvZp2,
                                    float l,
                                    double voxelSize,
                                    double binSize,
                                    float3 LBFPoint,
                                    float3 G,
                                    double ts_size,
                                    float expand_factor) {
    simParams->nvXp2 = nvXp2;
    simParams->nvYp2 = nvYp2;
    simParams->nvZp2 = nvZp2;
    simParams->l = l;
    simParams->voxelSize = voxelSize;
    simParams->binSize = binSize;
    simParams->LBFX = LBFPoint.x;
    simParams->LBFY = LBFPoint.y;
    simParams->LBFZ = LBFPoint.z;
    simParams->Gx = G.x;
    simParams->Gy = G.y;
    simParams->Gz = G.z;
    simParams->h = ts_size;
    simParams->beta = expand_factor;
    // Figure out how many bins there are in each direction
    simParams->nbX = (binID_t)(voxelSize * (double)((size_t)1 << nvXp2) / binSize) + 1;
    simParams->nbY = (binID_t)(voxelSize * (double)((size_t)1 << nvYp2) / binSize) + 1;
    simParams->nbZ = (binID_t)(voxelSize * (double)((size_t)1 << nvZp2) / binSize) + 1;
}

void DEMDynamicThread::allocateManagedArrays(size_t nClumpBodies,
                                             size_t nSpheresGM,
                                             unsigned int nClumpTopo,
                                             unsigned int nClumpComponents,
                                             unsigned int nMatTuples) {
    // Sizes of these arrays
    simParams->nSpheresGM = nSpheresGM;
    simParams->nClumpBodies = nClumpBodies;
    simParams->nDistinctClumpBodyTopologies = nClumpTopo;
    simParams->nDistinctClumpComponents = nClumpComponents;
    simParams->nMatTuples = nMatTuples;

    // Resize to the number of clumps
    TRACKED_VECTOR_RESIZE(inertiaPropOffsets, nClumpBodies, "inertiaPropOffsets", 0);
    TRACKED_VECTOR_RESIZE(familyID, nClumpBodies, "familyID", 0);
    TRACKED_VECTOR_RESIZE(voxelID, nClumpBodies, "voxelID", 0);
    TRACKED_VECTOR_RESIZE(locX, nClumpBodies, "locX", 0);
    TRACKED_VECTOR_RESIZE(locY, nClumpBodies, "locY", 0);
    TRACKED_VECTOR_RESIZE(locZ, nClumpBodies, "locZ", 0);
    TRACKED_VECTOR_RESIZE(oriQ0, nClumpBodies, "oriQ0", 1);
    TRACKED_VECTOR_RESIZE(oriQ1, nClumpBodies, "oriQ1", 0);
    TRACKED_VECTOR_RESIZE(oriQ2, nClumpBodies, "oriQ2", 0);
    TRACKED_VECTOR_RESIZE(oriQ3, nClumpBodies, "oriQ3", 0);
    TRACKED_VECTOR_RESIZE(hvX, nClumpBodies, "hvX", 0);
    TRACKED_VECTOR_RESIZE(hvY, nClumpBodies, "hvY", 0);
    TRACKED_VECTOR_RESIZE(hvZ, nClumpBodies, "hvZ", 0);
    TRACKED_VECTOR_RESIZE(hOmgBarX, nClumpBodies, "hOmgBarX", 0);
    TRACKED_VECTOR_RESIZE(hOmgBarY, nClumpBodies, "hOmgBarY", 0);
    TRACKED_VECTOR_RESIZE(hOmgBarZ, nClumpBodies, "hOmgBarZ", 0);
    TRACKED_VECTOR_RESIZE(h2aX, nClumpBodies, "h2aX", 0);
    TRACKED_VECTOR_RESIZE(h2aY, nClumpBodies, "h2aY", 0);
    TRACKED_VECTOR_RESIZE(h2aZ, nClumpBodies, "h2aZ", 0);
    TRACKED_VECTOR_RESIZE(h2AlphaX, nClumpBodies, "h2AlphaX", 0);
    TRACKED_VECTOR_RESIZE(h2AlphaY, nClumpBodies, "h2AlphaY", 0);
    TRACKED_VECTOR_RESIZE(h2AlphaZ, nClumpBodies, "h2AlphaZ", 0);

    // Resize to the number of spheres
    TRACKED_VECTOR_RESIZE(ownerClumpBody, nSpheresGM, "ownerClumpBody", 0);
    TRACKED_VECTOR_RESIZE(clumpComponentOffset, nSpheresGM, "sphereRadiusOffset", 0);
    TRACKED_VECTOR_RESIZE(materialTupleOffset, nSpheresGM, "materialTupleOffset", 0);

    // Resize to the length of the clump templates
    TRACKED_VECTOR_RESIZE(massClumpBody, nClumpTopo, "massClumpBody", 0);
    TRACKED_VECTOR_RESIZE(mmiXX, nClumpTopo, "mmiXX", 0);
    TRACKED_VECTOR_RESIZE(mmiYY, nClumpTopo, "mmiYY", 0);
    TRACKED_VECTOR_RESIZE(mmiZZ, nClumpTopo, "mmiZZ", 0);
    TRACKED_VECTOR_RESIZE(radiiSphere, nClumpComponents, "radiiSphere", 0);
    TRACKED_VECTOR_RESIZE(relPosSphereX, nClumpComponents, "relPosSphereX", 0);
    TRACKED_VECTOR_RESIZE(relPosSphereY, nClumpComponents, "relPosSphereY", 0);
    TRACKED_VECTOR_RESIZE(relPosSphereZ, nClumpComponents, "relPosSphereZ", 0);
    TRACKED_VECTOR_RESIZE(EProxy, (1 + nMatTuples) * nMatTuples / 2, "EProxy", 0);
    TRACKED_VECTOR_RESIZE(GProxy, (1 + nMatTuples) * nMatTuples / 2, "GProxy", 0);
    TRACKED_VECTOR_RESIZE(CoRProxy, (1 + nMatTuples) * nMatTuples / 2, "CoRProxy", 0);

    // Arrays for contact info
    // The lengths of contact event-based arrays are just estimates. My estimate of total contact pairs is 4n, and I
    // think the max is 6n (although I can't prove it). Note the estimate should be large enough to decrease the number
    // of reallocations in the simulation, but not too large that eats too much memory.
    TRACKED_VECTOR_RESIZE(idGeometryA, nClumpBodies * 4, "idGeometryA", 0);
    TRACKED_VECTOR_RESIZE(idGeometryB, nClumpBodies * 4, "idGeometryB", 0);
    TRACKED_VECTOR_RESIZE(contactForces, nClumpBodies * 4, "contactForces", make_float3(0));
    TRACKED_VECTOR_RESIZE(contactPointGeometryA, nClumpBodies * 4, "contactPointGeometryA", make_float3(0));
    TRACKED_VECTOR_RESIZE(contactPointGeometryB, nClumpBodies * 4, "contactPointGeometryB", make_float3(0));

    // Transfer buffer arrays
    // The following several arrays will have variable sizes, so here we only used an estimate.
    TRACKED_VECTOR_RESIZE(idGeometryA_buffer, nClumpBodies * 4, "idGeometryA_buffer", 0);
    TRACKED_VECTOR_RESIZE(idGeometryB_buffer, nClumpBodies * 4, "idGeometryB_buffer", 0);
}

void DEMDynamicThread::populateManagedArrays(const std::vector<unsigned int>& input_clump_types,
                                             const std::vector<float3>& input_clump_xyz,
                                             const std::vector<float3>& input_clump_vel,
                                             const std::vector<unsigned int>& input_clump_family,
                                             const std::vector<std::vector<unsigned int>>& input_clumps_sp_mat_ids,
                                             const std::vector<float>& clumps_mass_types,
                                             const std::vector<float3>& clumps_moi_types,
                                             const std::vector<std::vector<float>>& clumps_sp_radii_types,
                                             const std::vector<std::vector<float3>>& clumps_sp_location_types,
                                             const std::vector<float>& mat_k,
                                             const std::vector<float>& mat_g,
                                             const std::vector<float>& mat_CoR) {
    // Use some temporary hacks to get the info in the managed mem

    // First, load in material property (upper-triangle) matrix
    for (unsigned int i = 0; i < mat_k.size(); i++) {
        EProxy.at(i) = mat_k.at(i);
        GProxy.at(i) = mat_g.at(i);
        CoRProxy.at(i) = mat_CoR.at(i);
    }

    // Then load in clump mass and MOI
    // Remember this part should be quite different in the final version (due to being jitified)
    for (unsigned int i = 0; i < clumps_mass_types.size(); i++) {
        massClumpBody.at(i) = clumps_mass_types.at(i);
        float3 this_moi = clumps_moi_types.at(i);
        mmiXX.at(i) = this_moi.x;
        mmiYY.at(i) = this_moi.y;
        mmiZZ.at(i) = this_moi.z;
    }

    // Then, load in clump type info
    // Remember this part should be quite different in the final version (due to being jitified)

    size_t k = 0;
    std::vector<unsigned int> prescans;

    prescans.push_back(0);
    for (auto elem : clumps_sp_radii_types) {
        for (auto radius : elem) {
            radiiSphere.at(k) = radius;
            k++;
        }
        prescans.push_back(k);
    }
    prescans.pop_back();
    k = 0;

    for (auto elem : clumps_sp_location_types) {
        for (auto loc : elem) {
            relPosSphereX.at(k) = loc.x;
            relPosSphereY.at(k) = loc.y;
            relPosSphereZ.at(k) = loc.z;
            k++;
        }
        // std::cout << "sphere location types: " << elem.x << ", " << elem.y << ", " << elem.z << std::endl;
    }
    k = 0;

    // Then, load in input clumps
    for (size_t i = 0; i < input_clump_types.size(); i++) {
        auto type_of_this_clump = input_clump_types.at(i);
        inertiaPropOffsets.at(i) = type_of_this_clump;
        float3 LBF;
        LBF.x = simParams->LBFX;
        LBF.y = simParams->LBFY;
        LBF.z = simParams->LBFZ;
        auto this_CoM_coord = input_clump_xyz.at(i) - LBF;
        // std::cout << "CoM position: " << this_CoM_coord.x << ", " << this_CoM_coord.y << ", " << this_CoM_coord.z <<
        // std::endl;
        auto this_clump_no_sp_radii = clumps_sp_radii_types.at(type_of_this_clump);
        auto this_clump_no_sp_relPos = clumps_sp_location_types.at(type_of_this_clump);
        auto this_clump_no_sp_mat_ids = input_clumps_sp_mat_ids.at(type_of_this_clump);

        for (size_t j = 0; j < this_clump_no_sp_radii.size(); j++) {
            materialTupleOffset.at(k) = this_clump_no_sp_mat_ids.at(j);
            clumpComponentOffset.at(k) = prescans.at(type_of_this_clump) + j;
            ownerClumpBody.at(k) = i;
            k++;
            // std::cout << "Sphere Rel Pos offset: " << this_clump_no_sp_loc_offsets.at(j) << std::endl;
        }

        voxelID_t voxelNumX = (double)this_CoM_coord.x / simParams->voxelSize;
        voxelID_t voxelNumY = (double)this_CoM_coord.y / simParams->voxelSize;
        voxelID_t voxelNumZ = (double)this_CoM_coord.z / simParams->voxelSize;
        locX.at(i) = ((double)this_CoM_coord.x - (double)voxelNumX * simParams->voxelSize) / simParams->l;
        locY.at(i) = ((double)this_CoM_coord.y - (double)voxelNumY * simParams->voxelSize) / simParams->l;
        locZ.at(i) = ((double)this_CoM_coord.z - (double)voxelNumZ * simParams->voxelSize) / simParams->l;
        // std::cout << "Clump voxel num: " << voxelNumX << ", " << voxelNumY << ", " << voxelNumZ << std::endl;

        voxelID.at(i) += voxelNumX;
        voxelID.at(i) += voxelNumY << simParams->nvXp2;
        voxelID.at(i) += voxelNumZ << (simParams->nvXp2 + simParams->nvYp2);
        // std::cout << "Computed voxel num: " << voxelID.at(i) << std::endl;

        // Set initial velocity
        auto vel_of_this_clump = input_clump_vel.at(i);
        hvX.at(i) = vel_of_this_clump.x * simParams->h / simParams->l;
        hvY.at(i) = vel_of_this_clump.y * simParams->h / simParams->l;
        hvZ.at(i) = vel_of_this_clump.z * simParams->h / simParams->l;

        // Set family code
        familyID.at(i) = input_clump_family.at(i);
    }
}

void DEMDynamicThread::WriteCsvAsSpheres(std::ofstream& ptFile) const {
    ParticleFormatWriter pw;
    // pw.write(ptFile, ParticleFormatWriter::CompressionType::NONE, mass);
    std::vector<float> posX(simParams->nSpheresGM, 0);
    std::vector<float> posY(simParams->nSpheresGM, 0);
    std::vector<float> posZ(simParams->nSpheresGM, 0);
    std::vector<float> spRadii(simParams->nSpheresGM, 0);
    for (unsigned int i = 0; i < simParams->nSpheresGM; i++) {
        auto this_owner = ownerClumpBody.at(i);
        voxelID_t voxelIDX =
            voxelID.at(this_owner) & (((voxelID_t)1 << simParams->nvXp2) - 1);  // & operation here equals modulo
        voxelID_t voxelIDY = (voxelID.at(this_owner) >> simParams->nvXp2) & (((voxelID_t)1 << simParams->nvYp2) - 1);
        voxelID_t voxelIDZ = (voxelID.at(this_owner)) >> (simParams->nvXp2 + simParams->nvYp2);
        // std::cout << "this owner: " << this_owner << std::endl;
        // std::cout << "Out voxel ID: " << voxelID.at(this_owner) << std::endl;
        // std::cout << "Out voxel ID XYZ: " << voxelIDX << ", " << voxelIDY << ", " << voxelIDZ << std::endl;

        float this_sp_deviation_x = relPosSphereX.at(clumpComponentOffset.at(i));
        float this_sp_deviation_y = relPosSphereY.at(clumpComponentOffset.at(i));
        float this_sp_deviation_z = relPosSphereZ.at(clumpComponentOffset.at(i));
        float this_sp_rot_0 = oriQ0.at(this_owner);
        float this_sp_rot_1 = oriQ1.at(this_owner);
        float this_sp_rot_2 = oriQ2.at(this_owner);
        float this_sp_rot_3 = oriQ3.at(this_owner);
        hostApplyOriQ2Vector3<float, float>(this_sp_deviation_x, this_sp_deviation_y, this_sp_deviation_z,
                                            this_sp_rot_0, this_sp_rot_1, this_sp_rot_2, this_sp_rot_3);
        posX.at(i) = voxelIDX * simParams->voxelSize + locX.at(this_owner) * simParams->l + this_sp_deviation_x +
                     simParams->LBFX;
        posY.at(i) = voxelIDY * simParams->voxelSize + locY.at(this_owner) * simParams->l + this_sp_deviation_y +
                     simParams->LBFY;
        posZ.at(i) = voxelIDZ * simParams->voxelSize + locZ.at(this_owner) * simParams->l + this_sp_deviation_z +
                     simParams->LBFZ;
        // std::cout << "Sphere Pos: " << posX.at(i) << ", " << posY.at(i) << ", " << posZ.at(i) << std::endl;

        spRadii.at(i) = radiiSphere.at(clumpComponentOffset.at(i));
    }
    pw.write(ptFile, ParticleFormatWriter::CompressionType::NONE, posX, posY, posZ, spRadii);
}

inline void DEMDynamicThread::contactEventArraysResize(size_t nContactPairs) {
    TRACKED_QUICK_VECTOR_RESIZE(idGeometryA, nContactPairs);
    TRACKED_QUICK_VECTOR_RESIZE(idGeometryB, nContactPairs);
    TRACKED_QUICK_VECTOR_RESIZE(contactForces, nContactPairs);
    TRACKED_QUICK_VECTOR_RESIZE(contactPointGeometryA, nContactPairs);
    TRACKED_QUICK_VECTOR_RESIZE(contactPointGeometryB, nContactPairs);

    // Re-pack pointers in case the arrays got reallocated
    granData->idGeometryA = idGeometryA.data();
    granData->idGeometryB = idGeometryB.data();
    granData->contactForces = contactForces.data();
    granData->contactPointGeometryA = contactPointGeometryA.data();
    granData->contactPointGeometryB = contactPointGeometryB.data();
}

inline void DEMDynamicThread::unpackMyBuffer() {
    GPU_CALL(cudaMemcpy(stateOfSolver_resources.getNumContactsPointer(), &(granData->nContactPairs_buffer),
                        sizeof(size_t), cudaMemcpyDeviceToDevice));

    // Need to resize those contact event-based arrays before usage
    if (stateOfSolver_resources.getNumContacts() > idGeometryA.size()) {
        contactEventArraysResize(stateOfSolver_resources.getNumContacts());
    }

    GPU_CALL(cudaMemcpy(granData->idGeometryA, granData->idGeometryA_buffer,
                        stateOfSolver_resources.getNumContacts() * sizeof(bodyID_t), cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->idGeometryB, granData->idGeometryB_buffer,
                        stateOfSolver_resources.getNumContacts() * sizeof(bodyID_t), cudaMemcpyDeviceToDevice));
}

inline void DEMDynamicThread::sendToTheirBuffer() {
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_voxelID, granData->voxelID,
                        simParams->nClumpBodies * sizeof(voxelID_t), cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_locX, granData->locX, simParams->nClumpBodies * sizeof(subVoxelPos_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_locY, granData->locY, simParams->nClumpBodies * sizeof(subVoxelPos_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_locZ, granData->locZ, simParams->nClumpBodies * sizeof(subVoxelPos_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ0, granData->oriQ0, simParams->nClumpBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ1, granData->oriQ1, simParams->nClumpBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ2, granData->oriQ2, simParams->nClumpBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ3, granData->oriQ3, simParams->nClumpBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));
}

inline void DEMDynamicThread::calculateForces() {
    // reset force (acceleration) arrays for this time step and apply gravity
    size_t threads_needed_for_prep = simParams->nClumpBodies > stateOfSolver_resources.getNumContacts()
                                         ? simParams->nClumpBodies
                                         : stateOfSolver_resources.getNumContacts();
    size_t blocks_needed_for_prep = (threads_needed_for_prep + NUM_BODIES_PER_BLOCK - 1) / NUM_BODIES_PER_BLOCK;

    auto prep_force =
        JitHelper::buildProgram("DEMPrepForceKernels", JitHelper::KERNEL_DIR / "DEMPrepForceKernels.cu",
                                std::vector<JitHelper::Header>(), {"-I" + (JitHelper::KERNEL_DIR / "..").string()});

    prep_force.kernel("prepareForceArrays")
        .instantiate()
        .configure(dim3(blocks_needed_for_prep), dim3(NUM_BODIES_PER_BLOCK), sizeof(float) * TEST_SHARED_SIZE * 4,
                   streamInfo.stream)
        .launch(simParams, granData, stateOfSolver_resources.getNumContacts(), granTemplates);
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));

    // TODO: is there a better way??? Like memset?
    // GPU_CALL(cudaMemset(granData->contactForces, zeros, stateOfSolver_resources.getNumContacts() * sizeof(float3)));
    // GPU_CALL(cudaMemset(granData->h2AlphaX, 0, simParams->nClumpBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->h2AlphaY, 0, simParams->nClumpBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->h2AlphaZ, 0, simParams->nClumpBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->h2aX,
    //                     (double)simParams->h * (double)simParams->h * (double)simParams->Gx / (double)simParams->l,
    //                     simParams->nClumpBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->h2aY,
    //                     (double)simParams->h * (double)simParams->h * (double)simParams->Gy / (double)simParams->l,
    //                     simParams->nClumpBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->h2aZ,
    //                     (double)simParams->h * (double)simParams->h * (double)simParams->Gz / (double)simParams->l,
    //                     simParams->nClumpBodies * sizeof(float)));

    size_t blocks_needed_for_contacts =
        (stateOfSolver_resources.getNumContacts() + NUM_BODIES_PER_BLOCK - 1) / NUM_BODIES_PER_BLOCK;
    auto cal_force =
        JitHelper::buildProgram("DEMFrictionlessForceKernels", JitHelper::KERNEL_DIR / "DEMFrictionlessForceKernels.cu",
                                std::vector<JitHelper::Header>(), {"-I" + (JitHelper::KERNEL_DIR / "..").string()});

    // a custom kernel to compute forces
    cal_force.kernel("calculateNormalContactForces")
        .instantiate()
        .configure(dim3(blocks_needed_for_contacts), dim3(NUM_BODIES_PER_BLOCK), sizeof(float) * TEST_SHARED_SIZE * 5,
                   streamInfo.stream)
        .launch(simParams, granData, stateOfSolver_resources.getNumContacts(), granTemplates);
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
    // displayFloat3(granData->contactForces, stateOfSolver_resources.getNumContacts());
    // std::cout << "===========================" << std::endl;

    // Reflect those body-wise forces on their owner clumps
    // hostCollectForces(granData->inertiaPropOffsets, granData->idGeometryA, granData->idGeometryB,
    //                   granData->contactForces, granData->h2aX, granData->h2aY, granData->h2aZ,
    //                   granData->ownerClumpBody, granTemplates->massClumpBody, simParams->h,
    //                   stateOfSolver_resources.getNumContacts(),simParams->l);
    cubCollectForces(granData->inertiaPropOffsets, granData->idGeometryA, granData->idGeometryB,
                     granData->contactForces, granData->contactPointGeometryA, granData->contactPointGeometryB,
                     granData->h2aX, granData->h2aY, granData->h2aZ, granData->h2AlphaX, granData->h2AlphaY,
                     granData->h2AlphaZ, granData->ownerClumpBody, granTemplates->massClumpBody, granTemplates->mmiXX,
                     granTemplates->mmiYY, granTemplates->mmiZZ, simParams->h, stateOfSolver_resources.getNumContacts(),
                     simParams->nClumpBodies, simParams->l, contactPairArr_isFresh, streamInfo.stream,
                     stateOfSolver_resources, simParams->nDistinctClumpBodyTopologies);
    // displayArray<float>(granData->h2aX, simParams->nClumpBodies);
    // displayFloat3(granData->contactForces, stateOfSolver_resources.getNumContacts());
    // std::cout << stateOfSolver_resources.getNumContacts() << std::endl;

    // Calculate the torque on owner clumps from those body-wise forces
    // hostCollectTorques(granData->inertiaPropOffsets, granData->idGeometryA, granData->idGeometryB,
    //                    granData->contactForces, granData->contactPointGeometryA, granData->contactPointGeometryB,
    //                    granData->h2AlphaX, granData->h2AlphaY, granData->h2AlphaZ, granData->ownerClumpBody,
    //                    granTemplates->mmiXX, granTemplates->mmiYY, granTemplates->mmiZZ, simParams->h,
    //                    stateOfSolver_resources.getNumContacts(), simParams->l);
    // displayArray<float>(granData->oriQ0, simParams->nClumpBodies);
    // displayArray<float>(granData->oriQ1, simParams->nClumpBodies);
    // displayArray<float>(granData->oriQ2, simParams->nClumpBodies);
    // displayArray<float>(granData->oriQ3, simParams->nClumpBodies);
    // displayArray<float>(granData->h2AlphaX, simParams->nClumpBodies);
    // displayArray<float>(granData->h2AlphaY, simParams->nClumpBodies);
    // displayArray<float>(granData->h2AlphaZ, simParams->nClumpBodies);
}

inline void DEMDynamicThread::integrateClumpMotions() {
    size_t blocks_needed_for_clumps = (simParams->nClumpBodies + NUM_BODIES_PER_BLOCK - 1) / NUM_BODIES_PER_BLOCK;
    auto integrator =
        JitHelper::buildProgram("DEMIntegrationKernels", JitHelper::KERNEL_DIR / "DEMIntegrationKernels.cu",
                                std::vector<JitHelper::Header>(), {"-I" + (JitHelper::KERNEL_DIR / "..").string()});
    integrator.kernel("integrateClumps")
        .instantiate()
        .configure(dim3(blocks_needed_for_clumps), dim3(NUM_BODIES_PER_BLOCK), 0, streamInfo.stream)
        .launch(simParams, granData, granTemplates);
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
}

void DEMDynamicThread::workerThread() {
    // Set the gpu for this thread
    cudaSetDevice(streamInfo.device);
    cudaStreamCreate(&streamInfo.stream);
    int totGPU;
    cudaGetDeviceCount(&totGPU);
    printf("Total device: %d\n", totGPU);

    while (!pSchedSupport->dynamicShouldJoin) {
        {
            std::unique_lock<std::mutex> lock(pSchedSupport->dynamicStartLock);
            while (!pSchedSupport->dynamicStarted) {
                pSchedSupport->cv_DynamicStartLock.wait(lock);
            }
            // Ensure that we wait for start signal on next iteration
            pSchedSupport->dynamicStarted = false;
            if (pSchedSupport->dynamicShouldJoin) {
                break;
            }
        }

        // At the beginning of each user call, send kT a work order, b/c dT need results from CD to proceed. After this
        // one instance, kT and dT may work in an async fashion.
        {
            std::lock_guard<std::mutex> lock(pSchedSupport->kinematicOwnedBuffer_AccessCoordination);
            sendToTheirBuffer();
        }
        pSchedSupport->kinematicOwned_Cons2ProdBuffer_isFresh = true;
        contactPairArr_isFresh = true;
        pSchedSupport->schedulingStats.nKinematicUpdates++;
        // Signal the kinematic that it has data for a new work order.
        pSchedSupport->cv_KinematicCanProceed.notify_all();
        // Then dT will wait for kT to finish one initial run
        {
            std::unique_lock<std::mutex> lock(pSchedSupport->dynamicCanProceed);
            while (!pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh) {
                // loop to avoid spurious wakeups
                pSchedSupport->cv_DynamicCanProceed.wait(lock);
            }
        }

        for (size_t cycle = 0; cycle < nDynamicCycles; cycle++) {
            // if the produce is fresh, use it
            if (pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh) {
                {
                    // acquire lock and use the content of the dynamic-owned transfer buffer
                    std::lock_guard<std::mutex> lock(pSchedSupport->dynamicOwnedBuffer_AccessCoordination);
                    // std::this_thread::sleep_for(std::chrono::milliseconds(WAIT_GRANULARITY_MS));
                    unpackMyBuffer();
                    contactPairArr_isFresh = true;
                }
                // dT got the produce, now mark its buffer to be no longer fresh
                pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh = false;
                pSchedSupport->stampLastUpdateOfDynamic = cycle;
            }

            calculateForces();

            integrateClumpMotions();

            // calculateForces is done, set it to false
            // will be set to true next time it receives an update from kT
            contactPairArr_isFresh = false;

            // if it's the case, it's important at this point to let the kinematic know that this is the last dynamic
            // cycle; this is important otherwise the kinematic will hang waiting for communication swith the dynamic
            if (cycle == (nDynamicCycles - 1))
                pSchedSupport->dynamicDone = true;

            // if the kinematic is idle, give it the opportunity to get busy again
            if (!pSchedSupport->kinematicOwned_Cons2ProdBuffer_isFresh) {
                // acquire lock and refresh the work order for the kinematic
                {
                    std::lock_guard<std::mutex> lock(pSchedSupport->kinematicOwnedBuffer_AccessCoordination);
                    sendToTheirBuffer();
                }
                pSchedSupport->kinematicOwned_Cons2ProdBuffer_isFresh = true;
                pSchedSupport->schedulingStats.nKinematicUpdates++;
                // signal the kinematic that it has data for a new work order
                pSchedSupport->cv_KinematicCanProceed.notify_all();
            }

            // std::cout << "dT Total contact pairs: " << granData->nContactPairs_buffer << std::endl;
            // std::cout << "Dynamic side values. Cycle: " << cycle << std::endl;

            // dynamic wrapped up one cycle
            pSchedSupport->currentStampOfDynamic++;

            // check if we need to wait; i.e., if dynamic drifted too much into future, then we must wait a bit before
            // the next cycle begins
            if (pSchedSupport->dynamicShouldWait()) {
                // wait for a signal from the kinematic to indicate that
                // the kinematic has caught up
                pSchedSupport->schedulingStats.nTimesDynamicHeldBack++;
                std::unique_lock<std::mutex> lock(pSchedSupport->dynamicCanProceed);
                while (!pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh) {
                    // loop to avoid spurious wakeups
                    pSchedSupport->cv_DynamicCanProceed.wait(lock);
                }
            }
        }

        // When getting here, dT has finished one user call (although perhaps not at the end of the user script).
        userCallDone = true;
    }
}

void DEMDynamicThread::startThread() {
    std::lock_guard<std::mutex> lock(pSchedSupport->dynamicStartLock);
    pSchedSupport->dynamicStarted = true;
    pSchedSupport->cv_DynamicStartLock.notify_one();
}

bool DEMDynamicThread::isUserCallDone() {
    // return true if done, false if not
    return userCallDone;
}

void DEMDynamicThread::resetUserCallStat() {
    userCallDone = false;
    // Reset last kT-side data receiving cycle time stamp.
    pSchedSupport->stampLastUpdateOfDynamic = -1;
    pSchedSupport->currentStampOfDynamic = 0;
    // Reset dT stats variables, making ready for next user call
    pSchedSupport->dynamicDone = false;
    pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh = false;
    contactPairArr_isFresh = true;
}

}  // namespace sgps
