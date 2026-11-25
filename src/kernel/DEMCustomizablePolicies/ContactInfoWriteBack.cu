granData->contactPointGeometryA[myPatchContactID] = locCPA;
granData->contactPointGeometryB[myPatchContactID] = locCPB;
granData->contactForcesA[myPatchContactID] = force + A_only_force;
granData->contactForcesB[myPatchContactID] = -force + B_only_force;
granData->contactTorquesA[myPatchContactID] = A_only_torque;
granData->contactTorquesB[myPatchContactID] = B_only_torque;
if (simParams->storeNormal) {
    granData->contactNormals[myPatchContactID] = B2A;
}
