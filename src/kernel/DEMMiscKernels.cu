// DEM misc. kernels
#include <DEM/Defines.h>

__global__ void markOwnerToChange(deme::notStupidBool_t* idBool,
                                  float* ownerFactors,
                                  deme::bodyID_t* dIDs,
                                  float* dFactors,
                                  size_t n) {
    size_t myID = blockIdx.x * blockDim.x + threadIdx.x;
    if (myID < n) {
        deme::bodyID_t myOwner = dIDs[myID];
        float myFactor = dFactors[myID];
        idBool[myOwner] = 1;
        ownerFactors[myOwner] = myFactor;
    }
}

__global__ void dTModifyComponents(deme::DEMDataDT* granData, deme::notStupidBool_t* idBool, float* factors, size_t n) {
    size_t sphereID = blockIdx.x * blockDim.x + threadIdx.x;
    if (sphereID < n) {
        // Get my owner ID
        deme::bodyID_t myOwner = granData->ownerClumpBody[sphereID];
        // If not marked, we have nothing to do
        if (idBool[myOwner]) {
            float factor = factors[myOwner];
            // Expand radius and relPos
            granData->relPosSphereX[sphereID] *= factor;
            granData->relPosSphereY[sphereID] *= factor;
            granData->relPosSphereZ[sphereID] *= factor;
            granData->radiiSphere[sphereID] *= factor;
        }
    }
}

// How to template it???
__global__ void kTModifyComponents(deme::DEMDataKT* granData, deme::notStupidBool_t* idBool, float* factors, size_t n) {
    size_t sphereID = blockIdx.x * blockDim.x + threadIdx.x;
    if (sphereID < n) {
        // Get my owner ID
        deme::bodyID_t myOwner = granData->ownerClumpBody[sphereID];
        // If not marked, we have nothing to do
        if (idBool[myOwner]) {
            float factor = factors[myOwner];
            // Expand radius and relPos
            granData->relPosSphereX[sphereID] *= factor;
            granData->relPosSphereY[sphereID] *= factor;
            granData->relPosSphereZ[sphereID] *= factor;
            granData->radiiSphere[sphereID] *= factor;
        }
    }
}

__global__ void computeMarginFromAbsv(deme::DEMSimParams* simParams, deme::DEMDataKT* granData, size_t n) {
    size_t ownerID = blockIdx.x * blockDim.x + threadIdx.x;
    if (ownerID < n) {
        float absv = granData->marginSize[ownerID];
        if (absv > simParams->approxMaxVel) {
            absv = simParams->approxMaxVel;
        }
        granData->marginSize[ownerID] =
            (absv * simParams->expSafetyMulti + simParams->expSafetyAdder) * (granData->ts_buffer * granData->maxDrift);
    }
}

__global__ void fillValues(float* arr, float val, size_t n) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        arr[i] = val;
    }
}
