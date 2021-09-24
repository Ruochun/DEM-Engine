// *----------------------------------------
// GPU - Testing kernels

__global__ void dynamicTestKernel() {
    printf("Dynamic run\n");
}
__global__ void kinematicTestKernel(sgps::voxelID_default_t* data) {
    if (threadIdx.x == 0) {
        printf("Kinematic run\n");
    }

    if (threadIdx.x < N_INPUT_ITEMS) {
        // data[threadIdx.x] = 2 * data[threadIdx.x] + 1;
        // printf("%d\n", data[threadIdx.x]);
    }
}
// END of GPU Testing kernels
// *----------------------------------------