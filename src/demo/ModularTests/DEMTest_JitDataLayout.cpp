// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#include "core/utils/JitHelper.h"
#include <chrono>
#include <iostream>
#include <stdexcept>

// Changing an included data header at the SAME path must produce a new executable kernel. A subsequent unchanged
// build must reuse its key. This catches stale layout caches without modifying any real runtime header or cache.
int main() {
    try {
        const auto root =
            std::filesystem::temp_directory_path() /
            ("deme-jit-layout-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root / "DEM");
        const auto source = root / "probe.cu";
        std::ofstream(source) << "#include <DEM/Defines.h>\n"
                                 "extern \"C\" __global__ void readValue(int* out) { *out = DEME_TEST_VALUE; }\n";
        std::ofstream(root / "DEM" / "VariableTypes.h") << "// fixture\n";
        JitHelper::KERNEL_INCLUDE_DIR = root;
        int* device = nullptr;
        if (cudaMalloc(reinterpret_cast<void**>(&device), sizeof(int)) != cudaSuccess)
            throw std::runtime_error("CUDA allocation failed");
        cudaStream_t stream = nullptr;
        std::string previous;
        int iteration = 0;
        for (int value : {11, 29, 29}) {
            std::ofstream(root / "DEM" / "Defines.h") << "#define DEME_TEST_VALUE " << value << '\n';
            const auto program = JitHelper::buildProgram("data_layout_probe", source);
            if ((iteration == 1 && program.key() == previous) || (iteration == 2 && program.key() != previous))
                throw std::runtime_error("Header content change/reuse did not update the JIT cache key correctly");
            program.kernel("readValue").instantiate().configure(dim3(1), dim3(1), 0, stream).launch(device);
            int actual = 0;
            if (cudaMemcpy(&actual, device, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess || actual != value)
                throw std::runtime_error("Cached kernel used stale included header contents");
            previous = program.key();
            ++iteration;
        }
        cudaFree(device);
        std::cout << "PASS: included data-header changes invalidate JIT kernels; unchanged headers reuse the cache.\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
