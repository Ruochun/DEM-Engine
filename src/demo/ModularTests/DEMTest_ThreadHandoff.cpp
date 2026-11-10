// Copyright (c) 2021, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause

#include "DEM/API.h"

#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <thread>

using namespace deme;

// Exercise real worker handoffs with short calls, forced drift waits, sync/restart, and destruction.
// This is a scheduling stress regression, not a deterministic reproduction of every lost-wakeup interleaving.
int main() {
    std::mutex completion_mutex;
    std::condition_variable completion_cv;
    bool completed = false;
    // Fail rather than hang indefinitely if either worker loses a notification, including during destruction.
    // Include enough startup time for runtime kernel compilation with a cold JIT cache.
    std::thread watchdog([&]() {
        std::unique_lock<std::mutex> lock(completion_mutex);
        if (!completion_cv.wait_for(lock, std::chrono::seconds(600), [&]() { return completed; })) {
            std::cerr << "FAIL: worker handoff timed out\n" << std::flush;
            std::_Exit(EXIT_FAILURE);
        }
    });

    bool passed = true;
    for (int drift : {0, 1, 4}) {
        std::cout << "Initializing handoff test with update frequency " << drift << std::endl;
        DEMSolver solver;
        solver.SetVerbosity("ERROR");
        solver.InstructBoxDomainDimension(2.0, 2.0, 2.0);
        solver.SetInitBinNumTarget(8);
        solver.SetTimeStepSize(1.e-4);
        solver.SetGravitationalAcceleration(make_float3(0.f, 0.f, 0.f));
        solver.SetCDUpdateFreq(drift);
        solver.DisableAdaptiveUpdateFreq();
        auto material = solver.LoadMaterial({{"E", 1.e7f}, {"nu", 0.3f}, {"CoR", 0.2f}, {"mu", 0.4f}, {"Crr", 0.f}});
        auto sphere = solver.LoadSphereType(1.f, 0.05f, material);
        solver.AddClumps(sphere, make_float3(0.f, 0.f, 0.f));
        solver.Initialize();
        std::cout << "Exercising handoffs" << std::endl;

        // Half a timestep requests exactly one step without sensitivity to float timestep rounding.
        for (int call = 0; call < 200; ++call) {
            const double before = solver.GetSimTime();
            solver.DoDynamics(5.e-5);
            const double advance = solver.GetSimTime() - before;
            passed = passed && std::isfinite(advance) && std::abs(advance - 1.e-4) < 1.e-9;
            if (call % 7 == 0) {
                solver.DoDynamicsThenSync(0.0);
            }
            if (call % 3 == 0) {
                std::this_thread::yield();
            }
        }
        // Leave the final call asynchronous to exercise destruction with a potentially outstanding kT order.
    }

    {
        std::lock_guard<std::mutex> lock(completion_mutex);
        completed = true;
    }
    completion_cv.notify_one();
    watchdog.join();
    std::cout << (passed ? "PASS" : "FAIL") << ": worker handoff and restart\n";
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
