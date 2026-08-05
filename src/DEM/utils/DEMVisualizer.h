// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DEME_VISUALIZER_H
#define DEME_VISUALIZER_H

#include <cstdint>
#include <memory>
#include <string>

#include <cuda_runtime.h>

#include "../VariableTypes.h"

namespace deme {

class DEMSolver;

/// RGBA color used by the visualizer without exposing its rendering backend in DEME's public API.
struct DEMVisualizerColor {
    constexpr DEMVisualizerColor(std::uint8_t red = 255,
                                 std::uint8_t green = 255,
                                 std::uint8_t blue = 255,
                                 std::uint8_t alpha = 255)
        : r(red), g(green), b(blue), a(alpha) {}

    std::uint8_t r;
    std::uint8_t g;
    std::uint8_t b;
    std::uint8_t a;
};

/// Step-wise interactive viewer for the current state of a DEMSolver.
///
/// The visualizer never advances the simulation. Each Render call synchronously captures and draws exactly one solver
/// state, allowing applications to choose their own simulation-to-render cadence.
class DEMVisualizer {
  public:
    explicit DEMVisualizer(DEMSolver& solver);
    ~DEMVisualizer();

    DEMVisualizer(const DEMVisualizer&) = delete;
    DEMVisualizer& operator=(const DEMVisualizer&) = delete;
    DEMVisualizer(DEMVisualizer&&) = delete;
    DEMVisualizer& operator=(DEMVisualizer&&) = delete;

    /// Create the native window and graphics resources. Calling this more than once has no effect.
    void Initialize();
    /// Return true while the initialized visualization window remains open.
    bool Run() const;
    /// Capture the current solver state and draw one frame. This does not advance the solver.
    void Render();
    /// Release the window and graphics resources. Calling this more than once has no effect.
    void Close();

    void SetWindowSize(int width, int height);
    void SetWindowTitle(const std::string& title);
    void SetTargetFPS(int fps);
    void SetCameraPosition(float3 position);
    void SetCameraTarget(float3 target);
    void SetBackgroundColor(DEMVisualizerColor color);
    void SetFamilyColor(family_t family, DEMVisualizerColor color);

    /// Enable or disable component-sphere rendering. Enabled by default.
    void SetRenderSpheres(bool render = true);
    /// Enable or disable triangle rendering. Enabled by default.
    void SetRenderTriangles(bool render = true);
    bool IsRenderingSpheres() const;
    bool IsRenderingTriangles() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

}  // namespace deme

#endif
