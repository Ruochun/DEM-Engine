// Copyright (c) 2021, SBEL GPU Development Team
// Copyright (c) 2021, University of Wisconsin - Madison
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DEM/utils/DEMVisualizer.h"

#include <algorithm>
#include <string>
#include <unordered_map>

#include "DEM/API.h"
#include "core/utils/Logger.hpp"

// raylib defines convenience macros such as PI, so keep it after DEME headers to avoid altering their declarations.
#include <raylib.h>

namespace deme {
namespace {

Vector3 toRayVector(float3 value) {
    return Vector3{value.x, value.y, value.z};
}

Color toRayColor(DEMVisualizerColor color) {
    return Color{color.r, color.g, color.b, color.a};
}

DEMVisualizerColor defaultFamilyColor(family_t family) {
    // The multiplicative hash keeps adjacent family IDs visually distinct without maintaining a fixed-size palette.
    const std::uint32_t hash = static_cast<std::uint32_t>(family) * 2654435761u;
    return DEMVisualizerColor{static_cast<std::uint8_t>(80u + (hash & 0x7fu)),
                              static_cast<std::uint8_t>(80u + ((hash >> 8u) & 0x7fu)),
                              static_cast<std::uint8_t>(80u + ((hash >> 16u) & 0x7fu)), 255};
}

}  // namespace

struct DEMVisualizer::Impl {
    explicit Impl(DEMSolver& solver_in) : solver(solver_in) {}

    DEMSolver& solver;
    int width = 1280;
    int height = 720;
    int target_fps = 60;
    std::string title = "DEM-Engine 3 Visualizer";
    bool initialized = false;
    bool render_spheres = true;
    bool render_triangles = true;
    DEMVisualizerColor background{245, 245, 245, 255};
    Camera3D camera{Vector3{4.0f, 4.0f, 3.0f}, Vector3{0.0f, 0.0f, 0.0f}, Vector3{0.0f, 0.0f, 1.0f}, 45.0f,
                    CAMERA_PERSPECTIVE};
    std::unordered_map<family_t, DEMVisualizerColor> family_colors;

    Color familyColor(family_t family) const {
        const auto found = family_colors.find(family);
        return toRayColor(found == family_colors.end() ? defaultFamilyColor(family) : found->second);
    }
};

DEMVisualizer::DEMVisualizer(DEMSolver& solver) : m_impl(std::make_unique<Impl>(solver)) {}

DEMVisualizer::~DEMVisualizer() {
    Close();
}

void DEMVisualizer::Initialize() {
    if (m_impl->initialized) {
        return;
    }
    if (IsWindowReady()) {
        DEME_ERROR("Only one DEMVisualizer window can be initialized in a process at a time.");
    }
    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(m_impl->width, m_impl->height, m_impl->title.c_str());
    if (!IsWindowReady()) {
        DEME_ERROR("DEMVisualizer could not create a graphics window. Check DISPLAY/WSLg and OpenGL availability.");
    }
    ::SetTargetFPS(m_impl->target_fps);
    m_impl->initialized = true;
}

bool DEMVisualizer::Run() const {
    return m_impl->initialized && !WindowShouldClose();
}

void DEMVisualizer::Render() {
    if (!m_impl->initialized) {
        DEME_ERROR("DEMVisualizer::Render requires Initialize() to be called first.");
    }

    const DEMVisualizationSnapshot snapshot =
        m_impl->solver.GetVisualizationSnapshot(m_impl->render_spheres, m_impl->render_triangles);
    UpdateCamera(&m_impl->camera, CAMERA_FREE);

    BeginDrawing();
    ClearBackground(toRayColor(m_impl->background));
    BeginMode3D(m_impl->camera);

    if (m_impl->render_spheres) {
        for (const auto& sphere : snapshot.spheres) {
            DrawSphereEx(toRayVector(sphere.position), sphere.radius, 8, 12, m_impl->familyColor(sphere.family));
        }
    }
    if (m_impl->render_triangles) {
        for (const auto& triangle : snapshot.triangles) {
            const Color color = m_impl->familyColor(triangle.family);
            DrawTriangle3D(toRayVector(triangle.a), toRayVector(triangle.b), toRayVector(triangle.c), color);
            // Draw both windings so a simple inspection viewer does not hide intentionally one-sided mesh facets.
            DrawTriangle3D(toRayVector(triangle.c), toRayVector(triangle.b), toRayVector(triangle.a), color);
        }
    }

    // DEME is Z-up, whereas raylib's convenience grid is Y-up. Draw the reference grid explicitly in the XY plane.
    const Color grid_color = Color{200, 200, 200, 255};
    for (int line = -10; line <= 10; ++line) {
        DrawLine3D(Vector3{static_cast<float>(line), -10.0f, 0.0f}, Vector3{static_cast<float>(line), 10.0f, 0.0f},
                   grid_color);
        DrawLine3D(Vector3{-10.0f, static_cast<float>(line), 0.0f}, Vector3{10.0f, static_cast<float>(line), 0.0f},
                   grid_color);
    }
    EndMode3D();
    DrawText(TextFormat("t = %.6g s", snapshot.simulation_time), 10, 10, 20, DARKGRAY);
    DrawFPS(10, 36);
    EndDrawing();
}

void DEMVisualizer::Close() {
    if (!m_impl || !m_impl->initialized) {
        return;
    }
    CloseWindow();
    m_impl->initialized = false;
}

void DEMVisualizer::SetWindowSize(int width, int height) {
    if (width <= 0 || height <= 0) {
        DEME_ERROR("DEMVisualizer window dimensions must be positive (got %d x %d).", width, height);
    }
    m_impl->width = width;
    m_impl->height = height;
    if (m_impl->initialized) {
        ::SetWindowSize(width, height);
    }
}

void DEMVisualizer::SetWindowTitle(const std::string& title) {
    m_impl->title = title;
    if (m_impl->initialized) {
        ::SetWindowTitle(title.c_str());
    }
}

void DEMVisualizer::SetTargetFPS(int fps) {
    if (fps <= 0) {
        DEME_ERROR("DEMVisualizer target FPS must be positive (got %d).", fps);
    }
    m_impl->target_fps = fps;
    if (m_impl->initialized) {
        ::SetTargetFPS(fps);
    }
}

void DEMVisualizer::SetCameraPosition(float3 position) {
    m_impl->camera.position = toRayVector(position);
}

void DEMVisualizer::SetCameraTarget(float3 target) {
    m_impl->camera.target = toRayVector(target);
}

void DEMVisualizer::SetBackgroundColor(DEMVisualizerColor color) {
    m_impl->background = color;
}

void DEMVisualizer::SetFamilyColor(family_t family, DEMVisualizerColor color) {
    m_impl->family_colors[family] = color;
}

void DEMVisualizer::SetRenderSpheres(bool render) {
    m_impl->render_spheres = render;
}

void DEMVisualizer::SetRenderTriangles(bool render) {
    m_impl->render_triangles = render;
}

bool DEMVisualizer::IsRenderingSpheres() const {
    return m_impl->render_spheres;
}

bool DEMVisualizer::IsRenderingTriangles() const {
    return m_impl->render_triangles;
}

}  // namespace deme
