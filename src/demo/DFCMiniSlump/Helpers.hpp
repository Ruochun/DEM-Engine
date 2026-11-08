// CPU-only argument parsing and supported-surface diagnostics for the MiniSlump demo.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace minislump {
struct Options {
    std::filesystem::path particles =
        std::filesystem::path(DFC_MINISLUMP_DATA_DIR) / "DFCgeo000-data-particles-mini-2x4-dist.dat";
    std::filesystem::path output = "OUT_MiniSlump_CPP";
    std::string out_name = "output.txt";
    bool settle = true, no_cone = false, diagnostic = false, release = false;
    double frame_time = 0, duration = 9.5, plate = 50, sensitivity = 0;
    unsigned int bins = 500000;
    long long frame_steps = 800;
};
// Preserve Python's flags, including the counterintuitive --settle meaning skip settling.
inline Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];
        if (key == "--help" || key == "-h") {
            std::cout << "DFC MiniSlump: --particle_file/-p FILE --out-dir/-d DIR --out_name/-o FILE\n"
                         "--settle/-s (skip settling) --no-cone --frame-steps N --frame-time SECONDS\n"
                         "--slump-duration SECONDS --top-plate-height MM --init-bin-num N\n"
                         "--sensitivity-dt {1e-6,5e-7} --initial-force-diagnostic --verify-cone-release\n";
            std::exit(0);
        }
        if (key == "--settle" || key == "-s") {
            o.settle = false;
            continue;
        }
        if (key == "--no-cone") {
            o.no_cone = true;
            continue;
        }
        if (key == "--initial-force-diagnostic") {
            o.diagnostic = true;
            continue;
        }
        if (key == "--verify-cone-release") {
            o.release = o.diagnostic = true;
            continue;
        }
        if (++i == argc)
            throw std::runtime_error("Missing value for " + key);
        std::string value = argv[i];
        if (key == "--particle_file" || key == "-p")
            o.particles = value;
        else if (key == "--out-dir" || key == "-d")
            o.output = value;
        else if (key == "--out_name" || key == "-o")
            o.out_name = value;
        else {
            size_t end;
            double n = std::stod(value, &end);
            if (end != value.size() || !std::isfinite(n))
                throw std::runtime_error("Invalid value for " + key);
            if (key == "--frame-time" && n > 0)
                o.frame_time = n;
            else if (key == "--slump-duration" && n >= 0)
                o.duration = n;
            else if (key == "--top-plate-height")
                o.plate = n;
            else if (key == "--sensitivity-dt" && (n == 1e-6 || n == 5e-7))
                o.sensitivity = n;
            else if (key == "--frame-steps" && n >= 1 && n <= 2147483647 && std::floor(n) == n)
                o.frame_steps = static_cast<long long>(n);
            else if (key == "--init-bin-num" && n >= 1 && n <= 2147483647 && std::floor(n) == n)
                o.bins = static_cast<unsigned int>(n);
            else
                throw std::runtime_error("Unknown option or invalid value: " + key);
        }
    }
    if ((o.diagnostic && o.sensitivity) || ((o.diagnostic || o.sensitivity) && (o.no_cone || !o.settle)))
        throw std::runtime_error("Diagnostics/sensitivity require settling and the cone, and are mutually exclusive");
    if (o.sensitivity)
        o.frame_time = .001;
    return o;
}
// Fail on output errors instead of silently losing evidence from a stalled run.
inline std::ofstream output(const std::filesystem::path& path) {
    std::ofstream f(path);
    f.exceptions(std::ios::badbit | std::ios::failbit);
    f.precision(17);
    return f;
}
inline std::string read(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Cannot read " + path.string());
    return {std::istreambuf_iterator<char>(f), {}};
}
struct Surface {
    double top = 0, average = std::numeric_limits<double>::quiet_NaN();
    double top_radius = std::numeric_limits<double>::quiet_NaN(), bottom_radius = top_radius;
    size_t top_id = 0, top_count = 0, bottom_count = 0;
};
// Descend/ascend bounded bands until at least five particles support a surface;
// higher outliers are excluded instead of widening the band, matching Python.
inline Surface surface(const std::vector<float3>& p, const std::vector<double>& r) {
    Surface s;
    for (bool upper : {true, false}) {
        std::vector<double> levels;
        for (size_t i = 0; i < p.size(); ++i)
            levels.push_back(p[i].y + (upper ? r[i] : -r[i]));
        if (upper) {
            s.top_id = std::max_element(levels.begin(), levels.end()) - levels.begin();
            s.top = levels[s.top_id];
        }
        auto sorted = levels;
        std::sort(sorted.begin(), sorted.end());
        for (size_t band = 0; band < sorted.size(); ++band) {
            double anchor = sorted[upper ? sorted.size() - 1 - band : band];
            double low = upper ? anchor - 2.15 : anchor, high = upper ? anchor : anchor + 2.15;
            auto lo = std::lower_bound(sorted.begin(), sorted.end(), low);
            auto hi = std::upper_bound(sorted.begin(), sorted.end(), high);
            size_t count = hi - lo;
            if (count < 5)
                continue;
            double sum = 0, radius = 0;
            for (size_t i = 0; i < p.size(); ++i)
                if (levels[i] >= low && levels[i] <= high) {
                    sum += levels[i];
                    radius = std::max(radius, std::hypot(double(p[i].x), double(p[i].z)) + r[i]);
                }
            if (upper) {
                s.average = sum / count;
                s.top_count = count;
                s.top_radius = radius;
            } else {
                s.bottom_count = count;
                s.bottom_radius = radius;
            }
            break;
        }
    }
    return s;
}
}  // namespace minislump
