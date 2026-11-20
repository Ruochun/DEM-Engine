// Copyright (c) 2026, SBEL GPU Development Team
// SPDX-License-Identifier: BSD-3-Clause
#ifndef DEME_FORCE_MODEL_HELPERS_HPP
#define DEME_FORCE_MODEL_HELPERS_HPP

#include <cctype>
#include <string>

#include "HostSideHelpers.hpp"

namespace deme {

// Scan source tokens rather than comments/string literals so migration notes and diagnostics remain legal.
// This is intentionally not a preprocessor: even an inactive branch using the removed output must be migrated.
inline void validateForceModelOutputs(const std::string& source) {
    for (size_t i = 0; i < source.size();) {
        if (source.compare(i, 2, "//") == 0) {
            const size_t end = source.find('\n', i + 2);
            i = end == std::string::npos ? source.size() : end + 1;
        } else if (source.compare(i, 2, "/*") == 0) {
            const size_t end = source.find("*/", i + 2);
            i = end == std::string::npos ? source.size() : end + 2;
        } else if (source.compare(i, 2, "R\"") == 0) {
            const size_t open = source.find('(', i + 2);
            if (open == std::string::npos)
                return;
            const std::string ending = ")" + source.substr(i + 2, open - i - 2) + "\"";
            const size_t end = source.find(ending, open + 1);
            i = end == std::string::npos ? source.size() : end + ending.size();
        } else if (source[i] == '"' || source[i] == '\'') {
            const char quote = source[i++];
            while (i < source.size()) {
                if (source[i] == '\\') {
                    i += (i + 1 < source.size()) ? 2 : 1;
                } else if (source[i++] == quote) {
                    break;
                }
            }
        } else if (std::isalpha(static_cast<unsigned char>(source[i])) || source[i] == '_') {
            const size_t begin = i++;
            while (i < source.size() && (std::isalnum(static_cast<unsigned char>(source[i])) || source[i] == '_'))
                ++i;
            if (source.compare(begin, i - begin, "torque_only_force") == 0) {
                DEME_ERROR(
                    "torque_only_force is deprecated and unsupported. Use A_only_torque and B_only_torque "
                    "(global direct torques, N*m), or A_only_force and B_only_force (global contact forces, N). "
                    "Do not copy a force-like value directly into a torque; migrate the physical moment law.");
            }
        } else {
            ++i;
        }
    }
}

}  // namespace deme
#endif
