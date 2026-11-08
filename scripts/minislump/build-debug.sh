#!/usr/bin/env bash
# Build host debug symbols without changing the user's existing build tree.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
python3 "$repo/scripts/minislump/verify-assets.py"
cmake -S "$repo" -B "$repo/build-minislump-debug" \
    -DCMAKE_BUILD_TYPE=Debug -DDEME_BUILD_DEMOS=ON -DDEME_BUILD_PYTHON=OFF \
    -DDEME_BUILD_VISUALIZER=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=ON "$@"
cmake --build "$repo/build-minislump-debug" --target DEMdemo_DFCMiniSlump --parallel 1
