#!/usr/bin/env bash
# Each iteration is a fresh process, matching the Python repeat.sh reproduction.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
runs=${RUNS:-100}
binary=${MINISLUMP_BINARY:-$repo/build-minislump-debug/bin/DEMdemo_DFCMiniSlump}
[[ $runs =~ ^[1-9][0-9]*$ ]] || { echo 'RUNS must be a positive integer' >&2; exit 2; }
[[ -x $binary ]] || { echo "Build first: $repo/scripts/minislump/build-debug.sh" >&2; exit 2; }
python3 "$repo/scripts/minislump/verify-assets.py"
mkdir -p "$repo/build-minislump-debug/runs"
run_root=$(mktemp -d "$repo/build-minislump-debug/runs/repeat-XXXXXXXX")
child=
trap 'if [[ -n $child ]]; then kill -TERM "$child" 2>/dev/null || true; wait "$child" || true; fi; exit 130' INT TERM
for ((i=1; i<=runs; ++i)); do
    dir="$run_root/run-$i"
    mkdir "$dir"
    "$binary" "$@" --out-dir "$dir/output" > "$dir/run.log" 2>&1 &
    child=$!
    echo "$child" > "$dir/pid"
    echo "Run $i/$runs: PID $child; log $dir/run.log"
    echo 'If stalled, attach with VS Code MiniSlump: attach, then Pause and inspect all threads.'
    if wait "$child"; then child=; else status=$?; child=; echo "Run $i failed ($status); evidence: $dir" >&2; exit "$status"; fi
done
echo "Completed $runs fresh-process runs: $run_root"
