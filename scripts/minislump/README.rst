C++ DFC MiniSlump stall reproduction
==================================

This ports ``/home/ruochunz/MiniSlumpDemo/pyDEME_DFCMiniSlump.py`` as supplied
on 2026-11-07. The particle file, double-precision DFC force model, OBJ/MTL,
and Chrono transform are copied verbatim into ``src/demo/DFCMiniSlump``.
The device prerequisites are extracted from Python; ``ConeTransform.hpp``
contains the JSON's reference-frame-to-world position and XYZW quaternion.
``verify-assets.py`` checks the pinned inputs, including the mesh SHA256 from
the transform. If changing the cone, regenerate its transform and C++ constants
and deliberately update the hashes together.

Files and prerequisites
-----------------------

* ``build-debug.sh`` validates the inputs, configures ``build-minislump-debug``,
  and builds only ``DEMdemo_DFCMiniSlump`` and its dependencies, with one build job.
* ``verify-assets.py`` checks the bundled input hashes without building or running
  the simulation. Invoke it with ``python3 scripts/minislump/verify-assets.py``.
* ``repeat.sh`` runs the executable repeatedly as fresh processes, without gdb.
* ``launch.json`` supplies the single-run gdb launch and existing-process attach
  configurations for VS Code.
* ``tasks.json`` supplies the pre-launch build task used by ``launch.json``.

Build and run
-------------

The command examples below assume the repository root as the terminal's current
directory. Build and run::

    scripts/minislump/build-debug.sh
    build-minislump-debug/bin/DEMdemo_DFCMiniSlump --out-dir build-minislump-debug/first-run
    RUNS=100 scripts/minislump/repeat.sh

The build script configures a separate Debug tree with host symbols for both
the demo and solver. Extra arguments go to CMake, for example
``-DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc``. Requires the normal
DEME CUDA/CUB dependencies, Python 3 for hash validation, gdb and VS Code's
Microsoft C/C++ extension. This uses host gdb to inspect dT/kT synchronization;
it does not provide CUDA device breakpoints.

A default run settles for 1.5 simulated seconds and then flows for 9.5 simulated
seconds. It can take substantially longer in wall-clock time. The first run in
a build tree also compiles runtime CUDA kernels; INFO messages identify these
compilations. An existing build can be rebuilt directly with::

    cmake --build build-minislump-debug --target DEMdemo_DFCMiniSlump --parallel 1

To use the existing Release build instead::

    cmake --build build --target DEMdemo_DFCMiniSlump --parallel 1
    MINISLUMP_BINARY="$PWD/build/bin/DEMdemo_DFCMiniSlump" RUNS=2 scripts/minislump/repeat.sh --settle --slump-duration 0.001

Use an absolute ``MINISLUMP_BINARY`` path when selecting another build.
The scripts locate the repository relative to their own location, so they can
also be invoked by absolute path from another directory. Relative paths supplied
as demo arguments are resolved from the process's working directory.

VS Code: one run under gdb
-------------------------

The tracked ``launch.json`` and ``tasks.json`` beside this document are ready
to copy into ``.vscode/``. If configurations already exist, merge their
``configurations``, ``inputs`` and ``tasks`` arrays instead of overwriting them.
For this workspace, the entries have already been merged while preserving the
existing launch entry. The repository ignores ``.vscode/``; the tracked copies
make the setup reproducible on another checkout.

Open the **repository root** as the VS Code workspace, for example with
``code .`` from the root. Do not open ``build`` or ``build-minislump-debug`` as
the workspace: the configurations resolve executable and script paths relative
to ``${workspaceFolder}``, and VS Code reads the root's ``.vscode`` directory.
This is separate from the program's working directory: ``cwd`` is set to
``${workspaceFolder}/build-minislump-debug``. The terminal's current directory
does not control the working directory when pressing F5.

Choose **MiniSlump: launch** and press F5. Enter a new output directory name
for each launch. The pre-launch task builds the Debug executable. Add demo
arguments to the configuration's ``args`` array as needed.

**F5 runs the demo once.** When the process exits, the debugging session ends;
it does not start another repetition. For another debug run, press F5 again
and enter a new output directory name, such as ``OUT_debug_2``. The existing
configuration does not automatically debug a sequence of fresh processes.

For a short flow run, edit the launch configuration's ``args`` array to::

    "args": [
        "--settle", "--slump-duration", "0.001", "--frame-steps", "100",
        "--out-dir", "${workspaceFolder}/build-minislump-debug/${input:outputDirectory}"
    ]

Remove these extra flags to return to the full default workload. In particular,
``--settle`` means **skip settling**, despite its name.

The launch configuration uses a local ``pipeTransport`` through ``/bin/sh -c``
to start gdb directly. Program output appears in the Debug Console; this demo
does not require interactive stdin. This bypasses the integrated-terminal
handoff and its temporary ``Microsoft-MIEngine-*`` named pipes. That handoff
failed in this workspace with ``AggregateException: ... Cannot access a closed
file``. The local pipe configuration was checked using the installed C/C++
debug adapter: it stopped at a source location and completed ``--help`` normally.
The exact cause of the original adapter exception has not been established.

If that adapter error appears again, stop the failed session and select
**MiniSlump: launch** before pressing F5. Use a fresh output directory. For
additional adapter diagnostics, temporarily add ``"logging": {"engineLogging":
true}`` to this configuration and retain the Debug Console output preceding
the fatal error. The existing generic **(gdb) Launch** entry is unchanged.

Repeated runs without gdb
-------------------------

For repeated runs, run ``repeat.sh`` in a terminal. Each iteration launches a
new process, prints its PID and writes its log and output in a unique directory
under ``build-minislump-debug/runs``. It stops on the first nonzero exit and
leaves all evidence intact. ``RUNS`` defaults to 100; ``MINISLUMP_BINARY`` can
select another build. Additional arguments are forwarded to the demo::

    RUNS=10 scripts/minislump/repeat.sh --verify-cone-release
    RUNS=100 scripts/minislump/repeat.sh --settle --slump-duration 0.001

The second command is a short flow smoke test; it does not exercise settling.
The default full run is needed to reproduce the original workload. The script
chooses each run's output directory itself; do not supply ``--out-dir`` to it.
It runs sequentially and waits for one process to exit before starting the next.
If one process stalls, the repetition counter will not advance.

Each printed run directory contains ``pid``, ``run.log`` (stdout and stderr),
and ``output/``. To watch the active run, substitute its printed log path::

    tail -f build-minislump-debug/runs/repeat-XXXXXXXX/run-1/run.log

The script retains every run's files, so full repetitions accumulate output.
``MINISLUMP_BINARY`` changes the executable but does not change this log root.
These runs are not automatically under gdb, and starting **MiniSlump: launch**
creates a separate process rather than attaching to the active repetition.

Inspecting a suspected stall
---------------------------

If progress stops, select **MiniSlump: attach**, choose the printed PID, then
press Pause. Inspect all threads in Call Stack. In the Debug Console run::

    -exec info threads
    -exec thread apply all bt full

Save a trace to a file with::

    -exec set logging file /tmp/minislump-gdb.txt
    -exec set logging enabled on
    -exec thread apply all bt full
    -exec set logging enabled off

Useful breakpoints include ``deme::DEMSolver::DoDynamics``,
``deme::DEMSolver::DoDynamicsThenSync``, and the demo's cone-release and final
solver-destruction lines. Pause first to capture the natural stall; frequent
breakpoints can change its timing. If Linux ptrace policy prevents sibling
attachment, launch under gdb via F5 instead; no global ptrace setting change is
needed. Attachment was rejected in the development environment with
``ptrace: Inappropriate ioctl for device``, while launching under gdb worked.
Thus attaching to a stalled repetition may not be available on this machine.
For that case, use repeated manual F5 launches; automatic repeated execution
under gdb is not implemented by these scripts.

Ctrl-C in the repeat-script terminal terminates its active child.
There is deliberately no timeout that kills a stalled process before inspection.

Command-line gdb fallback
-------------------------

To debug one run without the VS Code adapter, run from the repository root::

    gdb --args build-minislump-debug/bin/DEMdemo_DFCMiniSlump --out-dir build-minislump-debug/OUT_gdb_1

At the gdb prompt use the following commands. In command-line gdb, omit the
``-exec`` prefix used by the VS Code Debug Console::

    break main
    run
    continue

If progress stops, press Ctrl-C inside gdb, then capture all thread stacks::

    info threads
    set logging file /tmp/minislump-gdb.txt
    set logging enabled on
    thread apply all bt full
    set logging enabled off

Use ``continue`` to resume or ``quit`` to end the debugging session. This also
runs only one simulation. For another run in the same gdb session, first use
``set args --out-dir build-minislump-debug/OUT_gdb_2`` and then ``run``; each run
needs a fresh output directory. Any additional demo flags must also be included
in ``set args`` because it replaces the previous argument list.

Diagnostic and smoke-test modes
-------------------------------

Each example below uses a different output directory. If repeating an example,
choose another name. List all accepted flags with::

    build-minislump-debug/bin/DEMdemo_DFCMiniSlump --help

Record initial owner state and first-step contact forces, then exit::

    build-minislump-debug/bin/DEMdemo_DFCMiniSlump --initial-force-diagnostic --out-dir build-minislump-debug/initial-forces-1

Also record forces on the first and refreshed steps after cone/plate release::

    build-minislump-debug/bin/DEMdemo_DFCMiniSlump --verify-cone-release --out-dir build-minislump-debug/cone-release-1

Exercise settling through the 1.5-second release boundary, with no flow phase::

    build-minislump-debug/bin/DEMdemo_DFCMiniSlump --slump-duration 0 --out-dir build-minislump-debug/settling-1

Run a fixed-step sensitivity test for 0.1 seconds while retaining confinement::

    build-minislump-debug/bin/DEMdemo_DFCMiniSlump --sensitivity-dt 1e-6 --out-dir build-minislump-debug/sensitivity-1

The other accepted sensitivity timestep is ``5e-7``. Sensitivity and force
diagnostics require the cone and settling to be enabled; do not combine them
with ``--no-cone`` or ``--settle``. Sensitivity and force diagnostics are also
mutually exclusive. ``--verify-cone-release`` implies the initial-force mode.
For normal runs, ``--frame-steps`` controls the output interval (default 800);
``--frame-time`` overrides it with a simulated-time interval.

Evidence from validation
------------------------

Release and Debug targets built successfully. Two short fresh-process Release
runs produced eight frames each with finite final positions and the expected
Python timestep schedule. The full settling check completed 33,072 steps and
42 frames, releasing at exactly 1.5 seconds. The cone-release diagnostic showed
237 active container contacts initially and zero after release. Gdb source
breakpoints and the direct-pipe cppdbg adapter launch were checked successfully.
The full 100-run default workload has not been validated; these checks do not
establish that the original intermittent stall is absent.

Matching behavior and deliberate differences
-------------------------------------------

The port retains coated-sphere mass normalization, -90 degree X rotation,
initial floor overlaps, fixed floor/container families, the top plane at 50 mm,
CUDA debug synchronization, disabled adaptive bins and JIT mass/template flags.
It calls DoDynamics once per requested integration step, with Python's settling
and flow timestep ramps and output synchronization points. Cone/plate contacts
are disabled at exactly 1.5 seconds, independent of kinetic energy. It retains
initial-force/cone-release and sensitivity diagnostics, supported-layer height
and radius measurements, translational/total energy, bottom reaction, frame
manifest, final positions and first-50-particle trajectories.

All Python CLI options are accepted (``--help`` lists them). In particular,
``--settle``/``-s`` means *skip* settling, just as in Python. Sensitivity mode
uses a fixed 1e-6 or 5e-7 step for 0.1 seconds, with 0.001-second frames.

Output defaults to ``OUT_MiniSlump_CPP`` and must be empty; unlike Python this
demo never recursively deletes an existing output directory. Input assets are
resolved from the source directory embedded at build time, independent of cwd.
``--particle_file`` overrides the bundled particle file. Run the build/repeat
scripts to validate asset hashes; direct executable invocation does not hash
inputs. Floating-point reductions and CSV formatting may differ slightly from
Python's ``math.fsum`` and formatting, so comparisons should use tolerances.
Console messages focus on initialization, frames, step progress and teardown;
the simulation CSV schema and timestep schedule are preserved.
