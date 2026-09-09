# Experimental Feng mesh geometry

The solver now offers two mesh–mesh geometry flavors. Both use the selected DEME force law, including the hopper's
default frictional Hertzian law. `default` uses the existing projection/voting geometry. `feng` uses intersection
boundary area, B2A normal and the point on the Feng contact line nearest the legacy point, while retaining legacy
penetration. The option defaults to `default` and must be set before `Initialize()`:

```cpp
solver.SetMeshUniversalContact(true);
solver.SetMeshMeshContactGeometry("feng");  // or "default"
```

Python exposes the same selector. This is the proposed geometry/legacy-depth hybrid, not Feng's original
contact-volume constitutive law; no exact energy conservation or speedup is claimed.

## Run the two hopper cases

Build the demo and focused checks:

```sh
cmake --build build --target DEMTest_JitDataLayout DEMTest_FengGeometry DEMTest_FengDiagnostics DEMTest_FengForces \
  DEMdemo_HopperSphereMeshedCylinder -- -j2
./build/bin/DEMTest_FengGeometry
./build/bin/DEMTest_FengDiagnostics
./build/bin/DEMTest_FengForces
./build/bin/DEMTest_JitDataLayout
```

Start with matching smoke runs, preserving sphere CSV and cylinder/hopper VTK frame sequences for ParaView:

```sh
./build/bin/DEMdemo_HopperSphereMeshedCylinder --smoke-test --comparison-report \
  --geometry=default --output-dir /tmp/hopper-default
./build/bin/DEMdemo_HopperSphereMeshedCylinder --smoke-test --comparison-report \
  --geometry=feng --output-dir /tmp/hopper-feng
python3 src/demo/ModularTests/compare_hopper_cases.py /tmp/hopper-default /tmp/hopper-feng \
  --output-dir /tmp/hopper-comparison
```

Remove `--smoke-test` from **both** commands for the original 5,250-cylinder / 3,500-sphere bed. Select fresh output
directories. Both runs use the same deterministic initial placement, timestep, materials, gate event and output
cadence. `--no-frames` suppresses visualization output if only CSV statistics are needed. `--fixed-cd` is optional;
use it on both runs or neither. It fixes maximum CD cadence at 10 steps and disables adaptive cadence/bin sizing,
but does not make the asynchronous solver bitwise deterministic.

The comparison script writes `comparison.json` and `differences.csv`, verifies that configuration/initial states
match and that the Feng run actually used Feng geometry, then reports differences without imposing a physics pass
threshold. Add `--plot` in a Python environment with matplotlib to produce `comparison.png` and `comparison.svg`:
discharge counts/mass, mean height, kinetic energies and cylinder orientation statistics. Open each run's
`DEMdemo_mesh_*.vtk` and `DEMdemo_spheres_*.csv` frame series for particle-motion comparisons at matching frame numbers.

The [report guide](HopperFengComparison.md) explains the sampling conventions. Compare curves across repeat runs and
timestep refinements, not just individual trajectories: phase-2 default repeats already showed substantial particle
trajectory divergence. Matching appearance or discharge rate alone does not establish the same physics.

## Eligibility and whole-patch fallback

The selector is deliberately restricted to the hopper's class of small rigid convex solids. Setup validation checks
actual triangle geometry, rather than trusting `SetConvex(true)`:

- One mesh/material patch, 4–256 triangles, no shell thickness.
- Finite nondegenerate triangles, outward convex face half-spaces and positive volume.
- Every directed edge has exactly one reverse edge after welding exactly coincident coordinates, with sphere topology.
- Simple patch combination, so all candidates for an owner pair share a reduction group.

At each active legacy contact the original diagnostic gates must pass, followed by an endpoint audit requiring
exactly one predecessor/successor per segment and one cycle containing every segment. This rejects missing edges,
duplicate loops, cancelling open chains, branches, reversed edges and disconnected loops. Groups above 512 primitive
candidates fall back to bound the quadratic audit. The Feng normal must also lie in the same hemisphere as the legacy
B2A normal. Candidate completeness still depends on the existing broad-phase/CD infrastructure: this adds no
independent all-triangle-pairs search. Setup tolerances are conservative and may reject otherwise usable meshes.

Only then are area, normal and point replaced together. Otherwise **all three remain legacy**. Legacy penetration,
contact ordering, patch identities, friction history and the Hertz force model are unchanged. The Hertz model already
projects tangential displacement onto the current contact plane. Selection may switch between steps as intersections
appear or degenerate; there is no blending or smoothing, so sensitivity at transitions remains part of the requested
physical comparison. Margin-only contact, containment without a surface boundary, coplanar ambiguity, concave shapes,
multiple patches, shells and flooded islands use legacy geometry. Mesh node edits invalidate the solid validation
until a fresh initialization; subsequent contacts use fallback.

The geometry still performs legacy projection/voting, then adds FP64 segment reductions and the bounded endpoint
audit. Depth-only optimization, general concave/multiple-loop support and performance tuning remain future work.

## Actual use, fallback reasons and verification

`GetMeshMeshFengDiagnostics()` returns a snapshot in Feng force mode even if `SetMeshMeshFengDiagnostics(false)`.
The original diagnostic fields remain available. Additional fields distinguish diagnostic success from selection:

- `ownersValidated`, `boundaryValidated`: stronger solid/boundary checks.
- `fengEligible`: an active legacy patch passed every force-selection check.
- `usedFeng`: this force evaluation actually consumed Feng geometry.
- `fallbackReason`: 0 eligible; 1 inactive legacy; 2 unsupported/invalidated solid; 3 unsupported grouping;
  4 diagnostic gate; 5 boundary audit/candidate limit; 6 inconsistent normal.

Hopper schema 2 adds `force_eligible`, `used_feng`, `used_elastic_load`, and the six fallback counters to `coverage.csv`.
`used_elastic_load` is the **legacy elastic-load proxy** on selected patches, not a measurement of the new force or
total force including damping/friction. The comparison JSON reports that coverage by category. Snapshots are taken
every 0.01 seconds, so these are sampled observations rather than every-step or time-integrated statistics.

Verification is staged:

1. CPU geometry/eligibility fixtures: analytic area/moment, refinement/diagonals, owner swap, topology counterexamples,
   inward/duplicated/concave solids and incomplete patches.
2. GPU diagnostics: existing mixed contact-key indexing, lifecycle and unchanged diagnostic-only motion.
3. GPU forces: positive Hertz response at two penetrations, force/torque against the selected geometry (including a distinguishable lever arm), action/reaction,
   unchanged sphere contacts, tangential-history accumulation, and exact default-force fallback on ambiguous/partial
   patches. Explicit deformation checks certificate invalidation.
4. Hopper: matched report validation, finite state and mass accounting, disabled-gate suppression, measured nonzero
   Feng use, plus visual/statistical differences for user evaluation. Passing these checks does not assert equivalence.

The JIT cache key now also includes the contents of `DEM/Defines.h` and `DEM/VariableTypes.h`. Kernel source text alone
does not capture changes to these included data layouts. The cache regression changes a header in place, verifies
the newly executed kernel, then verifies reuse with unchanged contents. Other included helper-header changes still
require care with cache refresh; this change specifically protects the shared data-layout contract.

Run timing includes stepping, reporting and frame output after initialization. Run cases sequentially on an otherwise
idle machine for timing comparisons; concurrent verification/build activity makes smoke-run timings unsuitable as
performance benchmarks.

## Verified smoke result (2026-09-09)

The focused build and all four modular checks above passed, including the distinguishable torque lever arm,
two-step Hertz friction history, automatic force-mode snapshots and JIT header-cache invalidation. Python bindings
passed a C++ syntax check. The comparison script also rejected a default report relabelled as Feng without actual use.

Matched `--smoke-test --fixed-cd` runs completed with 24 cylinders, 16 spheres and 0.30 seconds of simulated time.
Both passed report consistency, finite state, quaternion/axis normalization, energy, population and mass checks.
Each produced 32 finite sphere CSV / mesh VTK frames; disabling the gate removed its 12 triangles from visualization
and left no active gate contacts in discharge samples. Initial states matched exactly.

The Feng report contained 68 active sampled mesh contacts, all selected: 2 cylinder–cylinder, 16 cylinder–gate and
50 cylinder–hopper observations. The other 921 sampled candidates used inactive-legacy fallback. These are sampled
observations, not counts of unique contacts or all time steps. At the final sample, both runs had discharged 9 spheres;
the default run had discharged 11 cylinders and Feng 12. This single short pair does not establish physical equivalence.

Artifacts for this session are under `/tmp/deme-hopper-feng-final/{default,feng,comparison}`. The comparison directory
contains CSV/JSON differences and PNG/SVG plots. Recorded stepping/reporting times were 483.7 seconds (default) and
2861.8 seconds (Feng). Runs and build activity overlapped, so these are not controlled benchmark timings. The added
geometry/reduction/boundary work is not performance optimized; use the smoke cases before committing to large runs.
