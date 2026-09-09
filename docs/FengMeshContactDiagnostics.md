# Mesh–mesh intersection-boundary diagnostics

This documents the original diagnostic phase. The subsequent [experimental force selector](FengMeshContactGeometry.md)
adds stronger eligibility checks and optional whole-patch geometry replacement; the diagnostic switch alone still
does not change force selection.

The proposed geometry is applicable as an **experimental diagnostic**, but is not a safe unconditional replacement
for DEME patch voting. This implements Stage 1 of the supplied `need.md`: evaluate both geometries on the same
candidate/contact topology without changing forces, penetration, friction history, broad phase, or patch assignment.
The option defaults to off. It is not a performance optimization: enabled runs retain all legacy work and add two
CUDA kernels, four FP64 CUB reductions, and temporary storage.

## Usage

```cpp
solver.SetMeshMeshFengDiagnostics(true);  // before Initialize(), or while dynamics is idle
solver.Initialize();
solver.DoDynamics(step);
for (const auto& row : solver.GetMeshMeshFengDiagnostics()) {
    // Compare row.legacyArea with row.fengArea only after examining the diagnostic gates below.
    // row.vectorArea / geometricMoment / boundaryResidual are available even if the gates fail.
}
```

Python exposes the same methods; the getter returns a list of dictionaries with the C++ field names. Vectors are
three-element tuples of Python floats, preserving the double-precision geometry.

The getter synchronizes the dT stream and copies the latest mesh–mesh force-evaluation snapshot. It does not evaluate
new contacts. Coordinates describe the force evaluation **before integration**, rather than the post-step owner pose.
Call while dynamics is idle. Disabling/re-enabling clears the visible snapshot; a force step with no mesh–mesh
candidates also clears it. After moving a body, asynchronous contact detection may retain old candidate records for
one force evaluation; those records still get newly computed geometry at the current pose. Each `patchContact` is the global, transient patch-contact index for that evaluation,
not a stable contact-history identifier. Owner IDs identify the ordered A/B bodies.

## What is calculated

Each existing mesh–mesh primitive candidate is intersected using its actual, zero-thickness surface triangles.
Normals come from vertex winding. Segment direction follows `nA × nB`. Prism thickness, shell thickness and family
margins do not inflate these segments. Legacy projection rejection does not suppress diagnostic candidates: doing
so could discard part of a closed curve. No separate all-pairs search is performed, so completeness still depends
on the existing candidate infrastructure and contact grouping.

Vertices are rotated in FP64 and expressed relative to owner A's center, using the same origin for every primitive
in an owner pair. Four `double3` arrays hold vector area, moment, endpoint displacement, and statistics. Existing
CUB sum-by-key operations reduce them with the same global `geomToPatchMap` keys used by legacy voting. No atomic
accumulation, JIT kernel changes, or new host readback is added to the force loop. Scratch allocation follows the
existing reusable scratch-pool pattern; diagnostic snapshot buffers grow when needed.

The snapshot includes:

- `referenceOrigin`, `vectorArea`, `geometricMoment`, `boundaryResidual`;
- `boundaryLength`, `segmentCount`, `ambiguousPairCount`;
- `legacyArea`, `legacyNormal` (B2A), `legacyPenetration`, `legacyContactPoint`;
- `ownersWatertight`, `closurePassed`, `hasContactLine`;
- gated `fengArea`, `fengNormal` (B2A), `fengContactPoint`.

`closurePassed` checks at least three segments and
`|sum(end - start)| <= 1e-8 * total_segment_length`.
`hasContactLine` additionally requires watertight owners, no ambiguous candidates, finite line construction,
and `|S| > 1e-12 * total_segment_length^2`. Failed gates leave candidate line outputs zero; raw sums remain visible.
Neither flag certifies that the intersection boundary is complete, consistently wound, or free of duplicates.
For example, two unrelated open chains can have cancelling endpoint displacements. A duplicated closed curve also
passes the displacement check. The later force selector requires additional checks; neither of these flags alone authorizes replacement.

Degenerate triangles, coplanar/near-parallel candidates, edge-on-plane ambiguities, point touches, and very short
segments are excluded and counted as ambiguous. Separated planes/intervals contribute nothing. This conservative
policy can suppress otherwise usable diagnostic lines in aligned meshes. Existing geometry remains the force input
for **all** cases, including shell/margin-only contact and containment without surface intersections.

## Corrections needed to the proposal

### Existing patches do not guarantee closed boundaries

`DEMMesh::SetPatchIDs` permits arbitrary surface subsets, and `SetEachTriangleAsPatch` explicitly permits open
one-triangle patches. Flooded contact islands also need not be complete surface-intersection loops. Additivity alone
does not make an arbitrary subset a valid contact boundary.

For a shift of reference origin by `c`, the segment sum transforms as

```text
S' = S - 1/2 c × sum(end - start).
```

An open subset can therefore have a nonzero area that changes with reference origin. Validating individual triangle
pairs or checking only `|S| > tolerance` cannot solve this. Before a force-producing mode, the implementation needs a
validated boundary-completeness/multiplicity policy and a whole-patch fallback policy. Mixing projected scalar areas
and boundary vector-area contributions within one patch has no established geometric interpretation.

### The supplied moment coefficient fails a planar centroid test

For the definition in the attachment, `G = integral r × dS`, Stokes' theorem gives

```text
G = -1/2 integral_boundary |r|^2 dr
Gi = -1/2 [dot(p,q) + dot(q-p,q-p)/3] (q-p).
```

This is an independent derivation for that definition; the attachment's `-1/3` outer coefficient gives `2/3` of G.
For a unit square centered at `(2.5,3.5,5)` with vector area `(0,0,1)`, the correct moment is `(3.5,-2.5,0)`.
The attachment's coefficient puts the transverse line location at `(5/3,7/3)` instead of `(2.5,3.5)`.
The shared host/device helper uses `-1/2`; the analytic test checks against the surface moment, independently of
segment integration.

### Origin and normal conventions must stay coupled

Both S and G use local coordinates. With `n = S/|S|`, compute the local line base `cross(n,G)/|S|`, project the legacy
contact point onto that line, then add the reference origin. This implements the attachment's proposed nearest-line
hybrid, rather than world `lambda=0`, which would change the tangential-force lever arm under a translation along n.
For outward-wound A/B solids the diagnostic B2A normal is **`-n`**. Flipping n for the DEME convention must not also
flip the line base; equivalently, both S and G would have to change sign together. Owner swapping is tested.
Watertightness does not certify outward winding, so reversed-wound meshes require further validation before use.
For a nonplanar boundary, inspect `dot(G,n)` as well: a single force along n reproduces only the perpendicular part
of the geometric moment. These diagnostics do not introduce a separate torsional couple or certify mechanical
equivalence for a general nonplanar contact.

The geometry proposal is motivated by [Feng's contact-volume model](https://doi.org/10.1016/j.cma.2020.113493), whose
published abstract describes additive geometry and parallel surface-mesh evaluation. This implementation does not
claim the original constitutive model's energy conservation for DEME's arbitrary penetration-based force laws.

## Validation and remaining work

Build the focused tests:

```sh
cmake --build build --target DEMTest_FengGeometry DEMTest_FengDiagnostics -- -j2
./build/bin/DEMTest_FengGeometry
./build/bin/DEMTest_FengDiagnostics
```

The geometry test runs on the CPU using the same host/device functions: analytic triangle cuts, endpoint orientation,
owner swap, direct surface moments, translated square contact points, open-patch origin dependence, degenerate and
margin-only cases, and closed cubes with 2/8/32 triangles per face and diagonal flips.

The GPU integration test compares enabled/disabled owner velocities, verifies mixed contact-type key offsets,
checks cube geometry and refinement, exercises open per-triangle patch rejection, and checks snapshot lifecycle.
The initially deeply intersecting cube cases have zero legacy area with the initial candidate margin; their boundary
areas are about 0.940601 for all three resolutions (15, 35, and 73 segments). The velocity comparison includes active
sphere forces, but does not establish active mesh-force/torque equivalence or mechanical improvement.

Before enabling replacement geometry in force calculations, further work includes topology and winding validation,
transition/fallback behavior, sphere refinement and oblique impact studies, disconnected concave contact regions,
friction/torque and energy-drift comparisons, and performance/occupancy measurements. Stage 1 adds work, so no speedup
or large-simulation non-regression is claimed. Depth-only optimization remains a later stage.

Local verification: both focused targets built with CUDA 12.8; the shared CPU geometry tests and GPU integration
checks passed. Python bindings passed a C++ syntax check (the full Python extension was not built). CUDA memory
instrumentation was unavailable on this WSL/WDDM system: Compute Sanitizer reported "Failed to initialize WDDM
debugger interface" and "Device not supported". Its application run passed, but this is not a successful memcheck.
