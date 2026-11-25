# Changelog

## 3.2.0

### Breaking changes

- Force models now use global `A_only_force`, `B_only_force`, `A_only_torque`,
  and `B_only_torque` additive outputs. References to `torque_only_force` in
  executable model code or prerequisites produce a setup-time deprecation error.
- Rolling resistance now applies a dissipative, equal-and-opposite direct torque
  based on relative angular velocity. Recalibrate `Crr` when migrating from the
  former force-based rolling law.
- Contact force CSV output now names both contact sides explicitly: `f_a_x`, `f_a_y`,
  `f_a_z`, `f_b_x`, `f_b_y`, `f_b_z`, and `torque_a_x`, `torque_a_y`, `torque_a_z`,
  `torque_b_x`, `torque_b_y`, `torque_b_z`. These replace `f_x/y/z`, `torque_x/y/z`,
  `ForceB_X/Y/Z`, and `TorqueB_X/Y/Z`; update downstream column lookups. Legacy CSV
  aliases are not emitted or reserved by the wildcard reader. Torque columns store
  direct extra torques, excluding moments from forces.

### API clarification

- The frictional Hertzian model supports dissipative spinning friction through
  `Cspin` (default zero). Rolling and spinning have no rotational spring history.
- Rolling length bounds default to the intrinsic effective radius with multiplier
  one. `SetRollingResistanceLengthLimit` configures this multiplier and an optional
  absolute cap; `SetRollingResistanceLengthBoundsEnabled(false)` disables both.
  Nonfinite lengths still error. Apply runtime changes with a full `Update()`.
- Mesh `SetContactRadius`, `GetContactRadius`, and `ClearContactRadius` control the
  intrinsic size used by rolling bounds. Universal contact derives an unset radius
  from closed-mesh volume at initialization/full update; legacy contact treats
  meshes as walls. Open/failed meshes retain zero curvature. Explicit radii take
  precedence, and lightweight deformation does not recompute them.
- Detailed contact output exposes `GetForceA/B` and `GetTorqueA/B`. Existing
  `GetForce` and `GetTorque` methods continue to select side A.

### Contact processing fixes

- Complete primitive geometry, patch aggregation, and force evaluation for each
  contact type before advancing to the next type. This prevents compacted force
  writes from overwriting unconsumed mesh primitive geometry in mixed scenes.
- Reuse B-force storage for temporary primitive contact points, removing one
  `float3` scratch array per allocated contact slot. Built-in force laws are
  unchanged; custom models must not rely on the former cross-type ordering of
  writes to shared owner wildcards.

## 3.1.0

Changes since 3.0.19, the final 3.0.x release.

### Synchronized stepping

- Added `SetSynchronizedMode()` in C++ and Python, combining
  `SetCDUpdateFreq(0)` with `DisableAdaptiveUpdateFreq()`.
- Enforced same-state contact detection before every force evaluation in
  synchronized mode, including across consecutive `DoDynamics()` calls. Both
  workers remain active, but dynamics no longer advances using stale contacts.
- Documented setup and runtime mode changes. The mode setter does not wait for
  outstanding work or replace `DoDynamicsThenSync()` for host synchronization.

### Fixes and validation

- Fixed lost wakeups in dynamics/contact-detection worker handoffs that could
  leave simulations waiting indefinitely.
- Corrected contact-detection lag statistics to measure the age of valid results
  when consumed; synchronized mode reports zero lag.
- Added thread-handoff regressions covering repeated runs and zero-drift stepping.

### Build, examples, and documentation

- Disabled the native visualizer by default; enable it explicitly with
  `-DDEME_BUILD_VISUALIZER=ON`. Updated build and visualization instructions.
- Configured the meshed-cylinder hopper demo for synchronized stepping and
  finer, 32-sided cylinders, and revised its material settings.
- Fixed the MixerCubes demo's contact-detection update frequency and shortened
  its default simulation duration.
- Expanded synchronized-mode documentation, refreshed the Python API reference,
  and moved documentation maintenance guidance into the documentation README.
