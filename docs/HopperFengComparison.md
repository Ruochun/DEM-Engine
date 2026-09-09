# Hopper comparison reports

For the implemented force selector and two-flavor commands, see [Feng mesh geometry](FengMeshContactGeometry.md).
The following phase-2 record describes the diagnostic-only baseline; the current demo also accepts `--geometry=feng`.

This phase adds a reproducible measurement harness to `DEMdemo_HopperSphereMeshedCylinder`. It compares **default
forces with diagnostics off/on** and measures how often the current Feng diagnostic gates pass in a real hopper.
That phase did not implement or validate Feng forces. The force selector was added in the subsequent phase.

## Run and verify

```sh
cmake --build build --target DEMdemo_HopperSphereMeshedCylinder -- -j2
python3 src/demo/ModularTests/test_hopper_comparison.py \
  build/bin/DEMdemo_HopperSphereMeshedCylinder --output-root /tmp/hopper-phase2-check
```

Choose a new output root for each execution. The verifier runs two default smoke cases and one diagnostics-on smoke
case, each with 24 cylinders, 16 spheres, 0.16 seconds settling and 0.14 seconds discharge. It uses `--no-frames` and
`--fixed-cd` in all three runs. This selects `SetCDUpdateFreq(10)` and disables adaptive update frequency and bin sizing;
the timestep remains 5 microseconds. This is an opt-in verification setting, not a change to the normal demo.
Each child has a 900-second timeout. GPU access is required; a solver failure is a test failure, not a skipped pass.

The output root contains each run's CSVs, process logs, and `verification.json`. Existing comparison files are not
overwritten. To check saved outputs without rerunning the GPU simulations:

```sh
python3 src/demo/ModularTests/test_hopper_comparison.py \
  build/bin/DEMdemo_HopperSphereMeshedCylinder --output-root /tmp/hopper-phase2-check --validate-only
```

The phase passes only when:

1. Both targets of measurement share identical initial owner states and recorded configuration, apart from the
   diagnostic switch. Counts, mass accounting, finite states, energy signs and cylinder-axis normalization pass.
2. Every advanced sample has exactly four diagnostic-category rows; rejection reasons partition candidate patches.
   Force timestamps precede state timestamps by one timestep. Sampled pair runs account for all observed owner pairs.
3. The diagnostics-on hopper actually exercises mesh candidates, true surface intersections and positive legacy
   Hertzian elastic loads. **Zero eligible Feng contacts is a valid measurement**, not a manufactured pass for a
   force replacement; that result identifies work for the eligibility phase.
4. Diagnostics-on final positions, velocities, local angular velocities and quaternions agree within an explicit
   absolute floor or three times the measured default repeat difference, whichever is larger. The verifier also
   checks all sampled bulk metrics, excluding wall time, using a uniform envelope per species/observable derived
   from the maximum repeat difference over its whole curve. This avoids demanding exact equality where two noisy
   reference curves happen to cross. These envelopes test non-interference in this smoke case;
   they are not acceptance tolerances for a future Feng force model. Both repeat errors and diagnostic deviations
   are retained in `verification.json`.

## Individual and larger runs

```sh
# Default-force baseline, normal demo scheduler; preserve CSV/VTK visualization frames.
./build/bin/DEMdemo_HopperSphereMeshedCylinder --smoke-test --comparison-report \
  --geometry=default --output-dir /tmp/hopper-default

# Same settings, with Feng diagnostic computation and coverage reporting.
./build/bin/DEMdemo_HopperSphereMeshedCylinder --smoke-test --feng-diagnostics \
  --geometry=default --output-dir /tmp/hopper-diagnostics
```

Omit `--smoke-test` for the original 5,250-cylinder / 3,500-sphere setup and 0.70/7.50-second phases. Both runs must use
the same flags other than diagnostics and output directory. `--no-frames` suppresses visualization files without
changing the output/advance cadence. With no new flags the existing visualization demo remains available.
Automatic report directories end in `_default` or `_feng_diagnostics`, plus the existing `_smoke` suffix when used.

## Report interpretation

| File | Meaning |
| --- | --- |
| `run.csv` | Schema, actual force flavor, diagnostic switch, timestep/cadence, phases, counts, outlet plane and elastic material values. |
| `initial.csv`, `final.csv` | Every mobile owner's position, velocity, local angular velocity, quaternion, mass and principal MOI. |
| `metrics.csv` | Per-species sampled discharge/remaining mass, mean height, translational/rotational kinetic energy, maximum speed; cylinder world-axis second moments. |
| `coverage.csv` | Candidate/gate counts and elastic-load coverage, grouped by contact category at each advanced sample. |
| `pair_runs.csv` | Consecutive sampled candidate-owner-pair observations, observed span and end-of-run censoring. |

State metrics are sampled at the existing 0.01-second output cadence, plus initial, gate-opening and final states.
Discharge counts the first **observed particle center** below `z = -0.04 m`. Remaining mass means mass not yet observed
crossing that plane; it is not a geometric bed-volume estimate. A particle crossing and returning between samples
can be missed. `below_plane` separately measures instantaneous observed occupancy below the plane. Exact crossing
events, fine-grained contact lifetimes and time-integrated contact loads require finer or in-loop instrumentation.

Coverage samples use the last force-evaluation snapshot of each advance, before its integration. The initial and
gate-opening output events do not re-sample stale geometry. Gate classification uses owner identity, so changing its
family does not relabel earlier gate contacts. Categories are cylinder–cylinder, cylinder–hopper, cylinder–gate, and
other. Rejection precedence is non-watertight, ambiguous pair(s), no segment, open boundary, degenerate line, eligible.
`eligible` only means the current diagnostic gates passed; it is **not certified boundary completeness**.

An active patch here has positive legacy area and penetration. Its weighting is the elastic normal Hertzian term

```text
elastic_load = (4/3) E_effective sqrt(legacy_area / pi) legacy_penetration
```

The effective modulus uses the cylinder/flume materials supplied by this demo. Weighting excludes normal damping,
tangential force and rolling resistance. `eligible_elastic_load / legacy_elastic_load` therefore measures sampled
elastic-load coverage, not total-force coverage. No division is made when the denominator is zero; JSON reports
`null` for an undefined fraction. Summed loads across samples are not a force integral. Area ratios, normal cosines
and contact-point distances are summed only for active patches that also have a Feng line; divide by `comparable`
to obtain their sample mean.

Pair runs include **candidate** presence, including margin-only candidates. They aggregate all patches for an ordered
owner pair; they do not track friction-history islands. The observed span is the time between first/last sightings,
not a physical contact lifetime. A pair missing and returning between samples cannot be distinguished from continuous
presence. Runs still observed at the final sample are right-censored.

Elapsed wall time excludes initialization/JIT prewarming and includes stepping, report extraction and any requested
visualization output. The additional GPU geometry is evaluated every force step, although only output-cadence
snapshots are read. No speedup is expected in this diagnostic phase.

## Recorded smoke verification

The phase-2 check on 2026-09-09 completed all three GPU runs and passed the final verifier, including disabled-gate
contact suppression and a negative check that deliberately corrupted mass accounting is rejected. Local artifacts
are in `/tmp/deme-hopper-phase2-fixed`; they are not committed test fixtures.

| Sampled category | Candidate patches | Gate-passing patches | Eligible elastic-load fraction |
| --- | ---: | ---: | ---: |
| Cylinder–cylinder | 448 | 4 | 100% |
| Cylinder–hopper | 345 | 58 | 100% |
| Cylinder–gate | 216 | 4 | 100% |

Across the 30 advanced snapshots there were 381 intersection-segment observations. Rejections comprised 90 ambiguous
cylinder–cylinder candidates and 853 candidates without segments. These are repeated sample observations, not unique
contacts. All sampled positive elastic loads passed the current gates; this does not certify topology, cover every
force step, or establish equivalence of Feng and default forces.

Final sampled discharged counts were 8 spheres / 11 cylinders in each default run and 7 spheres / 11 cylinders with
diagnostics. Stepping/report wall times were approximately 581, 592 and 596 seconds. This is not a performance benchmark.

Default-repeat maximum final position difference was 0.13065 m, versus 0.13006 m for diagnostics-on against default A.
Individual trajectories therefore have substantial variability even with fixed CD settings. Passing the threefold
repeat envelope is only a coarse smoke check; it cannot establish negligible disturbance. Subsequent physics
validation needs controlled isolated contacts, additional repeats and timestep sensitivity, with physics tolerances
chosen independently of these large trajectory differences.

The running verifier initially used a pointwise bulk envelope and rejected an early cylinder-speed sample. The
whole-curve envelope described above was implemented while the simulations were still running, after observing
baseline variability; the completed artifacts were then rechecked with `--validate-only`. Both baseline variation
and diagnostics-on deviations remain available in `verification.json`.

## Gates recorded at the end of phase 2

- **Eligibility and fallback:** validate complete oriented boundaries and multiplicity for these meshes; add fixtures
  with deliberately missing/duplicated segments and incomplete patch grouping, and show correct whole-patch fallback.
- **Feng force mode:** demonstrate positive Hertzian forces, force–penetration behavior, stable torques/friction history,
  and transition behavior in isolated cylinder contacts before exposing `--geometry=feng`.
- **Hopper physics comparison:** run smoke, intermediate and full beds; compare settling/discharge curves and motion
  statistics against repeatability and timestep studies. Report actual Feng usage alongside the physics, so agreement
  caused by universal fallback cannot be mistaken for validation.

These gates extend the [Stage 1 geometry diagnostics](FengMeshContactDiagnostics.md); that reporting phase was committed as `ea74eef`.
