"""Verify hopper reporting and diagnostic non-interference; retain CSVs/logs for review.

Usage: python3 test_hopper_comparison.py build/bin/DEMdemo_HopperSphereMeshedCylinder --output-root /tmp/hopper-phase2
Recheck existing results with --validate-only. This does not test Feng force equivalence.
"""
import argparse
import csv
import json
import math
import subprocess
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def finite(row, exclude):
    for key, value in row.items():
        if key not in exclude:
            require(math.isfinite(float(value)), f"Non-finite {key}: {row}")


def validate(directory, diagnostics):
    meta = rows(directory / "run.csv")
    require(len(meta) == 1, "Expected one run metadata row")
    meta = meta[0]
    require(meta["geometry"] == "default" and meta["force_model"] == "frictional_hertzian", "Mislabelled force mode")
    require(int(meta["diagnostics"]) == diagnostics and meta["fixed_cd"] == "1" and meta["fixed_cd_steps"] == "10", "Unexpected run configuration")
    initial, final = rows(directory / "initial.csv"), rows(directory / "final.csv")
    require(len(initial) == len(final) == 40, "Smoke owner count")
    for data in (initial, final):
        require(len({r["owner"] for r in data}) == 40, "Duplicate owner state")
        for row in data:
            finite(row, {"species"})
            require(float(row["mass"]) > 0, "Nonpositive mass")
    metrics = rows(directory / "metrics.csv")
    require(len(metrics) == 66, "Expected initial, 30 advanced samples, gate event and final, each with two species")
    previous_count = {"sphere": 0, "cylinder": 0}
    mass_by_species = {species: sum(float(r["mass"]) for r in initial if r["species"] == species)
                       for species in previous_count}
    for row in metrics:
        finite(row, {"phase", "species"})
        species = row["species"]
        expected_count = 16 if species == "sphere" else 24
        require(int(row["count"]) == expected_count, "Particle population changed")
        ever = int(row["ever_below_plane"])
        require(previous_count[species] <= ever <= expected_count, "Discharged count is not monotone/bounded")
        require(0 <= int(row["below_plane"]) <= ever, "Plane counts inconsistent")
        previous_count[species] = ever
        require(math.isclose(float(row["discharged_mass"]) + float(row["remaining_mass"]),
                             mass_by_species[species], rel_tol=1e-10), "Mass accounting failed")
        require(float(row["translational_ke"]) >= 0 and float(row["rotational_ke"]) >= 0, "Negative energy")
        if species == "cylinder":
            require(abs(sum(float(row[k]) for k in ("axis_x2", "axis_y2", "axis_z2")) - 1) < 1e-4,
                    "Cylinder axis second moments do not sum to one")
    require(math.isclose(float(metrics[-1]["time"]), .30, abs_tol=1e-8), "Wrong simulation duration")
    coverage = rows(directory / "coverage.csv")
    pair_runs = rows(directory / "pair_runs.csv")
    if not diagnostics:
        require(not coverage and not pair_runs, "Disabled diagnostics emitted contact observations")
        return meta, initial, final, metrics, {}
    expected_samples = {r["sample"]: r for r in metrics if r["species"] == "cylinder" and r["phase"] in {"settling", "discharge"}}
    require(len(coverage) == 4 * len(expected_samples), "Missing category/sample coverage rows")
    seen = set()
    reasons = ("nonwatertight", "ambiguous", "no_segment", "open_boundary", "degenerate", "eligible")
    totals = {}
    for row in coverage:
        finite(row, {"phase", "category"})
        require(row["category"] in {'cylinder_cylinder', 'cylinder_hopper', 'cylinder_gate', 'other'},
                "Unknown coverage category")
        key = (row["sample"], row["category"])
        require(key not in seen, "Duplicate coverage sample/category")
        seen.add(key)
        state = expected_samples[row["sample"]]
        require(row["phase"] == state["phase"], "Coverage assigned to wrong gate phase")
        require(math.isclose(float(row["force_time"]), float(state["time"]) - float(meta["dt"]), abs_tol=1e-10),
                "Force snapshot timestamp must precede integration")
        require(sum(int(row[k]) for k in reasons) == int(row["patches"]), "Reasons do not partition candidate patches")
        require(0 <= int(row["active_eligible"]) <= int(row["active_patches"]) <= int(row["patches"]), "Active counts inconsistent")
        if row["phase"] == 'discharge' and row["category"] == 'cylinder_gate':
            require(int(row["active_patches"]) == 0, "Disabled gate still has active mesh contacts")
        require(0 <= float(row["eligible_elastic_load"]) <= float(row["legacy_elastic_load"]) + 1e-12, "Load coverage out of range")
        require(0 <= int(row["new_owner_pairs"]) <= int(row["owner_pairs"]), "Pair counts inconsistent")
        category = totals.setdefault(row["category"], {k: 0.0 for k in (*reasons, "patches", "segments", "legacy_elastic_load", "eligible_elastic_load")})
        for name in category:
            category[name] += float(row[name])
    require(sum(r["patches"] for r in totals.values()) > 0, "Test never exercised mesh candidates")
    require(sum(r["segments"] for r in totals.values()) > 0, "Test never exercised surface intersections")
    require(sum(r["legacy_elastic_load"] for r in totals.values()) > 0, "Test never exercised active Hertzian mesh contacts")
    require(bool(pair_runs), "No sampled candidate pair runs")
    for row in pair_runs:
        finite(row, {"category"})
        require(row["right_censored"] in {'0', '1'}, "Invalid pair-run censor flag")
        require(int(row["observations"]) >= 1 and float(row["observed_span"]) >= 0, "Invalid sampled pair span")
        require(math.isclose(float(row["last_force_time"]) - float(row["first_force_time"]),
                             float(row["observed_span"]), abs_tol=1e-10), "Pair timestamps inconsistent")
    require(sum(int(r["observations"]) for r in pair_runs) == sum(int(r["owner_pairs"]) for r in coverage),
            "Pair-run observations do not account for sampled owner pairs")
    for total in totals.values():
        total['eligible_patch_fraction'] = total['eligible'] / total['patches'] if total['patches'] else None
        total['eligible_elastic_load_fraction'] = (total['eligible_elastic_load'] / total['legacy_elastic_load']
                                                  if total['legacy_elastic_load'] else None)
    return meta, initial, final, metrics, totals


def compare(results):
    baseline, repeat, diagnostic = results
    for other in (repeat, diagnostic):
        require({k: v for k, v in baseline[0].items() if k != 'diagnostics'} ==
                {k: v for k, v in other[0].items() if k != 'diagnostics'}, "Run configurations differ")
    require(baseline[1] == repeat[1] == diagnostic[1], "Initial state/material mass properties differ")
    # Report default repeatability and diagnostics-on deviation independently. Floors are physical absolute units;
    # the threefold repeatability envelope is a smoke test rule, not a physical equivalence tolerance.
    differences = {}
    for fields, floor, label in ((('x', 'y', 'z'), 5e-7, 'position_m'),
                                 (('vx', 'vy', 'vz'), 1e-5, 'velocity_m_s'),
                                 (('wx', 'wy', 'wz'), 1e-3, 'angular_velocity_rad_s'),
                                 (('qx', 'qy', 'qz', 'qw'), 1e-5, 'quaternion')):
        errors = []
        for other in (repeat, diagnostic):
            require([r['owner'] for r in baseline[2]] == [r['owner'] for r in other[2]], "Owner ordering changed")
            errors.append(max(math.sqrt(sum((float(a[k]) - float(b[k])) ** 2 for k in fields))
                              for a, b in zip(baseline[2], other[2])))
        tolerance = max(floor, 3 * errors[0])
        differences[label] = {'default_repeat_max': errors[0], 'diagnostic_max': errors[1], 'tolerance': tolerance}
        require(errors[1] <= tolerance, f"Diagnostic non-interference failed for {label}: {differences[label]}")
    # Use one repeatability envelope per species/observable across the whole trace. A pointwise envelope would
    # spuriously demand exact equality wherever two noisy reference trajectories happen to cross.
    bulk = {}
    for a, b, d in zip(baseline[3], repeat[3], diagnostic[3]):
        for key in a:
            if key == 'wall_seconds':
                continue
            if key in {'sample', 'time', 'phase', 'species', 'count'}:
                require(a[key] == b[key] == d[key], f"Mismatched sample identity: {key}")
                continue
            av, bv, dv = map(float, (a[key], b[key], d[key]))
            error = bulk.setdefault(f"{a['species']}/{key}", {'default_repeat_max': 0., 'diagnostic_max': 0., 'scale': 0.})
            error['default_repeat_max'] = max(error['default_repeat_max'], abs(av - bv))
            error['diagnostic_max'] = max(error['diagnostic_max'], abs(av - dv))
            error['scale'] = max(error['scale'], abs(av))
    for key, error in bulk.items():
        error['tolerance'] = max(1e-10, 1e-5 * error.pop('scale'), 3 * error['default_repeat_max'])
        require(error['diagnostic_max'] <= error['tolerance'], f"Bulk trajectory differs for {key}: {error}")
    return {'final_states': differences, 'bulk_curves': bulk}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    executable = args.executable.resolve()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    names = ('default_a', 'default_b', 'diagnostics')
    if not args.validate_only:
        for name in names:
            require(not (root / name).exists(), f"Refusing to overwrite {root / name}")
        # Invalid geometry names must fail before constructing a GPU solver.
        # Loading a freshly linked CUDA executable can be slow even before main on WSL/cold filesystem caches.
        unsupported = subprocess.run([str(executable), '--geometry=invalid'], capture_output=True, text=True, timeout=120)
        require(unsupported.returncode != 0 and 'Unknown' in unsupported.stderr, 'Invalid geometry mode accepted')
        for name in names:
            command = [str(executable), '--smoke-test', '--comparison-report', '--geometry=default', '--no-frames',
                       '--fixed-cd', '--output-dir', str(root / name)]
            if name == 'diagnostics':
                command.append('--feng-diagnostics')
            print(f"Running {name}; log: {root / (name + '.log')}", flush=True)
            with (root / (name + '.log')).open('w') as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=900)
            require(result.returncode == 0, f"{name} failed; see {root / (name + '.log')}")
    results = [validate(root / name, name == 'diagnostics') for name in names]
    summary = {'claim': 'default-force diagnostic repeatability screening, not proof of non-interference or Feng force equivalence',
               'state_comparison': compare(results), 'sampled_coverage': results[2][4],
               'run_wall_seconds': {name: float(result[3][-1]['wall_seconds']) for name, result in zip(names, results)}}
    (root / 'verification.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f"PASS: hopper reports, active mesh coverage, sampled pair accounting and diagnostic repeatability screening. {root / 'verification.json'}")


if __name__ == '__main__':
    main()
