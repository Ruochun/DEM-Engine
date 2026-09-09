"""Compare completed default/Feng hopper reports; report differences without asserting physical equivalence.

Usage: python3 compare_hopper_cases.py DEFAULT_DIR FENG_DIR --output-dir COMPARISON_DIR [--plot]
CSV/JSON output uses only the standard library. --plot additionally requires matplotlib.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def load(directory, geometry):
    """Reject partial, mislabelled or internally inconsistent outputs before comparing the two trajectories."""
    metadata = read(directory / 'run.csv')
    require(len(metadata) == 1, 'Expected one run metadata row')
    meta = metadata[0]
    require(meta['schema'] == '2' and meta['geometry'] == geometry, 'Wrong report schema or geometry flavor')
    require(meta['force_model'] == 'frictional_hertzian', 'Expected the default Hertzian force model')
    initial, final, metrics = [read(directory / name) for name in ('initial.csv', 'final.csv', 'metrics.csv')]
    expected = int(meta['spheres']) + int(meta['cylinders'])
    require(len(initial) == len(final) == expected, 'Incomplete particle snapshots')
    initial_by_id = {r['owner']: r for r in initial}
    require(len(initial_by_id) == len({r['owner'] for r in final}) == expected, 'Duplicate owner IDs')
    require(set(initial_by_id) == {r['owner'] for r in final}, 'Owner identities changed')
    for row in final:
        require(all(row[k] == initial_by_id[row['owner']][k] for k in ('species', 'mass', 'ix', 'iy', 'iz')),
                'Particle identity or mass properties changed')
    for row in initial + final:
        require(float(row['mass']) > 0, 'Nonpositive particle mass')
        require(math.isclose(sum(float(row[k])**2 for k in ('qx', 'qy', 'qz', 'qw')), 1., abs_tol=1e-4),
                'Particle quaternion is not normalized')
    for row in initial + final + metrics:
        for key, value in row.items():
            if key not in {'species', 'phase'}:
                require(math.isfinite(float(value)), f'Non-finite {key}')
    require(len(metrics) >= 4 and all(r['phase'] == 'final' for r in metrics[-2:]), 'Run has no final sample')
    expected_time = float(meta['settling_s']) + float(meta['discharge_s'])
    require(math.isclose(float(metrics[-1]['time']), expected_time, abs_tol=2*float(meta['dt'])), 'Incomplete duration')
    for species, count in [('sphere', int(meta['spheres'])), ('cylinder', int(meta['cylinders']))]:
        require(sum(r['species'] == species for r in initial) == count, 'Incorrect species population')
        mass = sum(float(r['mass']) for r in initial if r['species'] == species)
        previous = 0
        for row in (r for r in metrics if r['species'] == species):
            require(int(row['count']) == count, 'Population changed')
            require(all(float(row[k]) >= 0 for k in ('translational_ke', 'rotational_ke', 'max_speed')),
                    'Negative energy or speed')
            if species == 'cylinder':
                axes = [float(row[k]) for k in ('axis_x2', 'axis_y2', 'axis_z2')]
                require(all(-1e-6 <= x <= 1 + 1e-6 for x in axes) and math.isclose(sum(axes), 1., abs_tol=1e-4),
                        'Invalid cylinder orientation statistics')
            ever = int(row['ever_below_plane'])
            require(previous <= ever <= count and 0 <= int(row['below_plane']) <= ever, 'Invalid discharge counts')
            previous = ever
            require(math.isclose(float(row['discharged_mass']) + float(row['remaining_mass']), mass, rel_tol=1e-9),
                    'Mass accounting failed')
    coverage = read(directory / 'coverage.csv')
    totals = {}
    expected_samples = {r['sample'] for r in metrics if r['phase'] in {'settling', 'discharge'}}
    if meta['diagnostics'] == '1':
        require(len(coverage) == 4 * len(expected_samples), 'Incomplete coverage report')
        require(len({(r['sample'], r['category']) for r in coverage}) == len(coverage), 'Duplicate coverage rows')
    else:
        require(not coverage, 'Unexpected diagnostic records')
    for row in coverage:
        for key, value in row.items():
            if key not in {'phase', 'category'}:
                require(math.isfinite(float(value)), f'Non-finite coverage {key}')
        require(row['sample'] in expected_samples, 'Coverage without a state sample')
        active, eligible, used = [int(row[k]) for k in ('active_patches', 'force_eligible', 'used_feng')]
        require(0 <= used <= eligible <= active <= int(row['patches']), 'Invalid force selection counts')
        require(used == (eligible if geometry == 'feng' else 0), 'Actual geometry use differs from selection')
        require(0 <= float(row['used_elastic_load']) <= float(row['legacy_elastic_load']) + 1e-12,
                'Invalid selected elastic-load coverage')
        rejected = sum(int(value) for key, value in row.items() if key.startswith('fallback_'))
        require(rejected + eligible == int(row['patches']), 'Fallback reasons do not partition patches')
        if row['category'] == 'cylinder_gate' and row['phase'] == 'discharge':
            require(active == 0, 'Disabled gate still has active mesh contacts')
        total = totals.setdefault(row['category'], {k: 0. for k in
            ('patches', 'active_patches', 'force_eligible', 'used_feng', 'legacy_elastic_load', 'used_elastic_load')})
        for key in total:
            total[key] += float(row[key])
    for total in totals.values():
        load = total['legacy_elastic_load']
        total['used_elastic_load_fraction'] = total['used_elastic_load'] / load if load else None
    if geometry == 'feng':
        require(sum(t['used_feng'] for t in totals.values()) > 0, 'No sampled Feng use: this is not a two-flavor test')
    return meta, initial, metrics, totals


def plot_curves(a, b, output):
    """Plot bulk observables; particle CSV/VTK frame sequences remain available for ParaView motion comparisons."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 2, figsize=(12, 11), constrained_layout=True)
    for ax, key, label in zip(axes.flat,
                             ('ever_below_plane', 'discharged_mass', 'mean_z', 'translational_ke', 'rotational_ke', 'axis_z2'),
                             ('Discharged count', 'Discharged mass [kg]', 'Mean height [m]',
                              'Translational KE [J]', 'Rotational KE [J]', 'Cylinder axis <z²>')):
        for rows, flavor, style in ((a, 'default', '-'), (b, 'Feng', '--')):
            for species, color in (('sphere', 'tab:blue'), ('cylinder', 'tab:orange')):
                if key == 'axis_z2' and species == 'sphere':
                    continue
                selected = [r for r in rows if r['species'] == species]
                ax.plot([float(r['time']) for r in selected], [float(r[key]) for r in selected],
                        linestyle=style, color=color, label=f'{flavor} {species}')
        ax.set(xlabel='Time [s]', ylabel=label)
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.savefig(output / 'comparison.png', dpi=180)
    fig.savefig(output / 'comparison.svg')
    plt.close(fig)


def main():
    """Match configurations and sample identities, then export signed differences and optional bulk plots."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('default_dir', type=Path)
    parser.add_argument('feng_dir', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    a, b = load(args.default_dir, 'default'), load(args.feng_dir, 'feng')
    comparable = lambda m: {k: v for k, v in m.items() if k not in {'geometry', 'diagnostics'}}
    require(comparable(a[0]) == comparable(b[0]), 'Simulation configurations differ')
    require(a[1] == b[1], 'Initial owner states differ')
    require(len(a[2]) == len(b[2]), 'Sample counts differ')
    differences, bulk = [], {}
    for x, y in zip(a[2], b[2]):
        for key in ('sample', 'time', 'phase', 'species', 'count'):
            require(x[key] == y[key], f'Sample mismatch: {key}')
        row = {key: x[key] for key in ('sample', 'time', 'phase', 'species')}
        for key in x:
            if key in row or key in {'count', 'wall_seconds'}:
                continue
            delta = float(y[key]) - float(x[key])
            row[key + '_feng_minus_default'] = delta
            entry = bulk.setdefault(x['species'] + '/' + key, {'max_abs_difference': 0., 'sum_square': 0., 'samples': 0})
            entry['max_abs_difference'] = max(entry['max_abs_difference'], abs(delta))
            entry['sum_square'] += delta * delta
            entry['samples'] += 1
            entry['final_default'], entry['final_feng'] = float(x[key]), float(y[key])
        differences.append(row)
    for entry in bulk.values():
        entry['rms_difference'] = math.sqrt(entry.pop('sum_square') / entry['samples'])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in ('comparison.json', 'differences.csv', 'comparison.png', 'comparison.svg'):
        require(not (args.output_dir / name).exists(), 'Choose a fresh comparison output directory')
    summary = {'claim': 'Completed matched two-flavor runs with measured Feng use; physical equivalence is not asserted',
               'bulk_differences': bulk, 'feng_sampled_coverage': b[3],
               'wall_seconds': {flavor: float(data[2][-1]['wall_seconds']) for flavor, data in [('default', a), ('feng', b)]}}
    (args.output_dir / 'comparison.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (args.output_dir / 'differences.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(differences[0]))
        writer.writeheader()
        writer.writerows(differences)
    if args.plot:
        plot_curves(a[2], b[2], args.output_dir)
    print(f'PASS: matched reports and nonzero sampled Feng use. Differences: {args.output_dir}')


if __name__ == '__main__':
    main()
