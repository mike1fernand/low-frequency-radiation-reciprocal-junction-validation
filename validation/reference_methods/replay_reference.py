"""Reconstruct the reference taper from three independent sets of saved numerical inputs."""
from __future__ import annotations

import json
import math
from pathlib import Path
import os
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
SURFACE = HERE / 'analysis/resolution_execution/surface'
sys.path.insert(0, str(SURFACE))
sys.path.insert(0, str(HERE / 'analysis/math_resolution/tetra'))

import combined_tetra_diagnostic as tetra
from evaluate_refinement import transmit
from spectral_histogram import load_histogram


def require_close(value, target, tolerance):
    if not math.isclose(value, target, rel_tol=0., abs_tol=tolerance):
        raise ValueError(f'{value!r} disagrees with {target!r} at tolerance {tolerance}')


def main():
    if len(sys.argv) != 2 or sys.flags.optimize:
        raise SystemExit('Usage: python -B replay_reference.py NEW_OUTPUT_DIRECTORY (without optimization)')
    out = Path(sys.argv[1]).resolve()
    if out == HERE or HERE in out.parents:
        raise SystemExit('Choose an output outside the reference_methods directory.')
    out.mkdir(parents=True, exist_ok=False)

    pv = np.mean([
        load_histogram(SURFACE / 'histograms' / f'sobol_seed{s}_p23_nseg10_bins16384.npz')
        .principal_value_jackson(1.2 ** 2, order=256)
        for s in (20260920, 20260921)], axis=0)
    taper = [.1, .2, .2, .2, .2, .1]
    surface_ladder = {}
    tetra_ladder = {}
    for grid, expected_surface in [(128,.12567058),(160,.12566614),(192,.12565926)]:
        with np.load(SURFACE / 'results' / f'omega1.2_N{grid}.npz') as data:
            surface_result = transmit(1.2, pv, [data[f'Im{i}'] for i in range(3)])
        require_close(surface_result['qT'], expected_surface, 5.1e-9)
        require_close(surface_result['Rq'] + surface_result['Tline'] + surface_result['qL'] + surface_result['qT'], 1., 2.e-11)
        surface_ladder[grid] = surface_result
        with np.load(HERE / 'tetra_ladder' / f'N{grid}.npz') as data:
            parts = SimpleNamespace(nseg=6, G0=data['G0'], Im_longitudinal_projector=data['ImL'],
                                    Im_branches=[data[f'Im{i}'] for i in range(3)],
                                    Im_branch_longitudinal_projectors=[data[f'ImBL{i}'] for i in range(3)])
            result = tetra.independent_scatter(1.2, taper, tetra.Parameters(), parts)
        expected = json.loads((HERE / 'tetra_ladder' / f'N{grid}.json').read_text())['scattering']
        keys = ['Rq_amplitude', 'Tq', 'TT_modal', 'TL_modal', 'Tbulk']
        for key in keys:
            require_close(getattr(result, key), expected[key], 2.e-12)
        tetra_ladder[grid] = {key:getattr(result, key) for key in keys}

    expected = json.loads((HERE / 'bloch/raw/generalized_bloch_reference_independent_core_v23.json').read_text())['methods']
    # The historical program writes beside its inputs. Preserve the source by running a copy.
    dest = out / 'bloch'
    shutil.copytree(HERE / 'bloch', dest)
    resultfile = dest / 'raw/generalized_bloch_reference_independent_core_v23.json'
    resultfile.unlink()
    with (out / 'bloch.log').open('w') as log:
        subprocess.run([sys.executable, '-B', str(dest / 'code/reproduce_generalized_bloch_reference.py')],
                       check=True, stdout=log, env={**os.environ, 'PYTHONOPTIMIZE':'0'})
    methods = json.loads(resultfile.read_text())['methods']
    if methods.keys() != expected.keys():
        raise ValueError('Generalized-Bloch extrapolation keys differ')
    for method, row in expected.items():
        for key, value in row.items():
            require_close(methods[method][key], value, 2.e-12)
    result = {'scope':'Reference taper at omega 1.2, finite-dimensional replay from historical integrations.',
              'surface':surface_ladder,'tetrahedron':tetra_ladder,'generalized_bloch':methods}
    (out / 'reference_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status':'PASS','surface_N192':surface_ladder[192]['qT'],
                      'tetra_N192':tetra_ladder[192]['TT_modal'],
                      'GB_low_pair':methods['low_pair_linear_richardson']['TT_modal']},indent=2))


if __name__ == '__main__':
    main()
