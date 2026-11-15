"""Run preserved numerical fixtures on fresh copies, with outputs kept outside the inputs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('profile', choices=('rigidity','scattering','checkpoint','clock','continuum','radial'))
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('These fixtures require assertions. Disable Python optimization.')
    out = args.output.resolve()
    if HERE == out or HERE in out.parents:
        parser.error('Use a new external directory for output.')
    out.mkdir(parents=True, exist_ok=False)
    env = {**os.environ, 'PYTHONDONTWRITEBYTECODE':'1', 'PYTHONOPTIMIZE':'0',
           'OPENBLAS_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1'}
    def run(path, *argv):
        with (out / (path.stem+'.log')).open('w') as log:
            subprocess.run([sys.executable, '-B', str(path), *map(str,argv)], check=True,
                           cwd=out, stdout=log, env=env)
    def require_pass(path, key='status'):
        val=json.loads(path.read_text())[key]
        if val is not True and str(val).lower() not in ('passed','pass','pass: finite-fixture coordinate regression only; no production transmission claim'):
            raise ValueError(f'{path.name} returned {key}={val!r}')
    if args.profile == 'rigidity':
        target = out / 'rigidity'
        shutil.copytree(HERE/'rigidity', target)
        (target/'records').mkdir()
        run(target/'run/p02_independent_checks_path_relative.py')
        require_pass(target/'records/P02_INDEPENDENT_CHECKS_PASS16.json','all_pass')
    elif args.profile == 'scattering':
        target = out / 'scattering'
        shutil.copytree(HERE/'scattering', target)
        (target/'results/P05_WEAK_COUPLING.json').unlink()
        for task, result in [('core_scattering','CORE_SCATTERING'),('matching_plane','MATCHING_PLANE'),
                             ('current_degeneracy','CURRENT_AND_DEGENERACY'),('weak_coupling','WEAK_COUPLING')]:
            run(target/'run'/f'p05_{task}.py')
            require_pass(target/'results'/f'P05_{result}.json')
    elif args.profile == 'checkpoint':
        run(HERE/'dynamics/run/checkpoint_smoke_test.py')
        require_pass(out/'checkpoint_smoke_test.log')
    elif args.profile == 'clock':
        target = out/'validate_r7_nonaffine_fixture.py'
        shutil.copyfile(HERE/'perspectives'/target.name,target)
        run(target)
        require_pass(out/'R7_NONAFFINE_FIXTURE_RESULTS.json')
    elif args.profile == 'continuum':
        run(HERE/'continuum/reproduce_continuum_benchmark.py','--output',out/'continuum.json')
        require_pass(out/'continuum.json')
    elif args.profile == 'radial':
        run(HERE/'shell/check_radial_shell.py', '--bulk-source',
            HERE/'analysis/math_resolution/tetra/combined_tetra_diagnostic.py','--output',out/'radial.json')
        require_pass(out/'radial.json')
    print(json.dumps({'profile':args.profile,'status':'PASS'}))


if __name__ == '__main__':
    main()
