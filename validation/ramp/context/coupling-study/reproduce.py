"""Replay all 40 saved cases with the independent spring-Hessian checker.

Requires Python 3 and NumPy. No external Paper One directory is needed.
Outputs go to a new directory; frozen results are never overwritten.
"""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys

ROOT = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output-dir', required=True, type=Path)
args = ap.parse_args()
out = args.output_dir.resolve()
if out.exists():
    raise SystemExit('Choose a new output directory to preserve existing results.')
if out == ROOT or ROOT in out.parents:
    raise SystemExit('Keep replay outputs outside the frozen assessment package.')

inputs = ROOT/'numerical_probe/inputs'
spec = json.loads((ROOT/'numerical_probe/PREDECLARATION.json').read_text())
for name, expected in spec['input_sha256'].items():
    p = inputs/name
    if hashlib.sha256(p.read_bytes()).hexdigest() != expected:
        raise SystemExit('Input hash mismatch: '+name)
for rel in ['numerical_probe/SCALE_PROBE.json',
            'numerical_probe/profile_comparison/PROFILE_COMPARISON.json',
            'independent_review/check_probe.py']:
    target = out/rel
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ROOT/rel, target)
shutil.copytree(inputs, out/'numerical_probe/inputs')
run = subprocess.run([sys.executable, str(out/'independent_review/check_probe.py')],
                     capture_output=True, text=True)
(out/'replay.stdout').write_text(run.stdout+'\n'+run.stderr)
if run.returncode:
    raise SystemExit('Replay failed; inspect '+str(out/'replay.stdout'))
review = json.loads((out/'independent_review/REVIEW.json').read_text())
assert review['status']=='PASS' and review['independent_case_count']==40
assert review['maximum_partition_disagreement'] < 1e-12
print(json.dumps({'status':'PASS','independent_cases':40,
                  'maximum_partition_disagreement':review['maximum_partition_disagreement'],
                  'scope':'Independent physical-spring assembly; existing approximate bare kernels reused.',
                  'output_dir':str(out)}, indent=2))
