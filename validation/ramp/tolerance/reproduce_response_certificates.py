"""Replay the frozen new response certificates into a fresh directory.

Requires Python3 and NumPy. This replays the exact Fraction response bounds,
not the separate full-domain contour certificate. Only the staged kernel input
location is changed; the frozen scripts/receipts remain untouched.
"""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys

HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output-dir',type=Path,required=True)
args=ap.parse_args();out=args.output_dir.resolve()
if out.exists() or out==HERE or HERE in out.parents:
    raise SystemExit('Choose a NEW output directory outside the frozen certificate directory.')
kernel=HERE/'inputs/REFERENCE_KERNEL.npz'
assert hashlib.sha256(kernel.read_bytes()).hexdigest()=='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
out.mkdir(parents=True);(out/'inputs').mkdir()
shutil.copyfile(kernel,out/'inputs/REFERENCE_KERNEL.npz')
ptext=(HERE/'verify_point_rational.py').read_text()
old="SOURCE=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs/REFERENCE_KERNEL.npz')"
assert ptext.count(old)==1
(out/'verify_point_rational.py').write_text(ptext.replace(old,"SOURCE=HERE/'inputs/REFERENCE_KERNEL.npz'",1))
shutil.copyfile(HERE/'verify_box_rational.py',out/'verify_box_rational.py')
for filename in ['verify_point_rational.py','verify_box_rational.py']:
    run=subprocess.run([sys.executable,str(out/filename)],capture_output=True,text=True)
    (out/(filename+'.stdout')).write_text(run.stdout+'\n'+run.stderr)
    if run.returncode:raise SystemExit('Replay failed: '+filename+'\n'+run.stderr)
point=json.loads((out/'POINT_RATIONAL_CERTIFICATE.json').read_text())
box=json.loads((out/'BOX_RATIONAL_CERTIFICATE.json').read_text())
assert point['status']==box['status']=='PASS'
result={'status':'PASS','point_interval':point['tighter_conditional_interval'],
        'whole_plus_minus5_percent_box_lower':box['tight_kernel_box_lower'],
        'scope':'Exact response propagation replay only; upstream contour-enclosure replay is a separate linked dependency.',
        'output_directory':str(out)}
print(json.dumps(result,indent=2))
