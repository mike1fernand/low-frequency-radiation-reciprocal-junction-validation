"""Replay the frozen P3/P1 symbol checks with bundled input in a fresh folder.

Requires Python 3 and NumPy. The public replay adapts the model location and document-identity field.
The mathematical calculations and numerical cases retain their original implementation.
"""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output-dir',type=Path,required=True)
args=ap.parse_args();out=args.output_dir.resolve()
if out.exists() or out==ROOT or ROOT in out.parents:
    raise SystemExit('Choose a new output directory outside the frozen package.')
source=ROOT/'frequency_bound/check_frequency_bound.py'
text=source.read_text()
old="MODEL=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs/model_api_v1_1.py')"
new="MODEL=HERE.parent/'inputs/model_api_v1_1.py'"
assert text.count(old)==1
(out/'frequency_bound').mkdir(parents=True)
(out/'inputs').mkdir()
# The historical Paper 3 draft supplied a provenance hash only.
shutil.copyfile(ROOT/'inputs/model_api_v1_1.py',out/'inputs/model_api_v1_1.py')
staged=out/'frequency_bound/check_frequency_bound.py'
reference_field = "'P3_text_sha256':H(P3),"
if text.count(reference_field) != 1:
    raise RuntimeError("Historical provenance field changed.")
staged.write_text(text.replace(old,new,1).replace(reference_field,
    "'paper3_document_verification':'outside public symbol replay',",1))
run=subprocess.run([sys.executable,str(staged)],capture_output=True,text=True)
(out/'replay.stdout').write_text(run.stdout+'\n'+run.stderr)
if run.returncode:raise SystemExit(run.stderr or run.stdout)
receipt=json.loads((out/'frequency_bound/FREQUENCY_BOUND_CHECK.json').read_text())
assert receipt['status']=='PASS'
print(json.dumps({'status':'PASS','scope':'Same frozen symbol/energy calculations; bundled model location and document-provenance field adapted.',
                  'frozen_script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                  'staged_script_sha256':hashlib.sha256(staged.read_bytes()).hexdigest(),
                  'output_dir':str(out)},indent=2))
