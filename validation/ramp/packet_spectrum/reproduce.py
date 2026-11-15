"""Replay the frozen fourteen-node spectrum calculation in a fresh directory."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output-dir',required=True,type=Path)
args=ap.parse_args();out=args.output_dir.resolve()
if out.exists() or out==ROOT or ROOT in out.parents:
    raise SystemExit('Choose a new directory outside the frozen packet-spectrum package.')
out.mkdir(parents=True)
for n in ['INPUT_BINDING.json','run_packet_spectrum.py']:shutil.copyfile(ROOT/n,out/n)
shutil.copytree(ROOT/'inputs',out/'inputs')
env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
run=subprocess.run([sys.executable,str(out/'run_packet_spectrum.py')],capture_output=True,text=True,env=env)
(out/'replay.stdout').write_text(run.stdout+'\n'+run.stderr)
if run.returncode:raise SystemExit(run.stderr or run.stdout)
old=json.loads((ROOT/'PACKET_SPECTRUM.json').read_text());new=json.loads((out/'PACKET_SPECTRUM.json').read_text())
assert new['status']=='PASS' and new['reference_regressions']==new['new_ramp_states']==70
differences=[]
for a,b in zip(old['families'],new['families']):
    assert a['label']==b['label']
    differences.append(abs(a['new_ramp_full_energy_in_band_contribution']['qT']-b['new_ramp_full_energy_in_band_contribution']['qT']))
assert max(differences)<1e-12
print(json.dumps({'status':'PASS','reference_regressions':70,'new_ramp_states':70,
                  'maximum_packet_prediction_difference':max(differences),
                  'scope':'Replay of the same saved cases, not new kernels or new predictions.',
                  'output_dir':str(out)},indent=2))
