"""Run one explicitly selected box/timestep control; never choose physics values."""
from pathlib import Path
import argparse, hashlib, json, os, resource, subprocess, time
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--half-box',type=int,required=True)
ap.add_argument('--dt',type=float,required=True)
ap.add_argument('--end',type=float,required=True)
ap.add_argument('--name',required=True)
ap.add_argument('--frozen-window',nargs=2,type=float)
args=ap.parse_args()
out=ROOT/args.name;assert not out.exists();out.mkdir()
config=json.loads((ROOT/'PREDECLARED_PACKET_CONFIGURATION.json').read_text())
config.update(bulk_shape=[120,2*args.half_box,2*args.half_box,3],quarter_shape=[120,args.half_box,args.half_box,3],dt=args.dt,target_time=args.end,target_final_step=round(args.end/args.dt),numerical_control_name=args.name)
if args.frozen_window:config['frozen_window']=args.frozen_window
assert abs(config['target_final_step']*args.dt-args.end)<1e-12
(out/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')
g=config['profile'];cmd=[str(ROOT/'evolve_candidate'),str(out),str(config['target_final_step']),format(args.dt,'.17g'),str(args.half_box)]+[format(x,'.17g') for x in g]+[str(round(12/args.dt))]
env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
t=time.perf_counter();run=subprocess.run(cmd,input=(ROOT/'initial_qv.txt').read_text(),capture_output=True,text=True,env=env)
(out/'evolution.stdout').write_text(run.stdout+'\n'+run.stderr)
assert run.returncode==0,run.stderr+run.stdout
r=json.loads(run.stdout);r.update(wall_seconds=time.perf_counter()-t,peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,source_sha256=hashlib.sha256((ROOT/'evolve_candidate.cpp').read_bytes()).hexdigest())
(out/'EVOLUTION_CHECK.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
