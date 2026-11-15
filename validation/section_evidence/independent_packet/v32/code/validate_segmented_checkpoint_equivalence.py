#!/usr/bin/env python3
"""Verify checkpointed and uninterrupted target-blind evolution are identical."""
from __future__ import annotations
import argparse,json,subprocess,sys,tempfile
from pathlib import Path


def strip_nondeterministic(d):
    return {k:v for k,v in d.items() if k not in {'elapsed_seconds','segmented_execution'}}


def maxdiff(a,b):
    if isinstance(a,dict): return max([maxdiff(a[k],b[k]) for k in a] or [0.0])
    if isinstance(a,list): return max([maxdiff(x,y) for x,y in zip(a,b)] or [0.0])
    if isinstance(a,(int,float)) and isinstance(b,(int,float)): return abs(float(a)-float(b))
    return 0.0 if a==b else 1.0


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase12-root',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    worker=Path(__file__).with_name('run_target_blind_global_td.py')
    base=[sys.executable,str(worker),'--phase12-root',a.phase12_root,'--nx','40','--ny','12','--nz','12','--dt','0.04','--tmax','10','--sample-start','8','--sample-step','1','--trace-step','1','--slab-start','20','--x0','5','--sigma','3']
    with tempfile.TemporaryDirectory(prefix='p1s_checkpoint_test_') as td:
        td=Path(td);full=td/'full.json';seg=td/'seg.json';c1=td/'c1.npz';c2=td/'c2.npz'
        subprocess.run(base+['--output',str(full)],check=True,capture_output=True,text=True)
        subprocess.run(base+['--output',str(seg),'--segment-end','4','--checkpoint-out',str(c1)],check=True,capture_output=True,text=True)
        subprocess.run(base+['--output',str(seg),'--resume-from',str(c1),'--segment-end','8','--checkpoint-out',str(c2)],check=True,capture_output=True,text=True)
        subprocess.run(base+['--output',str(seg),'--resume-from',str(c2)],check=True,capture_output=True,text=True)
        f=strip_nondeterministic(json.loads(full.read_text()));s=strip_nondeterministic(json.loads(seg.read_text()))
        out={'geometry':[40,12,12],'dt':0.04,'tmax':10.0,'full_equals_segmented':f==s,'maximum_absolute_numeric_difference':maxdiff(f,s),'snapshots':len(f['snapshots']),'trace_points':len(f['trace'])}
    Path(a.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,indent=2))
    if not out['full_equals_segmented']: raise SystemExit(1)
if __name__=='__main__':main()
