import sys,time,json,hashlib
from pathlib import Path
import numpy as np
import p1solver as p
ROOT=Path(__file__).resolve().parent
raw=ROOT/'reference_inputs'
out=[]
for label,eps,wr in [('squared_frequency',.01,1.2),('frequency',.024,np.sqrt(1.2**2-.01**2)),('linear_frequency_approx',.024,1.2)]:
 t=time.time()
 x=p.green_blocks(wr,6,256,eps=eps,onshell=False,quadrant=True)
 G=p.assemble(x['G'],6)
 G0=np.load(raw/'gb_eta010_g256_mean.npz')['G']
 np.savez_compressed(ROOT/f'candidate_{label}.npz',G=G,eps=eps,wr=wr)
 v=dict(label=label,eps=eps,wr=wr,seconds=time.time()-t,max_abs_diff=float(np.max(abs(G-G0))),rel_fro=float(np.linalg.norm(G-G0)/np.linalg.norm(G0)),first=str(G[0,0]),stats={k:(v if not isinstance(v,float) or np.isfinite(v) else None) for k,v in x['stats'].items()})
 out.append(v);print(json.dumps(v),flush=True)
(ROOT/'candidate_check.json').write_text(json.dumps(out,indent=2)+'\n')
