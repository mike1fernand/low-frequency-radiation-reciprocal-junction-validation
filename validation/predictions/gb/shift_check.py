import time,json
from pathlib import Path
import numpy as np
import p1solver as p
ROOT=Path(__file__).resolve().parent
raw=ROOT/'reference_inputs'
G0=np.load(raw/'gb_eta010_g256_mean.npz')['G']
out=[];Gs=[]
for shift in [(.25,.25),(.375,.375),(.25,.5),(.125,.375),(.21132486540518713,.21132486540518713),(.5,.5)]:
 t=time.time(); x=p.green_blocks(np.sqrt(1.2**2-.01**2),6,256,shift=shift,eps=.024,onshell=False)
 G=p.assemble(x['G'],6); Gs.append(G)
 v=dict(shift=shift,seconds=time.time()-t,max_abs_diff=float(np.max(abs(G-G0))),rel_fro=float(np.linalg.norm(G-G0)/np.linalg.norm(G0)),first=str(G[0,0]))
 out.append(v);print(json.dumps(v),flush=True)
np.savez_compressed(ROOT/'shift_candidates.npz',G=np.stack(Gs),shifts=np.array([v['shift'] for v in out]))
(ROOT/'shift_check.json').write_text(json.dumps(out,indent=2)+'\n')
