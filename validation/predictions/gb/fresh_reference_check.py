import json,time,sys,hashlib
from pathlib import Path
import numpy as np
import p1solver as p
ROOT=Path(__file__).resolve().parent
from independent_core_scattering import evaluate
raw=ROOT/'reference_inputs'
def clean_stats(stats):
 return {k:(v if not isinstance(v,float) or np.isfinite(v) else None) for k,v in stats.items()}
results=[]
for N in [384,768,1536]:
 gs=[];rec={'Nperp':N,'omega':1.2,'eta':[.01,.005,.0025],'shift':[.5,.5],'regulator':'z=(omega+i*eta)^2','kernels':[]}
 for eta in rec['eta']:
  t=time.time();o=p.green_blocks(np.sqrt(1.2**2-eta**2),6,N,eps=2*1.2*eta,onshell=False,quadrant=True)
  G=p.assemble(o['G'],6);gs.append(G)
  name=f'fresh_G_N{N}_eta{eta}.npz';np.savez_compressed(ROOT/name,G=G,blocks=o['G'],Nperp=N,omega=1.2,eta=eta,shift=[.5,.5])
  r={'file':name,'seconds':time.time()-t,'stats':clean_stats(o['stats'])};rec['kernels'].append(r);print(json.dumps({'N':N,'eta':eta,**r}),flush=True)
 t=time.time();o=p.green_blocks(1.2,6,N,eps=0,onshell=True,quadrant=True)
 IL=p.assemble_im(o['ImL'],6,6);IT=p.assemble_im(o['ImT'],6,6)
 name=f'fresh_residue_N{N}.npz';np.savez_compressed(ROOT/name,G=p.assemble(o['G'],6),ImL=IL,ImT=IT,Nperp=N,omega=1.2,shift=[.5,.5])
 rec['residue']={'file':name,'seconds':time.time()-t,'stats':clean_stats(o['stats'])}
 G10,G5,G25=gs
 fits={'low_pair_linear_richardson':2*G25-G5,'three_point_quadratic':G10/3-2*G5+8*G25/3}
 rec['results']={}
 for label,G in fits.items():
  I=(G-G.conj().T)/(2j)
  r=evaluate(1.2,[.1,.2,.2,.2,.2,.1],G,{'qL':IL,'qT':I-IL})
  rec['results'][label]={k:r[k] for k in ['Rq','Tq','Tbulk','closure','partials']}
 print(json.dumps(rec),flush=True);results.append(rec)
 (ROOT/'fresh_reference_results.json').write_text(json.dumps(results,indent=2)+'\n')
