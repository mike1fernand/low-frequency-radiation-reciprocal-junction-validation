"""Replay the supplied bandwidth solver and audit saved-kernel consistency.
No Green integration is rerun. Uses existing comparison rows without extending
frequencies or coupling profiles; records every evaluated variant explicitly.
"""
from pathlib import Path
import sys, json, hashlib, platform, subprocess
import numpy as np
import scipy
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
CHECK=BASE/'portable_checks'
sys.path.insert(0,str(CHECK/'rev'))
from alpha_solver import solve_alpha,load_g0

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 fs_path=CHECK/'comparison_inputs/FINAL_SCREEN.json'
 fs=json.loads(fs_path.read_text())
 out={'verification_scope':'Saved-kernel replay; no new Green integration; no continuous band certificate',
      'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
      'inputs':{'FINAL_SCREEN_sha256':digest(fs_path)},'counts':{},'kernels':[],'comparisons':[],'probe_runs':[]}
 for key in ('surface_rows','independent_rows'):
  rows=fs[key]
  keys=[(r['omega'],r['setting']) for r in rows]
  assert len(keys)==len(set(keys)),('duplicate keys',key)
  out['counts'][key]={'total':len(rows),'frequency_counts':{str(w):sum(r['omega']==w for r in rows) for w in sorted(set(r['omega'] for r in rows))}}
 surf={(r['omega'],r['setting']):r for r in fs['surface_rows']}
 tet={(r['omega'],r['setting']):r for r in fs['independent_rows']}
 for p in sorted((CHECK/'rev').glob('g0_*.npz')):
  with np.load(p) as d:
   w=float(d['omega']); assert w==float(p.stem[3:])
   n=int(d['nmax']); assert n>=6
   km={'filename':p.name,'sha256':digest(p),'omega':w,'nmax':n,'checks':[]}
   frequencies=sorted({r['omega'] for r in fs['surface_rows']})
   km['comparison_rows_available']=sum(abs(r['omega']-w)<1e-9 for r in fs['surface_rows'])
   km['shell_160_240_max_absolute_difference']=float(np.max(np.abs(d['ImT_160']+d['ImL_160']-d['ImT_240']-d['ImL_240'])))
   for N in (2048,4096):
    if f'scan_{N}' not in d: continue
    for deg in (2,3):
     V=np.vander(1j*d['eps'],deg+1,increasing=True)
     scan=d[f'scan_{N}']
     rebuilt=np.linalg.lstsq(V,scan.reshape(len(d['eps']),-1),rcond=None)[0][0].reshape(scan.shape[1:])
     saved=d[f'c0_{N}_deg{deg}']
     km['checks'].append({'Nperp':N,'degree':deg,'refit_c0_max_absolute_difference':float(np.max(np.abs(rebuilt-saved))),
       'extrapolated_im_vs_onshell_max_absolute_difference':float(np.max(np.abs(saved.imag-(d['ImT_240']+d['ImL_240'])[n-1:].real)))})
     rows=[r for r in fs['surface_rows'] if abs(r['omega']-w)<1e-9]
     if not rows: continue
     assert len(rows)==193
     assert sum(r['setting']=='nominal' for r in rows)==1
     assert sum(r['setting'].startswith('corner_') for r in rows)==64
     for nth in (160,240):
      g0,ImT,ImL=load_g0(w,Np=N,deg=deg,nth=nth)
      rec=[]
      for r in rows:
       match=tet.get((r['omega'],r['setting']))
       assert match is not None and np.array_equal(r['g'],match['g'])
       vals=solve_alpha(np.array(r['g']),w,g0,ImT,ImL,1.69,alpha=1.0)
       rec.append({'setting':r['setting'],'TqT':float(vals['TqT']),'surface_TqT':r['TqT'],'tetra_TqT':match['TT_modal'],'closure':float(vals['closure'])})
      mine=np.array([r['TqT'] for r in rec]); sur=np.array([r['surface_TqT'] for r in rec]); te=np.array([r['tetra_TqT'] for r in rec])
      mini=rec[int(np.argmin(mine))]
      out['comparisons'].append({'omega':w,'Nperp':N,'degree':deg,'shell_nth':nth,'count':len(rec),
         'minimum':mini,'maximum_absolute_difference_surface':float(np.max(np.abs(mine-sur))),
         'maximum_absolute_difference_tetra':float(np.max(np.abs(mine-te))),
         'maximum_closure_residual':max(abs(r['closure']) for r in rec),'rows':rec})
   out['kernels'].append(km)
 for arg in (['1.1125','1.2375','1.3'],['1.29']):
  proc=subprocess.run([sys.executable,str(CHECK/'rev/bandwidth_check.py'),*arg],capture_output=True,text=True)
  out['probe_runs'].append({'arguments':arg,'exitcode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr})
 (HERE/'bandwidth_saved_kernel_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 short={k:v for k,v in out.items() if k not in ('comparisons','probe_runs')}
 short['comparisons']=[{k:v for k,v in r.items() if k!='rows'} for r in out['comparisons']]
 short['probe_exitcodes']=[{'arguments':r['arguments'],'exitcode':r['exitcode']} for r in out['probe_runs']]
 (HERE/'bandwidth_summary.json').write_text(json.dumps(short,indent=2)+'\n')
 print(json.dumps(short,indent=2))
if __name__=='__main__': main()
