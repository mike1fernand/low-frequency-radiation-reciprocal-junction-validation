from pathlib import Path
import json,hashlib,shutil
import numpy as np
from scipy.optimize import root
from core import symbol,roots,states,source,E,R

P=Path(__file__).resolve().parent
source_dir=P.parents[1]/'review_r8/input/paper1-validation-integrated-r8/legacy_rc4/validation/current_numerics/analysis/resolution_execution/selected'
(P/'inputs').mkdir(exist_ok=True)
for name in ['REFERENCE_KERNEL.npz','REFERENCE_KERNEL.json']:
 destination=P/'inputs'/name
 if not destination.exists():
  if not (source_dir/name).exists():
   raise FileNotFoundError(f'Missing bundled input {destination}; original archive fallback is also unavailable')
  shutil.copyfile(source_dir/name,destination)
meta=json.loads((P/'inputs/REFERENCE_KERNEL.json').read_text())
assert hashlib.sha256((P/'inputs/REFERENCE_KERNEL.npz').read_bytes()).hexdigest()==meta['kernel_sha256']
eta,reference=states()
checks=[]
for direction in [1,-1]:
 for branch,kx in [(0,direction*2*np.arcsin(.6)),(2,direction*2*np.arcsin(1.2/np.sqrt(12)))]:
    k=np.array([kx,0.,0.]);val,vec=np.linalg.eigh(symbol(k));PT=vec[:,:2]@vec[:,:2].T
    checks.append({'k':k.tolist(),'shell_branch':branch,'shell_residual':float(abs(val[branch]-E)),
       'states':{name:{'qT_projector_weight':float(np.vdot(source(k,t),PT@source(k,t)).real),'longitudinal_weight':float(abs(source(k,t)[0])**2)} for name,t in eta.items()}})
# Dense full angular sample including every radial sheet by the archived analytic census.
mu,w=np.polynomial.legendre.leggauss(48);phi=(np.arange(64)+.5)*2*np.pi/64
m,p=np.meshgrid(mu,phi,indexing='ij');n=np.stack([m,np.sqrt(1-m*m)*np.cos(p),np.sqrt(1-m*m)*np.sin(p)],axis=-1).reshape(-1,3)
global_stats=[];seeds=[]
for branch in range(3):
 rr=roots(n,branch);gv=rr['grad'];angle=np.arctan2(np.linalg.norm(gv[:,1:],axis=1),abs(gv[:,0]))
 # Retain likely axial ray preimages and nonaxial low-angle local candidates.
 ix=np.argsort(angle)[:80]
 for j in ix:seeds.append((branch,rr['k'][j]))
 global_stats.append({'branch':branch,'min_angle_degrees':float(np.rad2deg(angle.min())),
   'min_radial_derivative':float(rr['radial'].min()),'max_shell_residual':float(np.max(abs(rr['values'][:,branch]-E))),
   'min_vg_component_times_k':np.min(gv*rr['k'],axis=0).tolist()})
sols=[];attempts=[]
for branch,seed in seeds:
 def f(k):
    d,j=symbol(k,True);v,e=np.linalg.eigh(d);u=e[:,branch];gr=np.einsum('i,aij,j->a',u,j,u)
    return np.array([v[branch]-E,gr[1],gr[2]])
 fit=root(f,seed,tol=1e-10);res=float(np.linalg.norm(f(fit.x)))
 if res<1e-7 and np.max(abs(fit.x))<R:
    k=fit.x;val=np.linalg.eigvalsh(symbol(k));gap=float(np.min(abs(np.delete(val,branch)-val[branch])))
    if not any(branch==s['branch'] and np.linalg.norm(k-np.array(s['k']))<1e-5 for s in sols):
     sols.append({'branch':branch,'k':k.tolist(),'residual':res,'gap':gap})
 attempts.append({'branch':branch,'success':bool(fit.success),'residual':res})
out={'reference':reference,'axial_cluster_checks':checks,'global_sample':global_stats,'axial_ray_candidates':sols,
 'search_attempts':len(attempts),'search_scope':'Numerical multi-start, not completeness proof; lower-plane degeneracies require separate treatment.'}
(P/'pilot.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
