"""Fresh sampled axial-root completeness check against a 70-digit cosine cubic.
Forty difficult transverse points are selected by largest A1 condition number on
the Nperp=1024 midpoint quadrant; forty random points use a recorded fixed seed.
This is a new sample, not recovery of the unspecified historical sample.
"""
from pathlib import Path
import sys,json,time
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'replay_sources'))
from p1solver import bulk_blocks,D_full,dD_dkx
mp.mp.dps=70; w=mp.mpf('2.9');t0=time.perf_counter();N=1024
k=(np.arange(N//2)+.5)*2*np.pi/N;ky,kz=np.meshgrid(k,k,indexing='ij');ky=ky.ravel();kz=kz.ravel()
A0,A1=bulk_blocks(ky,kz);cond=np.linalg.cond(A1);ids=np.argsort(cond)[-40:][::-1]
points=[('difficult',float(ky[i]),float(kz[i]),float(cond[i])) for i in ids]
seed=20260928;rng=np.random.default_rng(seed)
points.extend(('random',float(y),float(z),None) for y,z in rng.uniform(0,np.pi,(40,2)))
def add(a,b):
 out=[mp.mpf(0)]*max(len(a),len(b))
 for i,c in enumerate(a):out[i]+=c
 for i,c in enumerate(b):out[i]+=c
 return out
def scale(a,s):return [s*x for x in a]
def mul(a,b):
 out=[mp.mpf(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
records=[];det_resid=0
for typ,y,z,condition in points:
 Y=mp.mpf(y);Z=mp.mpf(z);cy,cz,sy,sz=mp.cos(Y),mp.cos(Z),mp.sin(Y),mp.sin(Z)
 a=[6-w*w,-2*(1+cy+cz)];b=[6-2*cy-2*cy*cz-w*w,-2*cy];c=[6-2*cz-2*cy*cz-w*w,-2*cz]
 p=mul(mul(a,b),c)
 p=add(p,scale([1,0,-1],16*sy*sy*sz*sz))
 p=add(p,scale(a,-4*sy*sy*sz*sz))
 p=add(p,scale(mul(b,[1,0,-1]),-4*sz*sz))
 p=add(p,scale(mul(c,[1,0,-1]),-4*sy*sy))
 croots=mp.polyroots(list(reversed(p)),maxsteps=1000,error=False,extraprec=100)
 exact=[]
 for r in croots:
  if abs(mp.im(r))<mp.mpf('1e-50') and -1<mp.re(r)<1:
   angle=float(mp.acos(mp.re(r)));exact.extend([-angle,angle])
 exact=np.array(sorted(exact))
 B0,B1=bulk_blocks(np.array([y]),np.array([z]));B0=B0[0];B1=B1[0]
 C=np.block([[np.zeros((3,3)),np.eye(3)],[-np.linalg.solve(B1,np.conj(B1)),-np.linalg.solve(B1,B0-float(w*w)*np.eye(3))]])
 roots=np.linalg.eigvals(C);unit=abs(abs(roots)-1)<=1e-6;got=np.sort(np.angle(roots[unit]));same=len(exact)==len(got)
 diff=float(np.max(abs(exact-got))) if same and len(got) else (0.0 if same else None)
 for xx in np.linspace(-2.9,2.9,9):
  val=sum(float(p[i])*np.cos(xx)**i for i in range(4));direct=np.linalg.det(D_full(np.array([xx,y,z]))-float(w*w)*np.eye(3));det_resid=max(det_resid,float(abs(val-direct)))
 err=float(np.max(np.min(abs(np.linalg.eigvalsh(D_full(np.column_stack([got,np.full(len(got),y),np.full(len(got),z)])))-float(w*w)),axis=1))) if len(got) else 0.
 records.append({'sample_type':typ,'ky':y,'kz':z,'A1_condition_number':condition or float(np.linalg.cond(B1)),'cosine_cubic_coefficients_low_to_high':[str(x) for x in p],'cosine_cubic_roots':[str(x) for x in croots],'high_precision_real_kx_roots':exact.tolist(),'companion_real_kx_roots':got.tolist(),'same_root_count':same,'maximum_kx_difference':diff,'maximum_shell_eigenvalue_residual':err})
result={'omega':2.9,'Nperp_difficult_selection':N,'difficult_selection':'40 largest A1 2-norm condition numbers on Nperp1024 midpoint quadrant','random_seed':seed,'random_count':40,'difficult_count':40,'precision_decimal_digits':70,'all_root_counts_match':all(r['same_root_count'] for r in records),'maximum_kx_difference':max(r['maximum_kx_difference'] or 0 for r in records),'maximum_shell_eigenvalue_residual':max(r['maximum_shell_eigenvalue_residual'] for r in records),'polynomial_vs_literal_nine_bond_determinant_max_residual':det_resid,'records':records,'elapsed_seconds':time.perf_counter()-t0,'scope':'New sampled completeness comparison. Does not recover historical sample identities, certify unsampled points, or bound integrated power error.'}
(ROOT/'root_checks_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
