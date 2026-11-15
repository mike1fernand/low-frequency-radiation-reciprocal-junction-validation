#!/usr/bin/env python3
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,hashlib
from pathlib import Path
import numpy as np
import low_frequency as m
ROOT=Path(__file__).resolve().parent

def load(tag,w):
 z=np.load(ROOT/tag/f'kernel_{w:.8f}.npz');return z['PV'],[z[f'Q{j}'] for j in range(3)]

def contour(w,nang=24,nrad=96,height=.15):
 # No eigendecomposition, shell roots, residue subtraction or finite eta here.
 # Deform the near-origin radial segment into lower complex t half-plane.
 u,wu=m.gl(nang);t,wt=m.gl(nrad);uu,vv=np.meshgrid(u,u,indexing='ij');aw=(wu[:,None]*wu[None,:]).ravel()
 base=np.pi*np.stack([np.ones(uu.size),uu.ravel(),vv.ravel()],axis=-1)
 out=np.zeros((6,3),complex);tc=4*w/np.pi
 z=tc*(t-1j*height*np.sin(np.pi*t));dz=tc*(1-1j*height*np.pi*np.cos(np.pi*t))
 far=tc+(1-tc)*t
 for p in [base,base[:,[1,0,2]],base[:,[1,2,0]]]:
  for a in range(0,len(p),32):
   pp=p[a:a+32]
   for arg,jac in [(z,dz),(far,np.full_like(t,1-tc))]:
    k=pp[:,None,:]*arg[None,:,None];iv=np.linalg.inv(m.dyn_bonds(k)-w*w*np.eye(3))
    v=arg[None,:,None,None]**2*m.projected(iv,k)
    out+=np.einsum('r,t,t,rtdc->dc',aw[a:a+32],wt,jac,v)
 return m.assemble(out.real)+1j*m.assemble(out.imag)

bT=1/(6*np.pi);bL=1/(12*np.pi*3**1.5);C=4*m.CQ*np.array([bT,bL]);ks=8/3
pred={'fixed_s_coefficients':{'qT':float(C[0]),'qL':float(C[1]),'bulk':float(C.sum())},'share':float(C[0]/C.sum()),'weak_coefficients':[float(bT*ks*ks/m.CQ),float(bL*ks*ks/m.CQ)],'kSigma':ks,'crossover_omega_over_s':ks/(2*m.CQ)}
checks={};ref=[]
for w in [0,.2,.05,.0125,.003125,.00078125]:
 P,bs=load('validation',w);P0,b0=load('refined',w)
 r={'omega':w,'PV_relative_grid_change':float(np.linalg.norm(P-P0)/np.linalg.norm(P))}
 if w:
  r.update(qT_relative_grid_change=float(np.linalg.norm(bs[0]+bs[1]-b0[0]-b0[1])/np.linalg.norm(bs[0]+bs[1])),qL_relative_grid_change=float(np.linalg.norm(bs[2]-b0[2])/np.linalg.norm(bs[2])))
  r['transmission_relative_changes']={str(s):abs(m.scatter(w,s,P,bs)['TqT']/m.scatter(w,s,P0,b0)['TqT']-1) for s in [1,.1,.01,.001,w]}
 ref.append(r)
checks['grid_refinement']=ref
# Compare direct retarded contour integral, with different D builder, to real-axis PV/shell method.
co=[]
for w in [.2,.05,.003125]:
 P,bs=load('validation',w);G=P+1j*sum(bs)
 for h in [.12,.24]:
  start=time.perf_counter();gc=contour(w,height=h)
  co.append({'omega':w,'height':h,'n_face':24,'n_radial':96,'PV_relative_difference':float(np.linalg.norm(gc.real-P)/np.linalg.norm(P)),'spectral_relative_difference':float(np.linalg.norm(gc.imag-G.imag)/np.linalg.norm(G.imag)),'seconds':time.perf_counter()-start})
checks['independent_contour']=co
# Different algebraic implementation of the coupled line/core solve.
aug=[]
for w in [.2,.003125]:
 P,bs=load('validation',w);G=P+1j*sum(bs)
 for s in [1,.01]:
  U,V,Ks,L=m.matrices(s*m.GREF);kap=2*np.arcsin(w/(2*m.CQ));gam=2*m.KQ*np.sin(kap)
  line=w*w*np.eye(7)-L-Ks;line=line.astype(complex);line[0,0]+=m.KQ*np.exp(1j*kap);line[-1,-1]+=m.KQ*np.exp(1j*kap)
  full=np.block([[line,-V.T],[V,np.linalg.inv(G)+U]]);rhs=np.zeros(19,complex);rhs[0]=np.sqrt(gam)
  sol=np.linalg.solve(full,rhs);eta=-np.linalg.solve(G,sol[7:]);tt=float(2*np.real(eta.conj()@(bs[0]+bs[1])@eta));row=m.scatter(w,s,P,bs)
  aug.append({'omega':w,'s':s,'augmented_qT_relative_difference':abs(tt/row['TqT']-1),'augmented_residual':float(np.linalg.norm(full@sol-rhs)/np.linalg.norm(rhs))})
checks['augmented_solve']=aug
# Held-out profile/segment controls. These are not input to deriving the universal coefficient.
profiles={'uniform6':[1/6]*6,'edge6':[.2,.15,.15,.15,.15,.2],'asymmetric6':[.05,.1,.15,.2,.25,.25],'single':[.3],'two':[.1,.1]}
held=[]
for name,g in profiles.items():
 n=len(g)
 for w in [.05,.003125,.00078125]:
  P,bs=load('validation',w);P=P[:2*n,:2*n];bs=[b[:2*n,:2*n] for b in bs]
  row=m.scatter(w,1,P,bs,np.array(g));row.update(profile_name=name,profile=g,qT_scaled_relative_error=row['TqT']/w**2/C[0]-1,qL_scaled_relative_error=row['TqL']/w**2/C[1]-1);held.append(row)
checks['held_out_profiles']=held
# Independent D bond formula and derivative finite difference.
rng=np.random.default_rng(6173);k=rng.uniform(-np.pi,np.pi,(137,3));p=rng.normal(size=(137,3));D,dp=m.dyn(k,p)
fd=(m.dyn_bonds(k+1e-5*p)-m.dyn_bonds(k-1e-5*p))/(2e-5)
checks['model_regressions']={'dynamical_matrix_relative':float(np.linalg.norm(D-m.dyn_bonds(k))/np.linalg.norm(D)),'radial_derivative_finite_difference_relative':float(np.linalg.norm(dp-fd)/np.linalg.norm(dp))}
# Leading rank-one spectral measure and actual asymptotic/crossover rows.
fixed=[];joint=[];spectral=[];v=np.tile([2.,0.],6)
rows=json.loads((ROOT/'refined/rows.json').read_text())
for row in rows:
 w=row['omega'];s=row['s']
 if s in [1.,.1,.01,.001]:
  fixed.append({**row,'qT_scaled_relative_error':row['TqT']/w**2/C[0]-1,'qL_scaled_relative_error':row['TqL']/w**2/C[1]-1})
 for lam in [.1,1,10]:
  if abs(s-lam*w)<1e-14:
   F=(lam*ks)**2/((lam*ks)**2+(2*m.CQ)**2)
   joint.append({**row,'lambda':lam,'predicted_R_limit':F,'R_absolute_error':row['R']-F,'qT_scaled_relative_error':row['TqT']/w**2/(C[0]*F)-1})
for w in [.2,.05,.0125,.003125,.00078125]:
 P,bs=load('validation',w)
 spectral.append({'omega':w,'qT_leading_matrix_relative_error':float(np.linalg.norm((bs[0]+bs[1])/w-bT*np.outer(v,v))/np.linalg.norm(bT*np.outer(v,v))),'qL_leading_matrix_relative_error':float(np.linalg.norm(bs[2]/w-bL*np.outer(v,v))/np.linalg.norm(bL*np.outer(v,v)))})
checks['spectral_asymptotic']=spectral
for name,obj in [('predictions.json',pred),('verification.json',checks),('fixed_coupling_rows.json',fixed),('joint_limit_rows.json',joint)]:m.dump(ROOT/name,obj)
print(json.dumps({'predictions':pred,'grid_refinement':ref,'contour':co,'model':checks['model_regressions'],'max_closure_all_rows':max(abs(r['closure']) for r in rows)},indent=2))
