"""Targeted audit; imports supplied portable solver, not original angular integrator."""
from pathlib import Path
import sys,json,hashlib,csv
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
sys.path[:0]=[str(BASE/'portable_checks'),str(BASE/'portable_checks/rev')]
from p1solver import D_full,gradD_dir,phi_vectors,reference_profile
from alpha_solver import solve_alpha,load_g0_ref12,coupling_mats_alpha
rng=np.random.default_rng(292709)
k=rng.uniform(-1.1,1.1,(37,3)); h=1e-5
analytic=np.stack([gradD_dir(k,np.broadcast_to(e,k.shape)) for e in np.eye(3)],axis=1)
fd=np.stack([(D_full(k+h*e)-D_full(k-h*e))/(2*h) for e in np.eye(3)],axis=1)
G,QT,QL=load_g0_ref12()
rows={}; symmetry=0.; source_difference=0.; coupled_difference=0.; cartesian_difference=0.
Nc=6;Nq=7;nb=72
for name,a,x in [('full',1.,1.),('X-only',0.,1.),('R-only',1.,0.)]:
 r=solve_alpha(reference_profile(Nc),1.2,G,QT,QL,1.69,alpha=a,xfac=x)
 U,V,Ks=coupling_mats_alpha(reference_profile(Nc),a,x)
 # Literal Cartesian quadratic energies and orthonormal source patterns.
 B=np.zeros((nb,12));H=np.zeros((Nq+nb,Nq+nb));sites=[]
 for s,g in enumerate(reference_profile(Nc)):
  for ci,(sy,sz) in enumerate([(-1,-1),(-1,1),(1,-1),(1,1)]):
   p=12*s+3*ci; sites.append([s,sy/2,sz/2]);B[p,2*s]=.5
   B[p+1,2*s+1]=sy/(2*np.sqrt(2));B[p+2,2*s+1]=sz/(2*np.sqrt(2))
   A=np.zeros(Nq+nb);A[Nq+p]=1;A[s]=A[s+1]=-.5
   R=np.zeros(Nq+nb);R[Nq+p+1]=sy;R[Nq+p+2]=sz;R[s]=-.5;R[s+1]=.5
   H+=(2*g/3)*(x*np.outer(A,A)+a*np.outer(R,R))
 cartesian_difference=max(cartesian_difference,float(np.max(abs(U-B.T@H[Nq:,Nq:]@B))),float(np.max(abs(V-B.T@H[Nq:,:Nq]))),float(np.max(abs(Ks-H[:Nq,:Nq]))))
 gamma=2*1.69*np.sin(r['kappa']); KD=1.69*(2*np.eye(Nq)-np.eye(Nq,k=1)-np.eye(Nq,k=-1))
 lead=np.zeros((Nq,Nq),complex);lead[0,0]=lead[-1,-1]=-1.69*np.exp(1j*r['kappa'])
 A=1.2**2*np.eye(Nq)-KD-Ks-lead
 joint=np.block([[A,-V.T],[G@V,np.eye(12)+G@U]])
 rhs=np.zeros(Nq+12,complex);rhs[0]=np.sqrt(gamma)
 sol=np.linalg.solve(joint,rhs);eta=V@sol[:Nq]+U@sol[Nq:]
 coupled_difference=max(coupled_difference,float(np.max(abs(eta-r['eta']))))
 eta=r['eta']; literal=np.einsum('ps,sc->pc',np.exp(-1j*k@np.asarray(sites).T),(B@eta).reshape(24,3))
 reduced=np.einsum('ps,pcb,sb->pc',np.exp(-1j*k[:,0,None]*np.arange(Nc)),phi_vectors(k[:,1],k[:,2]),eta.reshape(Nc,2))
 source_difference=max(source_difference,float(np.max(abs(literal-reduced))))
 def amplitudes(kk):
  f=np.einsum('ps,pcb,sb->pc',np.exp(-1j*kk[:,0,None]*np.arange(Nc)),phi_vectors(kk[:,1],kk[:,2]),eta.reshape(Nc,2))
  lam,vec=np.linalg.eigh(D_full(kk));return abs(np.einsum('pcb,pc->pb',vec,f))**2
 amp=amplitudes(k)
 for sy,sz,swap in [(sy,sz,swap) for sy in [-1,1] for sz in [-1,1] for swap in [False,True]]:
  kp=k.copy();kp[:,1]*=sy;kp[:,2]*=sz
  if swap:kp=kp[:,[0,2,1]]
  symmetry=max(symmetry,float(np.max(abs(amplitudes(kp)-amp))))
 rows[name]={kk:float(r[kk]) for kk in ['R','Tq','TqT','TqL','closure']}
res={'finite_difference_h':h,'gradient_max_abs_difference':float(np.max(abs(analytic-fd))),
'cartesian_coupling_max_abs_difference':cartesian_difference,
'simultaneous_vs_schur_eta_max_abs_difference':coupled_difference,
'literal_source_vs_projected_max_abs_difference':source_difference,
'C4v_branch_intensity_max_abs_difference':symmetry,
'new_kernel_states':rows,
'limitations':['Deterministic 37-point implementation checks, not global interval proofs.','All source cross-solves here share the supplied independent numerical Green blocks.','No fresh PV integration or global angular-preimage certification in this script.']}
(HERE/'angular_implementation_results.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
