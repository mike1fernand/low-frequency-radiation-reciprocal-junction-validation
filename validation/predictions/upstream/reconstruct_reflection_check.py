"""Fresh numerical check of the stated spatial-reflection condition; same kernels."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'replay_sources'))
from p1solver import assemble,assemble_im,coupling_mats,reference_profile
r=np.load(ROOT/'realaxis/w2.90_N1536_inputs.npz');w=float(r['omega']);Nc=int(r['Nc']);Kq=float(r['Kq'])
Im=r['Im'];g0=assemble(r['G'].real+1j*Im.sum((0,1))[Nc-1:].real,Nc)
Qp=sum(assemble_im(Im[b,0],Nc,Nc) for b in range(3));Qm=sum(assemble_im(Im[b,1],Nc,Nc) for b in range(3))
def solve(prof,side):
 U,V,Ks=coupling_mats(prof);nd=Nc+1;k=2*np.arcsin(w/(2*np.sqrt(Kq)));gam=2*Kq*np.sin(k);sig=-Kq*np.exp(1j*k)
 KD=2*Kq*np.eye(nd)-Kq*(np.eye(nd,k=1)+np.eye(nd,k=-1));A=(w*w*np.eye(nd)-KD-Ks).astype(complex);A[0,0]-=sig;A[-1,-1]-=sig
 H=np.block([[A,-V.T],[g0@V,np.eye(2*Nc)+g0@U]])
 rhs=np.zeros(nd+2*Nc,complex);rhs[0 if side=='left' else nd-1]=np.sqrt(gam)
 x=np.linalg.solve(H,rhs);eta=V@x[:nd]+U@x[nd:]
 plus=float(2*np.real(eta.conj()@Qp@eta));minus=float(2*np.real(eta.conj()@Qm@eta))
 return {'T_plus':plus,'T_minus':minus,'T_bulk':plus+minus,'upstream_share':(minus if side=='left' else plus)/(plus+minus)}
p=reference_profile(Nc);a=p*np.linspace(.8,1.2,Nc)
sets=[]
for label,prof in [('symmetric',p),('asymmetric_profile_also_reflected',a)]:
 left=solve(prof,'left');right=solve(prof[::-1],'right')
 sets.append({'case':label,'left_profile':prof.tolist(),'right_profile':prof[::-1].tolist(),'left_incidence':left,'right_incidence':right,'laboratory_hemisphere_swap_max_difference':max(abs(left['T_plus']-right['T_minus']),abs(left['T_minus']-right['T_plus'])),'relative_upstream_share_difference':abs(left['upstream_share']-right['upstream_share'])})
out={'input':'realaxis/w2.90_N1536_inputs.npz','rows':sets,'scope':'Numerical corroboration of the reflection condition using shared kernels; does not independently certify integration or establish directionality for an asymmetric profile left unreflected.'}
(ROOT/'reflection_check_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
