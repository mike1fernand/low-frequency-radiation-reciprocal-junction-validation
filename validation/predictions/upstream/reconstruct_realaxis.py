"""Fresh reconstruction of specified SM7E real-axis tests; not historical record recovery."""
from pathlib import Path
import sys, json, time, hashlib, platform
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'replay_sources'))
from p1solver import green_blocks, assemble, assemble_im, solve_junction, reference_profile, coupling_mats
from p1split import onshell_split
w=float(sys.argv[1]); N=int(sys.argv[2]); Nc=24; prof=reference_profile(Nc); Kq=2.25
OUT=ROOT/'realaxis';OUT.mkdir(exist_ok=True)
tag=f'w{w:.2f}_N{N}';start=time.perf_counter()
r=green_blocks(w,Nc,N,eps=0.0,onshell=False,quadrant=True)
print(json.dumps({'stage':'realaxis_green','omega':w,'Nperp':N,'stats':r['stats'],'seconds':time.perf_counter()-start}),flush=True)
Im,minv=onshell_split(w,Nc,N)
np.savez_compressed(OUT/f'{tag}_inputs.npz',G=r['G'],Im=Im,omega=w,Nc=Nc,Nperp=N,profile=prof,Kq=Kq)
M={(b,d):assemble_im(Im[b,d],Nc,Nc) for b in range(3) for d in range(2)}
Qp=sum(M[b,0] for b in range(3));Qm=sum(M[b,1] for b in range(3));Qt=Qp+Qm
g0=assemble(r['G'].real+1j*Im.sum((0,1))[Nc-1:].real,Nc)
sol=solve_junction(prof,w,g0,sum(M[b,d] for b in (0,1) for d in (0,1)),sum(M[2,d] for d in (0,1)),Kq=Kq)
eta=sol['eta'];minus=float(2*np.real(eta.conj()@Qm@eta));plus=float(2*np.real(eta.conj()@Qp@eta))
# Independent assembly of the simultaneous line/core equations, using the same kernels.
U,V,Ks=coupling_mats(prof);nd=Nc+1;kap=2*np.arcsin(w/(2*np.sqrt(Kq)));sig=-Kq*np.exp(1j*kap);gam=2*Kq*np.sin(kap)
KD=2*Kq*np.eye(nd)-Kq*(np.eye(nd,k=1)+np.eye(nd,k=-1))
A=w*w*np.eye(nd)-KD-Ks;A=A.astype(complex);A[0,0]-=sig;A[-1,-1]-=sig
H=np.block([[A,-V.T],[g0@V,np.eye(2*Nc)+g0@U]])
rhs=np.zeros(nd+2*Nc,complex);rhs[0]=np.sqrt(gam)
x=np.linalg.solve(H,rhs);chi2=x[:nd];u2=x[nd:];eta2=V@chi2+U@u2
minus2=float(2*np.real(eta2.conj()@Qm@eta2));plus2=float(2*np.real(eta2.conj()@Qp@eta2))
row={'omega':w,'Nperp':N,'T_minus':minus,'T_plus':plus,'B_minus':minus/(minus+plus),'Tbulk':float(sol['Tbulk']),'R':float(sol['R']),'Tq':float(sol['Tq']),'closure':float(sol['closure']),'simultaneous_T_minus':minus2,'simultaneous_T_plus':plus2,'simultaneous_chi_max_difference':float(np.max(abs(chi2-sol['chi']))),'simultaneous_eta_max_difference':float(np.max(abs(eta2-eta))),'simultaneous_equation_residual':float(np.max(abs(H@x-rhs))),'g0_imH_vs_shell_max_difference':float(np.max(abs((g0-g0.conj().T)/2j-Qt))),'raw_G_imag_vs_shell_max_difference':float(np.max(abs(r['G'].imag-Im.sum((0,1))[Nc-1:].real))),'min_abs_dlambda_dkx':float(minv),'green_stats':r['stats'],'elapsed_seconds':time.perf_counter()-start,'profile':prof.tolist(),'Kq':Kq,'a_mb_mq_Kaxis_Kdiag':[1,1,1,1,1],'numpy':np.__version__,'python':platform.python_version(),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'replay_sources').glob('*.py')},'scope':'Fresh real-axis analytic axial Green integration and directional shell at stated transverse grid. Imaginary part assembled from directional shell; closure and simultaneous-system agreement test algebra, not integration accuracy.'}
(OUT/f'{tag}_result.json').write_text(json.dumps(row,indent=2)+'\n')
print(json.dumps(row),flush=True)
