"""Fresh directional-shell checks with the available historical 2.6 real-Green fit.
Reuses saved complex Green inputs; does not regenerate their quadrature.
"""
from pathlib import Path
import sys,json,time,hashlib,ast
import numpy as np
ROOT=Path(__file__).resolve().parents[2]; SRC=Path(__file__).resolve().parent/'replay_sources';sys.path.insert(0,str(SRC))
from p1solver import assemble,assemble_im,solve_junction,reference_profile
from p1split import onshell_split
OUT=Path(__file__).resolve().parent/'dressed_26';OUT.mkdir(exist_ok=True)
N=int(sys.argv[1]); Nc=24;w=2.6;eps=np.array([.006,.008,.012,.016,.024])
saved={ast.literal_eval(k):v for k,v in np.load(SRC/'r3_eps_scan.npy',allow_pickle=True).item().items()}
Gs=np.array([saved[(w,float(e))] for e in eps]);c0=np.linalg.lstsq(np.vander(1j*eps,3,increasing=True),Gs.reshape(len(eps),-1),rcond=None)[0][0].reshape(Nc,2,2)
st=time.perf_counter();I,mv=onshell_split(w,Nc,N)
np.savez_compressed(OUT/f'shell_N{N}.npz',Im=I,omega=w,Nc=Nc,Nperp=N,min_abs_dlambda_dkx=mv)
M={(b,d):assemble_im(I[b,d],Nc,Nc) for b in range(3) for d in range(2)}
g0=assemble(c0.real+1j*I.sum((0,1))[Nc-1:].real,Nc)
r=solve_junction(reference_profile(Nc),w,g0,sum(M[b,d] for b in (0,1) for d in (0,1)),sum(M[2,d] for d in (0,1)),Kq=2.25)
eta=r['eta'];minus=float(2*np.real(eta.conj()@sum(M[b,1] for b in range(3))@eta))
res={'omega':w,'Nc':Nc,'cq':1.5,'s':1,'profile':reference_profile(Nc).tolist(),'Kaxis':1,'Kdiag':1,'mb':1,'mq':1,'a':1,'N_shell':N,'N_real_Green_historical':2048,'real_green_fit_degree':2,'real_green_epsilon':eps.tolist(),'z_convention':'omega**2+1j*epsilon','historical_Green_source':'audit_two_reports/evidence/paper1_round2_backward_emission/r3_eps_scan.npy','historical_Green_sha256':hashlib.sha256((SRC/'r3_eps_scan.npy').read_bytes()).hexdigest(),'T_minus':minus,'B_minus':minus/float(r['Tbulk']),'Tbulk':float(r['Tbulk']),'R':float(r['R']),'Tq':float(r['Tq']),'closure':float(r['closure']),'elapsed_seconds':time.perf_counter()-st,'scope':'Fresh directional shell; retained historical real-Green inputs. Not a fresh PV quadrature check or rigorous uncertainty enclosure.'}
(OUT/f'results_N{N}.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
