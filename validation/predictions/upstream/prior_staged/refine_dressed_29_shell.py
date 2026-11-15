"""Refine only the 2.9 shell to1536; reuse freshly generated N1024 PV inputs."""
from pathlib import Path
import sys,time,json,numpy as np
ROOT=Path(__file__).resolve().parents[2]; SRC=Path(__file__).resolve().parent/'replay_sources';sys.path.insert(0,str(SRC))
from p1solver import assemble,assemble_im,solve_junction,reference_profile
from p1split import onshell_split
OUT=Path(__file__).resolve().parent/'dressed_29';Nc=24;w=2.9;N=1536;st=time.perf_counter()
I,mv=onshell_split(w,Nc,N);np.savez_compressed(OUT/f'shell_N{N}.npz',Im=I,omega=w,Nc=Nc,Nperp=N,min_abs_dlambda_dkx=mv)
data=np.load(OUT/'green_N1024.npz');Gs=data['G'];eps=data['epsilon'];M={(b,d):assemble_im(I[b,d],Nc,Nc) for b in range(3) for d in range(2)}
rows=[]
for deg in (2,3):
    c0=np.linalg.lstsq(np.vander(1j*eps,deg+1,increasing=True),Gs.reshape(len(eps),-1),rcond=None)[0][0].reshape(Nc,2,2)
    g0=assemble(c0.real+1j*I.sum((0,1))[Nc-1:].real,Nc)
    r=solve_junction(reference_profile(Nc),w,g0,sum(M[b,d] for b in (0,1) for d in (0,1)),sum(M[2,d] for d in (0,1)),Kq=2.25)
    eta=r['eta'];minus=float(2*np.real(eta.conj()@sum(M[b,1] for b in range(3))@eta))
    rows.append({'N_shell':N,'N_green':1024,'fit_degree':deg,'T_minus':minus,'B_minus':minus/float(r['Tbulk']),'Tbulk':float(r['Tbulk']),'R':float(r['R']),'Tq':float(r['Tq']),'closure':float(r['closure'])})
res={'omega':w,'Nc':Nc,'cq':1.5,'s':1,'profile':reference_profile(Nc).tolist(),'z_convention':'omega**2+1j*epsilon','epsilons':eps.tolist(),'minimum_abs_dlambda_dkx':float(mv),'rows':rows,'seconds':time.perf_counter()-st,'scope':'Fresh1536 shell and previously generated fresh1024 complex Green inputs. Replays the specified SM7E alternative PV estimates, not the separate historical real-axis central estimate.'}
(OUT/'results_shell1536_green1024.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
