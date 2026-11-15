"""Fresh bounded verification of the stated 24-motif omega=2.9 result.
Uses inspected archived producers, without changing them. Saves generated inputs.
Run with a transverse size, e.g. `python recheck_dressed_29.py 256`.
"""
from pathlib import Path
import sys, json, hashlib, time, platform
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path(__file__).resolve().parent/'replay_sources'
sys.path.insert(0,str(SOURCE))
from p1solver import green_blocks,assemble,assemble_im,solve_junction,reference_profile
from p1split import onshell_split
OUT=Path(__file__).resolve().parent/'dressed_29'; OUT.mkdir(exist_ok=True)
N=int(sys.argv[1]); w=2.9; Nc=24; epss=np.array([.02,.03,.04,.06]); prof=reference_profile(Nc)
start=time.perf_counter()
Im,mv=onshell_split(w,Nc,N)
np.savez_compressed(OUT/f'shell_N{N}.npz',Im=Im,omega=w,Nc=Nc,Nperp=N,min_abs_dlambda_dkx=mv)
print(json.dumps({'stage':'shell','N':N,'seconds':time.perf_counter()-start,'minimum_abs_dlambda_dkx':float(mv)}),flush=True)
Gs=[]; stats=[]
for eps in epss:
    st=time.perf_counter(); r=green_blocks(w,Nc,N,eps=float(eps),onshell=False,quadrant=True)
    Gs.append(r['G']); stats.append(r['stats'])
    print(json.dumps({'stage':'complex_green','N':N,'epsilon':float(eps),'seconds':time.perf_counter()-st,'max_root_residual':float(r['stats']['max_resid'])}),flush=True)
Gs=np.array(Gs)
np.savez_compressed(OUT/f'green_N{N}.npz',G=Gs,epsilon=epss,omega=w,Nc=Nc,Nperp=N,profile=prof,Kq=2.25)
M={(b,d):assemble_im(Im[b,d],Nc,Nc) for b in range(3) for d in range(2)}
Qp=sum(M[b,0] for b in range(3)); Qm=sum(M[b,1] for b in range(3)); Qt=Qp+Qm
rows=[]
for deg in (2,3):
    V=np.vander(1j*epss,deg+1,increasing=True)
    c=np.linalg.lstsq(V,Gs.reshape(4,-1),rcond=None)[0][0].reshape(Nc,2,2)
    g0=assemble(c.real+1j*Im.sum((0,1))[Nc-1:].real,Nc)
    r=solve_junction(prof,w,g0,sum(M[b,d] for b in (0,1) for d in (0,1)),sum(M[2,d] for d in (0,1)),Kq=2.25)
    eta=r['eta']; minus=float(2*np.real(eta.conj()@Qm@eta)); plus=float(2*np.real(eta.conj()@Qp@eta))
    row={'N_shell':N,'N_green':N,'fit_degree':deg,'T_minus':minus,'T_plus':plus,'B_minus':minus/(minus+plus),'R':float(r['R']),'Tq':float(r['Tq']),'Tbulk':float(r['Tbulk']),'closure':float(r['closure']),'qminus_hermitian_residual':float(np.max(np.abs(Qm-Qm.conj().T))),'qminus_min_eigenvalue':float(np.linalg.eigvalsh((Qm+Qm.conj().T)/2).min()),'g0_imH_vs_Q_residual':float(np.max(np.abs((g0-g0.conj().T)/2j-Qt))),'fit_imag_vs_shell_max_entry_mismatch':float(np.max(np.abs(c.imag-Im.sum((0,1))[Nc-1:].real)))}
    rows.append(row); print(json.dumps({'stage':'dressed_result',**row}),flush=True)
meta={'omega':w,'Nc':Nc,'Kq':2.25,'cq':1.5,'s':1,'a_mb_mq_Kaxis_Kdiag':[1,1,1,1,1],'profile':prof.tolist(),'z_convention':'omega**2+1j*epsilon','epsilons':epss.tolist(),'Nperp':N,'quadrant_midpoint':True,'source_files_sha256':{n:hashlib.sha256((SOURCE/n).read_bytes()).hexdigest() for n in ('p1solver.py','p1split.py')},'numpy':np.__version__,'python':platform.python_version(),'rows':rows,'elapsed_seconds':time.perf_counter()-start,'interpretation':'Fresh generated Green blocks and fully dressed junction solve. Internal closure uses the same shell as the source Green imaginary part; it is not independent quadrature certification.'}
(OUT/f'results_N{N}.json').write_text(json.dumps(meta,indent=2)+'\n')
