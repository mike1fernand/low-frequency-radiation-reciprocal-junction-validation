"""Staged Cartesian midpoint reconstruction of SM7E broadened upstream diagnostics.

Uses literal nine-bond sums and four-corner transforms; no Bloch roots.
Each chunk is separately saved and combined before fitting powers in epsilon.
Usage: python reconstruct_cartesian.py N start stop [full]
       python reconstruct_cartesian.py N combine [full]
"""
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
PRODUCER_HASH=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
sys.path.insert(0,str(ROOT/'replay_sources'))
from p1solver import assemble,assemble_im,solve_junction,reference_profile
N=int(sys.argv[1]); mode=sys.argv[2]; full='full' in sys.argv
OUT=ROOT/('cartesian_full' if full else 'cartesian_quadrant')/f'N{N}';OUT.mkdir(parents=True,exist_ok=True)
w=2.9;Nc=24;eps=np.array([.2,.15,.1,.075,.05]);ns=np.arange(-Nc+1,Nc)
kxgrid=-np.pi+(np.arange(N)+.5)*2*np.pi/N
perp=kxgrid if full else (np.arange(N//2)+.5)*2*np.pi/N
if mode=='combine':
    files=sorted(OUT.glob('chunk_*.npz')); G=np.zeros((5,47,2,2),complex); Qm=G.copy(); Qp=G.copy(); indices=[]
    for p in files:
        r=np.load(p)
        assert int(r['N'])==N and float(r['omega'])==w and int(r['Nc'])==Nc
        assert np.array_equal(r['epsilon'],eps) and bool(r['full_grid'])==full
        assert str(r['producer_sha256'])==PRODUCER_HASH, 'Stale chunk: producer differs'
        G+=r['G'];Qm+=r['Qm'];Qp+=r['Qp'];indices.extend(r['indices'].tolist())
    assert sorted(indices)==list(range(N)), 'Chunks must cover each axial grid index exactly once'
    G/=N*len(perp)**2;Qm/=N*len(perp)**2;Qp/=N*len(perp)**2
    rows=[]
    for ei in range(5):
        g0=assemble(G[ei,Nc-1:],Nc);qm=assemble_im(Qm[ei],Nc,Nc);qp=assemble_im(Qp[ei],Nc,Nc)
        r=solve_junction(reference_profile(Nc),w,g0,qm+qp,np.zeros_like(qm),Kq=2.25)
        tm=float(2*np.real(r['eta'].conj()@qm@r['eta']));tp=float(2*np.real(r['eta'].conj()@qp@r['eta']))
        rows.append({'epsilon':float(eps[ei]),'T_minus_broadened_diagnostic':tm,'T_plus_broadened_diagnostic':tp,'closure':float(r['closure']),'ImG_vs_Q_max_residual':float(np.max(abs((g0-g0.conj().T)/2j-qm-qp)))})
    fits={str(d):float(np.polynomial.polynomial.polyfit(eps,[r['T_minus_broadened_diagnostic'] for r in rows],d)[0]) for d in (2,3,4)}
    result={'N':N,'omega':w,'Nc':Nc,'Kq':2.25,'epsilons':eps.tolist(),'grid':'full midpoint' if full else 'midpoint, using exact transverse reflection symmetries','rows':rows,'zero_width_fits_in_epsilon':fits,'scope':'New Cartesian integration reconstruction; finite-width powers and zero-width polynomial fits are numerical diagnostics, not physical broadening or certified bounds. Existing line solver and block assembly are reused.'}
    np.savez_compressed(OUT/'combined_inputs.npz',G=G,Qm=Qm,Qp=Qp,epsilon=eps,omega=w,N=N,Nc=Nc)
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    raise SystemExit
start=int(mode);stop=int(sys.argv[3]); t0=time.perf_counter()
ky,kz=np.meshgrid(perp,perp,indexing='ij');ky=ky.ravel();kz=kz.ravel();P=len(ky)
# Literal orthonormal centered X and R patterns on the four corners.
corners=np.array([[0,sy*.5,sz*.5] for sy in (-1,1) for sz in (-1,1)])
patterns=np.zeros((4,3,2));patterns[:,0,0]=.5
patterns[:,1,1]=corners[:,1]/np.sqrt(2);patterns[:,2,1]=corners[:,2]/np.sqrt(2)
phase=np.exp(-1j*(ky[:,None]*corners[:,1]+kz[:,None]*corners[:,2]))
Phi=np.einsum('pc,cva->pva',phase,patterns)
bonds=np.array([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)],float)
unit=bonds/np.linalg.norm(bonds,axis=1)[:,None]
G=np.zeros((5,47,2,2),complex);Qm=G.copy();Qp=G.copy()
for ix in range(start,stop):
    kx=kxgrid[ix]; D=np.zeros((P,3,3));dD=D.copy()
    for d,e in zip(bonds,unit):
        angle=kx*d[0]+ky*d[1]+kz*d[2];ee=np.outer(e,e)
        D+=2*(1-np.cos(angle))[:,None,None]*ee
        dD+=2*d[0]*np.sin(angle)[:,None,None]*ee
    lam,evec=np.linalg.eigh(D)
    # Resolve the velocity operator within tied eigenspaces. Arbitrary eigenvectors
    # at a degeneracy must not determine the direction mask of a broadened diagnostic.
    tie01=np.abs(lam[:,1]-lam[:,0])<1e-10;tie12=np.abs(lam[:,2]-lam[:,1])<1e-10
    for mask,inds in ((tie01 & ~tie12,[0,1]),(tie12 & ~tie01,[1,2]),(tie01 & tie12,[0,1,2])):
        sel=np.flatnonzero(mask)
        if len(sel):
            basis=evec[sel][:,:,inds]
            vel=np.einsum('pva,pvw,pwb->pab',basis,dD[sel],basis)
            _,rot=np.linalg.eigh(vel);newbasis=basis@rot
            for j,b in enumerate(inds):evec[sel,:,b]=newbasis[:,:,j]
    vx=np.einsum('pvb,pvw,pwb->pb',evec,dD,evec)
    pe=np.einsum('pva,pvb->pab',np.conj(Phi),evec)
    ee=np.einsum('pab,pcb->pbac',pe,np.conj(pe))
    delta=lam-w*w
    for ei,ep in enumerate(eps):
        green=1/(delta-1j*ep);im=green.imag
        bg=np.einsum('pb,pbac->ac',green,ee)
        bm=np.einsum('pb,pbac->ac',im*(vx<0),ee)
        bp=np.einsum('pb,pbac->ac',im*(vx>0),ee)
        phases=np.exp(1j*kx*ns)[:,None,None]
        G[ei]+=phases*bg;Qm[ei]+=phases*bm;Qp[ei]+=phases*bp
np.savez_compressed(OUT/f'chunk_{start:04d}_{stop:04d}.npz',G=G,Qm=Qm,Qp=Qp,indices=np.arange(start,stop),omega=w,N=N,Nc=Nc,epsilon=eps,full_grid=full,producer_sha256=PRODUCER_HASH)
print(json.dumps({'N':N,'start':start,'stop':stop,'full':full,'seconds':time.perf_counter()-t0}),flush=True)
