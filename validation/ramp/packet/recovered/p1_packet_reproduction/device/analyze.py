"""Independent direct finite-device dynamics, all cadence samples, modal degeneracies and prep control."""
from pathlib import Path
import sys,json
import numpy as np
path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
K=1.69;A=np.array([.1,.2,.2,.2,.2,.1]);dt=.04;N=120;M=56;h=M//2
line=np.memmap(path/'line.bin',dtype=float,mode='r').reshape(-1,2,N)
raw=np.memmap(path/'snap.bin',dtype=float,mode='r').reshape(46,2,N+3*N*h*h)
origin,initial_v=line[0]
E0=(np.sum(initial_v**2)+K*np.sum(np.diff(origin)**2))/2
link0=(2/3)*np.sum(A*(origin[60:66]**2+origin[61:67]**2))
def unfold(row):
    pos=row[N:].reshape(N,h,h,3)
    b=np.concatenate((pos[:,:,::-1]*[1,1,-1],pos),axis=2)
    return np.concatenate((b[:,::-1]*[1,-1,1],b),axis=1)
def masks(r):
    lo=tuple(slice(max(-d,0),s-max(d,0)) for d,s in zip(r,(N,M,M)))
    return lo
directions=[np.array(r) for r in ((1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1))]
ks=np.meshgrid(2*np.pi*np.fft.fftfreq(N),2*np.pi*np.fft.fftfreq(M),2*np.pi*np.fft.rfftfreq(M),indexing='ij',sparse=True)
H=np.zeros((N,M,M//2+1,3,3))
for r in directions:
    phase=sum(d*k for d,k in zip(r,ks) if d)
    H+=((2-2*np.cos(phase))/(r@r))[...,None,None]*np.outer(r,r)
lam,vec=np.linalg.eigh(H);lam=np.maximum(0,lam);n3=N*M*M
spatial_weights=np.broadcast_to(np.r_[1,np.ones(h-1)*2,1],lam.shape[:-1])
near_top=(lam[...,2]-lam[...,1])<=1e-12*np.maximum(1,lam[...,2])
triple=near_top & ((lam[...,2]-lam[...,0])<=1e-12*np.maximum(1,lam[...,2]))
def branch(x,v):
    coeff=[np.einsum('...ia,...i->...a',vec,np.fft.rfftn(m,axes=(0,1,2))) for m in [x,v]]
    u,p=coeff; w=spatial_weights
    e=(abs(p)**2+lam*abs(u)**2)*w[...,None]/(2*n3*E0)
    fixed=e[~near_top,:2].sum()+e[near_top&~triple,0].sum()
    lo=hi=fixed
    for start,rank,msk in [(1,1,near_top&~triple),(0,2,triple)]:
        px,ux=p[msk,start:],u[msk,start:];eig=lam[msk,start:].mean(axis=-1)
        gram=(px[:,:,None]*px[:,None,:].conj()+eig[:,None,None]*ux[:,:,None]*ux[:,None,:].conj())*w[msk,None,None]/(2*n3*E0)
        eigen=np.linalg.eigvalsh(gram);lo+=np.sum(eigen[:,:rank]);hi+=np.sum(eigen[:,-rank:])
    return float(lo),float(hi),float(e.sum())
times=np.arange(240,286)/5
rows=[];branch_rows=[]
for i,t in enumerate(times):
    fx,fv=raw[i];x=unfold(fx);v=unfold(fv)
    pbulk=0.;periodic=0.
    for r in directions:
        span=sum(x[...,j]*d for j,d in enumerate(r) if d)
        delta=span-np.roll(span,-r,(0,1,2));dbl=delta**2/(2*(r@r));pbulk+=dbl[masks(r)].sum();periodic+=dbl.sum()
    kin=(v*v).sum()/2; bulk=kin+pbulk
    shell=(kin-(v[:,2:-2,2:-2]**2).sum()/2)/kin
    q=fx[:N];p=fv[:N];e_line=(p@p+K*(np.diff(q)@np.diff(q)))/2
    c=fx[N:].reshape(N,h,h,3)[60:66,0,0]
    cv=fv[N:].reshape(N,h,h,3)[60:66,0,0]
    J=link=0
    for e,sx in enumerate((1,-1)):
        ext=sx*(c[:,0]-q[60+e:66+e])+c[:,1]+c[:,2]
        link+=np.sum(4*A*ext**2/6)
        J-=np.sum(4*A*ext*(sx*(p[60+e:66+e]+cv[:,0])+cv[:,1]+cv[:,2])/6)
    drift=abs((e_line+bulk+link)-(E0+link0))/(E0+link0)
    diff=(periodic-pbulk)/bulk
    passed=bool(link/E0<=.0025 and abs(J)/E0<=.0005 and shell<=.002 and diff<=.0025 and drift<=1e-5)
    rows.append(dict(t=float(t),open_bulk=float(bulk/E0),link=float(link/E0),J=float(J/E0),shell=float(shell),periodic_difference=float(diff),drift=float(drift),passed=passed))
    if 52.8<=t<=55.4:
        lo,hi,total=branch(x,v)
        assert abs(total-(periodic+kin)/E0)<1e-12
        branch_rows.append([float(t),lo,hi,total])
passing=[r['t'] for r in rows if r['passed']]
runlist=[];last=None
for t in passing:
    if last is None or t-last>.21:runlist.append([])
    runlist[-1].append(t);last=t
longest=max([g for g in runlist if g[-1]-g[0]>=2],key=len)
vel=np.fft.fft(line[:,1],axis=-1)
freq=2.6*np.abs(np.sin(np.pi*np.fft.fftfreq(N)));sense=np.sign(np.fft.fftfreq(N));sense[N//2]=0
u0=np.fft.fft(origin);v0=np.fft.fft(initial_v);all_times=np.arange(len(line))*dt
vf=v0*np.cos(all_times[:,None]*freq)-u0*freq*np.sin(all_times[:,None]*freq)
rms={}
for T in (5,10,20):
    p=slice(0,round(T/dt)+1); moving=sense!=0
    def ratio(er,ref):
        return float(np.sqrt(np.trapezoid(np.sum(np.abs(er)**2,axis=1),dx=dt)/np.trapezoid(np.sum(np.abs(ref)**2,axis=1),dx=dt)))
    sig=vf[p];err=vel[p]-sig
    rms[str(T)]=dict(running=ratio(err[:,moving],sig[:,moving]),all=ratio(err,sig))
corners=np.memmap(path/'corner.bin',dtype=float,mode='r').reshape(-1,6,2,3)
assert len(corners)==len(line)
J=np.zeros(len(line))
for e,sx in enumerate((1,-1)):
    c=corners[:,:,0,:]; cv=corners[:,:,1,:]
    ext=sx*(c[:,:,0]-line[:,0,60+e:66+e])+c[:,:,1]+c[:,:,2]
    J-=4*np.sum(A*ext*(sx*(line[:,1,60+e:66+e]+cv[:,:,0])+cv[:,:,1]+cv[:,:,2]),axis=1)/6
cum=np.r_[0,np.cumsum(dt*(J[1:]+J[:-1])/2)]/E0
valid=[r for r in rows if 52.8<=r['t']<=55.4]
assigned=np.array([r['open_bulk']+.5*(r['link']-link0/E0) for r in valid])
indices=np.array([round(r['t']/dt) for r in valid])
assert np.max(abs(J[indices]/E0-[r['J'] for r in valid]))<1e-14
residual_percent=100*np.max(np.abs(cum[indices]-assigned)/abs(assigned))
B=np.array(branch_rows)
def mean(B):
    return list(np.trapezoid(B[:,1:],x=B[:,0],axis=0)/(B[-1,0]-B[0,0]))
out=dict(open_line_E0=float(E0),initial_contact_E=float(link0),percent_current_balance_residual=float(residual_percent),passing_runs=[[g[0],g[-1]] for g in runlist],longest=[longest[0],longest[-1]],rows=rows,
         upper_middle_doublet_points=int(np.sum(spatial_weights[near_top&~triple])),
         triple_points=int(np.sum(spatial_weights[triple])),branch_samples=branch_rows,mean_low_high_total=mean(B),common_mean_low_high_total=mean(B[::3]),
         open_loaded_vs_free_periodic_velocity_relative_RMS=rms)
print(json.dumps(out,indent=2))
