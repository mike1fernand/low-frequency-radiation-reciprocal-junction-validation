"""Independent reference-symbol, source assembly, and radial-shell arithmetic."""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parent
OMEGA=1.2
E=OMEGA**2
R=np.arccos(.25)

def symbol(k, derivatives=False):
    k=np.asarray(k);c=np.cos(k);s=np.sin(k)
    D=2*s[..., :,None]*s[...,None,:]
    diag=6-2*c*(1+c.sum(axis=-1,keepdims=True)-c)
    for i in range(3):D[...,i,i]=diag[...,i]
    if not derivatives:return D
    J=np.zeros(k.shape[:-1]+(3,3,3)) # derivative coordinate, matrix i,j
    for a in range(3):
        for i in range(3):
            J[...,a,i,i]=2*s[...,a]*(1+c.sum(axis=-1)-c[...,a]) if a==i else 2*c[...,i]*s[...,a]
            for j in range(i):
                x=2*c[...,a]*s[...,j] if a==i else 2*c[...,a]*s[...,i] if a==j else 0
                J[...,a,i,j]=J[...,a,j,i]=x
    return D,J

def assembly():
    B=np.zeros((72,12)); HX=np.zeros((79,79)); HR=HX.copy(); Hp=HX.copy()
    profile=[.1,.2,.2,.2,.2,.1]
    for seg,g in enumerate(profile):
        for ci,(sy,sz) in enumerate([(-1,-1),(-1,1),(1,-1),(1,1)]):
            row=12*seg+3*ci;r=7+row
            B[row,2*seg]=.5;B[row+1,2*seg+1]=sy/(2*np.sqrt(2));B[row+2,2*seg+1]=sz/(2*np.sqrt(2))
            av=np.zeros(79);av[r]=1;av[seg]=av[seg+1]=-.5
            bv=np.zeros(79);bv[r+1]=sy;bv[r+2]=sz;bv[seg]=-.5;bv[seg+1]=.5
            HX+=2*g/3*np.outer(av,av);HR+=2*g/3*np.outer(bv,bv)
            for ep,sgn in [(seg,1),(seg+1,-1)]:
                n=np.array([sgn,sy,sz])/np.sqrt(3);v=np.zeros(79);v[ep]=-n[0];v[r:r+3]=n
                Hp+=g*np.outer(v,v)
    assert np.max(abs(Hp-HX-HR))<1e-15
    def blocks(H):return B.T@H[7:,7:]@B,B.T@H[7:,:7],H[:7,:7]
    return blocks(HX),blocks(HR),B

def states(kernel_path=None):
    if kernel_path is None:
        bound=json.loads((ROOT/'inputs/REFERENCE_KERNEL.json').read_text())
        assert hashlib.sha256((ROOT/'inputs/REFERENCE_KERNEL.npz').read_bytes()).hexdigest()==bound['kernel_sha256']
    z=np.load(kernel_path or ROOT/'inputs/REFERENCE_KERNEL.npz');G=z['G0'][:12,:12]
    Q=np.stack([z[f'Im_branch_{j}'][:12,:12] for j in range(3)])
    (Ux,Vx,Kx),(Ur,Vr,Kr),B=assembly()
    kap=2*np.arcsin(OMEGA/(2*np.sqrt(1.69)));gamma=2*1.69*np.sin(kap)
    KD=1.69*(2*np.eye(7)-np.eye(7,k=1)-np.eye(7,k=-1))
    lead=np.zeros((7,7),complex);lead[0,0]=lead[-1,-1]=-1.69*np.exp(1j*kap)
    out={};reference={}
    for name,a,b in [('full',1,1),('X_only',1,0),('R_only',0,1)]:
        U=a*Ux+b*Ur;V=a*Vx+b*Vr;Ks=a*Kx+b*Kr
        gc=np.linalg.solve(np.eye(12)+G@U,G);sig=-V.T@gc@V
        gd=np.linalg.inv(E*np.eye(7)-KD-Ks-lead-sig)
        chi=gd[:,0]*np.sqrt(gamma);A=np.linalg.inv(np.eye(12)+U@G)
        eta=A@V@chi;out[name]=eta
        powers=np.array([2*np.real(eta.conj()@q@eta) for q in Q])
        rq=abs(1-1j*gamma*gd[0,0])**2;tq=gamma**2*abs(gd[-1,0])**2
        reference[name]={'Rq':float(rq),'Tq':float(tq),'qT':float(powers[:2].sum()),'qL':float(powers[2]),'branch':powers.tolist(),'closure':float(abs(rq+tq+powers.sum()-1))}
        if name=='full':out['same_X']=A@Vx@chi;out['same_R']=A@Vr@chi
    assert np.max(abs(out['full']-out['same_X']-out['same_R']))<1e-15
    return out,reference

def source(k,eta):
    k=np.asarray(k);phase=np.exp(-1j*k[...,0,None]*np.arange(6))
    amp=np.einsum('...s,sa->...a',phase,np.asarray(eta).reshape(6,2))
    y=k[...,1]/2;z=k[...,2]/2
    f=np.zeros(k.shape,complex)
    f[...,0]=2*np.cos(y)*np.cos(z)*amp[...,0]
    f[...,1]=-np.sqrt(2)*1j*np.sin(y)*np.cos(z)*amp[...,1]
    f[...,2]=-np.sqrt(2)*1j*np.cos(y)*np.sin(z)*amp[...,1]
    return f

def roots(n,branch,iterations=44):
    n=np.asarray(n);lo=np.zeros(n.shape[:-1]);hi=R/np.max(abs(n),axis=-1)
    for _ in range(iterations):
        mid=(lo+hi)/2;values=np.linalg.eigvalsh(symbol(mid[...,None]*n))[...,branch]
        above=values>E;hi=np.where(above,mid,hi);lo=np.where(above,lo,mid)
    radius=(lo+hi)/2;k=radius[...,None]*n
    D,J=symbol(k,True);values,vectors=np.linalg.eigh(D);v=vectors[..., :,branch]
    grad=np.einsum('...i,...aij,...j->...a',v,J,v)
    return {'r':radius,'k':k,'values':values,'vectors':vectors,'v':v,'grad':grad,'radial':np.sum(grad*n,axis=-1)}

def density(root,eta):
    amp=np.sum(root['v'].conj()*source(root['k'],eta),axis=-1)
    return 2*np.pi/(2*np.pi)**3*root['r']**2/root['radial']*abs(amp)**2
