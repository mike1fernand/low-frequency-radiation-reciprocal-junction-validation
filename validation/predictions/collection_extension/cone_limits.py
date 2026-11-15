"""Low-frequency-first angular prediction, checked against the discrete model.

Uses bound low-frequency Green matrices, literal 24-point Fourier forces, and an
explicit nine-bond symbol. Finite numerical angular scans are not certificates.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,json,argparse,time,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'baseline/Prediction_Calculations_2026-09-26'
sys.path.insert(0,str(BASE/'low_frequency'))
import low_frequency as lf
B=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],float)
N=B/np.linalg.norm(B,axis=1)[:,None]
POS=np.array([[s,sy/2,sz/2] for s in range(6) for sy in [-1,1] for sz in [-1,1]])
def symbol(k):
    dot=np.einsum('...i,bi->...b',k,B)
    D=np.einsum('...b,bi,bj->...ij',4*np.sin(dot/2)**2,N,N)
    J=np.einsum('...b,ba,bi,bj->...aij',2*np.sin(dot),B,N,N)
    return D,J
def ray(theta,phi):
    t,p=np.broadcast_arrays(theta,phi)
    return np.stack([np.cos(t),np.sin(t)*np.cos(p),np.sin(t)*np.sin(p)],axis=-1)
def shell(direction,j,w):
    lo=np.zeros(direction.shape[:-1]);hi=np.full_like(lo,2*w)
    assert np.all(np.linalg.eigvalsh(symbol(hi[...,None]*direction)[0])[...,j]>w*w)
    for _ in range(43):
        mid=(lo+hi)/2;v=np.linalg.eigvalsh(symbol(mid[...,None]*direction)[0])[...,j]
        hi=np.where(v>w*w,mid,hi);lo=np.where(v<=w*w,mid,lo)
    r=(lo+hi)/2;k=r[...,None]*direction;D,J=symbol(k);ev,U=np.linalg.eigh(D);u=U[...,j]
    grad=np.einsum('...i,...aij,...j->...a',u,J,u)
    radial=np.sum(direction*grad,axis=-1)
    assert np.all(radial>0)
    return k,r,u,grad,radial
def sources(w,tag):
    path=BASE/f'low_frequency/{tag}/kernel_{w:.8f}.npz'
    z=np.load(path,allow_pickle=False);Q=[z[f'Q{j}'] for j in range(3)];G=z['PV']+1j*sum(Q)
    kap=2*np.arcsin(w/(2*lf.CQ));gamma=2*lf.KQ*np.sin(kap)
    results={}
    for kind in ['full','R']:
        U,V,Ks,K=lf.matrices(lf.GREF)
        if kind=='R':
            U[::2,:]=0;V[::2]=0;Ks[:]=0
            for j,g in enumerate(lf.GREF):Ks[j:j+2,j:j+2]+=2*g/3*np.array([[1.,-1.],[-1.,1.]])
        gc=np.linalg.solve(np.eye(12)+G@U,G)
        A=w*w*np.eye(7)-K-Ks+V.T@gc@V
        A[0,0]+=lf.KQ*np.exp(1j*kap);A[-1,-1]+=lf.KQ*np.exp(1j*kap)
        chi=np.linalg.solve(A,np.eye(7)[:,0]*np.sqrt(gamma))
        eta=np.linalg.solve(np.eye(12)+U@G,V@chi)
        force=np.array([[eta[2*s]/2,sy*eta[2*s+1]/(2*np.sqrt(2)),sz*eta[2*s+1]/(2*np.sqrt(2))] for s in range(6) for sy in [-1,1] for sz in [-1,1]])
        results[kind]={'force':force,'global':[float(2*np.real(eta.conj()@(Q[0]+Q[1])@eta)),float(2*np.real(eta.conj()@Q[2]@eta))]}
    return results,{'path':str(path.relative_to(HERE.parent)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
def density(out,force):
    k,r,u,grad,radial=out
    f=np.einsum('...p,pi->...i',np.exp(-1j*np.einsum('...i,pi->...p',k,POS)),force)
    return r*r/radial/(4*np.pi**2)*abs(np.sum(u.conj()*f,axis=-1))**2
def fractions(angle,kind):
    # Stable in d=1-c to avoid subtraction at very small angles.
    d=2*np.sin(np.deg2rad(angle)/2)**2
    if kind=='full':return np.array([d*d*(3-d)/4,d*(3-3*d+d*d)/2])
    return np.array([d*d*(15-25*d+15*d*d-3*d**3)/4,d**3*(20-15*d+3*d*d)/16])
def integrate(w,order,tag='validation'):
    start=time.monotonic();src,inp=sources(w,tag)
    x,wp=np.polynomial.legendre.leggauss(order);phi=(x+1)*np.pi/8;wp*=np.pi # C4v sector times8
    z,wz=np.polynomial.legendre.leggauss(order)
    angles=[10.,20.,30.];powers={a:{k:np.zeros(3) for k in src} for a in angles};scan_records=[]
    tt=np.linspace(0,np.pi,385)[:,None]
    for j in range(3):
        out=shell(ray(tt,phi[None,:]),j,w);cg=out[3][...,0]/np.linalg.norm(out[3],axis=-1)
        for a in angles:
            inside=cg>=np.cos(np.deg2rad(a));change=np.diff(inside.astype(int),axis=0)
            assert inside[0].all() and not inside[-1].any() and np.all((change==-1).sum(axis=0)==1) and not (change==1).any()
            ix=np.argmax(change==-1,axis=0);lo=tt[:,0][ix];hi=tt[:,0][ix+1]
            for _ in range(32):
                t=(lo+hi)/2;rr=shell(ray(t,phi),j,w);c=rr[3][...,0]/np.linalg.norm(rr[3],axis=-1)
                mask=c>=np.cos(np.deg2rad(a));lo=np.where(mask,t,lo);hi=np.where(mask,hi,t)
            theta=(lo+hi)/2;lower=np.cos(theta)
            mu=lower[None,:]+(z[:,None]+1)*(1-lower[None,:])/2
            wt=wz[:,None]*(1-lower[None,:])/2*wp[None,:]
            rr=shell(ray(np.arccos(mu),phi[None,:]),j,w)
            for kind in src:powers[a][kind][j]=float(np.sum(wt*density(rr,src[kind]['force'])))
            scan_records.append({'branch':j,'angle':a,'boundary_degrees':[float(np.rad2deg(theta.min())),float(np.rad2deg(theta.max()))]})
    cR=json.loads((BASE/'low_frequency/sector_predictions.json').read_text())
    coeff={'full':np.array([2*1.3/(3*np.pi),1.3/(3*np.pi*3**1.5)]),'R':np.array([cR['TqT_over_omega4'],cR['TqL_over_omega4']])}
    rows=[]
    for a in angles:
        for kind in src:
            actual=np.array([powers[a][kind][:2].sum(),powers[a][kind][2]])
            pred=coeff[kind]*fractions(a,kind)*w**(2 if kind=='full' else 4)
            rows.append({'angle':a,'kind':kind,'powers':actual.tolist(),'leading_prediction':pred.tolist(),'relative_corrections':(actual/pred-1).tolist(),'qT_share':float(actual[0]/actual.sum()),'predicted_share':float(pred[0]/pred.sum()),'actual_fractions_of_global':(actual/src[kind]['global']).tolist(),'limiting_fractions_of_global':fractions(a,kind).tolist()})
    return {'omega':w,'order':order,'source_tag':tag,'source_input':inp,'rows':rows,'scan':scan_records,'seconds':time.monotonic()-start,'scope':'Finite discrete-lattice samples approaching a low-frequency-first cone law. Polar support search numerical, not interval complete.'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--omega',type=float,required=True);ap.add_argument('--order',type=int,default=24);a=ap.parse_args()
    result=integrate(a.omega,a.order);name=f'cone_w{a.omega:.8f}_n{a.order}'.replace('.','p')+'.json'
    (HERE/name).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
