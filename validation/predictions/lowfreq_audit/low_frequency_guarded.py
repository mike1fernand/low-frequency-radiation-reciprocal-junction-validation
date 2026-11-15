#!/usr/bin/env python3
"""Independent exact-lattice ray quadrature, real-axis pole subtraction, low-w study.
No archived numerical kernels are imported. numpy only. All outputs are numerical,
not interval enclosures. Full Brillouin cube split into three octant face pyramids.
"""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
GREF=np.array([.1,.2,.2,.2,.2,.1]);CQ=1.3;KQ=CQ**2

def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')

def cached_kernel(path, omega, nang, nrad, script_sha256):
 """Reuse a kernel only when its computation identity matches this request."""
 meta=path.with_suffix('.json')
 if not path.exists() and not meta.exists():
  return None
 def refuse(reason):
  raise ValueError(f'Existing kernel cache {path}: {reason}; choose a fresh --tag to recompute without overwriting it.')
 if not path.exists() or not meta.exists():
  refuse('matrix or metadata is missing')
 try:
  dg=json.loads(meta.read_text())
 except (OSError,json.JSONDecodeError) as exc:
  refuse(f'metadata cannot be read: {exc}')
 if not isinstance(dg,dict):
  refuse('metadata is not a JSON object')
 expected={'omega':omega,'n_face':nang,'n_radial':nrad,'script_sha256':script_sha256}
 differences=[key for key,value in expected.items() if dg.get(key)!=value]
 if differences:
  refuse('computation identity mismatch: '+', '.join(differences))
 if dg.get('matrix_sha256')!=hashlib.sha256(path.read_bytes()).hexdigest():
  refuse('matrix SHA-256 mismatch')
 with np.load(path) as z:
  return z['PV'],[z[f'Q{j}'] for j in range(3)],dg

def gl(n,a=0.,b=1.):
 x,w=np.polynomial.legendre.leggauss(n);return a+(x+1)*(b-a)/2,w*(b-a)/2

def dyn(k,deriv=None):
 s=np.sin(k);z=2*np.sin(k/2)**2;c=1-z
 out=np.zeros(k.shape[:-1]+(3,3))
 for i in range(3):
  j,l=[v for v in range(3) if v!=i]
  out[...,i,i]=2*z[...,i]+2*(z[...,i]+z[...,j]-z[...,i]*z[...,j]+z[...,i]+z[...,l]-z[...,i]*z[...,l])
  for j in range(i):out[...,i,j]=out[...,j,i]=2*s[...,i]*s[...,j]
 if deriv is None:return out
 p=deriv;dp=np.zeros_like(out)
 for i in range(3):
  j,l=[v for v in range(3) if v!=i]
  dp[...,i,i]=2*s[...,i]*p[...,i]+2*(s[...,i]*p[...,i]*(c[...,j]+c[...,l])+c[...,i]*(s[...,j]*p[...,j]+s[...,l]*p[...,l]))
  for j in range(i):dp[...,i,j]=dp[...,j,i]=2*(c[...,i]*p[...,i]*s[...,j]+s[...,i]*c[...,j]*p[...,j])
 return out,dp

BONDS=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],float)
UNIT=BONDS/np.linalg.norm(BONDS,axis=1)[:,None]
def dyn_bonds(k,deriv=None):
 dot=np.einsum('...i,bi->...b',k,BONDS)
 out=np.einsum('...b,bi,bj->...ij',4*np.sin(dot/2)**2,UNIT,UNIT)
 if deriv is None:return out
 d=np.einsum('...i,bi->...b',deriv,BONDS)
 return out,np.einsum('...b,bi,bj->...ij',2*np.sin(dot)*d,UNIT,UNIT)

def projected(M,k,nseg=6):
 # Positive-octant parity-summed real blocks: XX(d), RR(d), XR(d), d>=0.
 cy=np.cos(k[...,1]/2);cz=np.cos(k[...,2]/2)
 ry=np.sin(k[...,1]/2)*cz;rz=cy*np.sin(k[...,2]/2)
 xx=4*cy**2*cz**2*M[...,0,0]
 rr=2*(ry**2*M[...,1,1]+rz**2*M[...,2,2]+2*ry*rz*M[...,1,2])
 xr=2*np.sqrt(2)*cy*cz*(ry*M[...,0,1]+rz*M[...,0,2])
 phase=k[...,0,None]*np.arange(nseg)
 return np.stack([xx[...,None]*np.cos(phase),rr[...,None]*np.cos(phase),xr[...,None]*np.sin(phase)],axis=-1)

def assemble(b):
 n=len(b);M=np.zeros((2*n,2*n))
 for i in range(n):
  for j in range(n):
   d=abs(i-j);sg=np.sign(i-j)
   M[2*i,2*j]=b[d,0];M[2*i+1,2*j+1]=b[d,1]
   M[2*i,2*j+1]=sg*b[d,2];M[2*i+1,2*j]=-sg*b[d,2]
 return M

def roots(p,omega,fun=dyn):
 # p has one direction per row; k=t*p. lambda_j(t) monotone in the root bracket
 # because all |k.d|<pi and every bond derivative is PSD.
 nr=len(p);p3=np.broadcast_to(p[:,None,:],(nr,3,3))
 lo=np.zeros((nr,3));hi=np.broadcast_to((2*omega/np.linalg.norm(p,axis=1))[:,None],(nr,3)).copy()
 vals=np.linalg.eigvalsh(fun(hi[...,None]*p3));v=np.take_along_axis(vals,np.arange(3)[None,:,None],axis=-1)[...,0]
 assert np.all(v>omega**2)
 for _ in range(48):
  t=(lo+hi)/2;ev=np.linalg.eigvalsh(fun(t[...,None]*p3));v=np.take_along_axis(ev,np.arange(3)[None,:,None],axis=-1)[...,0]
  hi=np.where(v>omega**2,t,hi);lo=np.where(v<=omega**2,t,lo)
 t=(lo+hi)/2;k=t[...,None]*p3;D,dD=fun(k,p3);ev,U=np.linalg.eigh(D)
 u=np.take_along_axis(U,np.broadcast_to(np.arange(3)[None,:,None,None],(nr,3,3,1)),axis=-1)[...,0]
 der=np.einsum('...i,...ij,...j->...',u,dD,u)
 assert np.all(der>0)
 residual=np.take_along_axis(ev,np.arange(3)[None,:,None],axis=-1)[...,0]-omega**2
 proj=projected(u[..., :,None]*u[...,None,:],k)
 return t,der,proj,{'max_root_residual':float(abs(residual).max()),'min_radial_derivative':float(der.min()),'max_k_bond_dot':float(abs(np.einsum('...i,bi->...b',k,BONDS)).max())}

def face_quadrature(omega,nang,nrad):
 t0=time.perf_counter();u,wu=gl(nang);t,wt=gl(nrad)
 uu,vv=np.meshgrid(u,u,indexing='ij');aw=(wu[:,None]*wu[None,:]).ravel()
 base=np.stack([np.ones(uu.size),uu.ravel(),vv.ravel()],axis=-1)*np.pi
 ps=[base,base[:,[1,0,2]],base[:,[1,2,0]]]
 P=np.zeros((6,3));I=np.zeros((3,6,3));maxres=0;mind=float('inf');maxdot=0
 for p in ps:
  if omega:
   tr,der,kproj,diag=roots(p,omega)
   residue=(tr**2/der)[...,None,None]*kproj
   subtraction=(2*tr)[...,None,None]*residue
   P+=np.einsum('r,rbdc,rb->dc',aw,residue,np.log((1-tr)/(1+tr)))
   I+=np.pi*np.einsum('r,rbdc->bdc',aw,residue)
   maxres=max(maxres,diag['max_root_residual']);mind=min(mind,diag['min_radial_derivative']);maxdot=max(maxdot,diag['max_k_bond_dot'])
  for a in range(0,len(p),32):
   pp=p[a:a+32];k=pp[:,None,:]*t[None,:,None]
   inv=np.linalg.inv(dyn(k)-omega**2*np.eye(3));v=t[None,:,None,None]**2*projected(inv,k)
   if omega:
    sub=subtraction[a:a+32];r=tr[a:a+32]
    for b in range(3):v-=sub[:,b,None,:,:]/(t[None,:,None,None]**2-r[:,b,None,None,None]**2)
   P+=np.einsum('r,t,rtdc->dc',aw[a:a+32],wt,v)
 G=assemble(P);bs=[assemble(x) for x in I]
 diag={'omega':omega,'n_face':nang,'n_radial':nrad,'seconds':time.perf_counter()-t0,'max_root_residual':maxres,'min_root_derivative':mind if omega else None,'max_root_bond_phase':maxdot,'PV_symmetry':float(np.linalg.norm(G-G.T)),'spectral_min_eigenvalues':[float(np.linalg.eigvalsh(x).min()) for x in bs]}
 return G,bs,diag

def independent_shell(omega,nmu,nphi):
 # Separate spherical-octant parameterization and explicit nine-bond D builder.
 mu,wm=gl(nmu);ph,wp=gl(nphi,0,np.pi/2);mu,ph=np.meshgrid(mu,ph,indexing='ij')
 st=np.sqrt(1-mu*mu);p=np.stack([mu,st*np.cos(ph),st*np.sin(ph)],axis=-1).reshape(-1,3)
 weight=(wm[:,None]*wp[None,:]).ravel()/np.pi**2
 r,der,proj,diag=roots(p,omega,dyn_bonds)
 bs=[assemble(b) for b in np.einsum('r,rb,rbdc->bdc',weight,r*r/der,proj)]
 return bs,diag

def matrices(g):
 n=len(g);U=np.zeros((2*n,2*n));V=np.zeros((2*n,n+1));Ks=np.zeros((n+1,n+1))
 for j in range(n):
  U[2*j:2*j+2,2*j:2*j+2]=2*g[j]/3*np.diag([1.,2.])
  V[2*j:2*j+2,j:j+2]=2*g[j]/3*np.array([[-1.,-1.],[-np.sqrt(2),np.sqrt(2)]])
  Ks[j:j+2,j:j+2]+=4*g[j]/3*np.eye(2)
 L=2*KQ*np.eye(n+1)-KQ*(np.eye(n+1,k=1)+np.eye(n+1,k=-1))
 return U,V,Ks,L

def scatter(omega,s,P,bs,g=GREF):
 U,V,Ks,L=matrices(s*g);n=len(g);Q=sum(bs);G=P+1j*Q
 kap=2*np.arcsin(omega/(2*CQ));gam=2*KQ*np.sin(kap)
 lead=np.zeros((n+1,n+1),complex);lead[0,0]=lead[-1,-1]=-KQ*np.exp(1j*kap)
 gc=np.linalg.solve(np.eye(2*n)+G@U,G)
 GD=np.linalg.inv(omega**2*np.eye(n+1)-L-Ks-lead+V.T@gc@V)
 chi=GD[:,0]*np.sqrt(gam);eta=np.linalg.solve(np.eye(2*n)+U@G,V@chi)
 powers=[float(2*np.real(eta.conj()@b@eta)) for b in bs]
 T=sum(powers);Tq=float(gam*abs(chi[-1])**2);R=float(abs(1-1j*np.sqrt(gam)*chi[0])**2)
 _,V1,_,_=matrices(g);chi0=-1j*np.exp(1j*kap*np.arange(n+1))/np.sqrt(gam);e1=V1@chi0
 AA=[float(2*np.real(e1.conj()@b@e1)) for b in bs]
 return {'omega':omega,'s':s,'TqT':powers[0]+powers[1],'TqL':powers[2],'Tbulk':T,'Tline':Tq,'R':R,'qT_share':(powers[0]+powers[1])/T if T else None,'closure':R+Tq+T-1,'Tbulk_over_omega2':T/omega**2,'AqT':AA[0]+AA[1],'AqL':AA[2],'weak_share':(AA[0]+AA[1])/sum(AA),'effective_resultant_over_sqrt_gamma':complex(2*eta[::2].sum()/np.sqrt(gam)).__str__()}

def static_checks(P,scales):
 rows=[];p=np.tile([2.,0.],6);e=np.ones(7);eL=np.eye(7)[:,0]
 for s in scales:
  U,V,Ks,L=matrices(s*GREF);L[0,0]-=KQ;L[-1,-1]-=KQ
  gc=np.linalg.solve(np.eye(12)+P@U,P);A=L+Ks-V.T@gc@V;Ds=np.linalg.solve(np.eye(12)+U@P,V)
  q=np.linalg.solve(A,eL)
  rows.append({'s':s,'min_static_device_eigenvalue':float(np.linalg.eigvalsh(A).min()),'condition_static_device':float(np.linalg.cond(A)),'force_resultant':float(p@Ds@q),'translation_Ve_Up':float(np.linalg.norm(V@e+U@p)),'translation_Kse_Vtp':float(np.linalg.norm(Ks@e+V.T@p)),'identity_residual':float(np.linalg.norm(e@A+p@Ds)),'Tline_leading_coefficient':float((2*CQ*q[-1])**2)})
 return rows

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--nang',type=int,default=16);ap.add_argument('--nrad',type=int,default=80);ap.add_argument('--frequencies',type=float,nargs='+',default=[0,.4,.2,.1,.05,.025,.0125,.00625]);ap.add_argument('--tag',required=True);ap.add_argument('--shellcheck',action='store_true');a=ap.parse_args()
 out=ROOT/a.tag;out.mkdir(exist_ok=True);status={'status':'running','args':vars(a),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'threads':1,'finished':[]}
 # Check the entire request before replacing any completed status or row files.
 names={}
 for w in a.frequencies:
  name=f'kernel_{w:.8f}.npz'
  if name in names and names[name]!=w:
   raise ValueError(f'Distinct requested frequencies round to {name}; run them under separate fresh --tag values.')
  names[name]=w
 cache={w:cached_kernel(out/f'kernel_{w:.8f}.npz',w,a.nang,a.nrad,status['script_sha256']) for w in a.frequencies}
 dump(out/'status.json',status)
 rows=[];checks=[]
 for w in a.frequencies:
  path=out/f'kernel_{w:.8f}.npz';meta=path.with_suffix('.json')
  cached=cache[w]
  if cached is None:
   P,bs,dg=face_quadrature(w,a.nang,a.nrad);np.savez_compressed(path,PV=P,**{f'Q{j}':b for j,b in enumerate(bs)});dg['matrix_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();dg['script_sha256']=status['script_sha256'];dump(meta,dg)
  else:
   P,bs,dg=cached
  if not w:
   static=static_checks(P,[1.,.1,.01,.001]);dump(out/'static_checks.json',{'min_Gstar_eigenvalue':float(np.linalg.eigvalsh(P).min()),'cases':static})
  else:
   scales=[1.,.1,.01,.001,w*.1,w,w*10]
   for s in scales:rows.append(scatter(w,s,P,bs))
   if a.shellcheck:
    independent,idg=independent_shell(w,a.nang*2,a.nang*2)
    err=[float(np.linalg.norm(x-y)/np.linalg.norm(x)) for x,y in zip(bs,independent)]
    # Individual lower branches can be degenerate; cluster is the physical check.
    checks.append({'omega':w,'branch_relative_differences':err,'qT_relative_difference':float(np.linalg.norm(bs[0]+bs[1]-independent[0]-independent[1])/np.linalg.norm(bs[0]+bs[1])),'qL_relative_difference':err[2],'independent_diagnostics':idg})
  status['finished'].append(w);dump(out/'status.json',status);dump(out/'rows.json',rows);dump(out/'independent_shell_checks.json',checks)
  print(json.dumps({'tag':a.tag,**dg}),flush=True)
 status['status']='completed';dump(out/'status.json',status)
if __name__=='__main__':main()
