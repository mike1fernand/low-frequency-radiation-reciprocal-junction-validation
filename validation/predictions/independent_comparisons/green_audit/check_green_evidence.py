"""Focused audit: recompute stored extrapolation diagnostics, and test axial-root
integration against a fresh full 3D resolvent using literal Cartesian form factors.
No epsilon-zero error certificate is inferred from the finite-width comparison.
"""
from pathlib import Path
import sys,json,ast,time,hashlib
import numpy as np
BASE=Path(__file__).resolve().parents[1]
INPUT=BASE/'input'
sys.path.insert(0,str(INPUT))
import p1solver as ps
out={'stored_fits':[], 'finite_width_comparison':[], 'source':{}}
for f in sorted((INPUT/'rev').glob('g0_*.npz')):
 d=np.load(f); nmax=int(d['nmax']); om=float(d['omega']); eps=d['eps']
 im=(d['ImT_240']+d['ImL_240'])[nmax-1:].real
 ent={'file':str(f.relative_to(BASE)), 'omega':om,'nmax':nmax,'epsilon':eps.tolist(),'fits':[],
      'shell_160_240_maxabs_total':float(np.max(np.abs(d['ImT_160']+d['ImL_160']-d['ImT_240']-d['ImL_240'])))}
 for k in sorted(d.files):
  if not k.startswith('scan_'):continue
  N=int(k.split('_')[1]); a=d[k]
  for degree in (2,3):
   V=np.vander(1j*eps,degree+1,increasing=True)
   c,*_=np.linalg.lstsq(V,a.reshape(len(eps),-1),rcond=None)
   c0=c[0].reshape(a.shape[1:]); stored=d[f'c0_{N}_deg{degree}']
   ent['fits'].append({'Nperp':N,'degree':degree,'refit_vs_saved_maxabs':float(np.max(np.abs(c0-stored))),
      'extrap_Im_vs_shell_maxabs':float(np.max(np.abs(c0.imag-im))),
      'fit_maxabs_residual':float(np.max(np.abs(V@c-a.reshape(len(eps),-1)))),
      'Re_deg2_vs_deg3':float(np.max(np.abs(d[f'c0_{N}_deg2'].real-d[f'c0_{N}_deg3'].real)))})
 if 'scan_4096' in d:
  ent['Re_2048_vs_4096_degree2_maxabs']=float(np.max(np.abs(d['c0_2048_deg2'].real-d['c0_4096_deg2'].real)))
  ent['Re_2048_vs_4096_degree3_maxabs']=float(np.max(np.abs(d['c0_2048_deg3'].real-d['c0_4096_deg3'].real)))
 out['stored_fits'].append(ent)
# Stored reference frequency: source metadata omissions are recorded separately.
st={ast.literal_eval(k):v for k,v in np.load(INPUT/'eps_scan_n48.npy',allow_pickle=True).item().items()}
sp=np.load(INPUT/'onshell_w1.2_n48.npy',allow_pickle=True).item(); nmax=int(sp['nmax'])
refs=[]
for N in (1536,2048):
 eps=np.array([.006,.008,.012,.016,.024]);a=np.array([st[(e,N)] for e in eps]);Im=(sp['T']+sp['L'])[nmax-1:].real
 for degree in (2,3):
  V=np.vander(1j*eps,degree+1,increasing=True);c,*_=np.linalg.lstsq(V,a.reshape(len(eps),-1),rcond=None);c0=c[0].reshape(a.shape[1:]);
  refs.append({'Nperp':N,'degree':degree,'nmax':nmax,'extrap_Im_vs_shell_maxabs_all48':float(np.max(np.abs(c0.imag-Im))),
      'extrap_Im_vs_shell_maxabs_first6':float(np.max(np.abs(c0[:6].imag-Im[:6])))})
out['reference_frequency_1p2']=refs
# Independent literal Cartesian construction. Each of four corners has a normalized
# axial displacement 1/2 or radial displacement (0, sy, sz)/(2 sqrt(2)).
def literal_D_and_phi(k):
 B=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],float)
 E=B/np.linalg.norm(B,axis=1)[:,None]
 D=np.einsum('...l,li,lj->...ij',2*(1-np.cos(k@B.T)),E,E)
 P=np.zeros(k.shape[:-1]+(3,2),complex)
 for sy in (-1,1):
  for sz in (-1,1):
   phase=np.exp(-.5j*(sy*k[...,1]+sz*k[...,2]))
   v=np.array([[.5,0],[0,sy/(2*np.sqrt(2))],[0,sz/(2*np.sqrt(2))]])
   P+=phase[...,None,None]*v
 return D,P
rng=np.random.default_rng(512);k=rng.uniform(-np.pi,np.pi,(300,3));D,P=literal_D_and_phi(k)
A0,A1=ps.bulk_blocks(k[:,1],k[:,2]);Db=A0+A1*np.exp(1j*k[:,0,None,None])+A1.conj()*np.exp(-1j*k[:,0,None,None])
out['source']['D_bond_vs_axial_blocks_maxabs']=float(np.max(np.abs(D-Db)))
out['source']['literal_phi_vs_solver_maxabs']=float(np.max(np.abs(P-ps.phi_vectors(k[:,1],k[:,2]))))
def direct(om,eps,nmax,N):
 ks=-np.pi+(np.arange(N)+.5)*2*np.pi/N;ky,kz=np.meshgrid(ks,ks,indexing='ij');G=np.zeros((nmax,2,2),complex)
 for kx in ks:
  k=np.stack([np.full_like(ky,kx),ky,kz],axis=-1);D,P=literal_D_and_phi(k)
  R=np.linalg.inv(D-(om**2+1j*eps)*np.eye(3));blk=np.einsum('yzca,yzcd,yzdb->ab',P.conj(),R,P)
  G+=np.exp(1j*kx*np.arange(nmax))[:,None,None]*blk
 return G/N**3
om=1.3;eps=.2;nmax=6
ar=ps.green_blocks(om,nmax,256,eps=eps,onshell=False,quadrant=True)
full=ps.green_blocks(om,nmax,128,eps=eps,onshell=False,quadrant=False)
quad=ps.green_blocks(om,nmax,128,eps=eps,onshell=False,quadrant=True)
out['source']['quadrant_full_N128_maxabs']=float(np.max(np.abs(full['G']-quad['G'])))
out['source']['axial_N128_256_maxabs']=float(np.max(np.abs(ar['G']-quad['G'])))
out['source']['axial_N256_stats']={k: (float(v) if np.isfinite(v) else None) for k,v in ar['stats'].items()}
last=None
for N in (64,96,128):
 t=time.time(); G=direct(om,eps,nmax,N)
 row={'omega':om,'epsilon_squared_frequency':eps,'N':N,'nmax':nmax,'maxabs_vs_axial_N256':float(np.max(np.abs(G-ar['G']))),'wall_seconds':time.time()-t}
 if last is not None:row['maxabs_vs_previous_direct']=float(np.max(np.abs(G-last)))
 last=G;out['finite_width_comparison'].append(row); print(row,flush=True)
out['notes']=['Saved fits recomputed from supplied arrays; no historical 2048/4096 scans rerun.',
 'Fresh full 3D comparison is only at z=1.3^2+0.2i and tests model/form-factor/contour agreement, not epsilon-zero extrapolation.',
 'Imaginary Green blocks in production assembly are set to the same shell used for branch powers; closure is not an independent quadrature certificate.',
 'No shell-root completeness proof or rigorous quadrature/error enclosure is furnished by this test.']
(BASE/'green_audit'/'green_evidence.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out['source'],indent=2),flush=True)
