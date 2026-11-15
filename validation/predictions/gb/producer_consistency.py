import json,time
from pathlib import Path
import numpy as np
import p1solver as p
root=Path(__file__).resolve().parent
rng=np.random.default_rng(1301);k=np.concatenate([rng.uniform(-np.pi,np.pi,(50,3)),np.array([[0,0,0],[np.pi,np.pi,0],[np.pi,np.pi,np.pi],[0,np.pi,0]])])
A0,A1=p.bulk_blocks(k[:,1],k[:,2]);D=A0+A1*np.exp(1j*k[:,0,None,None])+A1.conj()*np.exp(-1j*k[:,0,None,None]);Db=p.D_full(k)
rows=[]
# This is a smooth finite-regularization test, not a small-regulator error bound.
omega=1.2;eta=.1;wr=np.sqrt(omega**2-eta**2);eps=2*omega*eta
Ga=p.green_blocks(wr,6,256,eps=eps,onshell=False,quadrant=True)['G']
for N in [32,64,128]:
 t=time.time();Gd=p.green_blocks_direct(wr,eps,6,N)
 rows.append({'N3d':N,'max_abs_difference':float(abs(Ga-Gd).max()),'relative_frobenius':float(np.linalg.norm(Ga-Gd)/np.linalg.norm(Ga)),'seconds':time.time()-t})
a={'bulk_symbol_sample_count':len(k),'max_symbol_difference':float(abs(D-Db).max()),'smooth_regulator_test':{'omega':omega,'eta':eta,'regulator':'z=(omega+i*eta)^2','nseparations':6,'axial_Nperp':256,'shift':[.5,.5],'rows':rows},'limitations':'Same physical bulk and basis; independent direct 3D quadrature of inverse dynamical matrix versus axial-root residues. Smooth finite-regulator check is not an error enclosure for the small-regulator extrapolation.'}
(root/'producer_consistency.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2))
