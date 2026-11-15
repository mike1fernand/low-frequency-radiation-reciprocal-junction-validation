import json,hashlib
from pathlib import Path
import numpy as np
from core import roots,symbol,source,states,assembly,E,R

P=Path(__file__).resolve().parent
eta,reference=states();_,_,B=assembly()
axial=[]
for sgn in [1,-1]:
 for branch,kx in [(0,sgn*2*np.arcsin(.6)),(2,sgn*2*np.arcsin(1.2/np.sqrt(12)))]:
    k=np.array([kx,0.,0.]);val,vec=np.linalg.eigh(symbol(k));PT=vec[:,:2]@vec[:,:2].T
    row={'k':k.tolist(),'shell_branch':branch,'shell_residual':float(abs(val[branch]-E)),'states':{}}
    for name,t in eta.items():
        f=source(k,t);w=float(np.vdot(f,PT@f).real)
        assert abs(w)<1e-28
        row['states'][name]={'qT_projector_weight':w,'longitudinal_weight':float(abs(f[0])**2)}
    axial.append(row)
rng=np.random.default_rng(20260926)
# Independent literal Fourier sum on the 72 corner coordinates.
positions=np.array([(s,sy/2,sz/2) for s in range(6) for sy in [-1,1] for sz in [-1,1]])
source_error=0.;derivative_error=0.
for k in rng.uniform(-1.3,1.3,size=(40,3)):
    for t in list(eta.values())+[rng.normal(size=12)+1j*rng.normal(size=12)]:
        literal=(np.exp(-1j*positions@k)[:,None]*(B@t).reshape(24,3)).sum(axis=0)
        source_error=max(source_error,float(abs(literal-source(k,t)).max()))
    D,J=symbol(k,True)
    for axis in range(3):
        d=np.eye(3)[axis]*1e-6
        derivative_error=max(derivative_error,float(abs((symbol(k+d)-symbol(k-d))/(2e-6)-J[axis]).max()))
assert source_error<1e-13 and derivative_error<2e-9
mirrors=[]
for kx in [.3,.6,.9,1.1]:
    ky=np.arccos(2-E/2-np.cos(kx));k=np.array([kx,ky,0.]);vals,vec=np.linalg.eigh(symbol(k));T=vec[:,:2]
    angle=.731;O=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    row={'k':k.tolist(),'lower_shell_residual':float(abs(vals[:2]-E).max()),'states':{}}
    for name,t in eta.items():
        f=source(k,t);p=float(np.linalg.norm(T.conj().T@f)**2);pr=float(np.linalg.norm((T@O).conj().T@f)**2)
        row['states'][name]={'z_polarization_weight':float(abs(f[2])**2),'qT_cluster_weight':p,'rotation_residual':abs(pr-p)}
    mirrors.append(row)
# A second whole-domain cone-support check uses midpoint azimuths rather than
# Gaussian nodes, and explicitly tests opposite-kx ray contributions.
phi=(np.arange(192)+.5)*np.pi/(4*192);theta=np.linspace(0,np.pi/2,1537)
tt,pp=np.meshgrid(theta,phi,indexing='ij')
dirs=np.stack([np.cos(tt),np.sin(tt)*np.cos(pp),np.sin(tt)*np.sin(pp)],axis=-1)
preimages=[]
for branch in range(3):
    rr=roots(dirs,branch,iterations=40);g=rr['grad'];c=g[...,0]/np.linalg.norm(g,axis=-1)
    entry={'branch':branch,'minimum_forward_cosine':float(c.min()),'max_abs_shell_residual':float(abs(rr['values'][...,branch]-E).max()),'cones':[]}
    for deg in [.5,1,2,5,10,20,30]:
        threshold=np.cos(np.deg2rad(deg));inside=c>=threshold;delta=np.diff(inside.astype(int),axis=0)
        e={'angle':deg,'min_exits':int((delta==-1).sum(axis=0).min()),'max_exits':int((delta==-1).sum(axis=0).max()),
           'max_entries':int((delta==1).sum(axis=0).max()),'opposite_kx_preimage_nodes':int((c<=-threshold).sum())}
        assert e['min_exits']==e['max_exits']==1 and e['max_entries']==e['opposite_kx_preimage_nodes']==0
        entry['cones'].append(e)
    preimages.append(entry)
out={'status':'passed','literal_source_max_difference':source_error,'symbol_derivative_fd_max_difference':derivative_error,
 'axial_projector_checks':axial,'mirror_plane_checks':mirrors,'cone_support_check':{'theta_count':1537,'azimuth_count_sector':192,'branches':preimages},
 'limits':'All three radial sheets covered using the analytic census. Cone-support topology independently sampled over full hemispheres with midpoint azimuths; no interval proof excludes arbitrarily narrow unsampled angular islands. No exact all-preimage/dark-cone theorem is claimed.'}
(P/'SELECTION_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='mirror_plane_checks'},indent=2))
