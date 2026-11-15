"""Bounded multi-frequency reuse of the accepted corrected tetrahedron geometry.

This changes accumulation order and exploits exact inversion pairing. It does
not change scalar cell integrals, centroid residues, eigenvalues or parameters.
All matrices are derived from source outer products, without PSD projection.
"""
import json,hashlib,time,sys,os
from pathlib import Path
from dataclasses import asdict
import numpy as np
ROOT=Path(__file__).resolve().parent
ACCEPTED=ROOT
sys.path.insert(0,str(ACCEPTED))
import combined_tetra_diagnostic as accepted
from stable_scalar import stable_tetrahedron_resolvent
from centroid_weights import invariant_centroid_weights

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def source_hashes():
    return {str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ACCEPTED/'combined_tetra_diagnostic.py',ACCEPTED/'stable_scalar.py',ACCEPTED/'centroid_weights.py']}

def canonical(data):return json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False)

def prepare_geometry(directory,grid,shift=(0.,0.,0.),p=None):
    p=accepted.Parameters() if p is None else p;p.validate();shift=tuple(float(x) for x in shift)
    if int(grid)!=grid or grid<6 or len(shift)!=3 or not np.isfinite(shift).all():raise ValueError('invalid grid/shift')
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    settings={'grid':int(grid),'shift':list(shift),'parameters':asdict(p),'accepted_source_sha256':sha(ACCEPTED/'combined_tetra_diagnostic.py')}
    mf=directory/'geometry.json';ef=directory/'eigenvalues.npy'
    if mf.exists():
        meta=json.loads(mf.read_text())
        if meta['settings']!=settings or sha(ef)!=meta['eigenvalues_sha256']:raise ValueError('geometry cache identity mismatch')
        return np.load(ef,mmap_mode='r'),meta
    dk=2*np.pi/grid;axes=[-np.pi+dk*(np.arange(grid+1)+s) for s in shift]
    tmp=directory/'eigenvalues.pending.npy';ev=np.lib.format.open_memmap(tmp,mode='w+',dtype=float,shape=(grid+1,grid+1,grid+1,3));start=time.perf_counter()
    for begin in range(0,grid+1,16):
        kv=np.stack(np.meshgrid(axes[0][begin:begin+16],axes[1],axes[2],indexing='ij'),axis=-1)
        ev[begin:begin+16]=np.linalg.eigvalsh(accepted.bulk_stiffness(kv,p)/p.m_bulk)
    ev.flush();del ev;os.replace(tmp,ef)
    meta={'settings':settings,'eigenvalues_sha256':sha(ef),'seconds':time.perf_counter()-start}
    mf.write_text(json.dumps(meta,indent=2)+'\n');return np.load(ef,mmap_mode='r'),meta

def integrate_range(eigenvalues,omega,grid,start,stop,*,shift=(0.,0.,0.),p=None,nseg=6,regulator=1e-10,cube_chunk=1024):
    p=accepted.Parameters() if p is None else p;p.validate();omega=np.asarray(omega,float)
    if omega.ndim!=1 or len(omega)==0 or np.any(~np.isfinite(omega)) or np.any(omega<=0) or len(np.unique(omega))!=len(omega):raise ValueError('distinct finite positive frequencies required')
    if int(grid)!=grid or grid<6 or int(nseg)!=nseg or nseg<1 or int(cube_chunk)!=cube_chunk or cube_chunk<1:raise ValueError('invalid grid/nseg/cube_chunk')
    if np.asarray(eigenvalues).shape!=(grid+1,grid+1,grid+1,3):raise ValueError('geometry array shape mismatch')
    if np.asarray(shift).shape!=(3,) or not np.isfinite(shift).all():raise ValueError('invalid shift')
    if not (0<=start<stop<=grid**3) or int(start)!=start or int(stop)!=stop:raise ValueError('invalid cube range')
    if regulator<=0 or not np.isfinite(regulator):raise ValueError('invalid regulator')
    st=time.perf_counter();dk=2*np.pi/grid;norm=1/(grid**3*p.m_bulk);shift=np.asarray(shift,float)
    sums=np.zeros((len(omega),3,2,nseg,2,2),complex);mins=np.full(len(omega),np.inf);maxs=np.zeros(len(omega));clusters={};neg=np.zeros(len(omega),int)
    for begin in range(start,stop,cube_chunk):
        qq=np.arange(begin,min(begin+cube_chunk,stop));base=np.stack([qq//(grid*grid),(qq//grid)%grid,qq%grid],axis=-1);base_k=-np.pi+dk*(base+shift)
        for tet in accepted.TETRA:
            offsets=accepted.CUBE_VERTICES[tet];idx=base[:,None,:]+offsets[None,:,:];ev=eigenvalues[idx[...,0],idx[...,1],idx[...,2],:]
            centroid=base_k+dk*offsets.mean(axis=0);vals,vec=np.linalg.eigh(accepted.bulk_stiffness(centroid,p)/p.m_bulk);F=accepted.active_form_factor(centroid)
            principal=(centroid+np.pi)%(2*np.pi)-np.pi;kn=np.linalg.norm(principal,axis=1);kh=np.divide(principal,kn[:,None],out=np.zeros_like(principal),where=kn[:,None]>1e-14)
            W,WL,cc=invariant_centroid_weights(vals,vec,F,kh)
            for k,v in cc.items():clusters[k]=clusters.get(k,0)+v
            phase=np.exp(1j*centroid[:,0,None]*np.arange(nseg)[None,:])
            for branch in range(3):
                # Exact reciprocity average = explicit ±dx integration followed
                # by transpose averaging, not discarding an imaginary density.
                residue=np.concatenate([(phase[:,:,None,None]*x[:,None,branch]).real.reshape(len(qq),-1) for x in [W,WL]],axis=1)
                factors=np.stack([stable_tetrahedron_resolvent(ev[:,:,branch],complex(w*w,regulator)) for w in omega])
                mins=np.minimum(mins,factors.imag.min(axis=1));maxs=np.maximum(maxs,abs(factors).max(axis=1));neg+=(factors.imag<0).sum(axis=1)
                if not np.isfinite(factors).all():raise FloatingPointError('nonfinite scalar factor')
                sums[:,branch]+=(norm*(factors@residue)).reshape(len(omega),2,nseg,2,2)
    return {'sums':sums,'minimum_scalar_imaginary':mins,'maximum_abs_scalar':maxs,'negative_scalar_count':neg,'centroid_clusters':clusters,'seconds':time.perf_counter()-st,'start':start,'stop':stop}

def assemble(blocks,nseg):
    out=np.empty((2*nseg,2*nseg),complex)
    for s in range(nseg):
        for t in range(nseg):
            dx=s-t;out[2*s:2*s+2,2*t:2*t+2]=blocks[dx] if dx>=0 else blocks[-dx].T
    return out

def result_for_frequency(omega,grid,regulator,sums,diagnostics,nseg=6):
    b=[assemble(sums[i,0],nseg) for i in range(3)];bl=[assemble(sums[i,1],nseg) for i in range(3)];G=sum(b);im=[x.imag for x in b];il=sum(x.imag for x in bl)
    mineig=lambda x:float(np.linalg.eigvalsh((x+x.conj().T)/2).min())
    diagnostics=dict(diagnostics)|{'spectral_min_eigenvalue':mineig(G.imag),'branch_min_eigenvalues':[mineig(x) for x in im],'geometric_long_min_eigenvalue':mineig(il),'geometric_transverse_min_eigenvalue':mineig(G.imag-il),'reciprocity_residual_constructed':float(np.linalg.norm(G-G.T)/np.linalg.norm(G))}
    if min(diagnostics['spectral_min_eigenvalue'],*diagnostics['branch_min_eigenvalues'],diagnostics['geometric_long_min_eigenvalue'],diagnostics['geometric_transverse_min_eigenvalue']) < -1e-12*max(1,np.linalg.norm(G)) :raise ArithmeticError('spectral positivity failed')
    return accepted.TetrahedronResult(omega=float(omega),energy_imaginary_regulator=float(regulator),grid=int(grid),nseg=nseg,G0=G,Im_branches=im,Im_longitudinal_projector=il,Im_branch_longitudinal_projectors=[x.imag for x in bl],symmetry_residual=diagnostics['reciprocity_residual_constructed'],min_imaginary_density_eigenvalue=diagnostics['spectral_min_eigenvalue'],seconds=diagnostics.get('seconds',0),diagnostics=diagnostics)
