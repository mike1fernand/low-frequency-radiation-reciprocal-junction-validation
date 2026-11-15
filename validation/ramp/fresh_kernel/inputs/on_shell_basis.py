#!/usr/bin/env python3
"""Zero-broadening matrix spectral measure by iso-frequency surface integration."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from scipy.ndimage import map_coordinates
from skimage.measure import marching_cubes
from model_api_v1_1 import ModelParameters
from reduced_core import bulk_stiffness_batch,cross_section_form_factor


def _extend(a):
    b=np.concatenate([a,a[:1]],axis=0);b=np.concatenate([b,b[:,:1]],axis=1);return np.concatenate([b,b[:,:,:1]],axis=2)

def _interp_complex(a,coords):
    ae=_extend(a)
    return map_coordinates(ae.real,coords,order=1,mode='wrap')+1j*map_coordinates(ae.imag,coords,order=1,mode='wrap')

def _assemble_raw(pos,nseg):
    out=np.empty((2*nseg,2*nseg),complex)
    for s in range(nseg):
        for t in range(nseg):
            dx=s-t;b=pos[dx] if dx>=0 else pos[-dx].conj().T
            out[2*s:2*s+2,2*t:2*t+2]=b
    return out

def _assemble(pos,nseg):
    out=_assemble_raw(pos,nseg)
    # Exact real-reciprocal harmonic sector.
    return np.real(.5*(out+out.T))


def _degenerate_spectral_weights(vals, vec, coeff, kh):
    """Basis-invariant nodal representatives at backward-error-level ties.

    64*eps*max(1,max|lambda|) identifies floating equality, not a physical
    band-merging scale. Each member receives 1/rank of the complete cluster
    weight, so the sum retains the full spectral projector exactly.
    """
    W=np.einsum('...mA,...mB->...mAB',coeff.conj(),coeff,optimize=True)
    projected_k=np.einsum('...am,...a->...m',vec.conj(),kh,optimize=True)
    WL=W*abs(projected_k)[...,None,None]**2
    tie_tol=64*np.finfo(float).eps*np.maximum(1.,np.max(abs(vals),axis=-1))
    tie01=abs(vals[...,0]-vals[...,1])<=tie_tol
    tie12=abs(vals[...,1]-vals[...,2])<=tie_tol
    for group,mask in [((0,1),tie01&~tie12),
                       ((1,2),tie12&~tie01),
                       ((0,1,2),tie01&tie12)]:
        if not np.any(mask):continue
        rank=len(group)
        cluster_weight=W[mask][:,group].sum(axis=1)/rank
        # Keep all coherences in F^dagger P_C kh kh^dagger P_C F.
        c=coeff[mask][:,group];d=projected_k[mask][:,group]
        longitudinal_coeff=np.einsum('nm,nmA->nA',d.conj(),c)
        cluster_long=np.einsum('nA,nB->nAB',longitudinal_coeff.conj(),longitudinal_coeff)/rank
        for m in group:
            W[...,m,:,:][mask]=cluster_weight
            WL[...,m,:,:][mask]=cluster_long
    return W,WL

@dataclass
class OnShellMatrices:
    omega:float
    surface_grid:int
    Im_total:np.ndarray
    Im_longitudinal:np.ndarray
    Im_transverse:np.ndarray
    Im_branches:list[np.ndarray]
    Im_branch_longitudinal:list[np.ndarray]
    triangles:list[int]
    min_density_eigenvalue:float
    symmetry_residual_before_projection:float


def integrate_on_shell(omega:float,p:ModelParameters,*,nseg:int=6,grid:int=64)->OnShellMatrices:
    if omega<=0 or grid<24:raise ValueError
    k=-math.pi+2*math.pi*np.arange(grid)/grid;dk=2*math.pi/grid
    kx=k[:,None,None];ky=k[None,:,None];kz=k[None,None,:]
    vals,vec=np.linalg.eigh(bulk_stiffness_batch(kx,ky,kz,p)/p.m_bulk)
    F=cross_section_form_factor(ky,kz)
    coeff=np.einsum('...am,...aA->...mA',vec.conj(),F,optimize=True)
    KX,KY,KZ=np.meshgrid(k,k,k,indexing='ij');kn=np.sqrt(KX*KX+KY*KY+KZ*KZ)
    kh=np.zeros((grid,grid,grid,3));mask=kn>1e-14
    kh[...,0][mask]=KX[mask]/kn[mask];kh[...,1][mask]=KY[mask]/kn[mask];kh[...,2][mask]=KZ[mask]/kn[mask]
    weights,long_weights=_degenerate_spectral_weights(vals,vec,coeff,kh)
    branch_pos=[];branch_long=[];tri=[]
    phase_dx=[np.exp(1j*kx*d) for d in range(nseg)]
    for m in range(3):
        f=vals[...,m]-omega**2
        grad=[(np.roll(f,-1,axis=ax)-np.roll(f,1,axis=ax))/(2*dk) for ax in range(3)]
        try: verts,faces,_,_=marching_cubes(_extend(f),level=0,spacing=(1,1,1),allow_degenerate=False)
        except ValueError:
            branch_pos.append(np.zeros((nseg,2,2),complex));branch_long.append(np.zeros((nseg,2,2),complex));tri.append(0);continue
        cent=(verts[faces[:,0]]+verts[faces[:,1]]+verts[faces[:,2]])/3
        v0=verts[faces[:,1]]-verts[faces[:,0]];v1=verts[faces[:,2]]-verts[faces[:,0]]
        area=.5*np.linalg.norm(np.cross(v0,v1),axis=1)*dk*dk
        ci=cent.T
        gc=np.sqrt(sum(map_coordinates(_extend(g),ci,order=1,mode='wrap')**2 for g in grad))
        valid=gc>1e-11;fac=area[valid]/gc[valid]/((2*math.pi)**3*p.m_bulk)
        W=weights[...,m,:,:];WL=long_weights[...,m,:,:]
        pos=np.zeros((nseg,2,2),complex);pl=np.zeros_like(pos)
        for d in range(nseg):
            WD=W*phase_dx[d][...,None,None]
            for a in range(2):
                for b in range(2):
                    q=_interp_complex(WD[...,a,b],ci)[valid]
                    ql=_interp_complex((WL*phase_dx[d][...,None,None])[...,a,b],ci)[valid]
                    pos[d,a,b]=np.sum(fac*q);pl[d,a,b]=np.sum(fac*ql)
        branch_pos.append(pos);branch_long.append(pl);tri.append(int(len(faces)))
    rho_br=[_assemble(x,nseg) for x in branch_pos]
    rho_long=_assemble(sum(branch_long),nseg);rho_total=sum(rho_br)
    preA=_assemble_raw(sum(branch_pos),nseg)
    resid=float(max(np.linalg.norm(np.imag(preA))/max(np.linalg.norm(np.real(preA)),1e-30),
                    np.linalg.norm(preA-preA.T.conj())/max(np.linalg.norm(preA),1e-30)))
    Im_total=math.pi*rho_total;Im_long=math.pi*rho_long;Im_T=Im_total-Im_long
    eig=float(np.min(np.linalg.eigvalsh(.5*(rho_total+rho_total.T))))
    rho_bl=[_assemble(x,nseg) for x in branch_long]
    return OnShellMatrices(float(omega),int(grid),Im_total,Im_long,Im_T,[math.pi*x for x in rho_br],[math.pi*x for x in rho_bl],tri,eig,resid)
