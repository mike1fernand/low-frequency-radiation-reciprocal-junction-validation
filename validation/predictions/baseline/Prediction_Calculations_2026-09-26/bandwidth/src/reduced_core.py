#!/usr/bin/env python3
"""Exact active-symmetry reduction and finite-slab scattering.

For each coupled segment the scalar spine reaches only two normalized
cross-sectional patterns of the 12 bulk coordinates: equal axial motion X and
symmetric radial breathing motion R.  The reduction 12*nseg -> 2*nseg is exact.
"""
from __future__ import annotations
from dataclasses import dataclass
import cmath, math
from typing import Iterable
import numpy as np
from scipy import linalg
from model_api_v1_1 import AXIAL_BONDS, FACE_DIAGONAL_BONDS, ModelParameters, _unit

CORE_TRANS=np.array([[-1,-1],[-1,0],[0,-1],[0,0]],int)
CORE_SIGNS=np.array([[-1.,-1.],[-1.,1.],[1.,-1.],[1.,1.]])

def active_basis(nseg:int)->np.ndarray:
    if nseg<=0: raise ValueError('nseg must be positive')
    B=np.zeros((12*nseg,2*nseg),complex); nr=2*math.sqrt(2)
    for s in range(nseg):
        for a,(sy,sz) in enumerate(CORE_SIGNS):
            r=12*s+3*a
            B[r,2*s]=.5
            B[r+1,2*s+1]=sy/nr
            B[r+2,2*s+1]=sz/nr
    return B

def reduced_static_matrices(g_segments:Iterable[float],p:ModelParameters):
    g=np.asarray(tuple(g_segments),float)
    if g.ndim!=1 or len(g)==0 or np.any(g<0): raise ValueError('invalid g_segments')
    n=len(g); U=np.zeros((2*n,2*n),complex);V=np.zeros((2*n,n+1),complex);Kself=np.zeros((n+1,n+1))
    for s,gs in enumerate(g):
        U[2*s,2*s]=2*gs/3;U[2*s+1,2*s+1]=4*gs/3
        V[2*s,s]=V[2*s,s+1]=-2*gs/3
        r=2*math.sqrt(2)*gs/3
        V[2*s+1,s]=-r;V[2*s+1,s+1]=r
        Kself[s,s]+=4*gs/3;Kself[s+1,s+1]+=4*gs/3
    return U,V,Kself

def cross_section_form_factor(ky,kz):
    ky,kz=np.broadcast_arrays(np.asarray(ky,float),np.asarray(kz,float))
    F=np.zeros(ky.shape+(3,2),complex);nr=2*math.sqrt(2)
    for (yy,zz),(sy,sz) in zip(CORE_TRANS,CORE_SIGNS):
        ph=np.exp(-1j*(ky*yy+kz*zz))
        F[...,0,0]+=.5*ph;F[...,1,1]+=sy/nr*ph;F[...,2,1]+=sz/nr*ph
    return F

def bulk_stiffness_batch(kx,ky,kz,p:ModelParameters):
    kx,ky,kz=np.broadcast_arrays(np.asarray(kx,float),np.asarray(ky,float),np.asarray(kz,float))
    K=np.zeros(kx.shape+(3,3),float)
    bonds=[(d,p.k_axis) for d in AXIAL_BONDS]+[(d,p.k_diag) for d in FACE_DIAGONAL_BONDS]
    for d,ks in bonds:
        n=_unit(d);ph=kx*d[0]+ky*d[1]+kz*d[2]
        K+=2*ks*(1-np.cos(ph))[...,None,None]*np.outer(n,n)
    return K

def reduced_bare_green_finite_eta(omega,eta,nseg,p:ModelParameters,*,grid=96,shift=(0.,0.,0.),kx_chunk=4):
    p.validate()
    if omega<=0 or eta<=0 or grid<8: raise ValueError
    sh=np.asarray(shift,float);ar=np.arange(grid)
    bx=-math.pi+2*math.pi*(ar+sh[0])/grid
    by=-math.pi+2*math.pi*(ar+sh[1])/grid
    bz=-math.pi+2*math.pi*(ar+sh[2])/grid
    ky=by[None,:,None];kz=bz[None,None,:];F=cross_section_form_factor(ky,kz)
    abar=np.zeros((grid,2,2),complex);z=complex(omega,eta);I=np.eye(3)
    for st in range(0,grid,kx_chunk):
        en=min(grid,st+kx_chunk);kx=bx[st:en,None,None]
        A=bulk_stiffness_batch(kx,ky,kz,p).astype(complex)-p.m_bulk*z*z*I
        G=np.linalg.inv(A);GF=np.einsum('...ab,...bA->...aA',G,F,optimize=True)
        red=np.einsum('...aA,...aB->...AB',F.conj(),GF,optimize=True)
        abar[st:en]=red.mean(axis=(1,2))
    blocks={}
    for dx in range(-(nseg-1),nseg):
        blocks[dx]=np.einsum('k,kAB->AB',np.exp(1j*bx*dx),abar)/grid
    out=np.empty((2*nseg,2*nseg),complex)
    for s in range(nseg):
        for t in range(nseg): out[2*s:2*s+2,2*t:2*t+2]=blocks[s-t]
    return out

def analytic_q_lead_self_energy(z:complex,p:ModelParameters)->complex:
    a0=p.m_spine*z*z-2*p.k_spine;root=cmath.sqrt(a0*a0-4*p.k_spine**2)
    c=[.5*(a0-root),.5*(a0+root)];r=[x for x in c if x.imag<=1e-14]
    return min(r,key=abs) if r else min(c,key=lambda x:x.imag)

@dataclass(frozen=True)
class ReducedScatteringResult:
    omega:float;Rq:float;Tq:float;Tbulk:float;closure_error:float;min_gamma_bulk_eigenvalue:float

def scatter_from_reduced_green(omega,g_segments,p:ModelParameters,G0,*,eta_lead=1e-12,return_internal=False):
    g=np.asarray(tuple(g_segments),float);n=len(g);nq=n+1
    U,V,Kself=reduced_static_matrices(g,p)
    Gc=np.linalg.solve(np.eye(2*n)+G0@U,G0)
    Sigma=-(V.conj().T@Gc@V);Gamma=1j*(Sigma-Sigma.conj().T)
    Kq=np.zeros((nq,nq));np.fill_diagonal(Kq,2*p.k_spine)
    for i in range(nq-1):Kq[i,i+1]=Kq[i+1,i]=-p.k_spine
    z=complex(omega,eta_lead);sl=analytic_q_lead_self_energy(z,p)
    SL=np.zeros((nq,nq),complex);SR=np.zeros_like(SL);SL[0,0]=sl;SR[-1,-1]=sl
    GL=1j*(SL-SL.conj().T);GR=1j*(SR-SR.conj().T);gin=float(GL[0,0].real)
    A=p.m_spine*z*z*np.eye(nq)-Kq-Kself-Sigma-SL-SR
    rhs=np.zeros(nq,complex);rhs[0]=1;q=linalg.solve(A,rhs,check_finite=False)
    Tq=float(gin*np.real(q.conj()@GR@q));Tb=float(gin*np.real(q.conj()@Gamma@q));Rq=float(abs(1-1j*gin*q[0])**2)
    r=ReducedScatteringResult(float(omega),Rq,Tq,Tb,float(abs(Rq+Tq+Tb-1)),float(np.min(np.linalg.eigvalsh((Gamma+Gamma.conj().T)/2)).real))
    if return_internal:
        t=np.linalg.solve(np.eye(2*n)+U@G0,V@q)
        return r,{'q':q,'t':t,'U':U,'V':V,'G0':G0,'Gc':Gc,'Sigma_bulk':Sigma,'Gamma_bulk':Gamma,'gamma_in':gin,'Kself':Kself}
    return r
