#!/usr/bin/env python3
"""Independent finite-core and scalar-lead scattering evaluator.

This file deliberately does not import Phase-12 reduced_core.py. It re-derives
and implements the active-core static matrices, scalar lead self-energy, Schur
complement, retarded finite-core solve, flux coefficients, and partial bulk
power formula. It is used to ensure that the generalized-Bloch cross-check is
not sharing the load-bearing finite-core implementation.
"""
from __future__ import annotations
from dataclasses import dataclass
import cmath, math
import numpy as np

@dataclass(frozen=True)
class Parameters:
    m_spine: float=1.0
    k_spine: float=1.69


def active_core_static(g_segments):
    g=np.asarray(g_segments,float)
    if g.ndim!=1 or len(g)==0 or np.any(g<0): raise ValueError('invalid g')
    n=len(g); U=np.zeros((2*n,2*n),complex); C=np.zeros((2*n,n+1),complex); Kq_self=np.zeros((n+1,n+1),float)
    for s,gs in enumerate(g):
        # Four corner nodes, normalized equal-axial X and radial R patterns.
        U[2*s,2*s]=2.0*gs/3.0
        U[2*s+1,2*s+1]=4.0*gs/3.0
        C[2*s,s]=C[2*s,s+1]=-2.0*gs/3.0
        cr=2.0*math.sqrt(2.0)*gs/3.0
        C[2*s+1,s]=-cr; C[2*s+1,s+1]=cr
        Kq_self[s,s]+=4.0*gs/3.0;Kq_self[s+1,s+1]+=4.0*gs/3.0
    return U,C,Kq_self


def retarded_scalar_lead_self_energy(omega,eta,k,m):
    z=complex(omega,eta);a=m*z*z-2.0*k
    disc=cmath.sqrt(a*a-4.0*k*k)
    roots=(0.5*(a-disc),0.5*(a+disc))
    # Retarded self-energy has non-positive imaginary part and is the root
    # continuously connected to zero as |z|->infinity.
    admissible=[x for x in roots if x.imag<=1e-13]
    return min(admissible,key=abs) if admissible else min(roots,key=lambda x:x.imag)


def evaluate(omega,g_segments,G0,Im_parts=None,p=Parameters(),eta_lead=1e-12):
    g=np.asarray(g_segments,float);n=len(g);nq=n+1
    U,C,Kself=active_core_static(g)
    # Dressed bulk response with local core stiffness U: (G0^{-1}+U)^{-1}.
    Gd=np.linalg.solve(np.eye(2*n)+G0@U,G0)
    Sigma_bulk=-(C.conj().T@Gd@C)
    Gamma_bulk=1j*(Sigma_bulk-Sigma_bulk.conj().T)
    Kline=np.zeros((nq,nq),float);np.fill_diagonal(Kline,2.0*p.k_spine)
    for i in range(nq-1):Kline[i,i+1]=Kline[i+1,i]=-p.k_spine
    z=complex(omega,eta_lead);sl=retarded_scalar_lead_self_energy(omega,eta_lead,p.k_spine,p.m_spine)
    SL=np.zeros((nq,nq),complex);SR=np.zeros_like(SL);SL[0,0]=sl;SR[-1,-1]=sl
    GL=1j*(SL-SL.conj().T);GR=1j*(SR-SR.conj().T);gamma=float(GL[0,0].real)
    A=p.m_spine*z*z*np.eye(nq)-Kline-Kself-Sigma_bulk-SL-SR
    source=np.zeros(nq,complex);source[0]=1.0
    q=np.linalg.solve(A,source)
    Rq=float(abs(1.0-1j*gamma*q[0])**2)
    Tq=float(gamma*np.real(q.conj()@GR@q))
    Tbulk=float(gamma*np.real(q.conj()@Gamma_bulk@q))
    # Bulk active-coordinate response. Sign is immaterial for quadratic power.
    t=np.linalg.solve(np.eye(2*n)+U@G0,C@q)
    partial={}
    if Im_parts:
        for name,M in Im_parts.items():partial[name]=float(2.0*gamma*np.real(t.conj()@M@t))
    return {'Rq':Rq,'Tq':Tq,'Tbulk':Tbulk,'closure':float(abs(Rq+Tq+Tbulk-1.0)),'q':q,'t':t,'gamma_in':gamma,
            'Sigma_bulk':Sigma_bulk,'Gamma_bulk':Gamma_bulk,'partials':partial}
