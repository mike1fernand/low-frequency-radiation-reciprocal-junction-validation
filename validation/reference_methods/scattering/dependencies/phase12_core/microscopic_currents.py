#!/usr/bin/env python3
"""Microscopic harmonic bond-energy current and Bloch normalization."""
from __future__ import annotations
import math
import numpy as np
from model_api_v1_1 import AXIAL_BONDS,FACE_DIAGONAL_BONDS,ModelParameters,_unit,bulk_dynamical_matrix

def bond_current(v_i,v_j,force_on_i):
    """Instantaneous energy current from i to j for a pair spring."""
    return -0.5*float(np.dot(np.asarray(v_i)+np.asarray(v_j),np.asarray(force_on_i)))

def mode_group_velocity(k,e,omega,p:ModelParameters):
    k=np.asarray(k,float);e=np.asarray(e,complex);g=np.zeros(3,float)
    for d,ks in ([(d,p.k_axis) for d in AXIAL_BONDS]+[(d,p.k_diag) for d in FACE_DIAGONAL_BONDS]):
        d=np.asarray(d,float);n=_unit(d);proj=abs(np.vdot(n,e))**2
        g+=2*ks*math.sin(float(k@d))*d*proj/p.m_bulk
    return g/(2*omega)

def plane_wave_current(k,e,omega,p:ModelParameters):
    """Time-averaged current vector from the microscopic pair-current formula."""
    k=np.asarray(k,float);e=np.asarray(e,complex);V0=-1j*omega*e;J=np.zeros(3,float)
    for d,ks in ([(d,p.k_axis) for d in AXIAL_BONDS]+[(d,p.k_diag) for d in FACE_DIAGONAL_BONDS]):
        d=np.asarray(d,float);n=_unit(d);phase=np.exp(1j*float(k@d));Uj=e*phase;Vj=V0*phase
        Fi=ks*n*np.dot(n,Uj-e)
        # time average of -1/2 (vi+vj).Fi for real harmonic fields
        jb=-.25*np.real(np.vdot(V0+Vj,Fi))
        J+=d*jb
    return J

def plane_wave_energy_per_cell(omega,e,p:ModelParameters):
    return .5*p.m_bulk*omega**2*float(np.vdot(e,e).real)

def validate_mode_current(k,p:ModelParameters):
    vals,vec=np.linalg.eigh(bulk_dynamical_matrix(np.asarray(k,float),p))
    out=[]
    for m,w2 in enumerate(vals):
        if w2<=1e-12:continue
        w=math.sqrt(float(w2));e=vec[:,m];J=plane_wave_current(k,e,w,p);target=plane_wave_energy_per_cell(w,e,p)*mode_group_velocity(k,e,w,p)
        out.append({'mode':m,'omega':w,'J':J.tolist(),'target':target.tolist(),'relative_error':float(np.linalg.norm(J-target)/max(np.linalg.norm(target),1e-30))})
    return out
