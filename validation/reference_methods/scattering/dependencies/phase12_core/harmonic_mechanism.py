#!/usr/bin/env python3
"""Local harmonic geometry for a spine segment and one corner bulk node."""
from __future__ import annotations
import math
import numpy as np


def coupling_energy_and_force(q0: float, q1: float, u, sy: float, sz: float, g: float):
    """Return exact quadratic coupling energy and force on the bulk node.

    The node is at transverse signs sy,sz and is joined by two unstressed
    diagonal central springs to the segment endpoints.  No independent void,
    density, or occupancy variable is present.
    """
    u=np.asarray(u,float)
    qbar=0.5*(q0+q1); dq=q1-q0
    w=sy*u[1]+sz*u[2]+0.5*dq
    energy=(g/3.0)*((u[0]-qbar)**2+w*w)
    force=np.array([-(2*g/3)*(u[0]-qbar),
                    -(2*g/3)*sy*w,
                    -(2*g/3)*sz*w])
    return float(energy),force


def initial_transverse_force(delta_q: float, sy: float, sz: float, g: float):
    return np.array([-(g/3.0)*sy*delta_q, -(g/3.0)*sz*delta_q])


def finite_difference_force(q0,q1,u,sy,sz,g,h=1e-7):
    u=np.asarray(u,float); out=np.empty(3)
    for c in range(3):
        up=u.copy();um=u.copy();up[c]+=h;um[c]-=h
        ep,_=coupling_energy_and_force(q0,q1,up,sy,sz,g)
        em,_=coupling_energy_and_force(q0,q1,um,sy,sz,g)
        out[c]=-(ep-em)/(2*h)
    return out
