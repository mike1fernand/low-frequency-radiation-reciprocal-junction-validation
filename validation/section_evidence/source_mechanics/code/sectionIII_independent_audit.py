#!/usr/bin/env python3
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import sympy as sp

SEED=20260718
rng=np.random.default_rng(SEED)
OUT=Path(__file__).resolve().parents[1]/'results'/'sectionIII_independent_checks_v1.json'

# ---------- symbolic identities ----------
g,a=sp.symbols('g a', positive=True, real=True)
q0,q1=sp.symbols('q0 q1', real=True)
ux,uy,uz=sp.symbols('ux uy uz', real=True)
sy,sz=sp.symbols('sy sz', integer=True)
qbar=(q0+q1)/2
dq=q1-q0
S=sy*uy+sz*uz
A=ux-qbar
B=S+dq/2
eL=(ux-q0+sy*uy+sz*uz)/sp.sqrt(3)
eR=(-ux+q1+sy*uy+sz*uz)/sp.sqrt(3)
Hpair=sp.expand(g/2*(eL**2+eR**2))
Hred=sp.expand(g/3*(A**2+B**2))
F_u=sp.Matrix([-sp.diff(Hred,v) for v in (ux,uy,uz)])
Q0=-sp.diff(Hred,q0)
Q1=-sp.diff(Hred,q1)

uxd,Sd,qbd,dqd=sp.symbols('uxd Sd qbd dqd', real=True)
work_X=sp.simplify(2*g*A*(uxd-qbd)/3 + (-2*g*A*uxd/3 + 2*g*A*qbd/3))
work_R=sp.simplify(2*g*B*(Sd+dqd/2)/3 + (-2*g*B*Sd/3 - g*B*dqd/3))

symbolic={
 'paired_energy_residual': str(sp.simplify(Hpair-Hred)),
 'eL_residual': str(sp.simplify(eL-(A+B)/sp.sqrt(3))),
 'eR_residual': str(sp.simplify(eR-(B-A)/sp.sqrt(3))),
 'HX_spring_residual': str(sp.simplify(g*A**2/3-g*(eL-eR)**2/4)),
 'HR_spring_residual': str(sp.simplify(g*B**2/3-g*(eL+eR)**2/4)),
 'corner_force_residual': [str(x) for x in sp.simplify(F_u+2*g*sp.Matrix([A,B*sy,B*sz])/3)],
 'left_endpoint_generalized_force_residual': str(sp.simplify(Q0-g*(A+B)/3)),
 'right_endpoint_generalized_force_residual': str(sp.simplify(Q1-g*(A-B)/3)),
 'work_X_residual': str(work_X),
 'work_R_residual': str(work_R),
}

corners=[(-1,-1),(-1,1),(1,-1),(1,1)]

def H_num(x,g,sy,sz):
    ux,uy,uz,q0,q1=x
    qb=(q0+q1)/2; d=q1-q0
    aa=ux-qb; bb=sy*uy+sz*uz+d/2
    return g*(aa*aa+bb*bb)/3

def F_num(x,g,sy,sz):
    ux,uy,uz,q0,q1=x
    qb=(q0+q1)/2; d=q1-q0
    aa=ux-qb; bb=sy*uy+sz*uz+d/2
    return np.array([-2*g*aa/3,-2*g*bb*sy/3,-2*g*bb*sz/3,g*(aa+bb)/3,g*(aa-bb)/3])

def fd5(x,h,g,sy,sz):
    grad=np.zeros(5)
    for j in range(5):
        xmm=x.copy(); xm=x.copy(); xp=x.copy(); xpp=x.copy()
        xmm[j]-=2*h; xm[j]-=h; xp[j]+=h; xpp[j]+=2*h
        grad[j]=(H_num(xmm,g,sy,sz)-8*H_num(xm,g,sy,sz)+8*H_num(xp,g,sy,sz)-H_num(xpp,g,sy,sz))/(12*h)
    return -grad

max_force_err=0.0
for _ in range(5000):
    x=rng.normal(scale=rng.uniform(1e-3,2.0),size=5)
    gv=10**rng.uniform(-2,1)
    syv,szv=rng.choice([-1,1],size=2)
    exact=F_num(x,gv,syv,szv)
    L=max(1.0,np.linalg.norm(x[:3]),abs(x[3]),abs(x[4]))
    approx=fd5(x,1e-2*L,gv,syv,szv)
    scale=max(gv*L,np.linalg.norm(exact))
    max_force_err=max(max_force_err,float(np.linalg.norm(approx-exact)/scale))

# ---------- force-source decomposition and balance ----------
def full_corner_force(q0,q1,u,sy,sz,g=1.0):
    qb=(q0+q1)/2; d=q1-q0
    aa=u[0]-qb; ss=sy*u[1]+sz*u[2]; bb=ss+d/2
    return -(2*g/3)*np.array([aa,bb*sy,bb*sz],float)

def source_R(q0,q1,sy,sz,g=1.0):
    return -(g/3)*(q1-q0)*np.array([0,sy,sz],float)

def source_X(q0,q1,g=1.0):
    return (2*g/3)*((q0+q1)/2)*np.array([1,0,0],float)

def total_R(q0,q1,u,sy,sz,g=1.0):
    ss=sy*u[1]+sz*u[2]
    return -(2*g/3)*(ss+(q1-q0)/2)*np.array([0,sy,sz],float)

def endpoint_Q(q0,q1,u,sy,sz,g=1.0):
    qb=(q0+q1)/2; d=q1-q0
    aa=u[0]-qb; bb=sy*u[1]+sz*u[2]+d/2
    return g*(aa+bb)/3,g*(aa-bb)/3

def sums(q0,q1,us,g=1.0,a=1.0):
    F=np.zeros(3); FR=np.zeros(3); Fsrc=np.zeros(3)
    tau=np.zeros(3); tauR=np.zeros(3); tauSrc=np.zeros(3)
    Msrc=np.zeros((3,3)); MR=np.zeros((3,3)); qsum=0.0
    for u,(sy,sz) in zip(us,corners):
        r=np.array([0,a*sy/2,a*sz/2],float)
        f=full_corner_force(q0,q1,u,sy,sz,g)
        fr=total_R(q0,q1,u,sy,sz,g)
        fs=source_R(q0,q1,sy,sz,g)
        ql,qr=endpoint_Q(q0,q1,u,sy,sz,g)
        F+=f; FR+=fr; Fsrc+=fs; qsum+=ql+qr
        tau+=np.cross(r,f); tauR+=np.cross(r,fr); tauSrc+=np.cross(r,fs)
        Msrc+=np.outer(fs,r); MR+=np.outer(fr,r)
    return {'F_bulk':F,'F_R_total':FR,'F_R_source':Fsrc,'tau_bulk':tau,'tau_R_total':tauR,'tau_R_source':tauSrc,'M_R_source':Msrc,'M_R_total':MR,'endpoint_axial_sum':qsum}

# exact simple counterexamples
us_R=[np.zeros(3) for _ in corners]
us_R[corners.index((1,1))]=np.array([0,1,0],float)
ce_R=sums(0,0,us_R)
us_X=[np.zeros(3) for _ in corners]
us_X[corners.index((1,1))]=np.array([1,0,0],float)
ce_X=sums(0,0,us_X)

# centered A1 state
q0v,q1v=0.2,-0.1; UX,UR=0.12,-0.08; gv=0.7
us_A1=[np.array([UX,UR*sy,UR*sz],float) for sy,sz in corners]
a1=sums(q0v,q1v,us_A1,gv)
a1_total_with_line=a1['F_bulk']+np.array([a1['endpoint_axial_sum'],0,0])

# random proof that total R torque is always zero numerically
max_total_R_torque=0.0
for _ in range(1000):
    us=[rng.normal(size=3) for _ in corners]
    ss=sums(*rng.normal(size=2),us,g=float(10**rng.uniform(-1,1)),a=float(rng.uniform(.2,2)))
    max_total_R_torque=max(max_total_R_torque,float(np.linalg.norm(ss['tau_R_total'])))

# ---------- endpoint embedding ----------
def endpoint_vector_forces(qLvec,qRvec,u,sy,sz,g=1.0):
    nL=np.array([1,sy,sz],float)/math.sqrt(3)
    nR=np.array([-1,sy,sz],float)/math.sqrt(3)
    eL=float(nL@(u-qLvec)); eR=float(nR@(u-qRvec))
    return g*eL*nL,g*eR*nR

def endpoint_transverse_norm(us,q0,q1,g=1.0):
    FL=np.zeros(3); FR=np.zeros(3)
    qL=np.array([q0,0,0],float); qR=np.array([q1,0,0],float)
    for u,(sy,sz) in zip(us,corners):
        fl,fr=endpoint_vector_forces(qL,qR,u,sy,sz,g)
        FL+=fl; FR+=fr
    return float(np.linalg.norm(FL[1:])),float(np.linalg.norm(FR[1:])),FL,FR

a1_emb=endpoint_transverse_norm(us_A1,q0v,q1v,gv)
off_us=[rng.normal(size=3) for _ in corners]
off_emb=endpoint_transverse_norm(off_us,q0v,q1v,gv)

# ---------- Fourier transforms ----------
def local_direct(g,dq,ky,kz,a=1.0):
    out=np.zeros(3,dtype=complex)
    for sy,sz in corners:
        f=-(g/3)*dq*np.array([0,sy,sz],float)
        out+=f*np.exp(-1j*(ky*sy+kz*sz)*a/2)
    return out

def local_formula(g,dq,ky,kz,a=1.0):
    return (4j*g*dq/3)*np.array([0,np.sin(ky*a/2)*np.cos(kz*a/2),np.cos(ky*a/2)*np.sin(kz*a/2)],complex)

max_local_fourier=0.0
for _ in range(1000):
    gv=10**rng.uniform(-1,1); dqv=rng.normal()+1j*rng.normal(); ky,kz=rng.uniform(-2,2,size=2); av=rng.uniform(.2,2)
    d=local_direct(gv,dqv,ky,kz,av); f=local_formula(gv,dqv,ky,kz,av)
    max_local_fourier=max(max_local_fourier,float(np.linalg.norm(d-f)/max(1,np.linalg.norm(d),np.linalg.norm(f))))

def slab_direct(Q,kappa,kx,ky,kz,gs,a=1.0,i0=0):
    out=np.zeros(3,dtype=complex)
    for n,gv in enumerate(gs):
        i=i0+n; qi=Q*np.exp(1j*kappa*i*a); qip=Q*np.exp(1j*kappa*(i+1)*a); dq=qip-qi
        for sy,sz in corners:
            f=-(gv/3)*dq*np.array([0,sy,sz],float)
            out+=f*np.exp(-1j*(kx*(i+.5)*a+ky*sy*a/2+kz*sz*a/2))
    return out

def slab_formula(Q,kappa,kx,ky,kz,gs,a=1.0,i0=0):
    ap=sum(gv*np.exp(1j*(kappa-kx)*(i0+n+.5)*a) for n,gv in enumerate(gs))
    v=np.array([0,np.sin(ky*a/2)*np.cos(kz*a/2),np.cos(ky*a/2)*np.sin(kz*a/2)],complex)
    return -(8*Q/3)*np.sin(kappa*a/2)*ap*v

max_slab_fourier=0.0
for _ in range(500):
    Q=rng.normal()+1j*rng.normal(); kappa,kx,ky,kz=rng.uniform(-2,2,size=4); gs=rng.uniform(0,1,size=int(rng.integers(1,10))); av=rng.uniform(.2,1.5); i0=int(rng.integers(-5,5))
    d=slab_direct(Q,kappa,kx,ky,kz,gs,av,i0); f=slab_formula(Q,kappa,kx,ky,kz,gs,av,i0)
    max_slab_fourier=max(max_slab_fourier,float(np.linalg.norm(d-f)/max(1,np.linalg.norm(d),np.linalg.norm(f))))

# ---------- long-wave slopes ----------
def phase_removed(Q,g,kappa,ky,kz,a=1.0):
    dq=2j*Q*np.exp(1j*kappa*a/2)*np.sin(kappa*a/2)
    return np.exp(-1j*kappa*a/2)*local_formula(g,dq,ky,kz,a)

def leading(Q,g,kappa,ky,kz,a=1.0):
    return -(2*g*a*a*Q/3)*np.array([0,kappa*ky,kappa*kz],complex)

eps=np.logspace(-4,-1,16); direction=np.array([.8,.5,.3]); amps=[]; rems=[]
for e in eps:
    kappa,ky,kz=e*direction
    f=phase_removed(1+.2j,.7,kappa,ky,kz); l=leading(1+.2j,.7,kappa,ky,kz)
    amps.append(np.linalg.norm(f)); rems.append(np.linalg.norm(f-l))
slope_amp=float(np.polyfit(np.log(eps),np.log(amps),1)[0])
slope_rem=float(np.polyfit(np.log(eps[3:]),np.log(rems[3:]),1)[0])

# continuum fixed-direction SV/SH slopes
def pol(kx,ky,kz):
    k=np.array([kx,ky,kz],float); kn=np.linalg.norm(k); kp=math.hypot(ky,kz)
    return np.array([kp,-kx*ky/kp,-kx*kz/kp])/kn, np.array([0,-kz/kp,ky/kp])
vals_X=[]; vals_RSV=[]; vals_RSH=[]
for e in eps:
    kx,ky,kz=e*np.array([.7,.5,.3]); kappa=.8*e; esv,esh=pol(kx,ky,kz); R=phase_removed(1,1,kappa,ky,kz).real
    vals_X.append(abs(esv[0])); vals_RSV.append(abs(R@esv)); vals_RSH.append(abs(R@esh))
slopes={'X_to_SV':float(np.polyfit(np.log(eps),np.log(vals_X),1)[0]),'R_to_SV':float(np.polyfit(np.log(eps),np.log(vals_RSV),1)[0]),'R_to_SH':float(np.polyfit(np.log(eps[3:]),np.log(vals_RSH[3:]),1)[0])}

result={
 'seed':SEED,
 'symbolic':symbolic,
 'finite_difference':{'cases':5000,'five_point_step_factor':1e-2,'max_scale_aware_error':max_force_err},
 'source_vs_total_R':{
   'simple_total_R_force_counterexample':ce_R['F_R_total'].tolist(),
   'same_counterexample_total_R_torque':ce_R['tau_R_total'].tolist(),
   'simple_full_force_torque_counterexample':ce_X['tau_bulk'].tolist(),
   'max_random_total_R_torque_norm':max_total_R_torque,
   'note':'The total R force need not balance off symmetry, but every total-R corner force is radial, so total R torque about the segment midpoint is identically zero. Full X+R bulk torque can be nonzero off symmetry.'
 },
 'A1_balance':{
   'bulk_force':a1['F_bulk'].tolist(),
   'endpoint_axial_sum':a1['endpoint_axial_sum'],
   'bulk_plus_line_force':a1_total_with_line.tolist(),
   'bulk_torque':a1['tau_bulk'].tolist(),
 },
 'radial_source_moment':{'computed':a1['M_R_source'].tolist(),'expected_coefficient':float(-2*0.7*(q1v-q0v)/3)},
 'endpoint_embedding':{
   'A1_left_transverse_norm':a1_emb[0],
   'A1_right_transverse_norm':a1_emb[1],
   'off_symmetry_left_transverse_norm':off_emb[0],
   'off_symmetry_right_transverse_norm':off_emb[1],
 },
 'fourier':{'local_max_scaled_error':max_local_fourier,'full_slab_kappa_not_kx_max_scaled_error':max_slab_fourier},
 'long_wave':{'source_amplitude_slope':slope_amp,'post_leading_remainder_slope':slope_rem,'continuum_projection_slopes':slopes},
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2))
print(OUT)
print(json.dumps(result,indent=2)[:5000])
