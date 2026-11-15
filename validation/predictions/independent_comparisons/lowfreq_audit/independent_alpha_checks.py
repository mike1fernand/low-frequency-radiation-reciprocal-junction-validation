"""Audit supplied alpha solver from complete segment energy and joint block solve.

No new bulk integration: uses the supplied finite-frequency projected kernels.
Run from any directory. Outputs JSON beside this file.
"""
from pathlib import Path
import json, sys, platform
import numpy as np
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'portable_checks' / 'rev'))
from alpha_solver import coupling_mats_alpha, solve_alpha, load_g0, load_g0_ref12
from p1solver import reference_profile

def independent_hessian(g, alpha, xfac):
    # Eq. s5-segment-energy, global coordinate order q then (xi,zeta).
    nc = len(g); nd = nc + 1; h = np.zeros((nd+2*nc, nd+2*nc))
    for s, gs in enumerate(g):
        ax = np.zeros(len(h)); ar = np.zeros(len(h))
        ax[s:s+2] = -1; ax[nd+2*s] = 1
        ar[s:s+2] = [-1,1]; ar[nd+2*s+1] = np.sqrt(2)
        h += (2*gs/3)*(xfac*np.outer(ax,ax)+alpha*np.outer(ar,ar))
    return h

def direct_joint(g, w, cq, gg, imt, iml, alpha, xfac):
    nc=len(g); nd=nc+1; h=independent_hessian(g,alpha,xfac)
    ks=h[:nd,:nd]; v=h[nd:,:nd]; u=h[nd:,nd:]
    kq=cq*cq; kap=2*np.arcsin(w/(2*cq)); gam=2*kq*np.sin(kap)
    kd=2*kq*np.eye(nd)-kq*(np.eye(nd,k=1)+np.eye(nd,k=-1))
    lead=np.zeros((nd,nd),complex)
    lead[0,0]=lead[-1,-1]=-kq*np.exp(1j*kap)
    # Physical coupled equations: core displacement u_c=-G eta.
    j=np.block([[w*w*np.eye(nd)-kd-ks-lead,-v.T],
                [v,np.linalg.inv(gg)+u]])
    rhs=np.zeros(nd+2*nc,complex); rhs[0]=np.sqrt(gam)
    sol=np.linalg.solve(j,rhs); chi=sol[:nd]; uc=sol[nd:]
    eta=v@chi+u@uc
    q=lambda a,b: float(2*np.real(a.conj()@b@a))
    return dict(chi=chi,eta=eta,Tq=gam*abs(chi[-1])**2,
                R=abs(1-1j*np.sqrt(gam)*chi[0])**2,
                TqT=q(eta,imt),TqL=q(eta,iml))

g=reference_profile(6); records=[]; hessian_records=[]
for alpha,xfac in [(0,1),(.35,1),(1,1),(1,0),(0,0)]:
    u,v,ks=coupling_mats_alpha(g,alpha,xfac)
    h=independent_hessian(g,alpha,xfac); nd=7
    assembled=np.block([[ks,v.T],[v,u]])
    p=np.tile([2.,0.],6); e=np.ones(nd)
    hessian_records.append(dict(alpha=alpha,xfac=xfac,
        max_hessian_difference=float(np.max(abs(h-assembled))),
        minimum_hessian_eigenvalue=float(np.linalg.eigvalsh(h).min()),
        uniform_translation_residual=float(np.max(abs(h@np.r_[e,p])))))

for w,cq in [(1.2,1.3),(.8,.8),(.8,1.3)]:
    gg,imt,iml=load_g0_ref12() if w==1.2 else load_g0(.8)
    sols={}
    for name,alpha,xfac in [('full',1,1),('X',0,1),('R',1,0),('transparent',0,0)]:
        supplied=solve_alpha(g,w,gg,imt,iml,cq*cq,alpha=alpha,xfac=xfac)
        direct=direct_joint(g,w,cq,gg,imt,iml,alpha,xfac)
        sols[name]=supplied
        records.append(dict(omega=w,cq=cq,control=name,
            values={key:float(supplied[key]) for key in ['TqT','TqL','R','Tq','closure','XX','RR','XR']},
            max_eta_difference=float(np.max(abs(direct['eta']-supplied['eta']))),
            max_chi_difference=float(np.max(abs(direct['chi']-supplied['chi']))),
            max_power_difference=float(max(abs(direct[k]-supplied[k]) for k in ['TqT','TqL','R','Tq'])),
            same_state_sum_residual=float(supplied['XX']+supplied['RR']+supplied['XR']-supplied['TqT']),
            min_eig_ImT=float(np.linalg.eigvalsh(imt).min()),
            min_eig_ImL=float(np.linalg.eigvalsh(iml).min())))
    vals=[solve_alpha(g,w,gg,imt,iml,cq*cq,alpha=a)['TqT'] for a in np.linspace(0,1,21)]
    records.append(dict(omega=w,cq=cq,control='alpha_scan',alpha=np.linspace(0,1,21).tolist(),TqT=vals,
        min_successive_difference=float(np.min(np.diff(vals))),
        max_successive_difference=float(np.max(np.diff(vals))),
        full_minus_X=float(sols['full']['TqT']-sols['X']['TqT']),
        same_full_state_XR=float(sols['full']['XR']),
        X_self_loading_change=float(sols['full']['XX']-sols['X']['TqT'])))

kernel_variants=[]
for deg in [2,3]:
    for nth in [160,240]:
        gg,imt,iml=load_g0(.8,deg=deg,nth=nth)
        for cq in [.8,1.3]:
            full=solve_alpha(g,.8,gg,imt,iml,cq*cq)
            xx=solve_alpha(g,.8,gg,imt,iml,cq*cq,alpha=0)
            kernel_variants.append(dict(omega=.8,cq=cq,degree=deg,shell_order=nth,
                full_qT=full['TqT'],X_only_qT=xx['TqT'],
                full_minus_X=full['TqT']-xx['TqT'],same_full_state_XR=full['XR']))

Z=1.3; ct=1; cl=np.sqrt(3); rho=1
coefT=2*Z/(3*np.pi*rho*ct**3); coefL=Z/(3*np.pi*rho*cl**3)
result=dict(description=__doc__,python=platform.python_version(),numpy=np.__version__,
    hessian_checks=hessian_records,finite_frequency_checks=records,kernel_variants=kernel_variants,
    coefficient_arithmetic=dict(TqT_w2=coefT,TqL_w2=coefL,
        shear_share=coefT/(coefT+coefL),weak_qT=1/(6*np.pi)*(8/3)**2/Z,
        omega_0p4_qT_w2=.0518080040/.16,
        omega_0p4_relative_difference_from_limit=(.0518080040/.16)/coefT-1))
dest=Path(__file__).with_suffix('.json'); dest.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
