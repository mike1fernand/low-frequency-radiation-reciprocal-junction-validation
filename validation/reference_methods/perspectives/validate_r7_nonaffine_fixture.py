#!/usr/bin/env python3
"""Small actual-R7 geometry regression; not production scattering validation."""
from __future__ import annotations
import hashlib
import json
import argparse
from pathlib import Path
import platform
import sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
EXPECTED_SOURCE_HASHES = {
    'paper1_main.tex': '0a35330186d45d2345a943f6510c2ce5e9dc27320a95c60c4b974f6c9cee4b0d',
    'paper1_supplement.tex': 'fccc86e9a20159310b51112345dfee31a23f963b65c29f4327de4eb5f91a93ea',
}
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-dir',type=Path,
                    help='Optional directory of the exact reviewed R2 TeX files; verify their recorded hashes.')
args=parser.parse_args()
if args.source_dir is not None:
    actual_source_hashes={name:hashlib.sha256((args.source_dir/name).read_bytes()).hexdigest()
                          for name in EXPECTED_SOURCE_HASHES}
    if actual_source_hashes != EXPECTED_SOURCE_HASHES:
        parser.error('Source hashes differ from the reviewed R2 provenance; no new source audit is implied.')
    source_status='Reviewed R2 source hashes verified on this run.'
else:
    actual_source_hashes=None
    source_status='Source files not supplied; recorded reviewed-R2 provenance retained, not rechecked on this run.'
NQ = 4
bulk = [(i,j,k) for i in range(3) for j in (-1,0) for k in (-1,0)]
bulk_index = {r:n for n,r in enumerate(bulk)}
ND = NQ + 3*len(bulk)
mass = np.ones(ND)
Kq, Kaxis, Kdiag = 1.69, 1., 1.
g = np.array([.1,.2,.1])
BETA, TAUMAX = .15, 12.

def f(tau):
    return tau + .5*BETA*tau*tau

def inverse(t):
    # Stable form of (sqrt(1+2 beta t)-1)/beta.
    return 2*t/(np.sqrt(1+2*BETA*t)+1)

def clock_rate(tau):
    return 1+BETA*tau

K_parts = {kind: np.zeros((ND, ND)) for kind in ('line','bulk','coupling')}
bonds = []
coupling_links = []

def dofs(r):
    return np.arange(NQ+3*bulk_index[r], NQ+3*bulk_index[r]+3)

def add_bond(kind, stiffness, extension_vector, metadata):
    K_parts[kind] += stiffness*np.outer(extension_vector,extension_vector)
    bonds.append((kind,stiffness,extension_vector.copy(),metadata))

for i in range(NQ-1):
    b=np.zeros(ND); b[i]=-1; b[i+1]=1
    add_bond('line',Kq,b,{'endpoints':[i,i+1]})

positive_directions = [(1,0,0),(0,1,0),(0,0,1),
                       (1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
for r in bulk:
    for d in positive_directions:
        s=tuple(r[a]+d[a] for a in range(3))
        if s not in bulk_index:
            continue
        n=np.asarray(d,dtype=float); n/=np.linalg.norm(n)
        b=np.zeros(ND); b[dofs(r)]=-n; b[dofs(s)]=n
        stiffness=Kaxis if sum(abs(v) for v in d)==1 else Kdiag
        add_bond('bulk',stiffness,b,{'endpoints':[r,s],'direction':list(d)})

for i in range(3):
    for sy in (-1,1):
        for sz in (-1,1):
            r=(i,(sy-1)//2,(sz-1)//2)
            for line_site, sign in ((i,1),(i+1,-1)):
                n=np.array([sign,sy,sz],dtype=float)/np.sqrt(3)
                b=np.zeros(ND); b[line_site]=-n[0]; b[dofs(r)]=n
                add_bond('coupling',g[i],b,{'segment':i,'line_site':line_site,'bulk_site':r})
                coupling_links.append((line_site,dofs(r),n,g[i],b))

K = sum(K_parts.values())

def energy_components(x,p):
    eq = .5*np.dot(p[:NQ],p[:NQ]/mass[:NQ]) + .5*x@K_parts['line']@x
    eb = .5*np.dot(p[NQ:],p[NQ:]/mass[NQ:]) + .5*x@K_parts['bulk']@x
    ec = .5*x@K_parts['coupling']@x
    return np.array([eq,eb,ec])

def current(x,p):
    # Direct physical paired-spring force and velocity formula, rather than dE_b/dtau.
    velocity=p/mass
    ans=0.
    for iq,ib,n,stiffness,b in coupling_links:
        extension=float(b@x)
        endpoint_sum=velocity[ib].copy()
        endpoint_sum[0]+=velocity[iq]
        ans += -.5*stiffness*extension*float(n@endpoint_sum)
    return ans

def rhs_tau(tau,y):
    x,p=y[:ND],y[ND:2*ND]
    return np.r_[p/mass,-K@x,current(x,p)]

def rhs_time(t,y):
    F=clock_rate(inverse(t))
    return rhs_tau(inverse(t),y)/F

def rhs_lagrangian(t,y,omit_rate_term=False):
    x,velocity_t=y[:ND],y[ND:2*ND]
    F=clock_rate(inverse(t))
    acc=-(K@x)/mass/F**2
    if not omit_rate_term:
        acc -= BETA*velocity_t/F**2
    # p=∂L_t/∂x_t=M F x_t, L_t=F*T(x_t)-V/F.
    p=mass*F*velocity_t
    return np.r_[velocity_t,acc,current(x,p)/F]

x0=np.zeros(ND); p0=np.zeros(ND)
x0[:NQ]=[.01,.02,-.01,.005]
p0[:NQ]=[.01,-.015,.007,0.]
y0=np.r_[x0,p0,0.]
Eparts0=energy_components(x0,p0)
H0=Eparts0.sum()
tau=np.linspace(0,TAUMAX,1201)
time=f(tau)
F=clock_rate(tau)

# An independent spectral evolution of the constant quadratic Hamiltonian.
sqrtm=np.sqrt(mass)
eig,vec=eigh(K/sqrtm[:,None]/sqrtm[None,:])
scale=np.linalg.norm(K,2)
assert eig.min()>-2e-13*scale
eig_for_evolution=np.where(abs(eig)<1e-12*scale,0.,eig)
freq=np.sqrt(np.maximum(eig_for_evolution,0))
q0=vec.T@(sqrtm*x0)
v0=vec.T@(p0/sqrtm)
phase=freq[:,None]*tau
cos=np.cos(phase)
sin=np.sin(phase)
sinoverw=tau[None,:]*np.sinc(phase/np.pi)
modeq=q0[:,None]*cos+v0[:,None]*sinoverw
modev=-freq[:,None]*q0[:,None]*sin+v0[:,None]*cos
x_exact=(vec@modeq)/sqrtm[:,None]
p_exact=sqrtm[:,None]*(vec@modev)
state_exact=np.vstack([x_exact,p_exact])
state_scale=np.max(np.abs(state_exact))
E_exact=np.array([energy_components(x_exact[:,j],p_exact[:,j]) for j in range(len(tau))])
I_exact=E_exact[:,1]+.5*E_exact[:,2]-(Eparts0[1]+.5*Eparts0[2])

rows=[]
final_solutions={}
for tol in (1e-7,1e-9,1e-11):
    sol_tau=solve_ivp(rhs_tau,(0,TAUMAX),y0,method='DOP853',t_eval=tau,rtol=tol,atol=tol*1e-2)
    sol_time=solve_ivp(rhs_time,(0,f(TAUMAX)),y0,method='DOP853',t_eval=time,rtol=tol,atol=tol*1e-2)
    sol_lag=solve_ivp(rhs_lagrangian,(0,f(TAUMAX)),y0,method='DOP853',t_eval=time,rtol=tol,atol=tol*1e-2)
    assert sol_tau.success and sol_time.success and sol_lag.success
    lag_canonical=sol_lag.y.copy()
    lag_canonical[ND:2*ND]*=mass[:,None]*F[None,:]
    energies_tau=np.array([energy_components(sol_tau.y[:ND,j],sol_tau.y[ND:2*ND,j]) for j in range(len(tau))])
    energies_t=np.array([energy_components(sol_time.y[:ND,j],sol_time.y[ND:2*ND,j]) for j in range(len(tau))])
    energy_balance_t=energies_t[:,1]+.5*energies_t[:,2]-(Eparts0[1]+.5*Eparts0[2])
    row={
        'rtol':tol,'atol':tol*1e-2,
        'tau_nfev':sol_tau.nfev,'time_nfev':sol_time.nfev,'lagrangian_nfev':sol_lag.nfev,
        'tau_state_vs_modal_max_scaled':float(np.max(abs(sol_tau.y[:2*ND]-state_exact))/state_scale),
        'time_state_vs_modal_max_scaled':float(np.max(abs(sol_time.y[:2*ND]-state_exact))/state_scale),
        'same_event_tau_vs_time_state_max_scaled':float(np.max(abs(sol_tau.y[:2*ND]-sol_time.y[:2*ND]))/state_scale),
        'lagrangian_state_vs_modal_max_scaled':float(np.max(abs(lag_canonical[:2*ND]-state_exact))/state_scale),
        'tau_energy_drift_over_H0':float(np.max(abs(energies_tau.sum(axis=1)-H0))/H0),
        'time_energy_drift_over_H0':float(np.max(abs(energies_t.sum(axis=1)-H0))/H0),
        'time_integrated_current_vs_modal_transfer_over_H0':float(np.max(abs(sol_time.y[-1]-I_exact))/H0),
        'time_current_energy_balance_over_H0':float(np.max(abs(sol_time.y[-1]-energy_balance_t))/H0),
        'same_event_tau_vs_time_integrated_current_over_H0':float(np.max(abs(sol_tau.y[-1]-sol_time.y[-1]))/H0),
    }
    rows.append(row)
    final_solutions={'tau':sol_tau,'time':sol_time,'lagrangian':sol_lag}

# Independent reduced corner expression tests the construction, not only the ODE mapping.
rng=np.random.default_rng(20260923)
testx=rng.normal(size=ND)
reduced=0.
for i in range(3):
    qbar=(testx[i]+testx[i+1])/2
    delta=testx[i+1]-testx[i]
    for sy in (-1,1):
        for sz in (-1,1):
            u=testx[dofs((i,(sy-1)//2,(sz-1)//2))]
            reduced+=g[i]/3*((u[0]-qbar)**2+(sy*u[1]+sz*u[2]+delta/2)**2)
corner_energy_error=abs(reduced-.5*testx@K_parts['coupling']@testx)/max(abs(reduced),1.)

# Current independently integrated by both solvers above; sampled quadrature convergence below.
J_exact=np.array([current(x_exact[:,j],p_exact[:,j]) for j in range(len(tau))])
quadrature=[]
for stride in (16,8,4,2,1):
    # 1200 is not divisible by 16; retain endpoints explicitly.
    idx=np.unique(np.r_[np.arange(0,len(tau),stride),len(tau)-1])
    ts,taus,Fs,Js=time[idx],tau[idx],F[idx],J_exact[idx]
    Jt=Js/Fs
    old=np.trapezoid(Js,taus)
    transported=np.sum(np.diff(taus)*.5*(Fs[:-1]*Jt[:-1]+Fs[1:]*Jt[1:]))
    fresh=np.trapezoid(Jt,ts)
    quadrature.append({
        'nominal_delta_tau':.01*stride,'panels':len(idx)-1,
        'structural_trapezoid_error_over_H0':float(abs(old-I_exact[-1])/H0),
        'fresh_time_trapezoid_error_over_H0':float(abs(fresh-I_exact[-1])/H0),
        'transported_archive_identity_residual_over_H0':float(abs(old-transported)/H0),
        'fresh_vs_transported_trapezoid_over_H0':float(abs(fresh-transported)/H0),
    })

velocity_t=p_exact/mass[:,None]/F[None,:]
wrong_H=[]
for j in range(len(tau)):
    wrong_H.append(.5*np.dot(mass,velocity_t[:,j]**2)+.5*x_exact[:,j]@K@x_exact[:,j])
wrong_H=np.asarray(wrong_H)
wrong_lag=solve_ivp(lambda t,y:rhs_lagrangian(t,y,True),(0,f(TAUMAX)),y0,
                   method='DOP853',t_eval=time,rtol=1e-11,atol=1e-13)
assert wrong_lag.success
wrong_state=wrong_lag.y[:2*ND].copy(); wrong_state[ND:]*=mass[:,None]*F[None,:]

def modal_state_at(tau_value):
    ph=freq*tau_value
    a=q0*np.cos(ph)+v0*tau_value*np.sinc(ph/np.pi)
    adot=-freq*q0*np.sin(ph)+v0*np.cos(ph)
    return (vec@a)/sqrtm, sqrtm*(vec@adot)

def wrong_generator_work_integrand(tau_value):
    xx,pp=modal_state_at(tau_value)
    # A generator force Kx/F used in retained-H work introduces a spurious second 1/F.
    return current(xx,pp)/clock_rate(tau_value)

wrong_generator_transfer,wrong_generator_quad_error=quad(wrong_generator_work_integrand,0,TAUMAX,
                                                         epsabs=1e-14,epsrel=1e-11,limit=200)
negative={
    'naive_m_times_velocity_t_squared_energy_max_relative_error':float(np.max(abs(wrong_H-H0))/H0),
    'omit_nonaffine_rate_term_state_max_scaled_error':float(np.max(abs(wrong_state-state_exact))/state_scale),
    'omit_half_coupling_energy_current_balance_max_over_H0':float(np.max(abs(I_exact-(E_exact[:,1]-Eparts0[1])))/H0),
    'generator_force_used_in_retained_energy_work_endpoint_error_over_H0':float(abs(wrong_generator_transfer-I_exact[-1])/H0),
}

# Separate constant-clock finite resolvent identity; η is mapped with the real frequency.
F0=1.7
z_time=.8+.03j
z_tau=F0*z_time
M=np.diag(mass)
G_tau=np.linalg.inv(K-z_tau**2*M)
G_time=np.linalg.inv(K/F0-z_time**2*(F0*M))
force=np.zeros(ND); force[0]=1.
u_tau=G_tau@force
u_time=G_time@(force/F0)
resolvent={
    'F0':F0,'z_time':[z_time.real,z_time.imag],'z_tau':[z_tau.real,z_tau.imag],
    'relative_G_time_minus_F0_G_tau':float(np.linalg.norm(G_time-F0*G_tau)/np.linalg.norm(F0*G_tau)),
    'relative_response_error_with_force_divided_by_F0':float(np.linalg.norm(u_time-u_tau)/np.linalg.norm(u_tau)),
    'relative_response_error_if_generator_force_not_scaled':float(np.linalg.norm(G_time@force-u_tau)/np.linalg.norm(u_tau)),
    'scope':'Finite open fixture complex-frequency algebra, not the infinite retarded scattering calculation.'}

assert rows[-1]['same_event_tau_vs_time_state_max_scaled']<2e-8
assert rows[-1]['lagrangian_state_vs_modal_max_scaled']<2e-8
assert rows[-1]['time_current_energy_balance_over_H0']<2e-8
assert rows[-1]['time_integrated_current_vs_modal_transfer_over_H0']<2e-8
assert rows[-1]['time_energy_drift_over_H0']<2e-8
assert rows[-1]['same_event_tau_vs_time_state_max_scaled'] < rows[0]['same_event_tau_vs_time_state_max_scaled']/100
assert rows[-1]['time_current_energy_balance_over_H0'] < rows[0]['time_current_energy_balance_over_H0']/100
assert quadrature[-1]['structural_trapezoid_error_over_H0'] < quadrature[0]['structural_trapezoid_error_over_H0']/100
assert quadrature[-1]['fresh_time_trapezoid_error_over_H0'] < quadrature[0]['fresh_time_trapezoid_error_over_H0']/100
assert max(q['transported_archive_identity_residual_over_H0'] for q in quadrature)<1e-14
assert all(negative[k]>1e-3 for k in negative)
assert corner_energy_error<2e-14
assert resolvent['relative_G_time_minus_F0_G_tau']<2e-12
assert resolvent['relative_response_error_with_force_divided_by_F0']<2e-12
assert resolvent['relative_response_error_if_generator_force_not_scaled']>.5

result={
    'status':'PASS: finite-fixture coordinate regression only; no production transmission claim',
    'source_files':EXPECTED_SOURCE_HASHES,
    'source_verification':{'status':source_status,'actual_hashes_if_checked':actual_source_hashes},
    'source_equation_anchors':['eq:s2-Dax','eq:s2-Dfd','eq:s2-Hq','eq:s2-Hb',
                               'eq:s2-source-directions','eq:s2-Hc-corner','eq:s2-Hc-total',
                               'eq:s2-Hc-reduced','eq:bondcurrent'],
    'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
    'geometry':{'bulk_shape':[3,2,2],'bulk_nodes':len(bulk),'scalar_line_sites':NQ,'degrees_of_freedom':ND,
                'boundaries':'open/free; all external bonds omitted','coupling_segments':[0,1,2],
                'g':g.tolist(),'Kq':Kq,'Kaxis':Kaxis,'Kdiag':Kdiag,'masses':'all 1','a':1,
                'bond_counts':{kind:sum(b[0]==kind for b in bonds) for kind in K_parts}},
    'initial_state':{'line_q':x0[:NQ].tolist(),'line_p':p0[:NQ].tolist(),'bulk_q_and_p':'all zero',
                     'E_q_E_b_E_c':Eparts0.tolist(),'H0':H0},
    'clock':{'F':'1 + beta*tau','t':'tau + beta*tau^2/2','beta':BETA,'tau_interval':[0,TAUMAX],
             'time_interval':[0,float(f(TAUMAX))],'sample_events':len(tau)},
    'spectral_checks':{'min_eigenvalue':float(eig.min()),'max_eigenvalue':float(eig.max()),
                       'zero_modes_under_declared_tolerance':int(np.sum(eig_for_evolution==0)),
                       'zero_mode_tolerance':float(1e-12*scale),'stiffness_symmetry_max':float(np.max(abs(K-K.T))),
                       'reduced_vs_paired_corner_energy_relative_error':float(corner_energy_error)},
    'normalizations':{'state_scale_max_abs_modal_state':float(state_scale),'energy_scale':'H0'},
    'convergence':rows,'sampled_current_quadrature':quadrature,'negative_controls':negative,
    'wrong_generator_force_work_quadrature_error_estimate':float(wrong_generator_quad_error),
    'finite_affine_resolvent':resolvent,
    'limitations':['Prescribed state-independent common F only; not a physical clock model.',
                   'Finite fixture with reflected waves and three source segments; not the production six-segment infinite-bulk scattering geometry.',
                   'No claim about the production 12.6 percent conversion or 10.3 percent robustness bound.',
                   'Both coordinate integrations retain the same structural energy, masses, stiffnesses and node labels.'],
}
(HERE/'R7_NONAFFINE_FIXTURE_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
np.savez_compressed(HERE/'R7_NONAFFINE_FIXTURE_TRACES.npz',tau=tau,time=time,F=F,
                    exact_x=x_exact,exact_p=p_exact,exact_energy_components=E_exact,
                    exact_interface_transfer=I_exact,exact_current_tau=J_exact,
                    tau_integrated_state=final_solutions['tau'].y,time_integrated_state=final_solutions['time'].y,
                    lagrangian_integrated_state=final_solutions['lagrangian'].y)
print(json.dumps({'status':result['status'],'geometry':result['geometry'],'spectral_checks':result['spectral_checks'],
                  'finest':rows[-1],'negative_controls':negative,'quadrature':quadrature},indent=2))
