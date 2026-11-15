"""Reconstruct the paper's scalar input and evolve it three independent ways."""
import json
from math import pi
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import solve_ivp

n=120; c=1.3; k0=2*np.arcsin(6/13); j=np.arange(n); t=np.array([0.,5.,10.,20.,54.])
k=2*pi*np.fft.fftfreq(n); w=2*c*abs(np.sin(k/2)); sign=np.sign(k); sign[n//2]=0
u=.002*np.exp(-(j-18)**2/98)*np.cos(k0*(j-18)); U=np.fft.fft(u); V=-1j*sign*w*U; v=np.fft.ifft(V).real
# Explicit sum of real oscillator contributions, including the stationary Nyquist.
phase=np.exp(1j*k[:,None]*j)
direct_v=np.sum(V[:,None]*phase,axis=0).real/n
assert max(abs(direct_v-v))<1.e-16
W=(abs(V)**2+w**2*abs(U)**2)/(2*n)
eo=(v@v+c*c*sum(np.diff(u)**2))/2; edft=sum(W)
s=np.arange(60,66);g=np.array([.1,.2,.2,.2,.2,.1])
order=np.argsort(w); cov=(w>=.8)&(w<=1.6)
d=np.diff(np.eye(n),axis=0)
Hopen=c*c*d.T@d; closing=np.zeros(n);closing[-1]=1;closing[0]=-1
rows=[]
for bc,H in [('open',Hopen),('periodic',Hopen+c*c*np.outer(closing,closing))]:
    lam,Q=eigh(H); assert abs(lam[0])<1e-13;om=np.sqrt(np.maximum(0,lam));om[0]=0.
    up=Q.T@u;vp=Q.T@v
    uu=Q@(np.cos(om[:,None]*t)*up[:,None]+(t*np.sinc(om[:,None]*t/pi))*vp[:,None])
    vv=Q@(-(om*up)[:,None]*np.sin(om[:,None]*t)+vp[:,None]*np.cos(om[:,None]*t))
    if bc=='open':
        f=lambda _,y:np.r_[y[n:],c*c*np.diff(np.r_[0,np.diff(y[:n]),0])]
    else:
        f=lambda _,y:np.r_[y[n:],c*c*(np.roll(y[:n],1)-2*y[:n]+np.roll(y[:n],-1))]
    sol=solve_ivp(f,[0,t[-1]],np.r_[u,v],t_eval=t,method='DOP853',rtol=1e-12,atol=1e-16)
    assert sol.success
    ode_error=np.max(abs(sol.y-np.vstack((uu,vv))),axis=0)
    energy=(np.sum(vv**2,axis=0)+np.sum(uu*(H@uu),axis=0))/2
    assert max(ode_error)<1e-14
    answer=dict(model=bc,maximum_node_u_or_v_error=ode_error.tolist(),maximum_relative_energy_drift=float(max(abs(energy/energy[0]-1))))
    if bc=='periodic':
        uf=np.fft.ifft(U[:,None]*np.cos(w[:,None]*t)+(t*np.sinc(w[:,None]*t/pi))*V[:,None],axis=0).real
        vf=np.fft.ifft(V[:,None]*np.cos(w[:,None]*t)-(w*U)[:,None]*np.sin(w[:,None]*t),axis=0).real
        answer['max_DFT_vs_real_eigensystem_error']=float(max(np.max(abs(uf-uu)),np.max(abs(vf-vv))))
        assert answer['max_DFT_vs_real_eigensystem_error']<2e-15
    rows.append(answer)
print(json.dumps(dict(carrier=float(k0),group_speed=float(c*np.cos(k0/2)),periodic_energy=float(edft),open_energy=float(eo),relative_closing_energy=float((edft-eo)/eo),independent_closing_formula=float(c*c*(u[-1]-u[0])**2/(2*eo)),right_fraction=float(1-W[n//2]/edft),nyquist_fraction=float(W[n//2]/edft),inband_coverage=float(sum(W[cov])/edft),distinct_band_positive_frequencies=int(sum(cov[:n//2])),weighted_percentiles=w[order[np.searchsorted(W[order].cumsum()/edft,[.01,.99])]].tolist(),coupling_over_open_energy=float(2/3*sum(g*(u[s]**2+u[s+1]**2))/eo),times=t.tolist(),evolution=rows),indent=2))
