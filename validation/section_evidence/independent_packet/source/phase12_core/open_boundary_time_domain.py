#!/usr/bin/env python3
"""Independent open-boundary real-time solver for the original harmonic Hamiltonian."""
from __future__ import annotations
import math
import numpy as np
from numba import njit,prange
from model_api_v1_1 import ModelParameters,AXIAL_BONDS,FACE_DIAGONAL_BONDS,_unit
from reduced_core import bulk_stiffness_batch

BONDS=np.array(list(AXIAL_BONDS)+list(FACE_DIAGONAL_BONDS),np.int64)
NVECS=np.array([_unit(d) for d in BONDS],float)

def carrier_k(omega,p:ModelParameters):return 2*math.asin(omega/(2*math.sqrt(p.k_spine/p.m_spine)))
def exact_right_velocity(q,p:ModelParameters):
    n=len(q);k=2*np.pi*np.fft.fftfreq(n);om=2*np.sqrt(p.k_spine/p.m_spine)*np.abs(np.sin(k/2));sg=np.sign(k)
    if n%2==0:sg[n//2]=0
    return np.fft.ifft(-1j*sg*om*np.fft.fft(q)).real

def damping_profile_1d(n,width,gamma):
    out=np.zeros(n);w=int(width)
    if w<=0:return out
    for i in range(w):
        x=(w-i)/w;v=gamma*x*x
        out[i]=max(out[i],v);out[n-1-i]=max(out[n-1-i],v)
    return out

@njit(cache=True,fastmath=True)
def forces_open(q,u,gprof,Kq,ka,kd,fq,fu):
    nx,ny,nz,_=u.shape
    for i in range(nx):
        f=0.
        if i>0:f+=Kq*(q[i-1]-q[i])
        if i+1<nx:f+=Kq*(q[i+1]-q[i])
        fq[i]=f
    for i in range(nx):
      for j in range(ny):
       for k in range(nz):
        ux=u[i,j,k,0];uy=u[i,j,k,1];uz=u[i,j,k,2];f0=f1=f2=0.
        for b in range(BONDS.shape[0]):
          dx,dy,dz=BONDS[b];nv0,nv1,nv2=NVECS[b];ks=ka if b<3 else kd
          ip=i+dx;jp=j+dy;kp=k+dz
          if ip<nx and jp>=0 and jp<ny and kp>=0 and kp<nz:
            e=nv0*(u[ip,jp,kp,0]-ux)+nv1*(u[ip,jp,kp,1]-uy)+nv2*(u[ip,jp,kp,2]-uz)
            f0+=ks*e*nv0;f1+=ks*e*nv1;f2+=ks*e*nv2
          im=i-dx;jm=j-dy;km=k-dz
          if im>=0 and jm>=0 and jm<ny and km>=0 and km<nz:
            e=nv0*(u[im,jm,km,0]-ux)+nv1*(u[im,jm,km,1]-uy)+nv2*(u[im,jm,km,2]-uz)
            f0+=ks*e*nv0;f1+=ks*e*nv1;f2+=ks*e*nv2
        fu[i,j,k,0]=f0;fu[i,j,k,1]=f1;fu[i,j,k,2]=f2
    r3=math.sqrt(3.);y0=ny//2;z0=nz//2
    for i in range(nx-1):
      gi=gprof[i]
      if gi==0:continue
      for oy in (-1,0):
       sy=-1. if oy==-1 else 1.;j=y0+oy
       for oz in (-1,0):
        sz=-1. if oz==-1 else 1.;k=z0+oz;ux=u[i,j,k,0];uy=u[i,j,k,1];uz=u[i,j,k,2]
        eL=(ux+sy*uy+sz*uz-q[i])/r3;eR=(-ux+sy*uy+sz*uz+q[i+1])/r3
        fu[i,j,k,0]+=-gi*(eL-eR)/r3;fu[i,j,k,1]+=-gi*sy*(eL+eR)/r3;fu[i,j,k,2]+=-gi*sz*(eL+eR)/r3
        fq[i]+=gi*eL/r3;fq[i+1]+=-gi*eR/r3

@njit(cache=True)
def coupling_current(q,vq,u,vu,gprof):
    """Current from q to bulk with half of each coupling spring on each side."""
    nx,ny,nz,_=u.shape;r3=math.sqrt(3.);y0=ny//2;z0=nz//2;J=0.
    for i in range(nx-1):
      gi=gprof[i]
      if gi==0:continue
      for oy in (-1,0):
       sy=-1. if oy==-1 else 1.;j=y0+oy
       for oz in (-1,0):
        sz=-1. if oz==-1 else 1.;k=z0+oz;ub=u[i,j,k];vb=vu[i,j,k]
        # Left spring n=(+1,sy,sz)/sqrt(3).
        eL=(ub[0]+sy*ub[1]+sz*ub[2]-q[i])/r3
        fbx=-gi*eL/r3;fby=-gi*eL*sy/r3;fbz=-gi*eL*sz/r3
        J+=.5*(vb[0]*fbx+vb[1]*fby+vb[2]*fbz+vq[i]*fbx)
        # Right spring n=(-1,sy,sz)/sqrt(3).
        eR=(-ub[0]+sy*ub[1]+sz*ub[2]+q[i+1])/r3
        fbx=gi*eR/r3;fby=-gi*eR*sy/r3;fbz=-gi*eR*sz/r3
        J+=.5*(vb[0]*fbx+vb[1]*fby+vb[2]*fbz+vq[i+1]*fbx)
    return J

@njit(cache=True)
def subsystem_energies(q,vq,u,vu,gprof,M,m,Kq,ka,kd):
    nx,ny,nz,_=u.shape;Eq=.5*M*np.sum(vq*vq);Eb=.5*m*np.sum(vu*vu);Ec=0.
    for i in range(nx-1):Eq+=.5*Kq*(q[i+1]-q[i])**2
    for i in range(nx):
     for j in range(ny):
      for k in range(nz):
       ui=u[i,j,k]
       for b in range(BONDS.shape[0]):
        dx,dy,dz=BONDS[b];ip=i+dx;jp=j+dy;kp=k+dz
        if ip<nx and jp>=0 and jp<ny and kp>=0 and kp<nz:
         nv=NVECS[b];e=nv[0]*(u[ip,jp,kp,0]-ui[0])+nv[1]*(u[ip,jp,kp,1]-ui[1])+nv[2]*(u[ip,jp,kp,2]-ui[2]);ks=ka if b<3 else kd;Eb+=.5*ks*e*e
    r3=math.sqrt(3.);y0=ny//2;z0=nz//2
    for i in range(nx-1):
     gi=gprof[i]
     if gi==0:continue
     for oy in (-1,0):
      sy=-1. if oy==-1 else 1.;j=y0+oy
      for oz in (-1,0):
       sz=-1. if oz==-1 else 1.;k=z0+oz;ub=u[i,j,k]
       eL=(ub[0]+sy*ub[1]+sz*ub[2]-q[i])/r3;eR=(-ub[0]+sy*ub[1]+sz*ub[2]+q[i+1])/r3;Ec+=.5*gi*(eL*eL+eR*eR)
    return Eq,Eb,Ec

def modal_bulk_energy(u,vu,p:ModelParameters):
    nx,ny,nz,_=u.shape;N=nx*ny*nz;uh=np.fft.fftn(u,axes=(0,1,2));vh=np.fft.fftn(vu,axes=(0,1,2))
    kx=2*np.pi*np.fft.fftfreq(nx)[:,None,None];ky=2*np.pi*np.fft.fftfreq(ny)[None,:,None];kz=2*np.pi*np.fft.fftfreq(nz)[None,None,:]
    vals,vec=np.linalg.eigh(bulk_stiffness_batch(kx,ky,kz,p)/p.m_bulk)
    au=np.einsum('...am,...a->...m',vec.conj(),uh,optimize=True);av=np.einsum('...am,...a->...m',vec.conj(),vh,optimize=True)
    E=.5*p.m_bulk*(np.abs(av)**2+vals*np.abs(au)**2)/N
    KX,KY,KZ=np.meshgrid(2*np.pi*np.fft.fftfreq(nx),2*np.pi*np.fft.fftfreq(ny),2*np.pi*np.fft.fftfreq(nz),indexing='ij');kn=np.sqrt(KX*KX+KY*KY+KZ*KZ);kh=np.zeros((nx,ny,nz,3));msk=kn>1e-14
    kh[...,0][msk]=KX[msk]/kn[msk];kh[...,1][msk]=KY[msk]/kn[msk];kh[...,2][msk]=KZ[msk]/kn[msk]
    lf=np.abs(np.einsum('...a,...am->...m',kh,vec,optimize=True))**2
    total=float(E.sum());long=float((E*lf).sum())
    return {'total':total,'longitudinal':long,'transverse':total-long,'transverse_fraction':(total-long)/total if total else 0.}

def incident_spectrum(q,vq,p:ModelParameters):
    n=len(q);qh=np.fft.fft(q);vh=np.fft.fft(vq);k=2*np.pi*np.fft.fftfreq(n);om=2*np.sqrt(p.k_spine/p.m_spine)*np.abs(np.sin(k/2))
    e=.5*p.m_spine*(np.abs(vh)**2+om**2*np.abs(qh)**2)/n
    return om,e

def run(*,nx=180,ny=30,nz=30,dt=.04,tmax=95.,omega0=1.2,sigma=10.,x0=35.,slab_start=88,
        g_segments=(.1,.2,.2,.2,.2,.1),amplitude=2e-3,sponge_width=8,sponge_gamma=.55,p=None,snapshot_time=70.,compute_snapshot=True,verbose=False,snapshot_method='kinetic'):
    if p is None:p=ModelParameters()
    x=np.arange(nx);k0=carrier_k(omega0,p);q=amplitude*np.exp(-.5*((x-x0)/sigma)**2)*np.cos(k0*(x-x0));vq=exact_right_velocity(q,p)
    u=np.zeros((nx,ny,nz,3));vu=np.zeros_like(u);gprof=np.zeros(nx-1)
    gprof[slab_start:slab_start+len(g_segments)]=np.asarray(g_segments)
    dq=damping_profile_1d(nx,sponge_width,sponge_gamma);dbx=dq[:,None,None];dby=damping_profile_1d(ny,sponge_width,sponge_gamma)[None,:,None];dbz=damping_profile_1d(nz,sponge_width,sponge_gamma)[None,None,:];db=np.maximum(np.maximum(dbx,dby),dbz)
    fq=np.zeros_like(q);fu=np.zeros_like(u);forces_open(q,u,gprof,p.k_spine,p.k_axis,p.k_diag,fq,fu)
    E0=subsystem_energies(q,vq,u,vu,gprof,p.m_spine,p.m_bulk,p.k_spine,p.k_axis,p.k_diag)[0]
    om_spec,e_spec=incident_spectrum(q,vq,p);steps=int(round(tmax/dt));hist=[];cum=0.;Jold=coupling_current(q,vq,u,vu,gprof);snap=None
    dampq=np.exp(-.5*dq*dt);dampb=np.exp(-.5*db*dt)
    for n in range(steps):
        vq*=dampq;vu*=dampb[...,None]
        vq+=.5*dt*fq/p.m_spine;vu+=.5*dt*fu/p.m_bulk;q+=dt*vq;u+=dt*vu
        forces_open(q,u,gprof,p.k_spine,p.k_axis,p.k_diag,fq,fu)
        vq+=.5*dt*fq/p.m_spine;vu+=.5*dt*fu/p.m_bulk
        vq*=dampq;vu*=dampb[...,None]
        J=coupling_current(q,vq,u,vu,gprof);cum+=.5*(Jold+J)*dt;Jold=J;t=(n+1)*dt
        if compute_snapshot and snap is None and t>=snapshot_time:
            if snapshot_method=='kinetic':
                _,eb,ec=subsystem_energies(q,vq,u,vu,gprof,p.m_spine,p.m_bulk,p.k_spine,p.k_axis,p.k_diag);ft=kinetic_transverse_fraction(vu);tot=eb+.5*ec;snap={'total':tot,'longitudinal':tot*(1-ft),'transverse':tot*ft,'transverse_fraction':ft,'time':t,'cross_energy':0.}
            else:snap=projector_bulk_energy(u.copy(),vu.copy(),p)|{'time':t}
        if n%100==0 or n==steps-1:
            if verbose: print('td step',n+1,'of',steps,'time',t,'T',cum/E0,flush=True)
            Eq,Eb,Ec=subsystem_energies(q,vq,u,vu,gprof,p.m_spine,p.m_bulk,p.k_spine,p.k_axis,p.k_diag)
            hist.append({'time':t,'current_integral':cum,'T_current':cum/E0,'Eq':Eq,'Ebulk':Eb,'Ecoupling':Ec})
    if snap is None:
        
        if compute_snapshot and snapshot_method=='kinetic':
            _,eb,ec=subsystem_energies(q,vq,u,vu,gprof,p.m_spine,p.m_bulk,p.k_spine,p.k_axis,p.k_diag);ft=kinetic_transverse_fraction(vu);tot=eb+.5*ec;snap={'total':tot,'longitudinal':tot*(1-ft),'transverse':tot*ft,'transverse_fraction':ft,'time':tmax,'cross_energy':0.}
        elif compute_snapshot:snap=projector_bulk_energy(u,vu,p)|{'time':tmax}
        else:snap={'total':0.,'longitudinal':0.,'transverse':0.,'transverse_fraction':0.,'time':tmax,'cross_energy':0.}
    return {'parameters':{'nx':nx,'ny':ny,'nz':nz,'dt':dt,'tmax':tmax,'omega0':omega0,'sigma':sigma,'x0':x0,'slab_start':slab_start,'g_segments':list(g_segments),'amplitude':amplitude,'sponge_width':sponge_width,'sponge_gamma':sponge_gamma},
            'initial_incident_energy':E0,'Tbulk_current':cum/E0,'snapshot':snap|{'Tbulk':snap['total']/E0,'TL':snap['longitudinal']/E0,'TT':snap['transverse']/E0},
            'history':hist,'spectrum_omega':om_spec.tolist(),'spectrum_energy':e_spec.tolist()}

def projector_bulk_energy(u,vu,p:ModelParameters):
    """Fast basis-invariant Helmholtz energy partition without branch labelling.

    The exact lattice quadratic energy is partitioned symmetrically between
    k-parallel and k-orthogonal fields; the L/T cross term is split equally so
    the two reported components sum exactly to the periodic quadratic energy.
    """
    nx,ny,nz,_=u.shape;N=nx*ny*nz
    uh=np.fft.fftn(u,axes=(0,1,2));vh=np.fft.fftn(vu,axes=(0,1,2))
    kx=2*np.pi*np.fft.fftfreq(nx)[:,None,None];ky=2*np.pi*np.fft.fftfreq(ny)[None,:,None];kz=2*np.pi*np.fft.fftfreq(nz)[None,None,:]
    kn2=kx*kx+ky*ky+kz*kz;safe=np.where(kn2>1e-28,kn2,1.)
    du=kx*uh[...,0]+ky*uh[...,1]+kz*uh[...,2];dv=kx*vh[...,0]+ky*vh[...,1]+kz*vh[...,2]
    uL=np.stack([kx*du/safe,ky*du/safe,kz*du/safe],axis=-1);vL=np.stack([kx*dv/safe,ky*dv/safe,kz*dv/safe],axis=-1)
    uT=uh-uL;vT=vh-vL
    K=bulk_stiffness_batch(kx,ky,kz,p)
    KuL=np.einsum('...ab,...b->...a',K,uL,optimize=True);KuT=np.einsum('...ab,...b->...a',K,uT,optimize=True)
    kinL=.5*p.m_bulk*np.sum(np.abs(vL)**2)/N;kinT=.5*p.m_bulk*np.sum(np.abs(vT)**2)/N
    potL=.5*np.real(np.sum(uL.conj()*KuL))/N;potT=.5*np.real(np.sum(uT.conj()*KuT))/N
    cross=np.real(np.sum(uL.conj()*KuT))/N
    EL=float(kinL+potL+.5*cross);ET=float(kinT+potT+.5*cross);tot=EL+ET
    return {'total':tot,'longitudinal':EL,'transverse':ET,'transverse_fraction':ET/tot if tot else 0.,'cross_energy':float(cross)}

def kinetic_transverse_fraction(vu):
    """Fast basis-invariant kinematic L/T split of a real-space bulk snapshot."""
    nx,ny,nz,_=vu.shape;vh=np.fft.fftn(vu,axes=(0,1,2))
    kx=2*np.pi*np.fft.fftfreq(nx)[:,None,None];ky=2*np.pi*np.fft.fftfreq(ny)[None,:,None];kz=2*np.pi*np.fft.fftfreq(nz)[None,None,:]
    k2=kx*kx+ky*ky+kz*kz;safe=np.where(k2>1e-28,k2,1.)
    dot=kx*vh[...,0]+ky*vh[...,1]+kz*vh[...,2]
    l2=np.sum(np.abs(dot)**2/safe);tot=np.sum(np.abs(vh)**2)
    return float(1-l2/tot) if tot else 0.
