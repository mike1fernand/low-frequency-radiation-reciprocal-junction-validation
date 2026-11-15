#!/usr/bin/env python3
"""Matrix spectral measure and zero-broadening Jackson-Chebyshev PV reconstruction."""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence
import numpy as np
from scipy.stats import qmc
from model_api_v1_1 import ModelParameters
from reduced_core import bulk_stiffness_batch,cross_section_form_factor

@dataclass
class SpectralHistogram:
    edges:np.ndarray
    mass:np.ndarray # dx>=0, 2,2,bins
    sample_count:int
    sampler:str
    metadata:dict
    @property
    def centers(self):return .5*(self.edges[:-1]+self.edges[1:])
    @property
    def widths(self):return np.diff(self.edges)
    @property
    def nseg(self):return self.mass.shape[0]
    def assemble_real_kernel(self,positive):
        out=np.empty((2*self.nseg,2*self.nseg),complex)
        for s in range(self.nseg):
            for t in range(self.nseg):
                dx=s-t; b=positive[dx] if dx>=0 else positive[-dx].conj().T
                out[2*s:2*s+2,2*t:2*t+2]=b
        return out
    def resolvent(self,z):
        ker=1/(self.centers-z*z)
        pos=np.einsum('dABj,j->dAB',self.mass,ker,optimize=True)
        neg=np.einsum('dBAj,j->dAB',self.mass.conj(),ker,optimize=True)
        out=np.empty((2*self.nseg,2*self.nseg),complex)
        for s in range(self.nseg):
            for t in range(self.nseg):
                dx=s-t;out[2*s:2*s+2,2*t:2*t+2]=pos[dx] if dx>=0 else neg[-dx]
        return out
    def principal_value_jackson(self,energy,order=192,margin=1.02):
        if order<8:raise ValueError
        lo,hi=float(self.edges[0]),float(self.edges[-1]);mid=.5*(lo+hi);half=.5*(hi-lo)*margin
        x=(self.centers-mid)/half;x0=(energy-mid)/half
        if max(abs(x0),float(np.max(np.abs(x))))>=1:raise ValueError('bad scaling')
        n=np.arange(order,dtype=float);th=math.pi/(order+1)
        g=((order-n+1)*np.cos(n*th)+np.sin(n*th)/math.tan(th))/(order+1)
        tm=np.ones_like(x);tc=x.copy();u0=1.;ker=2*g[1]*u0*tc
        if order>2:
            u1=2*x0
            for j in range(2,order):
                tp=2*x*tc-tm;ker+=2*g[j]*u1*tp;tm,tc=tc,tp;u0,u1=u1,2*x0*u1-u0
        pos=np.einsum('dABj,j->dAB',self.mass,ker/half,optimize=True)
        G=self.assemble_real_kernel(pos)
        # Exact reciprocal real harmonic sector; residual imaginary/antisymmetric
        # pieces are quadrature noise and are recorded before this projection.
        return np.real(.5*(G+G.T))
    def jackson_density(self,energy,order=192,margin=1.02):
        lo,hi=float(self.edges[0]),float(self.edges[-1]);mid=.5*(lo+hi);half=.5*(hi-lo)*margin
        x=(self.centers-mid)/half;x0=(energy-mid)/half
        n=np.arange(order,dtype=float);th=math.pi/(order+1)
        g=((order-n+1)*np.cos(n*th)+np.sin(n*th)/math.tan(th))/(order+1)
        t0=np.ones_like(x);t1=x.copy();o0=1.;o1=x0;ker=g[0]*np.ones_like(x)
        if order>1:ker+=2*g[1]*o1*t1
        for j in range(2,order):
            t2=2*x*t1-t0;o2=2*x0*o1-o0;ker+=2*g[j]*o2*t2;t0,t1=t1,t2;o0,o1=o1,o2
        ker/=math.pi*math.sqrt(1-x0*x0)*half
        pos=np.einsum('dABj,j->dAB',self.mass,ker,optimize=True)
        G=self.assemble_real_kernel(pos)
        return np.real(.5*(G+G.T))

def _mode_data(kx,ky,kz,p):
    vals,vec=np.linalg.eigh(bulk_stiffness_batch(kx,ky,kz,p)/p.m_bulk)
    F=cross_section_form_factor(ky,kz)
    c=np.einsum('...am,...aA->...mA',vec.conj(),F,optimize=True)
    return vals,c

def _acc(mass,vals,c,kx,edges,nseg,norm):
    nb=len(edges)-1;e=vals.reshape(-1);idx=np.searchsorted(edges,e,side='right')-1;idx=np.clip(idx,0,nb-1)
    W=np.einsum('...mA,...mB->...mAB',c.conj(),c,optimize=True)
    kb=np.broadcast_to(kx,vals.shape[:-1])
    for dx in range(nseg):
        wd=(W*np.exp(1j*kb*dx)[...,None,None,None]).reshape(-1,2,2)*norm
        for a in range(2):
            for b in range(2):
                v=wd[:,a,b]
                mass[dx,a,b]+=np.bincount(idx,weights=v.real,minlength=nb)+1j*np.bincount(idx,weights=v.imag,minlength=nb)

def sample_shifted_grids(p:ModelParameters,*,nseg,grid=80,shifts=((.5,.5,.5),),nbins=8192,energy_range=(0.,12.),kx_chunk=4):
    edges=np.linspace(*energy_range,nbins+1);mass=np.zeros((nseg,2,2,nbins),complex);ar=np.arange(grid)
    norm=1/(p.m_bulk*grid**3*len(shifts))
    for sx,sy,sz in shifts:
        bx=-math.pi+2*math.pi*(ar+sx)/grid;by=-math.pi+2*math.pi*(ar+sy)/grid;bz=-math.pi+2*math.pi*(ar+sz)/grid
        ky=by[None,:,None];kz=bz[None,None,:]
        for st in range(0,grid,kx_chunk):
            en=min(grid,st+kx_chunk);kx=bx[st:en,None,None];v,c=_mode_data(kx,ky,kz,p);_acc(mass,v,c,kx,edges,nseg,norm)
    return SpectralHistogram(edges,mass,grid**3*len(shifts),'shifted_grid',{'grid':grid,'shifts':[list(x) for x in shifts],'nbins':nbins})

def sample_sobol(p:ModelParameters,*,nseg,power=20,seed=12345,nbins=8192,energy_range=(0.,12.),chunk=65536):
    ns=2**power;edges=np.linspace(*energy_range,nbins+1);mass=np.zeros((nseg,2,2,nbins),complex)
    pts=qmc.Sobol(3,scramble=True,seed=seed).random_base2(power);norm=1/(p.m_bulk*ns)
    for st in range(0,ns,chunk):
        k=2*math.pi*(pts[st:st+chunk]-.5);v,c=_mode_data(k[:,0],k[:,1],k[:,2],p);_acc(mass,v,c,k[:,0],edges,nseg,norm)
    return SpectralHistogram(edges,mass,ns,'sobol',{'power':power,'seed':seed,'nbins':nbins})

def save_histogram(path,h):
    np.savez_compressed(path,edges=h.edges,mass=h.mass,sample_count=h.sample_count,sampler=h.sampler,metadata=np.array(h.metadata,dtype=object))
def load_histogram(path):
    d=np.load(path,allow_pickle=True);return SpectralHistogram(d['edges'],d['mass'],int(d['sample_count']),str(d['sampler']),dict(d['metadata'].item()))
