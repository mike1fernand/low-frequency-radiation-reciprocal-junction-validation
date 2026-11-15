#!/usr/bin/env python3
"""Full 12-coordinate/segment reference used only to verify the exact reduction."""
from __future__ import annotations
import math
import numpy as np
from model_api_v1_1 import ModelParameters,AXIAL_BONDS,FACE_DIAGONAL_BONDS,_unit
CORE_TRANS=np.array([[-1,-1],[-1,0],[0,-1],[0,0]],int)
CORE_SIGNS=np.array([[-1.,-1.],[-1.,1.],[1.,-1.],[1.,1.]])
def finite_core_matrices(g_segments,p):
 g=np.asarray(g_segments,float);n=len(g);coords=np.empty((4*n,3),int);U=np.zeros((12*n,12*n),complex);V=np.zeros((12*n,n+1),complex);Ks=np.zeros((n+1,n+1))
 for s,gs in enumerate(g):
  for a,((yy,zz),(sy,sz)) in enumerate(zip(CORE_TRANS,CORE_SIGNS)):
   site=4*s+a;coords[site]=[s,yy,zz];sl=slice(3*site,3*site+3);nm=np.array([1.,sy,sz])/math.sqrt(3);np_=np.array([-1.,sy,sz])/math.sqrt(3)
   U[sl,sl]+=gs*(np.outer(nm,nm)+np.outer(np_,np_));V[sl,s]+=gs*nm*(-1/math.sqrt(3));V[sl,s+1]+=gs*np_*(1/math.sqrt(3));Ks[s,s]+=gs/3;Ks[s+1,s+1]+=gs/3
 return coords,U,V,Ks
def bare_green(omega,eta,g_segments,p,grid=32):
 n=len(g_segments);coords,_,_,_=finite_core_matrices(g_segments,p);N=grid;k=2*math.pi*np.fft.fftfreq(N);kx=k[:,None,None];ky=k[None,:,None];kz=k[None,None,:]
 K=np.zeros((N,N,N,3,3))
 for d,ks in ([(d,p.k_axis) for d in AXIAL_BONDS]+[(d,p.k_diag) for d in FACE_DIAGONAL_BONDS]):
  nh=_unit(d);ph=kx*d[0]+ky*d[1]+kz*d[2];K+=2*ks*(1-np.cos(ph))[...,None,None]*np.outer(nh,nh)
 Gk=np.linalg.inv(K.astype(complex)-p.m_bulk*complex(omega,eta)**2*np.eye(3));gr=np.empty_like(Gk)
 for a in range(3):
  for b in range(3):gr[...,a,b]=np.fft.ifftn(Gk[...,a,b])
 G=np.empty((12*n,12*n),complex);sh=np.array(gr.shape[:3])
 for i,ri in enumerate(coords):
  for j,rj in enumerate(coords):
   d=(ri-rj)%sh;G[3*i:3*i+3,3*j:3*j+3]=gr[int(d[0]),int(d[1]),int(d[2])]
 return G
