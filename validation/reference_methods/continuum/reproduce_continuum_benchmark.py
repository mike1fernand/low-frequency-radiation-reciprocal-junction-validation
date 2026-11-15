#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
omega=1.2;cT=1.;cL=math.sqrt(3);kappa=0.9594572753515053;g=np.array([.1,.2,.2,.2,.2,.1]);mid=np.arange(6)+.5

def source_components(kx,ky,kz,mode):
 A=np.sum(g[:,None]*np.exp(1j*(kappa-kx.ravel())[None,:]*mid[:,None]),axis=0).reshape(kx.shape)
 cy=np.cos(ky/2);sy=np.sin(ky/2);cz=np.cos(kz/2);sz=np.sin(kz/2);out=np.zeros(kx.shape+(3,),complex)
 if mode in ('X','full'):out[...,0]=(8/3)*math.cos(kappa/2)*A*cy*cz
 if mode in ('R','full'):
  out[...,1]=-(8/3)*math.sin(kappa/2)*A*sy*cz;out[...,2]=-(8/3)*math.sin(kappa/2)*A*cy*sz
 return out

def integrate(nmu,nphi,mode):
 mu,w=np.polynomial.legendre.leggauss(nmu);phi=(np.arange(nphi)+.5)*2*np.pi/nphi;MU,PH=np.meshgrid(mu,phi,indexing='ij');st=np.sqrt(np.maximum(0,1-MU*MU));n=np.stack([MU,st*np.cos(PH),st*np.sin(PH)],axis=-1);wt=w[:,None]*(2*np.pi/nphi);kL=(omega/cL)*n;kT=(omega/cT)*n;FL=source_components(kL[...,0],kL[...,1],kL[...,2],mode);FT=source_components(kT[...,0],kT[...,1],kT[...,2],mode);P=float(np.sum(wt*np.abs(np.sum(n*FL,axis=-1))**2)/cL**3);S=float(np.sum(wt*(np.sum(np.abs(FT)**2,axis=-1)-np.abs(np.sum(n*FT,axis=-1))**2))/cT**3);return {'P_weight':P,'S_weight':S,'S_to_P_ratio':S/P,'S_fraction':S/(S+P)}
conv=[];last=None
for nmu,nphi in [(60,120),(120,360),(240,720)]:
 r={m:integrate(nmu,nphi,m) for m in ('X','R','full')};conv.append({'n_mu':nmu,'n_phi':nphi,'results':r});last=r
point=2*(cL/cT)**3/(1+2*(cL/cT)**3)
out={'schema':'paper1.p06.continuum-first-born.v1','parameters':{'omega':omega,'c_T':cT,'c_L':cL,'kappa':kappa,'profile':g.tolist()},'formula':'P weight = integral |n dot F(k_L n)|^2/c_L^3; S weight = integral (|F(k_T n)|^2-|n dot F|^2)/c_T^3','point_force_S_fraction':point,'quadrature_convergence':conv,'finite_source':last,'status':'passed' if max(x['results']['full']['S_fraction'] for x in conv)-min(x['results']['full']['S_fraction'] for x in conv)<5e-14 else 'failed'}
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
