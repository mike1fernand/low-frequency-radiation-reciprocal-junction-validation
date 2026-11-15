#!/usr/bin/env python3
"""Check all changed-medium samples retain tied modes in a common energy bin."""
from pathlib import Path
import sys,json,hashlib,time
sys.dont_write_bytecode=True
import numpy as np
from scipy.stats import qmc
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[2];sys.path.insert(0,str(HERE/'src'))
import spectral_histogram as sh
from model_api_v1_1 import ModelParameters

def need(x,msg):
 if not x:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rows=[];t0=time.perf_counter()
for kd in (.75,.9,1.,1.1,1.25):
 p=ModelParameters(k_diag=kd);edges=np.linspace(0.,4+8*kd,6145);pairs=split=outside=points=0
 def add(kx,ky,kz):
  global pairs,split,outside,points
  vals,_=sh._mode_data(kx,ky,kz,p)
  tol=64*np.finfo(float).eps*np.maximum(1.,np.max(abs(vals),axis=-1))
  tied=np.diff(vals,axis=-1)<=tol[...,None]
  ids=np.clip(np.searchsorted(edges,vals,side='right')-1,0,len(edges)-2)
  pairs+=int(tied.sum());split+=int((tied & (np.diff(ids,axis=-1)!=0)).sum())
  outside+=int(((vals<edges[0])|(vals>edges[-1])).sum());points+=vals.size//3
 for sx,sy,sz in ((.17,.39,.61),(.67,.89,.11)):
  ar=np.arange(64);x=-np.pi+2*np.pi*(ar+sx)/64;y=-np.pi+2*np.pi*(ar+sy)/64;z=-np.pi+2*np.pi*(ar+sz)/64
  for st in range(0,64,4):add(x[st:st+4,None,None],y[None,:,None],z[None,None,:])
 k=2*np.pi*(qmc.Sobol(3,scramble=True,seed=97531).random_base2(19)-.5)
 for st in range(0,len(k),65536):add(k[st:st+65536,0],k[st:st+65536,1],k[st:st+65536,2])
 need(points==2*64**3+2**19,'Stiffness sample coverage count differs')
 need(split==0,'Numerically tied eigenspace split between histogram bins')
 need(outside==0,'A sample is outside analytic spectral support')
 mass=[]
 for sampler in ('grid64x2','sobol19'):
  file=HERE/'stiffness'/f'kdiag{kd:g}_N80_{sampler}_hist.npz'
  with np.load(file,allow_pickle=False) as z:M=z['mass'].sum(axis=-1)
  expected=np.zeros_like(M);expected[0]=np.eye(2)
  residual=float(np.max(abs(M-expected)));need(residual<1e-4,'Histogram gross mass sum-rule failure')
  mass.append(dict(sampler=sampler,maximum_sum_rule_error=residual,file_sha256=digest(file)))
 rows.append(dict(k_diag=kd,sample_points=points,numerically_tied_pairs=pairs,tied_pairs_split_across_bins=split,outside_spectral_support=outside,histogram_sum_rule_checks=mass))
result=dict(status='all_changed_medium_sample_checks_passed',rows=rows,seconds=time.perf_counter()-t0,source_sha256={str(p.relative_to(BASE)):digest(p) for p in [Path(__file__),*sorted((HERE/'src').glob('*.py'))]},scope='Actual original sampler points, no omitted samples; ties use the same finite-precision 64eps rule as the nodal source audit. No quadrature convergence or exact-error conclusion.')
(HERE/'STIFFNESS_SAMPLE_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
