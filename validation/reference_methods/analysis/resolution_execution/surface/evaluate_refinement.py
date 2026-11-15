#!/usr/bin/env python3
"""Separate spectral-grid, sample-count, seed, PV-order and bin sensitivities.
Reported spreads and changes are observed diagnostics, never true-error bounds.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'src'))
import spectral_histogram as sh
from reduced_core import scatter_from_reduced_green
from model_api_v1_1 import ModelParameters

def transmit(w,pv,branches):
 res,x=scatter_from_reduced_green(w,[.1,.2,.2,.2,.2,.1],ModelParameters(),pv[:12,:12]+1j*sum(branches),return_internal=True)
 bs=[float(2*x['gamma_in']*np.vdot(x['t'],b@x['t']).real) for b in branches]
 return {'qT':sum(bs[:2]),'qL':bs[2],'Tbulk':res.Tbulk,'Rq':res.Rq,'Tline':res.Tq,'amplitude_closure':res.closure_error,'partition_residual':abs(sum(bs)-res.Tbulk)}

def evaluate():
 recs=[];means=[];spectral=[];pvmatrices={}
 available={}
 for w in [1.2,1.0]:
  rows=[]
  for p in sorted((HERE/'results').glob(f'omega{w:g}_N*.json')):
   r=json.loads(p.read_text());
   if r.get('spectral_support_segments',6)!=6:continue
   data=np.load(HERE/'results'/r['matrix_file']);rows.append((r['grid'],r,data));spectral.append(r)
  if not rows:continue
  N,record,data=max(rows,key=lambda x:x[0]);available[w]=(N,[data[f'Im{i}'] for i in range(3)])
 for power in [20,21,22,23]:
  hists=[]
  use10=all((HERE/'histograms'/f'sobol_seed{seed}_p{power}_nseg10_bins16384.npz').exists() for seed in [20260920,20260921])
  for seed in [20260920,20260921]:
   file=HERE/'histograms'/(f'sobol_seed{seed}_p{power}'+('_nseg10' if use10 else '')+'_bins16384.npz')
   if not file.exists():continue
   hists.append((seed,sh.load_histogram(file)))
  for bins in [8192,16384]:
   for order in [192,256,320,384]:
    for w,(N,branches) in available.items():
     pvs=[]
     for seed,h in hists:
      work=h if bins==16384 else sh.SpectralHistogram(h.edges[::2],h.mass.reshape(h.nseg,2,2,8192,2).sum(axis=-1),h.sample_count,h.sampler,{**h.metadata,'nbins':8192,'coarsened_from_bins':16384})
      pv=work.principal_value_jackson(w*w,order=order);pvs.append(pv)
      recs.append({'omega':w,'surface_grid':N,'power':power,'sample_count':h.sample_count,'seed':seed,'nseg':h.nseg,'bins':bins,'order':order,'coefficient':transmit(w,pv,branches),'PV_norm':float(np.linalg.norm(pv))})
     if len(pvs)==2:
      mean=np.mean(pvs,axis=0);means.append({'omega':w,'surface_grid':N,'power':power,'sample_count_per_seed':2**power,'seeds':[20260920,20260921],'nseg':h.nseg,'bins':bins,'order':order,'coefficient':transmit(w,mean,branches),'PV_norm':float(np.linalg.norm(mean))})
      pvmatrices[f'omega{w:g}_p{power}_bins{bins}_M{order}']=mean
 np.savez_compressed(HERE/'PV_MEAN_MATRICES.npz',**pvmatrices)
 summary=[]
 for w,(N,branches) in available.items():
  for power in [20,21,22,23]:
   for order in [192,256,320,384]:
    vals=[r['coefficient']['qT'] for r in recs if r['omega']==w and r['power']==power and r['bins']==16384 and r['order']==order]
    mm=[r for r in means if r['omega']==w and r['power']==power and r['bins']==16384 and r['order']==order]
    if len(vals)==2 and mm:
     mean=mm[0]['coefficient']['qT'];summary.append({'omega':w,'N':N,'power':power,'order':order,'mean_qT':mean,'seed_qT_min':min(vals),'seed_qT_max':max(vals),'seed_span_relative':(max(vals)-min(vals))/mean})
 out={'scope':'Observed sensitivities from fresh sample/count/bin/order/surface changes; no true-error bound or statistical confidence claim.','surface_records':spectral,'PV_individual_records':recs,'PV_mean_records':means,'summary':summary,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'REFINEMENT_RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':evaluate()
