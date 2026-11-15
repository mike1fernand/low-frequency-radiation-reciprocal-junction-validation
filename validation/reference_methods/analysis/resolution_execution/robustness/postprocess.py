#!/usr/bin/env python3
"""Build the active 85-row ledger after an explicit fresh-PV selection."""
from pathlib import Path
import argparse,hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[2];sys.path.insert(0,str(HERE/'src'))
import spectral_histogram as hist
from model_api_v1_1 import ModelParameters
from reduced_core import scatter_from_reduced_green
from worker import calc
BASE_G=np.array([.1,.2,.2,.2,.2,.1])

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bound_path(item):
 p=(BASE/item['path']).resolve();need(BASE in p.parents and p.is_file(),'Missing/outside input');need(digest(p)==item['sha256'],'Input hash mismatch: '+str(p));return p

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pv-contract',type=Path,required=True);a=ap.parse_args();contract=json.loads(a.pv_contract.read_text())
 need(contract.get('schema')=='paper1.robustness-pv-contract.v1' and contract.get('status')=='selected_for_declared_finite_robustness','No selected fresh PV contract')
 need(contract.get('physical_bulk')==dict(k_axis=1.,k_diag=1.,m_bulk=1.),'PV physical medium differs')
 order=contract['order'];bins=contract['bins'];need(isinstance(order,int) and order>=8 and isinstance(bins,int) and bins>0,'Invalid PV order/bins')
 need(len(contract['histograms'])==2,'This declared contract requires two equally weighted fresh histograms')
 paths=[bound_path(x) for x in contract['histograms']];hists=[hist.load_histogram(p) for p in paths]
 for h in hists:
  need(h.nseg==10 and h.mass.shape==(10,2,2,bins),'Fresh histogram must contain ten actual separation blocks at declared bins')
  need(np.isfinite(h.mass).all() and np.isfinite(h.edges).all(),'Nonfinite selected histogram')
  need(h.edges[0]==0 and h.edges[-1]==12 and np.all(np.diff(h.edges)>0),'Selected histogram spectral support differs')
 sample_audit=json.loads((HERE/'STIFFNESS_SAMPLE_AUDIT.json').read_text())
 need(sample_audit['status']=='all_changed_medium_sample_checks_passed' and len(sample_audit['rows'])==5,'Stiffness sample audit incomplete')
 need(all(r['outside_spectral_support']==0 and r['tied_pairs_split_across_bins']==0 for r in sample_audit['rows']),'Stiffness range/degenerate-bin audit failed')
 run=json.loads((HERE/'CHECKPOINT_RUN.json').read_text());need(run['status']=='completed_all_54_kernel_jobs' and run['completed']==54,'Kernel jobs incomplete')
 for record in run['records']:
  p=BASE/record['record'];need(digest(p)==record['sha256'],'Checkpoint JSON changed')
  r=json.loads(p.read_text());need(digest(p.parent/r['matrix_file'])==r['matrix_sha256'],'Checkpoint matrix changed')
  need(len(r['geometry'])==3 and all(g['omitted_gradient_cells']==0 and g['nonpositive_areas']==0 for g in r['geometry']),'Geometry omissions in checkpoint')
  if r['family']=='stiffness':
   need(r['spectral_range_audit']['roundoff_outside_support']==0,'Histogram samples occurred outside chosen support')
   for x in r['histogram_inputs']:need(digest(p.parent/x['file'])==x['sha256'],'Stiffness histogram changed')
 def load(w,n):
  stem=HERE/'kernels'/f'omega{w:.6f}_N80_nseg{n}';record=json.loads(Path(str(stem)+'.json').read_text())
  with np.load(Path(str(stem)+'.npz'),allow_pickle=False) as z:return [z[f'Im_branch_{i}'] for i in range(3)]
 def pv(w,n):
  full=np.mean([h.principal_value_jackson(w*w,order=order) for h in hists],axis=0)
  need(full.shape==(20,20) and np.isfinite(full).all() and np.array_equal(full,full.T),'Invalid full-support PV reconstruction')
  return full[:2*n,:2*n]
 refbs10=load(1.2,10);refpv10=pv(1.2,10);refbs=[m[:12,:12] for m in refbs10];refpv=refpv10[:12,:12]
 ref=calc(1.2,BASE_G,ModelParameters(),refpv,refbs)
 coupling=[dict(scale=s,**calc(1.2,BASE_G*s,ModelParameters(),refpv,refbs)) for s in (0.,.05,.075,.1,.15,.2,.25,.5,1.,1.5,2.)]
 profiles={'reference':BASE_G,'uniform':np.full(6,.2),'soft_ramp':np.array([.05,.12,.2,.2,.12,.05]),'short_core':np.array([0.,.1,.2,.2,.1,0.]),'edge_heavy':np.array([.2,.1,.1,.1,.1,.2])}
 prows=[dict(name=k,g=g.tolist(),**calc(1.2,g,ModelParameters(),refpv,refbs)) for k,g in profiles.items()]
 speeds=[dict(c_spine_longwave=c,**calc(1.2,BASE_G,ModelParameters(k_spine=c*c),refpv,refbs)) for c in (.65,.8,.9,1.,1.1,1.2,1.3,1.5,1.7,2.)]
 units=[dict(a=aa,**calc(1.2,BASE_G,ModelParameters(a=aa),refpv,refbs)) for aa in (.25,.5,1.,2.,4.)]
 lengths=[]
 for n in (2,3,4,6,8,10):
  g=np.full(n,.2);g[0]=g[-1]=.1
  lengths.append(dict(nseg=n,**calc(1.2,g,ModelParameters(),refpv10[:2*n,:2*n],[M[:2*n,:2*n] for M in refbs10])))
 frequencies=[]
 for j in range(49):
  w=round(.6+.025*j,6);n=10 if j==24 else 6
  frequencies.append(dict(index=j,omega=w,method=dict(surface_grid=80,PV_order=order,PV_bins=bins),**calc(w,BASE_G,ModelParameters(),pv(w,6),[M[:12,:12] for M in load(w,n)])))
 stiffness=[]
 for kd in (.75,.9,1.,1.1,1.25):
  file=HERE/'stiffness'/f'kdiag{kd:g}_N80.json';r=json.loads(file.read_text());stiffness.append(dict(k_diag_over_k_axis=kd,method=dict(surface_grid=80,PV_order=160,PV_bins=6144,energy_range=r['histogram_energy_range']),**r['coefficient']))
 families=dict(frequency=frequencies,coupling_scale=[r for r in coupling if r['scale']>0],length=lengths,profile=prows,source_speed=speeds,stiffness_ratio=stiffness)
 expected=dict(frequency=49,coupling_scale=10,length=6,profile=5,source_speed=10,stiffness_ratio=5)
 need({k:len(v) for k,v in families.items()}==expected,'Active family enumeration changed')
 allrows=[r for rows in families.values() for r in rows];need(len(allrows)==85,'85-row coverage failed')
 need(all(r['T_qT_modal']>0 for r in allrows),'A sampled qT coefficient is not positive')
 weak=[r for r in coupling if 0<r['scale']<=.2];slope=float(np.polyfit(np.log([r['scale'] for r in weak]),np.log([r['T_qT_modal'] for r in weak]),1)[0])
 ranges={k:dict(rows=len(rows),qT_min=min(r['T_qT_modal'] for r in rows),qT_max=max(r['T_qT_modal'] for r in rows),qT_fraction_min=min(r['qT_fraction_of_bulk'] for r in rows),qT_fraction_max=max(r['qT_fraction_of_bulk'] for r in rows),maximum_closure=max(r['flux_closure_abs'] for r in rows),maximum_partition_residual=max(r['modal_sum_residual_abs'] for r in rows)) for k,rows in families.items()}
 need(max(r['flux_closure_abs'] for r in allrows)<1e-8,'Global accounting threshold failed')
 zero=coupling[0];need(abs(zero['Tbulk'])<1e-14 and abs(zero['T_qT_modal'])<1e-14,'Zero-coupling control failed')
 need(max(abs(r['T_qT_modal']-ref['T_qT_modal']) for r in units)<1e-14,'Dimensional re-encoding control failed')
 result=dict(schema='paper1.corrected-active-robustness.v1',status='passed_for_85_declared_finite_samples',scope='Finite sampled positivity and numerical ranges at declared resolutions; not a uniform accuracy bound, nonlinear extremum certificate or exact-lattice guarantee.',reference_at_robustness_N80=ref,PV_contract=contract,PV_contract_sha256=digest(a.pv_contract),kernel_run_sha256=digest(HERE/'CHECKPOINT_RUN.json'),stiffness_sample_audit_sha256=digest(HERE/'STIFFNESS_SAMPLE_AUDIT.json'),families=families,ranges=ranges,coupling_including_zero=coupling,dimensional_unit_covariance=units,weak_coupling_log_slope_sampled=slope,reference_is_not_best_sampled_profile=max(r['T_qT_modal'] for r in prows)>ref['T_qT_modal'],active_row_count=len(allrows),minimum_active_qT=min(r['T_qT_modal'] for r in allrows),maximum_closure=max(r['flux_closure_abs'] for r in allrows),maximum_partition_residual=max(r['modal_sum_residual_abs'] for r in allrows),producer_sha256=digest(Path(__file__)),changes_from_historical=['Basis-invariant nodal spectral weights.','N80 shells throughout the finite robustness sample.','Selected fresh two-ensemble PV contract for the common physical bulk.','Exact model-dependent spectral support for stiffness ratios; no physical eigenvalue clipping.'])
 (HERE/'ACTIVE85.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ('status','active_row_count','minimum_active_qT','maximum_closure','maximum_partition_residual','weak_coupling_log_slope_sampled','ranges')},indent=2))
if __name__=='__main__':main()
