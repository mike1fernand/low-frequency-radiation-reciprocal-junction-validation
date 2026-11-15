#!/usr/bin/env python3
"""One restartable corrected finite-robustness kernel at a time."""
from pathlib import Path
import argparse,hashlib,json,sys,time,resource
sys.dont_write_bytecode=True
import numpy as np
from scipy.ndimage import map_coordinates
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[2];sys.path.insert(0,str(HERE/'src'))
from model_api_v1_1 import ModelParameters
from reduced_core import scatter_from_reduced_green
import on_shell_basis as shell
import spectral_histogram as hist

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def source_hashes():return {str(p.relative_to(BASE)):digest(p) for p in [Path(__file__),*sorted((HERE/'src').glob('*.py'))]}
def save(path,obj):path.write_text(json.dumps(obj,indent=2)+'\n')
geometry=[];original_mc=shell.marching_cubes

def strict_mc(volume,*args,**kwargs):
 need(np.isfinite(volume).all(),'Nonfinite shell eigenvalue grid')
 need(volume.min()<0<volume.max(),'Expected propagating branch lacks a bracket: do not silently skip')
 try:result=original_mc(volume,*args,**kwargs)
 except Exception as exc:raise RuntimeError('Marching exception: do not replace branch by zero') from exc
 vertices,faces=result[:2];cent=vertices[faces].mean(axis=1);N=volume.shape[0]-1;dk=2*np.pi/N;f=volume[:-1,:-1,:-1]
 grad=[(np.roll(f,-1,axis=i)-np.roll(f,1,axis=i))/(2*dk) for i in range(3)]
 gn=np.sqrt(sum(map_coordinates(shell._extend(v),cent.T,order=1,mode='wrap')**2 for v in grad))
 area=.5*np.linalg.norm(np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]]),axis=1)*dk*dk
 need(np.isfinite(gn).all() and np.all(gn>1e-11),'Unresolved-gradient shell cells would be omitted')
 need(np.isfinite(area).all() and np.all(area>0),'Nonpositive or nonfinite shell triangle area')
 geometry.append(dict(triangles=len(faces),minimum_gradient=float(gn.min()),minimum_area=float(area.min()),omitted_gradient_cells=0,nonpositive_areas=0))
 return result
shell.marching_cubes=strict_mc
range_audit={};original_acc=hist._acc

def strict_acc(mass,vals,c,kx,edges,nseg,norm):
 need(np.isfinite(vals).all() and np.isfinite(c).all(),'Nonfinite spectral sample')
 tol=128*np.finfo(float).eps*max(1.,abs(edges[-1]))
 lo=float(vals.min());hi=float(vals.max());need(lo>=edges[0]-tol and hi<=edges[-1]+tol,'Physical spectrum outside histogram support')
 range_audit['minimum_sampled_eigenvalue']=min(range_audit.get('minimum_sampled_eigenvalue',lo),lo)
 range_audit['maximum_sampled_eigenvalue']=max(range_audit.get('maximum_sampled_eigenvalue',hi),hi)
 range_audit['sampled_modes']=range_audit.get('sampled_modes',0)+vals.size
 range_audit['roundoff_outside_support']=range_audit.get('roundoff_outside_support',0)+int(np.count_nonzero((vals<edges[0])|(vals>edges[-1])))
 range_audit['modes_above_historical_edge12']=range_audit.get('modes_above_historical_edge12',0)+int(np.count_nonzero(vals>12))
 return original_acc(mass,vals,c,kx,edges,nseg,norm)
hist._acc=strict_acc

def arrays(on):return dict(Im_total=on.Im_total,Im_branch_0=on.Im_branches[0],Im_branch_1=on.Im_branches[1],Im_branch_2=on.Im_branches[2],Im_longitudinal=on.Im_longitudinal,Im_transverse=on.Im_transverse)
def check_shell(on):
 mins=[]
 for M in (on.Im_total,*on.Im_branches):
  need(np.isfinite(M).all() and np.array_equal(M,M.T),'Nonfinite or nonsymmetric shell matrix')
  mineig=float(np.linalg.eigvalsh(M).min());need(mineig>=-256*np.finfo(float).eps*max(1.,float(np.linalg.norm(M,2))),'Shell matrix PSD error exceeds roundoff')
  mins.append(mineig)
 need(len(geometry)==3,'Expected exactly three nonempty branch geometries')
 need(np.allclose(sum(on.Im_branches),on.Im_total,rtol=0,atol=1e-14),'Branch matrices do not sum to total')
 return mins

def calc(w,g,p,pv,bs):
 res,x=scatter_from_reduced_green(w,g,p,pv+1j*sum(bs),return_internal=True)
 powers=[float(2*x['gamma_in']*np.vdot(x['t'],M@x['t']).real) for M in bs]
 row=dict(Rq=res.Rq,Tq=res.Tq,Tbulk=res.Tbulk,T_qT_modal=sum(powers[:2]),T_qL_modal=powers[2],qT_fraction_of_bulk=sum(powers[:2])/res.Tbulk if res.Tbulk else 0.,flux_closure_abs=res.closure_error,modal_sum_residual_abs=abs(sum(powers)-res.Tbulk),min_bulk_linewidth_eigenvalue=res.min_gamma_bulk_eigenvalue)
 need(all(np.isfinite(x) for x in row.values()),'Nonfinite response diagnostic')
 need(row['flux_closure_abs']<1e-8 and row['modal_sum_residual_abs']<1e-10,'Scattering accounting check failed')
 need(min(row[k] for k in ('Rq','Tq','Tbulk','T_qT_modal','T_qL_modal'))>=-1e-12,'Negative passive response')
 return row

def run(a):
 t0=time.perf_counter();sources=source_hashes()
 if a.family=='shell':
  p=ModelParameters();tag=f'omega{a.omega:.6f}_N{a.grid}_nseg{a.nseg}';folder=HERE/'kernels';w=a.omega
 else:
  p=ModelParameters(k_diag=a.kdiag);tag=f'kdiag{a.kdiag:g}_N{a.grid}';folder=HERE/'stiffness';w=1.2;a.nseg=6
 output=folder/(tag+'.json');matrix=folder/(tag+'.npz')
 need(not output.exists() and not matrix.exists(),'Checkpoint already exists; inspect it before rerunning')
 on=shell.integrate_on_shell(w,p,nseg=a.nseg,grid=a.grid);mins=check_shell(on);payload=arrays(on);sampling=[]
 record=dict(schema='paper1.corrected-robustness-kernel.v1',family=a.family,omega=w,nseg=a.nseg,surface_grid=a.grid,model_parameters=p.__dict__,geometry=geometry,shell_min_eigenvalues=mins,raw_reciprocity_residual=on.symmetry_residual_before_projection,source_sha256=sources)
 if a.family=='stiffness':
  upper=(4*p.k_axis+8*p.k_diag)/p.m_bulk;edges=(0.,upper)
  sampling_start=time.perf_counter()
  h1=hist.sample_shifted_grids(p,nseg=6,grid=64,shifts=((.17,.39,.61),(.67,.89,.11)),nbins=6144,energy_range=edges)
  h2=hist.sample_sobol(p,nseg=6,power=19,seed=97531,nbins=6144,energy_range=edges)
  pv=.5*(h1.principal_value_jackson(w*w,160)+h2.principal_value_jackson(w*w,160))
  payload['PV']=pv
  for name,h in [('grid64x2',h1),('sobol19',h2)]:
   hp=folder/(tag+'_'+name+'_hist.npz');np.savez_compressed(hp,edges=h.edges,mass=h.mass,sample_count=h.sample_count,sampler=h.sampler,metadata_json=json.dumps(h.metadata))
   sampling.append(dict(file=hp.name,sha256=digest(hp),sample_count=h.sample_count,sampler=h.sampler,metadata=h.metadata))
  record.update(k_diag_over_k_axis=a.kdiag,PV_order=160,histogram_energy_range=list(edges),spectral_range_audit=range_audit,histogram_inputs=sampling,sampling_seconds=time.perf_counter()-sampling_start,coefficient=calc(w,[.1,.2,.2,.2,.2,.1],p,pv,on.Im_branches))
 np.savez_compressed(matrix,**payload)
 record.update(status='completed_for_declared_finite_sample_scope',matrix_file=matrix.name,matrix_sha256=digest(matrix),seconds=time.perf_counter()-t0,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='Corrected finite-resolution numerical sample; no uniform accuracy or exact-lattice error bound.')
 need(sources==source_hashes(),'Sources changed during checkpoint')
 save(output,record);print(json.dumps(dict(status=record['status'],family=a.family,output=str(output.relative_to(BASE)),seconds=record['seconds']),indent=2))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--family',choices=['shell','stiffness'],required=True);ap.add_argument('--omega',type=float,default=1.2);ap.add_argument('--grid',type=int,default=80);ap.add_argument('--nseg',type=int,default=6);ap.add_argument('--kdiag',type=float,default=1.);a=ap.parse_args()
 need(np.isfinite(a.omega) and a.omega>0 and a.grid>=24 and a.nseg>0 and np.isfinite(a.kdiag) and a.kdiag>0,'Invalid physical/numerical input')
 try:run(a)
 except Exception as exc:
  name=f'FAILED_{a.family}_{a.omega:g}_{a.kdiag:g}_{time.time_ns()}.json'
  save(HERE/name,dict(status='failed',args=vars(a),error_type=type(exc).__name__,error=str(exc),geometry_completed=geometry,source_sha256=source_hashes()))
  raise
if __name__=='__main__':main()
