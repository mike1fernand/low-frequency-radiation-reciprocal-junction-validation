#!/usr/bin/env python3
"""Fresh invariant-nodal spectral matrices, one frequency/grid checkpoint.
Strict wrappers reject hidden marching failures and omitted small-gradient cells.
PV in the summary is the OLD five-ensemble M192 mean, a fixed comparison baseline.
"""
from pathlib import Path
import sys,json,argparse,time,hashlib,resource
import numpy as np
from scipy.ndimage import map_coordinates
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(HERE/'src'))
import on_shell_basis as basis
import spectral_histogram as sh
from reduced_core import scatter_from_reduced_green
from model_api_v1_1 import ModelParameters
ap=argparse.ArgumentParser();ap.add_argument('--omega',required=True,type=float);ap.add_argument('--grid',required=True,type=int);ap.add_argument('--nseg',type=int,default=6);a=ap.parse_args()
if a.nseg not in [6,10]:raise ValueError('Only declared6/10 support allowed')
if not np.isfinite(a.omega) or not a.omega>0:raise ValueError('omega must be finite positive')
geo=[];original_mc=basis.marching_cubes

def strict_mc(volume,*args,**kwargs):
 if not np.all(np.isfinite(volume)):raise RuntimeError('Nonfinite spectral grid')
 if not volume.min()<0<volume.max():raise RuntimeError('This worker requires three resolved propagating branches')
 try:ans=original_mc(volume,*args,**kwargs)
 except Exception as exc:raise RuntimeError('Marching failure must not become zero branch') from exc
 verts,faces=ans[:2];cent=(verts[faces[:,0]]+verts[faces[:,1]]+verts[faces[:,2]])/3
 f=volume[:-1,:-1,:-1];N=f.shape[0];dk=2*np.pi/N
 grad=[(np.roll(f,-1,axis=i)-np.roll(f,1,axis=i))/(2*dk) for i in range(3)]
 gn=np.sqrt(sum(map_coordinates(basis._extend(v),cent.T,order=1,mode='wrap')**2 for v in grad))
 if not np.all(np.isfinite(gn)) or np.any(gn<=1e-11):raise RuntimeError('Would discard unresolved-gradient surface cells')
 geo.append({'triangles':len(faces),'min_gradient':float(gn.min()),'max_gradient':float(gn.max()),'skipped_gradient_cells':0})
 return ans
basis.marching_cubes=strict_mc
start=time.time();on=basis.integrate_on_shell(a.omega,ModelParameters(),nseg=a.nseg,grid=a.grid);seconds=time.time()-start
for matrix in [on.Im_total,*on.Im_branches,on.Im_longitudinal,on.Im_transverse]:
 if not np.all(np.isfinite(matrix)):raise RuntimeError('Nonfinite matrix')
 if np.linalg.eigvalsh(matrix).min() < -256*np.finfo(float).eps*max(1.,np.linalg.norm(matrix,2)):raise RuntimeError('Spectral PSD failure beyond roundoff')
oldraw=ROOT/'analysis/math_resolution/surface/inputs/raw'
histnames=['grid80x4','grid96x2','sobol20a','sobol20b','sobol20c']
pv=np.mean([sh.load_histogram(oldraw/(name+'.npz')).principal_value_jackson(a.omega*a.omega,192)[:2*a.nseg,:2*a.nseg] for name in histnames],axis=0)
res,x=scatter_from_reduced_green(a.omega,[.1,.2,.2,.2,.2,.1],ModelParameters(),pv[:12,:12]+1j*on.Im_total[:12,:12],return_internal=True)
bs=[float(2*x['gamma_in']*np.vdot(x['t'],M[:12,:12]@x['t']).real) for M in on.Im_branches]
tag=f'omega{a.omega:.15g}_N{a.grid}'+('' if a.nseg==6 else f'_nseg{a.nseg}');file=HERE/'results'/(tag+'.npz')
np.savez_compressed(file,Im_total=on.Im_total,Im0=on.Im_branches[0],Im1=on.Im_branches[1],Im2=on.Im_branches[2],Im_long=on.Im_longitudinal,Im_transverse=on.Im_transverse,PV_old_mean192=pv)
record={'schema':'paper1.surface-convergence.v1','omega':a.omega,'grid':a.grid,'spectral_support_segments':a.nseg,'model':{'a':1.,'m_bulk':1.,'m_spine':1.,'k_axis':1.,'k_diag':1.,'k_spine':1.69},'profile':[.1,.2,.2,.2,.2,.1], 'seconds_surface':seconds,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'geometry':geo,'raw_reciprocity_residual':on.symmetry_residual_before_projection,'branch_min_eigenvalues':[float(np.linalg.eigvalsh(M).min()) for M in on.Im_branches], 'fixed_old_PV_baseline':{'order':192,'ensembles':histnames,'Rq':res.Rq,'Tline':res.Tq,'Tbulk':res.Tbulk,'qT':sum(bs[:2]),'qL':bs[2],'closure':res.closure_error,'partition_residual':abs(sum(bs)-res.Tbulk)}, 'matrix_file':file.name,'matrix_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/'src/on_shell_basis.py').read_bytes()).hexdigest(),'scope':'Fresh basis-invariant spectral quadrature; coefficient baseline uses old PV; no true-error certificate.'}
(HERE/'results'/(tag+'.json')).write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
