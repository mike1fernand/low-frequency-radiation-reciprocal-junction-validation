"""Staged fresh quadrature at one frozen P1 candidate; no parameter search.

Pilot: surface N48 and two new scrambled Sobol2^18 PV histograms.
Refinement: surface N80,N128; two new Sobol2^20 PV histograms; M192,256,320.
The original source files and historical kernels are never modified.
"""
from pathlib import Path
import argparse, dataclasses, hashlib, json, math, os, platform, resource, shutil, sys, time
import numpy as np
from scipy.ndimage import map_coordinates

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PACKAGE=Path('/workspace/scratch/044aa7b17d80/p1_implemented_2026_11_08/revision/paper1-significance')
INPUT=HERE/'inputs';RESULTS=HERE/'results';INPUT.mkdir(parents=True,exist_ok=True);RESULTS.mkdir(exist_ok=True)
ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['pilot','refine'],required=True);args=ap.parse_args()
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p,d:p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
sourcefiles={
 'model_api_v1_1.py':'recheck/purity_inputs/model_api_v1_1.py',
 'reduced_core.py':'recheck/purity_inputs/reduced_core.py',
 'on_shell_basis.py':'tradeoff_validation/inputs/refined_surface_pv/surface/src/on_shell_basis.py',
 'spectral_histogram.py':'tradeoff_validation/inputs/refined_surface_pv/surface/src/spectral_histogram.py',
 'REFERENCE_KERNEL.npz':'tradeoff_validation/inputs/refined_surface_pv/selected/REFERENCE_KERNEL.npz',
 'REFERENCE_KERNEL.json':'tradeoff_validation/inputs/refined_surface_pv/selected/REFERENCE_KERNEL.json',
}
before={rel:H(PACKAGE/rel) for rel in sourcefiles.values()}
for name,rel in sourcefiles.items():
    target=INPUT/name
    if target.exists(): assert H(target)==H(PACKAGE/rel)
    else:shutil.copyfile(PACKAGE/rel,target)
assert H(INPUT/'REFERENCE_KERNEL.npz')=='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
sys.path.insert(0,str(INPUT))
import on_shell_basis as basis
import spectral_histogram as sh
from reduced_core import scatter_from_reduced_green
from model_api_v1_1 import ModelParameters
p=ModelParameters();omega=1.2;g=np.arange(1,7,dtype=float)/7
grids=[48] if args.phase=='pilot' else [80,128]
power=18 if args.phase=='pilot' else 20
seeds=[20261108,20261109];orders=[192,256,320]
assert np.isclose(g.sum(),3)
plan={'status':'predeclared','phase':args.phase,'omega':omega,'g':g.tolist(),
      'surface_grids':grids,'PV_power':power,'PV_seeds':seeds,'PV_orders':orders,
      'nseg':6,'nbins':16384,'new_parameter_search':False,
      'phase_finish':'Only these declared matrices and their candidate responses; no automatic next refinement.',
      'input_sha256':{name:H(INPUT/name) for name in sourcefiles},'script_sha256':H(Path(__file__))}
pre=HERE/(args.phase+'_PREDECLARATION.json');assert not pre.exists();J(pre,plan)

geom=[];original_mc=basis.marching_cubes
def strict_mc(volume,*a,**kw):
    assert np.isfinite(volume).all() and volume.min()<0<volume.max()
    ans=original_mc(volume,*a,**kw);verts,faces=ans[:2]
    cent=(verts[faces[:,0]]+verts[faces[:,1]]+verts[faces[:,2]])/3
    f=volume[:-1,:-1,:-1];N=f.shape[0];dk=2*np.pi/N
    grad=[(np.roll(f,-1,axis=i)-np.roll(f,1,axis=i))/(2*dk) for i in range(3)]
    gn=np.sqrt(sum(map_coordinates(basis._extend(v),cent.T,order=1,mode='wrap')**2 for v in grad))
    assert np.isfinite(gn).all() and np.all(gn>1e-11)
    geom.append({'triangles':len(faces),'minimum_gradient':float(gn.min()),'maximum_gradient':float(gn.max()),'skipped_cells':0})
    return ans
basis.marching_cubes=strict_mc
start=time.perf_counter();spectra={};timings=[]
for N in grids:
    t=time.perf_counter();geom.clear();on=basis.integrate_on_shell(omega,p,nseg=6,grid=N)
    matrices={'Im_total':on.Im_total,**{'Im_branch_'+str(i):M for i,M in enumerate(on.Im_branches)}}
    assert all(np.isfinite(M).all() for M in matrices.values())
    mins=[float(np.linalg.eigvalsh(M).min()) for M in matrices.values()]
    assert min(mins)>-1e-12
    file=RESULTS/f'surface_N{N}.npz';assert not file.exists();np.savez_compressed(file,**matrices)
    rec={'surface_grid':N,'seconds':time.perf_counter()-t,'geometry':geom.copy(),
         'minimum_eigenvalues':mins,'raw_reciprocity_residual':on.symmetry_residual_before_projection,
         'matrix_sha256':H(file)}
    J(RESULTS/f'surface_N{N}.json',rec);timings.append(rec);spectra[N]=matrices
    print(json.dumps({'stage':'surface','N':N,'seconds':rec['seconds']}),flush=True)

hists=[];hist_records=[]
for seed in seeds:
    t=time.perf_counter();hist=sh.sample_sobol(p,nseg=6,power=power,seed=seed,nbins=16384)
    file=RESULTS/f'PV_hist_p{power}_seed{seed}.npz';assert not file.exists();sh.save_histogram(file,hist)
    rec={'seed':seed,'power':power,'sample_count':hist.sample_count,'seconds':time.perf_counter()-t,'sha256':H(file)}
    hists.append(hist);hist_records.append(rec)
    print(json.dumps({'stage':'PV_histogram',**rec}),flush=True)

pvs={};states=[];saved={}
for M in orders:
    pcs=[h.principal_value_jackson(omega**2,M) for h in hists]
    for label,PV in [('seedA',pcs[0]),('seedB',pcs[1]),('mean',(pcs[0]+pcs[1])/2)]:
        assert np.isfinite(PV).all() and np.max(abs(PV-PV.T))<1e-12
        pvs[(M,label)]=PV
        for N,mats in spectra.items():
            G=PV+1j*mats['Im_total']
            result,inside=scatter_from_reduced_green(omega,g,p,G,eta_lead=1e-12,return_internal=True)
            BT=mats['Im_branch_0']+mats['Im_branch_1'];BL=mats['Im_branch_2'];t=inside['t'];q=inside['q'];gamma=inside['gamma_in']
            tqT=float(2*gamma*np.vdot(t,BT@t).real);tqL=float(2*gamma*np.vdot(t,BL@t).real)
            rec={'surface_N':N,'PV_power':power,'PV_order':M,'PV_selection':label,
                 'Rq':result.Rq,'Tline':result.Tq,'TqT':tqT,'TqL':tqL,
                 'Tbulk':tqT+tqL,'bulk_qT_share':tqT/(tqT+tqL),
                 'closure':abs(result.Rq+result.Tq+tqT+tqL-1),
                 'spectral_selfenergy_residual':abs(tqT+tqL-result.Tbulk)}
            assert min(rec[k] for k in ['Rq','Tline','TqT','TqL'])>=-1e-12
            assert rec['closure']<1e-8 and rec['spectral_selfenergy_residual']<1e-10
            key=f'N{N}_M{M}_{label}'
            for name,v in [('PV',PV),('G',G),('BT',BT),('BL',BL),('q',q),('t',t),*inside.items()]:
                if isinstance(v,(np.ndarray,int,float,complex,np.number)):saved[key+'__'+name]=np.asarray(v)
            states.append(rec)
with np.load(INPUT/'REFERENCE_KERNEL.npz') as z:
    Gref=z['G0'][:12,:12];BTref=(z['Im_branch_0']+z['Im_branch_1'])[:12,:12];BLref=z['Im_branch_2'][:12,:12]
rr,ii=scatter_from_reduced_green(omega,g,p,Gref,eta_lead=1e-12,return_internal=True)
reference=float(2*ii['gamma_in']*np.vdot(ii['t'],BTref@ii['t']).real)
assert abs(reference-.7358742883948081)<1e-12
statefile=RESULTS/(args.phase+'_INTERNAL_STATES.npz');np.savez_compressed(statefile,**saved)
assert all(H(PACKAGE/rel)==digest for rel,digest in before.items())
record={'status':'PASS','phase':args.phase,'omega':omega,'g':g.tolist(),
        'model_parameters':dataclasses.asdict(p),'spectral_timings':timings,'PV_histograms':hist_records,
        'states':states,'coupled_states':len(states),'new_surface_matrices':len(grids),'new_PV_histograms':2,
        'cached_selected_reference_qT':reference,'seconds':time.perf_counter()-start,
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'source_files_unchanged':True,'input_sha256':plan['input_sha256'],
        'internal_states_sha256':H(statefile),'runtime':{'python':platform.python_version(),'numpy':np.__version__},
        'accuracy_scope':'Fresh independent samples/grids under the archived same quadrature algorithms; convergence/sensitivity evidence, not a rigorous exact-lattice enclosure or independent-method proof.'}
J(RESULTS/(args.phase+'_RESULT.json'),record)
print(json.dumps({'status':'PASS','phase':args.phase,'seconds':record['seconds'],'states':states},indent=2),flush=True)
