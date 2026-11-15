"""Frozen Stage3 comparison: exactly two responses, no Green-kernel generation.
Run only from a fresh exported copy to reproduce; never overwrite frozen results.
"""
from pathlib import Path
import dataclasses, hashlib, json, math, platform, sys, time
import numpy as np
import scipy

ROOT=Path(__file__).resolve().parent
INPUT=ROOT/'inputs/refined_surface_pv/selected'
MODULES=ROOT.parent/'recheck/purity_inputs'
EXPECTED={
 'model_api_v1_1.py':'aa4128940dc7d39682af53c87cf04ae2532de4af880d25393ef4a4fe0b067980',
 'reduced_core.py':'0c6967e11c84d6c2fffeadc250c9d02f39374f30b1eade90d30660aa7144ed57',
}
KERNEL_SHA='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
for name,digest in EXPECTED.items():assert sha(MODULES/name)==digest,name
assert sha(INPUT/'REFERENCE_KERNEL.npz')==KERNEL_SHA
meta=json.loads((INPUT/'REFERENCE_KERNEL.json').read_text())
assert meta['omega']==1.2 and meta['energy']==1.44
assert meta['surface_grid']==192 and meta['PV_selection']['Jackson_order']==256
assert meta['PV_selection']['sample_count_per_seed']==2**23
assert meta['PV_selection']['seeds']==[20260920,20260921]
assert meta['matrix_shape']==[20,20] and meta['support_segments']==10
assert meta['profile']==[.1,.2,.2,.2,.2,.1]
sys.path.insert(0,str(MODULES))
from model_api_v1_1 import ModelParameters
from reduced_core import scatter_from_reduced_green
p=ModelParameters(**meta['model_parameters']);p.validate()
assert dataclasses.asdict(p)=={'a':1.,'m_bulk':1.,'m_spine':1.,'k_axis':1.,'k_diag':1.,'k_spine':1.69,'k_couple':.2}
with np.load(INPUT/'REFERENCE_KERNEL.npz',allow_pickle=False) as z:
 arrays={k:z[k].copy() for k in z.files}
for k,a in arrays.items():assert a.shape==(20,20) and np.isfinite(a).all(),k
G=arrays['G0'][:12,:12].copy()
B=[arrays['Im_branch_'+str(j)][:12,:12].copy() for j in range(3)]
BT=B[0]+B[1];BL=B[2];BS=BT+BL
MATRIX_TOL=1e-12;REGRESSION_TOL=1e-11;FLUX_TOL=1e-9
matrix={
 'G_reciprocity_max_abs':float(np.max(np.abs(G-G.T))),
 'PV_hermiticity_max_abs':float(np.max(np.abs(G.real-G.real.T))),
 'branch_hermiticity_max_abs':max(float(np.max(np.abs(b-b.conj().T))) for b in B),
 'spectral_branch_sum_max_abs':float(np.max(np.abs(BS-arrays['Im_total'][:12,:12]))),
 'G_imaginary_spectral_sum_max_abs':float(np.max(np.abs(G.imag-BS))),
 'qT_key_sum_max_abs':float(np.max(np.abs(BT-arrays['Im_qT'][:12,:12]))),
 'qL_key_max_abs':float(np.max(np.abs(BL-arrays['Im_qL'][:12,:12]))),
 'PV_key_max_abs':float(np.max(np.abs(G.real-arrays['PV'][:12,:12]))),
 'branch_min_eigenvalues':[float(np.linalg.eigvalsh((b+b.conj().T)/2).min()) for b in B],
}
assert all(v<=MATRIX_TOL for k,v in matrix.items() if k!='branch_min_eigenvalues'),matrix
assert min(matrix['branch_min_eigenvalues'])>=-MATRIX_TOL,matrix
profiles=[('reference',[.1,.2,.2,.2,.2,.1]),('endpoint_only',[.5,0,0,0,0,.5])]
for _,g in profiles:assert math.isclose(sum(g),1,rel_tol=0,abs_tol=1e-15)
receipt_path=ROOT/'TWO_STATE_EXECUTION.json';output_path=ROOT/'TWO_PROFILE_RESULT.json'
assert not output_path.exists(),'Frozen result exists; do not spend another response budget.'
receipt={'status':'started','maximum_coupled_states':2,'coupled_states_executed':0,'new_kernel_evaluations':0,'profiles':profiles,'source_sha256':sha(Path(__file__))}
with receipt_path.open('x') as f:json.dump(receipt,f,indent=2)
states=[];t0=time.perf_counter()
for name,g in profiles:
 # Each call below is one coupled state; all matrix-only and power checks reuse it.
 result,inside=scatter_from_reduced_green(1.2,g,p,G,eta_lead=1e-12,return_internal=True)
 receipt['coupled_states_executed']+=1
 receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
 gamma=inside['gamma_in'];f=inside['t'];q=inside['q']
 tqT=float(2*gamma*np.real(f.conj()@BT@f))
 tqL=float(2*gamma*np.real(f.conj()@BL@f));tb=tqT+tqL
 state={'name':name,'g':g,'sum_g':sum(g),'physical_total_spring_stiffness':8*sum(g),
        'active_motifs':sum(x>0 for x in g),'active_springs':8*sum(x>0 for x in g),
        'peak_g':max(g),'Rq':result.Rq,'Tline':result.Tq,'TqT':tqT,'TqL':tqL,'Tbulk':tb,
        'bulk_qT_share':tqT/tb,'flux_closure_abs':abs(result.Rq+result.Tq+tb-1),
        'source_vs_selfenergy_bulk_residual':abs(tb-result.Tbulk),
        'gamma_in':gamma,'P_inc_unit_displacement':1.2*gamma/4,
        'solver_result':dataclasses.asdict(result)}
 # gamma convention and unit amplitude reproduce P=k*omega*sin(qlead)/2.
 lead_k=2*math.asin(1.2/(2*math.sqrt(p.k_spine/p.m_spine)))
 state['unit_incident_power_vs_dispersion_abs']=abs(state['P_inc_unit_displacement']-.5*p.k_spine*1.2*math.sin(lead_k))
 assert state['flux_closure_abs']<FLUX_TOL,state
 assert state['source_vs_selfenergy_bulk_residual']<MATRIX_TOL,state
 assert min(result.Rq,result.Tq,tqT,tqL)>-MATRIX_TOL,state
 if name=='reference':
  ref=meta['scattering']
  state['reference_regression_abs']={k:abs(state[k]-ref[v]) for k,v in [('Rq','Rq'),('Tline','Tline'),('Tbulk','Tbulk'),('TqT','qT'),('TqL','qL')]}
  assert max(state['reference_regression_abs'].values())<REGRESSION_TOL,state
 states.append(state)
ref,end=states
result={'status':'PASS','omega':1.2,'method':'existing N192 basis-invariant marching-surface spectrum and M256 mean two scrambled Sobol23 PV; no new quadrature',
        'source_basis':'X1,R1,...,X6,R6; leading12 of archived20 for ten available segments',
        'model_parameters':dataclasses.asdict(p),'kernel_sha256':KERNEL_SHA,'kernel_metadata_sha256':sha(INPUT/'REFERENCE_KERNEL.json'),
        'module_sha256':EXPECTED,'tolerances':{'matrix_absolute':MATRIX_TOL,'reference_absolute':REGRESSION_TOL,'flux_absolute':FLUX_TOL},
        'matrix_checks':matrix,'states':states,
        'differences':{'TqT_endpoint_minus_reference':end['TqT']-ref['TqT'],'bulk_qT_share_endpoint_minus_reference':end['bulk_qT_share']-ref['bulk_qT_share'],
                       'TqT_relative_increase':end['TqT']/ref['TqT']-1,'TqL_relative_increase':end['TqL']/ref['TqL']-1},
        'direction_survives':end['TqT']>ref['TqT'] and end['bulk_qT_share']<ref['bulk_qT_share'],
        'accuracy_scope':'Approximate selected-kernel point results. Internal consistency and reference regression; no endpoint-specific exact-lattice error enclosure or statistical95percent confidence.',
        'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'seconds':time.perf_counter()-t0}}
with output_path.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
receipt.update(status='completed',output_sha256=sha(output_path),elapsed_seconds=time.perf_counter()-t0)
receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'PASS','states':states,'differences':result['differences'],'coupled_states_executed':2,'new_kernel_evaluations':0},indent=2))
