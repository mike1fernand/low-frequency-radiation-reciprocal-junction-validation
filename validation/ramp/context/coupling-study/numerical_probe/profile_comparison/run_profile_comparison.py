"""Twelve predeclared positive/zero complete-motif profiles; two cached kernels.

This is a fixed finite comparison after the scale probe, not an optimizer.
"""
from pathlib import Path
import dataclasses, hashlib, json, math, sys, time
import numpy as np

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
INPUT=HERE.parent/'inputs'
PACKAGE=Path('/workspace/scratch/044aa7b17d80/p1_implemented_2026_11_08/revision/paper1-significance')
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p,d:p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
assert not (HERE/'PROFILE_COMPARISON.json').exists(),'Do not overwrite frozen outputs.'
TOTALS=[2.5,3.,4.]
SHAPES={
    'uniform':np.ones(6)/6,
    'endpoint_only':np.array([.5,0,0,0,0,.5]),
    'increasing_linear':np.arange(1,7,dtype=float)/21,
    'decreasing_linear':np.arange(6,0,-1,dtype=float)/21,
}
for w in SHAPES.values():assert np.all(w>=0) and abs(w.sum()-1)<1e-15
before={str(p.relative_to(PACKAGE)):H(p) for p in PACKAGE.rglob('*') if p.is_file()}
meta=json.loads((INPUT/'REFERENCE_KERNEL.json').read_text())
prior=json.loads((HERE.parent/'SCALE_PROBE.json').read_text())
for n,h in prior['input_and_script_sha256'].items():
    if n!='run_scale_probe.py':assert H(INPUT/n)==h
assert prior['same_model_parameter_cache_states_executed']==16
sys.path.insert(0,str(INPUT))
from model_api_v1_1 import ModelParameters
from reduced_core import scatter_from_reduced_green,analytic_q_lead_self_energy
p=ModelParameters(**meta['model_parameters']);p.validate()
with np.load(INPUT/'REFERENCE_KERNEL.npz',allow_pickle=False) as z:ref={k:z[k].copy() for k in z.files}
with np.load(INPUT/'on_shell_N80.npz',allow_pickle=False) as z:old={k:z[k].copy() for k in z.files}
sl=np.s_[:12,:12]
BT192=ref['Im_branch_0'][sl]+ref['Im_branch_1'][sl];BL192=ref['Im_branch_2'][sl]
BT80=old['Im_branch_0'][sl]+old['Im_branch_1'][sl];BL80=old['Im_branch_2'][sl]
G80=np.load(INPUT/'pv_mean_M192.npy',allow_pickle=False)[sl]+1j*(BT80+BL80)
caches=[('N192_M256',ref['G0'][sl],BT192,BL192),('N80_M192',G80,BT80,BL80)]
pre={
    'date':'2026-11-08','protocol_fixed_before_response_loop':True,
    'omega':1.2,'total_sum_g':TOTALS,'normalized_shapes':{k:w.tolist() for k,w in SHAPES.items()},
    'same_model_profiles':12,'parameter_cache_states':24,'previous_scale_probe_states':16,'combined_state_budget':40,
    'new_green_evaluations':0,'optimizer':False,'adaptive_expansion':False,
    'model':'Unchanged actual P1 complete paired-spring contact on six available motif positions; gi are actual per-link spring stiffnesses.',
    'cache_scope':'N192/M256 selected kernel plus historical N80/M192 sensitivity cache. N80/M192 is not the published N80/M256 family.',
    'comparison':'At each fixed sum_g compare four shapes against existing reference-taper scale from the prior batch. Across totals physical contact stiffness changes.',
    'active_count_warning':'Endpoint-only uses16activecontactsprings instead of48; available span/line/bulk retained. This is not a fixed active-count comparison.',
    'finish_condition':'Exactly12profiles/24parameter-cache states, full saved q/t systems and flux/residual checks; stop even if a boundary profile wins.',
    'input_hashes':{n:h for n,h in prior['input_and_script_sha256'].items() if n!='run_scale_probe.py'},
    'prior_scale_result_sha256':H(HERE.parent/'SCALE_PROBE.json'),
}
J(HERE/'PREDECLARATION.json',pre)
rows=[];saved={};start=time.perf_counter()
for cache,G,BT,BL in caches:
    assert max(np.max(np.abs(B-B.conj().T)) for B in [BT,BL])<1e-12
    assert min(np.linalg.eigvalsh(B).min() for B in [BT,BL])>-1e-12
    for total in TOTALS:
        for name,w in SHAPES.items():
            g=total*w
            response,x=scatter_from_reduced_green(1.2,g,p,G,eta_lead=1e-12,return_internal=True)
            q=x['q'];t=x['t'];gin=x['gamma_in'];nq=len(q);nb=len(t);z=complex(1.2,1e-12)
            qT=float(2*gin*np.vdot(t,BT@t).real);qL=float(2*gin*np.vdot(t,BL@t).real)
            K=np.diag(np.full(nq,2*p.k_spine))
            for i in range(nq-1):K[i,i+1]=K[i+1,i]=-p.k_spine
            leads=np.zeros((nq,nq),complex);leads[0,0]=leads[-1,-1]=analytic_q_lead_self_energy(z,p)
            A0=p.m_spine*z*z*np.eye(nq)-K-x['Kself']-leads;A=A0-x['Sigma_bulk']
            rhs=np.zeros(nq,complex);rhs[0]=1
            F=np.block([[A0,x['V'].conj().T@G],[-x['V'],np.eye(nb)+x['U']@G]])
            frhs=np.r_[rhs,np.zeros(nb,complex)];fx=np.r_[q,t]
            row={
                'cache':cache,'shape':name,'sum_g':float(g.sum()),'g':g.tolist(),
                'physical_total_spring_stiffness':float(8*g.sum()),
                'available_motifs':6,'active_motifs':int(np.count_nonzero(g)),'active_springs':int(8*np.count_nonzero(g)),
                'peak_g':float(g.max()),'Rq':response.Rq,'Tline':response.Tq,'TqT':qT,'TqL':qL,'Tbulk':qT+qL,
                'bulk_qT_share':qT/(qT+qL),'four_channel_closure':abs(response.Rq+response.Tq+qT+qL-1),
                'source_vs_selfenergy_residual':abs(qT+qL-response.Tbulk),
                'line_system_relative_residual':float(np.linalg.norm(A@q-rhs)/(np.linalg.norm(A)*np.linalg.norm(q)+np.linalg.norm(rhs))),
                'full_block_relative_residual':float(np.linalg.norm(F@fx-frhs)/(np.linalg.norm(F)*np.linalg.norm(fx)+np.linalg.norm(frhs))),
                'line_2norm_condition':float(np.linalg.cond(A)),'full_block_2norm_condition':float(np.linalg.cond(F)),
                'minimum_loading_gamma_eigenvalue':response.min_gamma_bulk_eigenvalue,
            }
            assert row['four_channel_closure']<1e-9
            assert row['source_vs_selfenergy_residual']<1e-12
            assert row['line_system_relative_residual']<1e-12 and row['full_block_relative_residual']<1e-12
            assert min(response.Rq,response.Tq,qT,qL)>-1e-12
            assert response.min_gamma_bulk_eigenvalue>-1e-9
            prefix=f'{cache}_G{total:g}_{name}'
            for key,val in {'A':A,'rhs':rhs,'full_block':F,'full_rhs':frhs,**x}.items():saved[prefix+'_'+key]=np.asarray(val)
            prior_ref=next(v for v in prior['states'] if v['cache']==cache and v['scale']==total)
            row['matching_total_reference_taper_qT']=prior_ref['TqT']
            row['qT_minus_matching_total_reference_taper']=qT-prior_ref['TqT']
            rows.append(row)
np.savez_compressed(HERE/'INTERNAL_STATES.npz',**saved)
assert before=={str(p.relative_to(PACKAGE)):H(p) for p in PACKAGE.rglob('*') if p.is_file()}
assert len(rows)==24
high=[r for r in rows if r['cache']=='N192_M256']
comparison=[]
for r in high:
    l=next(v for v in rows if v['cache']=='N80_M192' and v['shape']==r['shape'] and abs(v['sum_g']-r['sum_g'])<1e-14)
    comparison.append({'sum_g':r['sum_g'],'shape':r['shape'],'refined_qT':r['TqT'],'older_qT':l['TqT'],'qT_absolute_difference':abs(r['TqT']-l['TqT'])})
result={
    'status':'PASS','date':'2026-11-08','omega':1.2,'model_parameters':dataclasses.asdict(p),
    'states':rows,'cache_comparison':comparison,'best_sampled_refined':max(high,key=lambda r:r['TqT']),
    'additional_parameter_cache_states':24,'prior_scale_states':16,'combined_state_count':40,
    'new_green_evaluations':0,'optimization_performed':False,'source_package_unchanged':True,
    'source_payloads_checked':len(before),'script_sha256':H(Path(__file__)),
    'predeclaration_sha256':H(HERE/'PREDECLARATION.json'),'internal_states_sha256':H(HERE/'INTERNAL_STATES.npz'),
    'elapsed_seconds':time.perf_counter()-start,
    'accuracy_scope':'Approximate point predictions on existing kernels, no exact-lattice error enclosure. Cache differences are sensitivity, not confidence intervals. Flux/residual checks are internal consistency. No optimum, bandwidth, robustness, independent waveform test or scientific-priority conclusion.',
}
J(HERE/'PROFILE_COMPARISON.json',result)
print(json.dumps({'status':'PASS','states':[{k:r[k] for k in ['cache','sum_g','shape','TqT','Rq','Tline','TqL','bulk_qT_share']} for r in rows],'best_sampled_refined':result['best_sampled_refined'],'combined_states':40,'seconds':result['elapsed_seconds']},indent=2))
