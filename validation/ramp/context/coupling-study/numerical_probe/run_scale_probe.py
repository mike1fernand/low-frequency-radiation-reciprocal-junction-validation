"""Bounded eight-scale diagnostic of the unchanged P1 six-motif contact.

No new Green kernel, optimization or trajectory. Frozen source package is read
only. Output is kept here, and source modules/caches are copied into inputs/.
"""
from pathlib import Path
import dataclasses, hashlib, json, math, platform, shutil, sys, time
import numpy as np
import scipy

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PACKAGE = Path('/workspace/scratch/044aa7b17d80/p1_implemented_2026_11_08/revision/paper1-significance')
INPUT = HERE/'inputs'
OUTPUT = HERE/'SCALE_PROBE.json'
SCALES = [1., 1.5, 2., 2.5, 3., 4., 6., 10.]
OMEGA = 1.2
ETA_LEAD = 1e-12
REFERENCE = np.array([.1, .2, .2, .2, .2, .1])
H = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
J = lambda p, x: p.write_text(json.dumps(x, indent=2, allow_nan=False)+'\n')
assert not OUTPUT.exists(), 'Choose a fresh output; do not overwrite results.'

package_before = {str(p.relative_to(PACKAGE)):H(p) for p in PACKAGE.rglob('*') if p.is_file()}
INPUT.mkdir(parents=True, exist_ok=True)
copy_files = {
    'model_api_v1_1.py': 'recheck/purity_inputs/model_api_v1_1.py',
    'reduced_core.py': 'recheck/purity_inputs/reduced_core.py',
    'on_shell_N80.npz': 'recheck/purity_inputs/on_shell_N80.npz',
    'pv_mean_M192.npy': 'recheck/purity_inputs/pv_mean_M192.npy',
    'N80_INPUT_PROVENANCE.json': 'recheck/purity_inputs/INPUT_PROVENANCE.json',
    'REFERENCE_KERNEL.npz': 'tradeoff_validation/inputs/refined_surface_pv/selected/REFERENCE_KERNEL.npz',
    'REFERENCE_KERNEL.json': 'tradeoff_validation/inputs/refined_surface_pv/selected/REFERENCE_KERNEL.json',
}
for dst, src in copy_files.items():
    shutil.copyfile(PACKAGE/src, INPUT/dst)
assert H(INPUT/'model_api_v1_1.py') == 'aa4128940dc7d39682af53c87cf04ae2532de4af880d25393ef4a4fe0b067980'
assert H(INPUT/'reduced_core.py') == '0c6967e11c84d6c2fffeadc250c9d02f39374f30b1eade90d30660aa7144ed57'
assert H(INPUT/'REFERENCE_KERNEL.npz') == '82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
prov80 = json.loads((INPUT/'N80_INPUT_PROVENANCE.json').read_text())
for item in prov80['files']:
    assert H(INPUT/item['path'])==item['sha256']
meta = json.loads((INPUT/'REFERENCE_KERNEL.json').read_text())
assert meta['omega']==OMEGA and meta['energy']==OMEGA**2
assert meta['surface_grid']==192 and meta['PV_selection']['Jackson_order']==256
assert meta['PV_selection']['seeds']==[20260920,20260921]
assert meta['PV_selection']['sample_count_per_seed']==2**23
assert meta['profile']==REFERENCE.tolist() and meta['matrix_shape']==[20,20]
sys.path.insert(0, str(INPUT))
from model_api_v1_1 import ModelParameters
from reduced_core import scatter_from_reduced_green, analytic_q_lead_self_energy
p=ModelParameters(**meta['model_parameters']);p.validate()
assert dataclasses.asdict(p)==dataclasses.asdict(ModelParameters())

with np.load(INPUT/'REFERENCE_KERNEL.npz',allow_pickle=False) as z:
    refined={k:z[k].copy() for k in z.files}
with np.load(INPUT/'on_shell_N80.npz',allow_pickle=False) as z:
    older={k:z[k].copy() for k in z.files}
pv80=np.load(INPUT/'pv_mean_M192.npy',allow_pickle=False)
sl=np.s_[:12,:12]
G192=refined['G0'][sl]
BT192=refined['Im_branch_0'][sl]+refined['Im_branch_1'][sl]
BL192=refined['Im_branch_2'][sl]
BT80=older['Im_branch_0'][sl]+older['Im_branch_1'][sl]
BL80=older['Im_branch_2'][sl]
G80=pv80[sl]+1j*(BT80+BL80)
assert np.max(np.abs(G192.imag-BT192-BL192))<1e-12
assert np.max(np.abs(BT80+BL80-older['Im_total'][sl]))<1e-12

predeclare={
    'date':'2026-11-08', 'scope':'Finite predeclared scale continuation of the actual complete paired-spring contact at fixed omega1.2, six motifs, unchanged bare bulk and line parameters.',
    'scales':SCALES,'reference_g':REFERENCE.tolist(),'parameter_cache_states':16,
    'new_green_evaluations':0,'new_optimization':False,
    'primary_cache':'N192 spectral / M256 PV selected kernel',
    'sensitivity_cache':'Historical N80 spectral / M192 PV; NOT the published N80/M256 family.',
    'input_sha256':{x:H(INPUT/x) for x in copy_files},
    'finish_condition':'Eight scales on each of the two existing caches, after accepted N192 reference regression, with vectors/residuals saved; no adaptive continuation or search expansion.',
}
J(HERE/'PREDECLARATION.json',predeclare)

rows=[];saved={};start=time.perf_counter()
for cache,G,BT,BL in [('N192_M256',G192,BT192,BL192),('N80_M192',G80,BT80,BL80)]:
    assert max(np.max(np.abs(M-M.conj().T)) for M in [BT,BL])<1e-12
    assert min(np.linalg.eigvalsh(M).min() for M in [BT,BL])>-1e-12
    for index,scale in enumerate(SCALES):
        g=scale*REFERENCE
        result,internal=scatter_from_reduced_green(OMEGA,g,p,G,eta_lead=ETA_LEAD,return_internal=True)
        q=internal['q'];t=internal['t'];gamma=internal['gamma_in']
        tqT=float(2*gamma*np.vdot(t,BT@t).real)
        tqL=float(2*gamma*np.vdot(t,BL@t).real)
        nq=len(q);nb=len(t);z=complex(OMEGA,ETA_LEAD)
        Kq=np.diag(np.full(nq,2*p.k_spine))
        for j in range(nq-1):Kq[j,j+1]=Kq[j+1,j]=-p.k_spine
        sigma_lead=analytic_q_lead_self_energy(z,p)
        leads=np.zeros((nq,nq),complex);leads[0,0]=leads[-1,-1]=sigma_lead
        Aq0=p.m_spine*z*z*np.eye(nq)-Kq-internal['Kself']-leads
        A=Aq0-internal['Sigma_bulk'];rhs=np.zeros(nq,complex);rhs[0]=1
        block=np.block([[Aq0,internal['V'].conj().T@G],[-internal['V'],np.eye(nb)+internal['U']@G]])
        full_rhs=np.r_[rhs,np.zeros(nb,complex)]
        full_x=np.r_[q,t]
        state={
            'cache':cache,'scale':scale,'g':g.tolist(),'omega':OMEGA,
            'physical_total_spring_stiffness':float(8*g.sum()),'active_springs':48,
            'Rq':result.Rq,'Tline':result.Tq,'TqT':tqT,'TqL':tqL,
            'Tbulk':tqT+tqL,'bulk_qT_share':tqT/(tqT+tqL),
            'four_channel_closure':abs(result.Rq+result.Tq+tqT+tqL-1),
            'source_vs_selfenergy_bulk_residual':abs(tqT+tqL-result.Tbulk),
            'line_system_relative_residual':float(np.linalg.norm(A@q-rhs)/(np.linalg.norm(A)*np.linalg.norm(q)+np.linalg.norm(rhs))),
            'full_block_relative_residual':float(np.linalg.norm(block@full_x-full_rhs)/(np.linalg.norm(block)*np.linalg.norm(full_x)+np.linalg.norm(full_rhs))),
            'line_system_2norm_condition':float(np.linalg.cond(A)),
            'full_block_2norm_condition':float(np.linalg.cond(block)),
            'bulk_loading_min_eigenvalue':result.min_gamma_bulk_eigenvalue,
            'gamma_in':gamma,'P_inc_unit_displacement':OMEGA*gamma/4,
            'solver_result':dataclasses.asdict(result),
        }
        assert state['four_channel_closure']<1e-9
        assert state['source_vs_selfenergy_bulk_residual']<1e-12
        assert state['line_system_relative_residual']<1e-12
        assert state['full_block_relative_residual']<1e-12
        assert min(result.Rq,result.Tq,tqT,tqL)>-1e-12
        assert result.min_gamma_bulk_eigenvalue>-1e-9
        if cache=='N192_M256' and scale==1:
            ref=meta['scattering']
            state['reference_regression_max_abs']=max(abs(state[k]-ref[v]) for k,v in [('Rq','Rq'),('Tline','Tline'),('Tbulk','Tbulk'),('TqT','qT'),('TqL','qL')])
            assert state['reference_regression_max_abs']<1e-11
        prefix=f'{cache}_s{index}'
        for key,value in {'q':q,'t':t,'A':A,'rhs':rhs,'full_block':block,'full_rhs':full_rhs,**internal}.items():
            saved[prefix+'_'+key]=np.asarray(value)
        rows.append(state)

np.savez_compressed(HERE/'INTERNAL_STATES.npz',**saved)
high=[r for r in rows if r['cache']=='N192_M256']
low=[r for r in rows if r['cache']=='N80_M192']
comparison=[{'scale':a['scale'],'N192_M256_qT':a['TqT'],'N80_M192_qT':b['TqT'],'qT_absolute_difference':abs(a['TqT']-b['TqT']),'bulk_share_absolute_difference':abs(a['bulk_qT_share']-b['bulk_qT_share'])} for a,b in zip(high,low)]
package_after={str(p.relative_to(PACKAGE)):H(p) for p in PACKAGE.rglob('*') if p.is_file()}
assert package_before==package_after,'Source package changed'
result={
    'status':'PASS','date':'2026-11-08','model_parameters':dataclasses.asdict(p),
    'eta_lead':ETA_LEAD,'basis':'X1,R1,...,X6,R6; leading12 of cached20',
    'predeclared_scales':SCALES,'states':rows,'cache_comparison':comparison,
    'best_sampled_refined':max(high,key=lambda x:x['TqT']),
    'same_model_parameter_cache_states_executed':len(rows),
    'new_green_evaluations':0,'optimization_performed':False,
    'input_and_script_sha256':{**{x:H(INPUT/x) for x in copy_files},'run_scale_probe.py':H(Path(__file__))},
    'internal_states_sha256':H(HERE/'INTERNAL_STATES.npz'),
    'source_package_unchanged':True,'source_package_payloads_checked':len(package_before),
    'accuracy_scope':'Approximate selected-kernel point results. Cache disagreement is numerical sensitivity, not a rigorous error enclosure. Flux closure/linear residuals establish internal consistency, not Green-kernel accuracy. Historical N80 cache uses M192 PV, not the published M256 family. No packet, robustness, bandwidth or global-optimum claim.',
    'elapsed_seconds':time.perf_counter()-start,
    'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
}
J(OUTPUT,result)
print(json.dumps({'status':result['status'],'states':[{k:r[k] for k in ['cache','scale','Rq','Tline','TqT','TqL','bulk_qT_share','four_channel_closure']} for r in rows],'best_sampled_refined':result['best_sampled_refined'],'states_executed':len(rows)},indent=2))
