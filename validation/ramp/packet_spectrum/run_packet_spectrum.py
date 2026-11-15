"""Recompute the new physical-profile finite-DFT packet prediction.

No curve interpolation, new kernel generation or time-domain data are used.
Requires Python3, NumPy and SciPy. Refuse to overwrite frozen outputs.
"""
from pathlib import Path
import dataclasses,hashlib,json,math,platform,sys,time
import numpy as np

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;INPUT=HERE/'inputs';OUT=HERE/'PACKET_SPECTRUM.json'
assert not OUT.exists(),'Use a fresh exported copy, preserving frozen results.'
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p,x:p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
binding=json.loads((HERE/'INPUT_BINDING.json').read_text())
for n,h in binding['copied_sources'].items():assert H(INPUT/n)==h,n
sys.path.insert(0,str(INPUT))
from model_api_v1_1 import ModelParameters
from reduced_core import scatter_from_reduced_green
from packet_definition import packet_spectrum
config=json.loads((INPUT/'preparation.json').read_text());contract=json.loads((INPUT/'FREQUENCY_CONTRACT.json').read_text())
spec,omega,W,inside,q,v=packet_spectrum(config)
n=config['n'];k=2*np.pi*np.fft.fftfreq(n,d=config['a']);prop=(np.abs(k)*config['a']>0)&(np.abs(k)*config['a']<np.pi)
retained=inside&prop
ids=[i for i in range(1,n//2) if retained[i]]
assert len(ids)==len(contract['central_nodes'])==14
assert abs(float(W[retained].sum()/W.sum())-contract['inside_coverage'])<2e-15
for i,node in zip(ids,contract['central_nodes']):
    assert i==node['positive_k_index'] and n-i==node['negative_k_index']
    assert abs(omega[i]-node['omega'])<2e-13
    assert abs((W[i]+W[-i])/W.sum()-node['weight_over_full_periodic_energy'])<2e-17
assert spec['parseval_relative_residual']<1e-14
g_ref=np.array([.1,.2,.2,.2,.2,.1]);g_new=np.arange(1,7,dtype=float)/7
p=ModelParameters();p.validate();assert dataclasses.asdict(p)=={'a':1.,'m_bulk':1.,'m_spine':1.,'k_axis':1.,'k_diag':1.,'k_spine':1.69,'k_couple':.2}
pre={'date':'2026-11-08','new_profile_g':g_new.tolist(),'reference_g':g_ref.tolist(),'new_profile_states':70,'reference_regression_states':70,
 'families':[f['label'] for f in binding['families']],'frequency_nodes':[r['omega'] for r in contract['central_nodes']],
 'normalization':'Propagating-only numerator; full signed120-bin energy denominator, both conjugate members retained.','new_kernels':0,'interpolation':False,
 'finish_condition':'Five fixed kernel families at exactly fourteen frozen packet nodes, with matching reference regressions and saved internal q/t; no adaptive extra points.'}
J(HERE/'PREDECLARATION.json',pre)
saved={};families=[];start=time.perf_counter();regression_max=0.;closure_max=0.;new_count=0;ref_count=0
for family in binding['families']:
    rows=[];S={name:0. for name in ['qT','qL','Tline','Rq']};Sref=0.
    for index,(node,inp) in enumerate(zip(contract['central_nodes'],family['nodes'])):
        kp=HERE/inp['kernel'];mp=HERE/inp['metadata'];assert H(kp)==inp['kernel_sha256'] and H(mp)==inp['metadata_sha256']
        meta=json.loads(mp.read_text());om=node['omega'];assert abs(meta['omega']-om)<2e-13
        with np.load(kp,allow_pickle=False) as z:
            G=z['G0'].copy()
            if family['method']=='surface':B=[z['Im_branch_'+str(j)].copy() for j in range(3)];reference=meta['coefficient'];mapkeys={'qT':'qT','qL':'qL','Tline':'Tline','Rq':'Rq'}
            else:B=[z['Im'+str(j)].copy() for j in range(3)];reference=meta['scattering'];mapkeys={'qT':'TT_modal','qL':'TL_modal','Tline':'Tq','Rq':'Rq_amplitude'}
        BT=B[0]+B[1];BL=B[2]
        assert G.shape==(12,12) and all(np.isfinite(x).all() for x in [G,BT,BL])
        assert np.max(abs(G.imag-BT-BL))<2e-12
        assert max(np.max(abs(x-x.conj().T)) for x in [BT,BL])<2e-12
        assert min(np.linalg.eigvalsh(x).min() for x in [BT,BL])>-1e-12
        states={}
        for name,g in [('reference',g_ref),('ramp',g_new)]:
            result,internal=scatter_from_reduced_green(om,g,p,G,eta_lead=1e-12,return_internal=True)
            f=internal['t'];gamma=internal['gamma_in'];qt=float(2*gamma*np.vdot(f,BT@f).real);ql=float(2*gamma*np.vdot(f,BL@f).real)
            part={'qT':qt,'qL':ql,'Tline':result.Tq,'Rq':result.Rq};closure=abs(sum(part.values())-1)
            assert min(part.values())>=-1e-12 and max(part.values())<=1+1e-10
            assert closure<1e-8 and abs(qt+ql-result.Tbulk)<2e-11
            closure_max=max(closure_max,closure)
            key=f"{family['label']}_{node['node_id']}_{name}"
            for field in ['q','t','U','V','Kself','Sigma_bulk']:
                if field in internal:saved[key+'_'+field]=internal[field]
            states[name]={**part,'Tbulk':qt+ql,'closure_abs':closure,'gamma_in':gamma,'input_displacement_power':om*gamma/4}
            if name=='reference':
                ref_count+=1;err=max(abs(part[k]-reference[mapkeys[k]]) for k in mapkeys);assert err<1e-10,(family['label'],node['node_id'],err)
                regression_max=max(regression_max,err);states[name]['stored_reference_regression_max_abs']=err
            else:new_count+=1
        weight=node['weight_over_full_periodic_energy']
        for name in S:S[name]+=weight*states['ramp'][name]
        Sref+=weight*states['reference']['qT']
        rows.append({'node_id':node['node_id'],'omega':om,'positive_index':node['positive_k_index'],'negative_index':node['negative_k_index'],'paired_weight_over_full_periodic_energy':weight,'states':states,'kernel_sha256':inp['kernel_sha256']})
    c=spec['coverage'];tail=1-c
    families.append({'label':family['label'],'method':family['method'],'grid':family['grid'],'bulk_regulator':family.get('bulk_regulator'),
       'new_ramp_full_energy_in_band_contribution':S,'new_ramp_qT_passivity_tail_interval':[S['qT'],S['qT']+tail],
       'new_ramp_qT_conditional_in_band_mean':S['qT']/c,'original_reference_qT_in_band_contribution':Sref,
       'weighted_partition_sum':sum(S.values()),'rows':rows})
assert new_count==ref_count==70
Ehc=float(2/3*np.sum(g_new*(q[60:66]**2+q[61:67]**2)))
Eopen=spec['open_line_energy'];Etotal=Eopen+Ehc;factor=spec['periodic_parseval_energy']/Etotal
for f in families:
    f['qT_in_band_contribution_per_initial_open_total_energy']=f['new_ramp_full_energy_in_band_contribution']['qT']*factor
    f['qT_tail_interval_per_initial_open_total_energy']=[x*factor for x in f['new_ramp_qT_passivity_tail_interval']]
pairwise=[]
for a in range(len(families)):
    for b in range(a+1,len(families)):
        aa=families[a];bb=families[b];delta=[x['paired_weight_over_full_periodic_energy']*(y['states']['ramp']['qT']-x['states']['ramp']['qT']) for x,y in zip(aa['rows'],bb['rows'])]
        pairwise.append({'a':aa['label'],'b':bb['label'],'qT_weighted_signed_difference':sum(delta),'qT_weighted_L1_difference':sum(abs(x) for x in delta),'scope':'Observed cached-method/grid sensitivity, not a true-error bound.'})
np.savez_compressed(HERE/'INTERNAL_STATES.npz',**saved)
out={'status':'PASS','date':'2026-11-08','new_profile_g':g_new.tolist(),'physical_total_spring_stiffness':8*sum(g_new),'active_springs':48,
 'model_parameters':dataclasses.asdict(p),'preparation':config,'spectrum_summary':spec,'nyquist_fraction':float(W[n//2]/W.sum()),'propagating_fraction':float(W[prop].sum()/W.sum()),
 'initial_ramp_contact_energy':Ehc,'initial_open_total_energy':Etotal,'periodic_to_initial_open_total_energy_factor':factor,
 'families':families,'pairwise':pairwise,'reference_regressions':ref_count,'new_ramp_states':new_count,'new_green_evaluations':0,
 'maximum_reference_regression_abs':regression_max,'maximum_flux_closure_abs':closure_max,'saved_internal_array_count':len(saved),
 'input_binding_sha256':H(HERE/'INPUT_BINDING.json'),'predeclaration_sha256':H(HERE/'PREDECLARATION.json'),'producer_sha256':H(Path(__file__)),'internal_states_sha256':H(HERE/'INTERNAL_STATES.npz'),
 'limits':['The central finiteDFT weighted prediction differs from the monochromatic carrier value and from a finite-time stored-energy average.',
 'Each tail interval bounds omitted spectral weight only; it is not a numerical quadrature, finite-box, timestep, cadence or finite-time estimator error enclosure.',
 'Existing kernels are reused for a new contact profile. The regressions and closure checks are internal consistency; no new kernel production or rigorous kernel accuracy certificate is supplied here.',
 'Both signed propagating conjugate bins describe one real right-moving preparation; qT already counts both lower branches once.'],
 'elapsed_seconds':time.perf_counter()-start,'runtime':{'python':platform.python_version(),'numpy':np.__version__}}
J(OUT,out)
print(json.dumps({'status':'PASS','reference_regressions':ref_count,'new_ramp_states':new_count,'maximum_reference_regression_abs':regression_max,
 'maximum_flux_closure_abs':closure_max,'coverage':spec['coverage'],'initial_contact_energy':Ehc,
 'families':[{k:f[k] for k in ['label','new_ramp_full_energy_in_band_contribution','new_ramp_qT_passivity_tail_interval','original_reference_qT_in_band_contribution']} for f in families]},indent=2))
