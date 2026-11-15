"""Bind independent archived tetrahedron kernels and evaluate the new profile.

No new tetrahedron quadrature. Every matching old reference row is regressed
before the new ramp response, using the independently implemented solver.
"""
from pathlib import Path
import dataclasses, hashlib, json, shutil, sys, time
import numpy as np
from scipy import linalg

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
R8=Path('/workspace/scratch/044aa7b17d80/p1_736_validation_2026_11_08/r8_selected/paper1-validation-integrated-r8/legacy_rc4/validation/current_numerics')
ACCEPTED=R8/'analysis/math_resolution/tetra'
TETRA=R8/'analysis/resolution_execution/tetra'
OUT=HERE/'independent_tetra';INPUT=OUT/'inputs';INPUT.mkdir(parents=True,exist_ok=True)
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p,d:p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
assert not (OUT/'RESULT.json').exists(),'Preserve frozen results.'
for name in ['combined_tetra_diagnostic.py','stable_scalar.py','centroid_weights.py']:
    shutil.copyfile(ACCEPTED/name,INPUT/name)
sys.path.insert(0,str(INPUT))
import combined_tetra_diagnostic as accepted
p=accepted.Parameters();p.validate()
assert dataclasses.asdict(p)=={'a':1.,'m_bulk':1.,'m_spine':1.,'k_axis':1.,'k_diag':1.,'k_spine':1.69}
omega=1.2;ramp=np.arange(1,7,dtype=float)/7;ref=np.array([.1,.2,.2,.2,.2,.1])
cases=[('reference_N96','omega_001'),('reference_N128','omega_001'),
       ('reference_N128_halfshift','omega_001'),('reference_N160','omega_001'),('reference_N192','omega_000')]
plan={'date':'2026-11-08','omega':omega,'candidate_g':ramp.tolist(),
      'cases':cases,'new_tetrahedron_quadratures':0,'new_candidate_states':5,
      'reference_regression_states':5,'finish':'Bind kernels to metadata/config/source, regress existing reference, then calculate ramp once per kernel.',
      'producer_modules':{n:H(INPUT/n) for n in ['combined_tetra_diagnostic.py','stable_scalar.py','centroid_weights.py']}}
J(OUT/'PREDECLARATION.json',plan)
rows=[];states={};start=time.perf_counter()
for case,stem in cases:
    directory=TETRA/case;meta=json.loads((directory/(stem+'.json')).read_text());conf=json.loads((directory/'CONFIG.json').read_text())
    assert meta['omega']==omega and meta['grid']==conf['config']['grid']
    assert conf['config']['profile']==ref.tolist() and conf['config']['nseg']==6
    assert conf['config']['parameters']==dataclasses.asdict(p)
    assert conf['config']['omega'][int(stem[-3:])]==omega
    assert H(directory/(stem+'.npz'))==meta['kernel_sha256']
    for name in plan['producer_modules']:
        key='analysis/math_resolution/tetra/'+name
        assert conf['config']['source_sha256'][key]==H(INPUT/name)
    for name in [stem+'.npz',stem+'.json','CONFIG.json']:
        dest=INPUT/case/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(directory/name,dest)
    with np.load(directory/(stem+'.npz'),allow_pickle=False) as z:
        mats={k:z[k].copy() for k in z.files}
    assert all(v.shape==(12,12) and np.isfinite(v).all() for v in mats.values())
    G=mats['G0'];branches=[mats['Im'+str(j)] for j in range(3)]
    assert np.max(abs(G.imag-sum(branches)))<1e-12
    assert min(np.linalg.eigvalsh(M).min() for M in branches)>-1e-12
    kernel=accepted.TetrahedronResult(omega=omega,energy_imaginary_regulator=meta['regulator'],
        grid=meta['grid'],nseg=6,G0=G,Im_branches=branches,
        Im_longitudinal_projector=mats['ImL'],
        Im_branch_longitudinal_projectors=[mats['ImBL'+str(j)] for j in range(3)],
        symmetry_residual=meta['diagnostics']['reciprocity_residual_constructed'],
        min_imaginary_density_eigenvalue=meta['diagnostics']['spectral_min_eigenvalue'],
        seconds=0.,diagnostics=meta['diagnostics'])
    rr=accepted.independent_scatter(omega,ref,p,kernel);published=meta['scattering']
    regression=max(abs(getattr(rr,key)-published[key]) for key in ['Rq_amplitude','Tq','Tbulk','TT_modal','TL_modal'])
    assert regression<1e-11,regression
    result=accepted.independent_scatter(omega,ramp,p,kernel)
    assert result.closure_amplitude<1e-8 and result.component_closure_modal<1e-12
    # Save internal candidate vectors and matrices using this solver's own
    # independently implemented construction, without calling reduced_core.
    U,V,Kself=accepted.static_line_core_matrices(ramp)
    Gc=linalg.solve(np.eye(12)+G@U,G);SigmaB=-V.T@Gc@V
    nq=7;Kq=np.diag(np.full(nq,2*p.k_spine))
    for j in range(nq-1):Kq[j,j+1]=Kq[j+1,j]=-p.k_spine
    eta=1e-12;z=complex(omega,eta);s=accepted.surface_self_energy_scalar_chain(omega,p,eta)
    leads=np.zeros((nq,nq),complex);leads[0,0]=leads[-1,-1]=s
    A=p.m_spine*z*z*np.eye(nq)-Kq-Kself-SigmaB-leads
    rhs=np.zeros(nq,complex);rhs[0]=1;q=linalg.solve(A,rhs)
    f=linalg.solve(np.eye(12)+U@G,V@q);gamma=float(-2*s.imag)
    tp=float(2*gamma*np.vdot(f,(branches[0]+branches[1])@f).real)
    assert abs(tp-result.TT_modal)<1e-12
    row={'case':case,'omega':omega,'grid':meta['grid'],'shift':meta['shift'],
         'regulator_energy_imaginary':meta['regulator'],'reference_regression_max_abs':regression,
         'candidate':dataclasses.asdict(result),'bulk_qT_share':result.TT_modal/result.Tbulk,
         'relative_line_residual':float(np.linalg.norm(A@q-rhs)/(np.linalg.norm(A)*np.linalg.norm(q)+np.linalg.norm(rhs))),
         'candidate_kernel_sha256':meta['kernel_sha256'],'config_sha256':H(directory/'CONFIG.json')}
    rows.append(row)
    for name,val in [('G',G),('BT',branches[0]+branches[1]),('BL',branches[2]),('U',U),('V',V),('Kself',Kself),('A',A),('q',q),('f',f),('SigmaB',SigmaB)]:states[case+'__'+name]=val
sf=OUT/'INTERNAL_STATES.npz';np.savez_compressed(sf,**states)
record={'status':'PASS','method':'Recovered independent corrected linear-tetrahedron bare kernels; new profile evaluated by independent tetrahedron solver.',
        'omega':omega,'g':ramp.tolist(),'model_parameters':dataclasses.asdict(p),'rows':rows,
        'new_candidate_states':5,'reference_regression_states':5,'new_quadrature_evaluations':0,
        'input_source_archive':'PAPER1_VALIDATION_INTEGRATED_R8(1).zip, selected extraction; source and kernel/config hashes checked per row.',
        'module_sha256':plan['producer_modules'],'internal_states_sha256':H(sf),'seconds':time.perf_counter()-start,
        'accuracy_scope':'Independent quadrature route and mesh/shift ladder. Numerical convergence/sensitivity only; no rigorous infinite-lattice error bound.'}
J(OUT/'RESULT.json',record)
print(json.dumps({'status':'PASS','seconds':record['seconds'],'rows':[{'case':x['case'],'qT':x['candidate']['TT_modal'],'qL':x['candidate']['TL_modal'],'regression':x['reference_regression_max_abs']} for x in rows]},indent=2))
