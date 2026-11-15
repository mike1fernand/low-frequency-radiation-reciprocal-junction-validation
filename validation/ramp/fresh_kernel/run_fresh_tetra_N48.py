"""One fresh corrected N48 tetrahedron kernel, retained as two checkpoints."""
from pathlib import Path
import dataclasses, hashlib, json, shutil, sys, time
import numpy as np
from scipy import linalg

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'fresh_tetra_N48';OUT.mkdir(exist_ok=True)
R8=Path('/workspace/scratch/044aa7b17d80/p1_736_validation_2026_11_08/r8_selected/paper1-validation-integrated-r8/legacy_rc4/validation/current_numerics')
SOURCE=OUT/'sources';H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
J=lambda p,d:p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
files=['analysis/resolution_execution/tetra/multifrequency_tetra.py',
       'analysis/math_resolution/tetra/combined_tetra_diagnostic.py',
       'analysis/math_resolution/tetra/stable_scalar.py',
       'analysis/math_resolution/tetra/centroid_weights.py']
for rel in files:
    dst=SOURCE/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(R8/rel,dst)
assert not (OUT/'RESULT.json').exists()
sys.path.insert(0,str(SOURCE/'analysis/resolution_execution/tetra'))
import multifrequency_tetra as m
N=48;omega=1.2;regulator=1e-10;g=np.arange(1,7,dtype=float)/7;p=m.accepted.Parameters()
plan={'omega':omega,'grid':N,'shift':[0,0,0],'regulator_energy_imaginary':regulator,
      'candidate_g':g.tolist(),'new_kernel_evaluations':1,'additional_frequencies':0,
      'ranges':[[0,65536],[65536,N**3]],'cube_chunk':1024,
      'model_parameters':dataclasses.asdict(p),'source_hashes':{rel:H(SOURCE/rel) for rel in files},
      'scope':'Fresh independent corrected tetrahedron integration; one fixed point, no search or rigorous error enclosure.'}
J(OUT/'PREDECLARATION.json',plan);start=time.perf_counter()
ev,geometry=m.prepare_geometry(OUT/'geometry',N,p=p)
parts=[]
for first,last in plan['ranges']:
    t=time.perf_counter();part=m.integrate_range(ev,[omega],N,first,last,p=p,nseg=6,regulator=regulator,cube_chunk=1024)
    assert np.all(part['negative_scalar_count']==0) and np.all(part['minimum_scalar_imaginary']>=0)
    assert np.isfinite(part['sums']).all()
    fp=OUT/f'part_{first}_{last}.npz';np.savez_compressed(fp,**{k:part[k] for k in ['sums','minimum_scalar_imaginary','maximum_abs_scalar','negative_scalar_count']})
    info={k:v for k,v in part.items() if not isinstance(v,np.ndarray)}|{'sha256':H(fp)}
    J(OUT/f'part_{first}_{last}.json',info);parts.append(part)
    print(json.dumps({'completed_range':[first,last],'seconds':time.perf_counter()-t}),flush=True)
sums=sum(x['sums'][0] for x in parts);minimum=min(float(x['minimum_scalar_imaginary'][0]) for x in parts)
diag={'seconds':sum(x['seconds'] for x in parts),'minimum_scalar_imaginary':minimum,'range_count':2,'negative_scalar_count':0}
kernel=m.result_for_frequency(omega,N,regulator,sums,diag,nseg=6)
candidate=m.accepted.independent_scatter(omega,g,p,kernel)
reference=m.accepted.independent_scatter(omega,[.1,.2,.2,.2,.2,.1],p,kernel)
assert candidate.closure_amplitude<1e-8 and candidate.component_closure_modal<1e-10
np.savez_compressed(OUT/'KERNEL.npz',G0=kernel.G0,**{'Im'+str(j):x for j,x in enumerate(kernel.Im_branches)},ImL=kernel.Im_longitudinal_projector)
G=kernel.G0;U,V,Kself=m.accepted.static_line_core_matrices(g);Gc=linalg.solve(np.eye(12)+G@U,G);Sigma=-V.T@Gc@V
z=complex(omega,1e-12);sl=m.accepted.surface_self_energy_scalar_chain(omega,p,1e-12)
Kq=np.diag(np.full(7,2*p.k_spine))
for j in range(6):Kq[j,j+1]=Kq[j+1,j]=-p.k_spine
leads=np.zeros((7,7),complex);leads[0,0]=leads[-1,-1]=sl
A=p.m_spine*z*z*np.eye(7)-Kq-Kself-Sigma-leads;rhs=np.zeros(7,complex);rhs[0]=1
q=linalg.solve(A,rhs);f=linalg.solve(np.eye(12)+U@G,V@q)
np.savez_compressed(OUT/'INTERNAL_STATES.npz',G0=G,U=U,V=V,Kself=Kself,Sigma=Sigma,A=A,q=q,f=f)
result={'status':'PASS','omega':omega,'grid':N,'candidate_g':g.tolist(),'geometry':geometry,
        'candidate':dataclasses.asdict(candidate),'reference_control_same_fresh_kernel':dataclasses.asdict(reference),
        'kernel_diagnostics':kernel.diagnostics,'fresh_quadrature_seconds':time.perf_counter()-start,
        'kernel_sha256':H(OUT/'KERNEL.npz'),'state_sha256':H(OUT/'INTERNAL_STATES.npz'),
        'source_hashes':plan['source_hashes'],'original_sources_unchanged':all(H(R8/r)==H(SOURCE/r) for r in files),
        'scope':plan['scope']}
J(OUT/'RESULT.json',result)
print(json.dumps({'status':'PASS','seconds':result['fresh_quadrature_seconds'],'candidate':result['candidate'],'reference_qT':reference.TT_modal},indent=2))
