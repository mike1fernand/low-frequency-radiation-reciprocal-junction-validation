#!/usr/bin/env python3
"""Independent audit of stored kernels; all writes remain in this audit folder."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import importlib.util, itertools, json, hashlib, subprocess, sys, shutil, tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent/'baseline/Prediction_Calculations_2026-09-26/low_frequency'
spec=importlib.util.spec_from_file_location('baseline_lowfreq',BASE/'low_frequency.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def direct_form(k):
    # Eq. s4-XR-basis and Eq. s5-form-factor: explicit 72 physical coordinates.
    F=np.zeros((12,3),complex)
    for j in range(6):
        for sy,sz in itertools.product([-1,1],repeat=2):
            phase=np.exp(1j*np.dot(k,[j,sy/2,sz/2]))
            F[2*j,0]+=phase/2
            F[2*j+1,1]+=sy*phase/(2*np.sqrt(2))
            F[2*j+1,2]+=sz*phase/(2*np.sqrt(2))
    return F

out={}
rng=np.random.default_rng(926)
errs=[]
for _ in range(31):
    k=rng.uniform(.01,np.pi,3)
    direct=np.zeros((12,12),complex)
    for signs in itertools.product([-1,1],repeat=3):
        kk=k*np.array(signs);F=direct_form(kk)
        D=np.zeros((3,3))
        for d in m.BONDS:
            unit=d/np.linalg.norm(d)
            D+=4*np.sin(np.dot(kk,d)/2)**2*np.outer(unit,unit)
        direct+=F@np.linalg.inv(D)@F.conj().T/8
    fast=m.assemble(m.projected(np.linalg.inv(m.dyn(k)),k))
    errs.append(float(np.linalg.norm(direct-fast)/np.linalg.norm(direct)))
out['explicit_72_coordinate_parity_max_relative_error']=max(errs)

# Separate spherical acoustic angular integration, with no reduced projection routine.
x,w=np.polynomial.legendre.leggauss(64)
phi=np.linspace(0,2*np.pi,128,endpoint=False)
n=np.stack(np.broadcast_arrays(x[:,None],np.sqrt(1-x[:,None]**2)*np.cos(phi),np.sqrt(1-x[:,None]**2)*np.sin(phi)),axis=-1)
weights=w[:,None]*np.ones((1,len(phi)))*2*np.pi/len(phi)
pl=n[...,0]**2
integ_L=float(np.sum(weights*pl));integ_T=float(np.sum(weights*(1-pl)))
b=np.array([integ_T,integ_L/3**1.5])/(16*np.pi**2)
out['independent_point_force_coefficients']={'angular_T':integ_T,'angular_L':integ_L,'T_coefficients':(4*m.CQ*b).tolist(),'relative_error':(b/np.array([1/(6*np.pi),1/(12*np.pi*3**1.5)])-1).tolist()}

z=np.load(BASE/'validation/kernel_0.00000000.npz');G=z['PV']
z=np.load(BASE/'validation/kernel_0.00078125.npz');P=z['PV'];bs=[z[f'Q{j}'] for j in range(3)]
profiles={'only_first':[.3,0,0,0,0,0],'only_last':[0,0,0,0,0,.3],'disconnected_active_segments':[0,.2,0,0,.1,0]}
rows=[]
e=np.ones(7);p=np.tile([2.,0.],6);el=np.eye(7)[:,0]
for name,g in profiles.items():
    g=np.array(g);U,V,Ks,L=m.matrices(g);L[0,0]-=m.KQ;L[-1,-1]-=m.KQ
    A=L+Ks-V.T@np.linalg.solve(np.eye(12)+G@U,G)@V
    D=np.linalg.solve(np.eye(12)+U@G,V)
    row=m.scatter(.00078125,1,P,bs,g)
    rows.append({'name':name,'profile':g.tolist(),'min_static_eigenvalue':float(np.linalg.eigvalsh(A).min()),'static_resultant':float(p@D@np.linalg.solve(A,el)),'force_balance_residual':float(np.linalg.norm(e@A+p@D)),'qT_coefficient':row['TqT']/.00078125**2,'qT_coefficient_relative_error':row['TqT']/.00078125**2/(4*m.CQ*b[0])-1,'closure':row['closure']})
out['nonnegative_sparse_profile_checks']=rows

# A new off-grid kernel, independently rerun through the numerical entry point.
Pnew,bnew,diag=m.face_quadrature(.017,18,96)
row=m.scatter(.017,1,Pnew,bnew)
out['fresh_offgrid']={**diag,'qT_coefficient':row['TqT']/.017**2,'qL_coefficient':row['TqL']/.017**2,'closure':row['closure']}

# Reproduce the original stale-cache bug and verify the guarded copy via its CLI.
copy=ROOT/'cache_original.py';shutil.copyfile(BASE/'low_frequency.py',copy)
def run(script,tag,nang):
    result=subprocess.run([sys.executable,str(script),'--tag',tag,'--nang',str(nang),'--nrad','12','--frequencies','0'],capture_output=True,text=True,check=True)
    return json.loads((ROOT/tag/'kernel_0.00000000.json').read_text())
before=run(copy,'cache_bug',4);after=run(copy,'cache_bug',5)
guard=ROOT/'low_frequency_guarded.py'
with tempfile.TemporaryDirectory(prefix='cache_strict_',dir=ROOT) as temp:
    tag=Path(temp).name
    gbefore=run(guard,tag,4)
    original_hash=gbefore['matrix_sha256']
    status_bytes=(Path(temp)/'status.json').read_bytes()
    failed=subprocess.run([sys.executable,str(guard),'--tag',tag,'--nang','5','--nrad','12','--frequencies','0'],capture_output=True,text=True)
    status_unchanged=status_bytes==(Path(temp)/'status.json').read_bytes()
    greuse=run(guard,tag,4)
    specg=importlib.util.spec_from_file_location('guarded_lowfreq',guard)
    mg=importlib.util.module_from_spec(specg);specg.loader.exec_module(mg)
    path=Path(temp)/'kernel_0.00000000.npz';meta=path.with_suffix('.json')
    checks={}
    for name,request in [('angular',(0,5,12,gbefore['script_sha256'])),('radial',(0,4,13,gbefore['script_sha256'])),('source',(0,4,12,'bad')),('rounded_frequency',(1e-12,4,12,gbefore['script_sha256']))]:
        try:mg.cached_kernel(path,*request);checks[name]=False
        except ValueError:checks[name]=True
    for name,changed in [('matrix_hash',{**gbefore,'matrix_sha256':'bad'}),('missing_source_hash',{k:v for k,v in gbefore.items() if k!='script_sha256'}),('invalid_json',None)]:
        meta.write_text(json.dumps(changed) if changed is not None else '{broken')
        try:mg.cached_kernel(path,0,4,12,gbefore['script_sha256']);checks[name]=False
        except ValueError:checks[name]=True
    meta.unlink()
    try:mg.cached_kernel(path,0,4,12,gbefore['script_sha256']);checks['missing_metadata']=False
    except ValueError:checks['missing_metadata']=True
    out['cache_reproduction']={'original_first_actual_n_face':before['n_face'],'original_second_requested_n_face':5,'original_second_actual_n_face':after['n_face'],'guarded_mismatch_refused':failed.returncode!=0,'guarded_refusal_message_mentions_fresh_tag':'fresh --tag' in failed.stderr,'guarded_matrix_unchanged':hashlib.sha256(path.read_bytes()).hexdigest()==original_hash,'guarded_status_unchanged_after_refusal':status_unchanged,'guarded_identical_request_reuses_matrix':gbefore['matrix_sha256']==greuse['matrix_sha256'],'guarded_identical_request_reuses_timing_metadata':gbefore['seconds']==greuse['seconds'],'strict_refusal_checks':checks}

(ROOT/'audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
