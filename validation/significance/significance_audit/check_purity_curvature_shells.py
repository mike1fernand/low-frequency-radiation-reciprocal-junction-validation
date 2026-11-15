#!/usr/bin/env python3
"""Small prescribed-source spectral-shell check, not a coupled junction solve."""
from pathlib import Path
import hashlib, json, sys, tempfile
import numpy as np

HERE=Path(__file__).resolve().parent
INPUTS=HERE.parent/'recheck/purity_inputs'
if not INPUTS.is_dir():
    INPUTS=HERE.parent/'paper1_next_correspondence/package/paper1-significance/recheck/purity_inputs'

def main():
    bonds=np.array([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)],float)
    unit=bonds/np.linalg.norm(bonds,axis=1)[:,None]
    projectors=np.einsum('ei,ej->eij',unit,unit)
    mu,wmu=np.polynomial.legendre.leggauss(12)
    phi=2*np.pi*np.arange(24)/24
    direction=np.array([[z,np.sqrt(1-z*z)*np.cos(p),np.sqrt(1-z*z)*np.sin(p)] for z in mu for p in phi])
    weight=np.repeat(wmu,24)*2*np.pi/24
    assert abs(weight.sum()-4*np.pi)<1e-13
    nb=direction@bonds.T
    eta0=6*np.sqrt(3)/(6*np.sqrt(3)+1)
    fixtures=[{'name':'one_motif_axial','X':[.5],'R':[0.],'xi':[0.]},
              {'name':'one_motif_with_radial','X':[.5],'R':[.6],'xi':[0.]},
              {'name':'two_motif_synthetic_force','X':[.1,.4],'R':[-.8,.9],'xi':[0.,3.]}]
    known={x['path']:x for x in json.loads((INPUTS/'INPUT_PROVENANCE.json').read_text())['files']}
    payload=(INPUTS/'model_api_v1_1.py').read_bytes()
    assert hashlib.sha256(payload).hexdigest()==known['model_api_v1_1.py']['sha256']
    max_root=0.;max_engine=0.;tie_evaluations=0;rows=[]
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='p1-curvature-shell-') as tmp:
        (Path(tmp)/'model_api_v1_1.py').write_bytes(payload);sys.path.insert(0,tmp)
        from model_api_v1_1 import ModelParameters,bulk_dynamical_matrix
        for case in fixtures:
            X,R,xi=(np.array(case[k]) for k in ['X','R','xi'])
            assert abs(2*X.sum()-1)<1e-14
            M=2*np.dot(xi,X);N=R.sum()/np.sqrt(2);P=np.dot(xi,R)/np.sqrt(2)
            beta=-2/21+2*(P-M*N)/3-N*N/3
            values=[]
            for x in [.08,.04,.02]:
                powers=[];saved=[]
                for branch,a0 in enumerate([1.,1.,3.]):
                    r=np.full(len(direction),x/np.sqrt(a0))
                    for _ in range(6):
                        phase=r[:,None]*nb
                        D=np.einsum('ne,eij->nij',4*np.sin(phase/2)**2,projectors)
                        vals,vec=np.linalg.eigh(D);mode=vec[:,:,branch]
                        dD=np.einsum('ne,eij->nij',2*np.sin(phase)*nb,projectors)
                        derivative=np.einsum('ni,nij,nj->n',mode,dD,mode)
                        r-=(vals[:,branch]-x*x)/derivative
                    phase=r[:,None]*nb
                    D=np.einsum('ne,eij->nij',4*np.sin(phase/2)**2,projectors)
                    vals,vec=np.linalg.eigh(D);mode=vec[:,:,branch]
                    derivative=np.einsum('ni,nij,nj->n',mode,np.einsum('ne,eij->nij',2*np.sin(phase)*nb,projectors),mode)
                    root=float(np.max(abs(vals[:,branch]-x*x))/(x*x));max_root=max(max_root,root)
                    assert root<2e-11
                    k=r[:,None]*direction
                    exp=np.exp(-1j*k[:,0,None]*xi[None,:])
                    Fx=2*(exp@X)*np.cos(k[:,1]/2)*np.cos(k[:,2]/2)
                    Fy=-1j*np.sqrt(2)*(exp@R)*np.sin(k[:,1]/2)*np.cos(k[:,2]/2)
                    Fz=-1j*np.sqrt(2)*(exp@R)*np.cos(k[:,1]/2)*np.sin(k[:,2]/2)
                    F=np.column_stack([Fx,Fy,Fz])
                    integrand=abs(np.sum(mode*F,axis=1))**2*r*r/derivative
                    if branch<2:
                        # At true transverse ties use the full radial residue.
                        # Each lower-sheet pass receives half the same rank-two
                        # contribution, avoiding two unrelated eigenbasis choices.
                        ties=abs(vals[:,1]-vals[:,0])<=128*np.finfo(float).eps*np.max(abs(vals),axis=1)
                        if np.any(ties):
                            VT=vec[ties,:,:2]
                            dD=np.einsum('ne,eij->nij',2*np.sin(phase[ties])*nb[ties],projectors)
                            slopes=np.einsum('nik,nij,njl->nkl',VT,dD,VT)
                            FT=np.einsum('nik,ni->nk',VT,F[ties])
                            residue=np.einsum('ni,ni->n',FT.conj(),np.linalg.solve(slopes,FT[...,None])[...,0]).real
                            integrand[ties]=r[ties]**2*residue/2
                            tie_evaluations+=int(np.count_nonzero(ties))
                    powers.append(float(weight@integrand))
                    saved.append({'branch':branch,'radial_roots':r.tolist(),'spectral_integrands':integrand.tolist()})
                    for idx in [0,61,145]:
                        err=float(np.max(abs(D[idx]-bulk_dynamical_matrix(k[idx],ModelParameters()))));max_engine=max(max_engine,err);assert err<2e-14
                eta=(powers[0]+powers[1])/sum(powers)
                curvature=(eta-eta0)/(eta0*(1-eta0)*x*x)
                values.append({'x':x,'branch_powers_common_normalization':powers,'eta_T':eta,'curvature_estimate':float(curvature),'shell_state':saved})
            extrapolated=(4*values[-1]['curvature_estimate']-values[-2]['curvature_estimate'])/3
            error=abs(extrapolated-beta)
            assert error<5e-5,(case['name'],beta,extrapolated)
            # Translating every xi only gives one global Fourier phase.
            h=7.;assert abs((P+h*N)-(M+h)*N-(P-M*N))<2e-15
            if len(X)==1:assert beta<=-2/21
            rows.append({**case,'M':float(M),'N':float(N),'P':float(P),'predicted_beta':float(beta),'extrapolated_beta':float(extrapolated),'absolute_extrapolation_error':float(error),'rows':values})
    result={'status':'PASS','scope':'Fixed prescribed complete-motif force patterns and actual H18 dispersion; no physical g profile solved, no Green kernel or production junction response.',
            'angular_orders':[12,24],'frequency_points':[.08,.04,.02],'tolerance_frozen_before_run':5e-5,
            'degenerate_cluster_evaluations':tie_evaluations,'max_relative_radial_root_residual':max_root,'max_actual_engine_symbol_difference':max_engine,
            'engine_sha256':hashlib.sha256(payload).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'cases':rows,'limits':['Synthetic force moments do not prove passive physical attainability','No finite-frequency validity window certified','Threshold and curvature proofs remain the basis of the theorem, not this floating-point check']}
    (HERE/'PURITY_CURVATURE_SHELL_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
    print(json.dumps([{k:v for k,v in c.items() if k!='rows'} for c in rows],indent=2))

if __name__=='__main__':main()
