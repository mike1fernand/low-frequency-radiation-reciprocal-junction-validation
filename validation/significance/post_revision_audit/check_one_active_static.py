#!/usr/bin/env python3
"""One-active-motif static identity: exact parity algebra + archived blocks.

Synthetic SPD compliance fixtures validate the reduction only. They are not
new physical static kernels, numerical predictions for a specific bulk, or
full scattering calculations.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,tempfile
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
INPUTS=ROOT/'recheck/purity_inputs'

def main():
    exact=[]
    fixtures=[(F(169,100),F(1,5),F(3,10),F(4,5)),
              (F(3,10),F(2),F(7,10),F(1,10)),
              (F(4),F(1,20),F(1,5),F(6,5)),
              (F(7,10),F(9,10),F(2,5),F(3,5))]
    for K,g,chiR,chiX in fixtures:
        b=2*g/3
        ae=2*b/(1+b*chiX);ao=2*K+2*b/(1+2*b*chiR)
        qL=(1/ae+1/ao)/2;qR=(1/ae-1/ao)/2
        X=b*(qL+qR)/(1+b*chiX)
        nu=b*(qL-qR)/(1+2*b*chiR)
        target=1/(2+3*K/g+4*K*chiR)
        assert X==F(1,2) and nu==target and 0<nu<F(1,2)
        beta=-F(2,21)-nu*nu/3
        assert -F(5,28)<beta<-F(2,21)
        exact.append({'Kq':str(K),'g_actual':str(g),'chi_R':str(chiR),'chi_X':str(chiX),'nu_R':str(nu),'beta_share':str(beta)})
    K=F(169,100);chi=F(3,10)
    gs=[F(1,20),F(1,5),F(2)]
    nus=[1/(2+3*K/g+4*K*chi) for g in gs]
    betas=[-F(2,21)-nu*nu/3 for nu in nus]
    assert nus[0]<nus[1]<nus[2] and betas[0]>betas[1]>betas[2]
    metadata=json.loads((INPUTS/'INPUT_PROVENANCE.json').read_text())
    known={x['path']:x for x in metadata['files']}
    modules={n:(INPUTS/n).read_bytes() for n in ['model_api_v1_1.py','reduced_core.py']}
    for n,data in modules.items():assert hashlib.sha256(data).hexdigest()==known[n]['sha256']
    rng=np.random.default_rng(20261107);records=[];maxerr=0.;adverse=None
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='one-active-static-') as td:
        for n,data in modules.items():(Path(td)/n).write_bytes(data)
        sys.path.insert(0,td)
        from model_api_v1_1 import ModelParameters
        from reduced_core import reduced_static_matrices
        for K,g,chiR,chiX in fixtures:
            par=ModelParameters(k_spine=float(K),m_spine=2.1,m_bulk=.8)
            for nseg in [1,4,6]:
                for active in sorted({0,nseg//2,nseg-1}):
                    profile=np.zeros(nseg);profile[active]=float(g)
                    U,V,Ks=reduced_static_matrices(profile,par)
                    d=rng.uniform(.5,1.5,2*nseg);d[2*active]=float(chiX);d[2*active+1]=float(chiR)
                    off=rng.normal(size=(2*nseg,2*nseg));off=(off+off.T)/2
                    np.fill_diagonal(off,0);off[2*active,2*active+1]=off[2*active+1,2*active]=0
                    off=.1*off/max(1,np.linalg.norm(off,2))
                    G=np.sqrt(d)[:,None]*(np.eye(2*nseg)+off)*np.sqrt(d)[None,:]
                    assert np.linalg.eigvalsh(G).min()>0
                    L=np.zeros((nseg+1,nseg+1))
                    for j in range(nseg):L[j,j]+=float(K);L[j+1,j+1]+=float(K);L[j,j+1]-=float(K);L[j+1,j]-=float(K)
                    A=L+Ks-V.conj().T@np.linalg.solve(np.linalg.inv(G)+U,V)
                    D=np.linalg.solve(np.eye(2*nseg)+U@G,V)
                    for side,idx in [('left',0),('right',nseg)]:
                        e=np.eye(nseg+1)[:,idx];q=np.linalg.solve(A,e);f=-D@q
                        sign=1 if side=='left' else -1
                        target=sign/float(2+3*K/g+4*K*chiR)
                        nu=float(f[2*active+1].real/np.sqrt(2));X=float(f[2*active].real)
                        inactive=np.delete(f,[2*active,2*active+1])
                        err=max(abs(X-.5),abs(nu-target),float(np.max(abs(inactive),initial=0)),float(np.linalg.norm(A@q-e)))
                        maxerr=max(maxerr,err);assert err<1e-11
                        records.append({'nseg':nseg,'active':active,'incidence':side,'nu_R':nu,'formula':target,'max_error':err})
                    if adverse is None and nseg==1:
                        bad=G.copy();bad[0,1]=bad[1,0]=.05*np.sqrt(float(chiR*chiX))
                        Ab=L+Ks-V.conj().T@np.linalg.solve(np.linalg.inv(bad)+U,V)
                        Db=np.linalg.solve(np.eye(2)+U@bad,V)
                        fb=-Db@np.linalg.solve(Ab,np.array([1.,0.]))
                        diff=abs(float(fb[1].real/np.sqrt(2))-1/float(2+3*K/g+4*K*chiR))
                        assert diff>1e-5
                        adverse={'broken_local_XR_reflection_difference':diff,'scope':'Rejects using the diagonal-local-compliance formula after breaking its symmetry premise.'}
    out={'status':'PASS','scope':'Exact one-active parity identity and finite actual-block/synthetic-SPD algebra checks only; no physical static kernel or scattering calculation.',
         'exact_rational_fixtures':exact,'matrix_cases':len(records),'max_matrix_error':maxerr,
         'matrix_fixtures':records,'adverse_control':adverse,'coupling_ordering_exact':True,
         'source_hash':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'engine_hashes':{n:hashlib.sha256(data).hexdigest() for n,data in modules.items()}}
    (HERE/'ONE_ACTIVE_STATIC_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['matrix_fixtures','exact_rational_fixtures']},indent=2))

if __name__=='__main__':main()
