"""Check normalization and static identities used in the P1 predecessor comparison.

No bulk quadrature, scattering production or novelty test is performed.
The static checks use exact archived contact blocks and synthetic SPD compliance
matrices to discriminate an algebraic identity, not to model a new material.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, json, sys, tempfile
import numpy as np


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=here.parent/'recheck'/'purity_inputs')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a new output filename; frozen reports are retained')

    # Integral over polar cosine: integral mu^2 dmu = 2/3.
    # Angular integrals divided by pi; the transverse projector includes both
    # transverse polarizations, rather than assigning an extra degeneracy later.
    angular_L = 2*Q(2, 3)
    angular_T = 4-angular_L
    # pi*delta(c^2 k^2-omega^2) with Fourier measure d^3k/(2pi)^3
    # gives b_C*pi*rho*c_C^3 = (angular_C/pi)/16.
    bL = angular_L/16
    bT = angular_T/16
    assert bL == Q(1, 12) and bT == Q(1, 6)
    # Pout=(omega/2)|F|^2 Im(chi), F=2 Z omega A, Pin=Z omega^2 |A|^2/2.
    # Cancelling produces T_C=4 Z b_C omega^2.
    assert 4*bT == Q(2, 3) and 4*bL == Q(1, 3)

    Z=Q(13, 10); ksum=Q(8, 3)
    crossover=[]
    for lam in [Q(0), Q(1, 10), Q(1), Q(10)]:
        scale=lam*ksum
        reflection=scale*scale/(scale*scale+(2*Z)**2)
        transmission=(2*Z)**2/(scale*scale+(2*Z)**2)
        assert reflection+transmission == 1
        crossover.append({'lambda':str(lam),'R_leading_exact':str(reflection),
                          'Tline_leading_exact':str(transmission)})

    provenance=json.loads((args.inputs/'INPUT_PROVENANCE.json').read_text())
    records={x['path']:x for x in provenance['files']}
    modules=['model_api_v1_1.py','reduced_core.py']
    payloads={name:(args.inputs/name).read_bytes() for name in modules}
    for name,data in payloads.items():
        assert hashlib.sha256(data).hexdigest() == records[name]['sha256']
    rng=np.random.default_rng(20261106)
    cases=[]
    sys.dont_write_bytecode=True
    with tempfile.TemporaryDirectory(prefix='p1-comparison-contact-') as td:
        stage=Path(td)
        for name,data in payloads.items():
            (stage/name).write_bytes(data)
        sys.path.insert(0,str(stage))
        from model_api_v1_1 import ModelParameters
        from reduced_core import reduced_static_matrices
        for n in [1,2,6]:
            profiles=[np.ones(n), np.linspace(.1,.6,n),
                      np.array([.4]+[0.]*(n-1))]
            for profile_index,g in enumerate(profiles):
                for compliance_scale in [.2, 3.]:
                    U,V,Kself=reduced_static_matrices(g,ModelParameters())
                    U=U.real;V=V.real
                    X=rng.normal(size=(2*n,2*n))
                    G=compliance_scale*(np.eye(2*n)+X@X.T/(2*n))
                    L=np.zeros((n+1,n+1))
                    for j in range(n):
                        v=np.zeros(n+1);v[j]=-1;v[j+1]=1
                        L+=(.8+.1*j)*np.outer(v,v)
                    H0=np.linalg.inv(G)+U
                    A0=L+Kself-V.T@np.linalg.solve(H0,V)
                    D0=np.linalg.solve(np.eye(2*n)+U@G,V)
                    e=np.ones(n+1);p=np.tile([2.,0.],n);eL=np.eye(n+1)[:,0]
                    errors=[np.max(np.abs(V@e+U@p)),
                            np.max(np.abs(Kself@e+V.T@p)),
                            np.max(np.abs(e@A0+p@D0)),
                            abs(p@D0@np.linalg.solve(A0,eL)+1)]
                    assert np.linalg.eigvalsh(A0).min()>0
                    assert max(errors)<1e-11,errors
                    cases.append({'motifs':n,'profile_index':profile_index,
                                  'compliance_scale':compliance_scale,
                                  'maximum_identity_error':float(max(errors))})
        # Adverse controls: the theorem cannot invert an uncoupled contact or
        # an unattached component of a disconnected line.
        n=6;g=np.array([.4,0.,0.,0.,0.,0.])
        U,V,Kself=reduced_static_matrices(g,ModelParameters());G=np.eye(2*n)
        L=np.zeros((n+1,n+1))
        for j in range(n):
            if j==3:continue
            v=np.zeros(n+1);v[j]=-1;v[j+1]=1;L+=np.outer(v,v)
        A=L+Kself-V.T@np.linalg.solve(np.linalg.inv(G)+U,V)
        assert abs(np.linalg.eigvalsh(A).min())<1e-12
        Lconnected=np.zeros_like(L)
        for j in range(n):
            v=np.zeros(n+1);v[j]=-1;v[j+1]=1;Lconnected+=np.outer(v,v)
        assert abs(np.linalg.eigvalsh(Lconnected).min())<1e-12

    out={'status':'PASS','scope':'Independent normalization/static-algebra checks for a literature comparison, not a new radiation result, full theorem replay or novelty proof.',
         'angular_integrals_divided_by_pi':{'T':str(angular_T),'L':str(angular_L)},
         'pi_rho_c_cubed_times_b':{'T':str(bT),'L':str(bL)},
         'pi_rho_c_cubed_over_Z_times_C':{'T':str(4*bT),'L':str(4*bL)},
         'synthetic_SPD_static_cases':cases,
         'adverse_controls':'PASS: all-zero contact and disconnected unattached line component have singular static matrices',
         'joint_limit_leading_line_cases':crossover,
         'input_modules_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in payloads.items()},
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({'status':out['status'],'static_cases':len(cases),
                      'maximum_static_error':max(x['maximum_identity_error'] for x in cases),
                      'output':str(args.output)}))


if __name__=='__main__':main()
