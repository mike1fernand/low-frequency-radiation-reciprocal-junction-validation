"""Independent exact Q-positivity and scalar point-bound review.

The scalar endpoints below are rigorous conditional deductions from declared
norm caps. This script does not certify those matrix norm caps or the external
Green-error premise; those require their separate validated receipts.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json
import numpy as np

ROOT=Path(__file__).resolve().parent
INPUT=Path('/workspace/scratch/044aa7b17d80/p1_coupling_assessment_2026_11_08/numerical_probe/inputs/REFERENCE_KERNEL.npz')
expected='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
assert hashlib.sha256(INPUT.read_bytes()).hexdigest()==expected
with np.load(INPUT,allow_pickle=False) as z:
    Q=z['Im_branch_0'][:12,:12]+z['Im_branch_1'][:12,:12]
    assert np.array_equal(Q,z['Im_qT'][:12,:12])
assert np.array_equal(Q,Q.conj().T) and np.count_nonzero(Q.imag)==0
q=[[F.from_float(float(x)) for x in row] for row in Q.real]
n=len(q);L=[[F(i==j) for j in range(n)] for i in range(n)];D=[]
for j in range(n):
    dj=q[j][j]-sum((L[j][k]**2*D[k] for k in range(j)),F(0))
    assert dj>0,('nonpositive exact LDL pivot',j)
    D.append(dj)
    for i in range(j+1,n):
        L[i][j]=(q[i][j]-sum((L[i][k]*D[k]*L[j][k] for k in range(j)),F(0)))/dj
assert all(sum((L[i][k]*D[k]*L[j][k] for k in range(n)),F(0))==q[i][j]
           for i in range(n) for j in range(n))

caps={'H_frobenius':F(947,1000),'J_frobenius':F(3,5),
      'eta_c_norm':F(569,1000),'a_lower':F(8578,10000),'a_upper':F(8579,10000)}
def enclosure(eG,eQ):
    Z=caps['eta_c_norm']/(1-caps['H_frobenius']*eG)
    d=caps['J_frobenius']*eG*Z
    assert 1-caps['H_frobenius']*eG>0
    assert caps['a_lower']>d
    lower=(caps['a_lower']-d)**2-2*eQ*Z**2
    upper=(caps['a_upper']+d)**2+2*eQ*Z**2
    return {'eG':str(eG),'eQ':str(eQ),'Z':float(Z),'d':float(d),
            'lower':float(lower),'upper':float(upper),
            'lower_rational':str(lower),'upper_rational':str(upper),
            'greater_than_72percent':lower>F(18,25),
            'greater_than_73percent':lower>F(73,100)}
out={
    'status':'PASS','date':'2026-11-08',
    'scope':'Exact dyadic Q SPD and exact rational scalar propagation conditional on matrix norm caps and Green-error premises. Not standalone kernel or matrix-norm certification.',
    'kernel_sha256':expected,'Q_branchsum_exactly_equals_Im_qT':True,
    'Q_exact_dyadic_SPD':True,'LDL_exact_reconstruction':True,
    'LDL_positive_pivots':[{'index':i,'approximation':float(x),
                           'numerator_bits':x.numerator.bit_length(),
                           'denominator_bits':x.denominator.bit_length()} for i,x in enumerate(D)],
    'norm_caps_assumed_not_certified_here':{k:str(v) for k,v in caps.items()},
    'advertised_kernel_caps':enclosure(F(1,50),F(1,200)),
    'decimal_ledger_approximation_not_used_for_certification':enclosure(F('0.003903451143041911'),F('0.0019669086816087237')),
    'safe_upward_caps':enclosure(F(1,250),F(1,500)),
    'derivation':[
        'A=Ac+Dg DeltaG Peta, R=Ac^-1, eta=Peta y, eta_c=Peta yc.',
        'eta=(I+H DeltaG)^-1 eta_c with H=Peta R Dg.',
        'Fc y=Fc yc-J DeltaG eta with J=Fc R Dg.',
        '||eta||<=Z=||eta_c||/(1-||H||eG), ||Fc(y-yc)||<=||J||eG Z.',
        'True T differs from ||Fc y||^2 by 2 eta^dagger DeltaQstar eta, bounded by 2eQ Z^2.',
        'Therefore (a_lower-d)^2-2eQ Z^2 <= T <= (a_upper+d)^2+2eQ Z^2.'
    ]}
ROOT.mkdir(parents=True,exist_ok=True)
(ROOT/'EXACT_Q_AND_SCALAR_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['LDL_positive_pivots','derivation']},indent=2))
