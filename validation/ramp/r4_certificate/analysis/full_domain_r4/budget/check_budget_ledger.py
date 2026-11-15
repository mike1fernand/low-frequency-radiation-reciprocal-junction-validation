"""Independent exact/algebraic controls for the R4 budget gate."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import budget_ledger as bl
import matrix_budget as mb


def run():
    checks = []
    def need(ok, name):
        if not ok: raise ArithmeticError(name)
        checks.append(name)
    def rejects(name, call):
        try: call()
        except (ValueError, TypeError, ArithmeticError): checks.append(name)
        else: raise ArithmeticError(name+' accepted invalid input')
    z = lambda re, im=0: (F(re), F(im))
    need(bl.exact_positive_definite([[z(2),z(0,1)],[z(0,-1),z(2)]]), 'Complex Hermitian eigenvalues1,3 accepted')
    need(not bl.exact_positive_definite([[z(1),z(2)],[z(2),z(1)]]), 'Indefinite eigenvalues-1,3 rejected')
    need(not bl.exact_positive_definite([[z(1),z(1)],[z(1),z(1)]]), 'Singular eigenvalues0,2 rejected')
    rejects('Non-Hermitian positivity input rejected',lambda:bl.exact_positive_definite([[z(1),z(1,1)],[z(1,1),z(2)]]))
    d=[[z(1),z(0,1)],[z(0,1),z(-1)]]
    u,_=bl.exact_operator_upper(d)
    need(F(2)<u<F(200000001,100000000), 'Nonnormal nilpotent norm2 bounded; spectral radius0 not used')
    d=[[z(1,1),z(0)],[z(0),z(1,-1)]]
    u,_=bl.exact_operator_upper(d)
    need(u*u>2 and u<F(141421357,100000000), 'Complex sqrt2 diagonal norm bounded exactly')
    original=np.linalg.norm
    try:
        np.linalg.norm=lambda *a,**k:0.0
        rejects('False floating norm proposal cannot pass exact LDL',lambda:bl.exact_operator_upper(d))
    finally:
        np.linalg.norm=original
    zero_target=[[[0,0] for _ in range(12)] for _ in range(12)]
    def uniform(real,imag):
        return [[[-real,real],[-imag,imag]] for _ in range(17)]
    good=bl.assess(uniform(F(1,4000),F(0)),zero_target,'Q')
    need(good['arithmetic_budget_pass'], 'Known Q row-bound0.0046475 passes0.005')
    need(good['scaled_operator_upper']==F(1859,400000), 'Exact11nonzeroentries and1.69scaling retained')
    bad=bl.assess(uniform(F(1,3000),F(0)),zero_target,'Q')
    need(not bad['arithmetic_budget_pass'], 'Known excessive Q sufficient bound rejected')
    both=bl.assess(uniform(F(1,1100),F(1,1100)),zero_target,'g')
    only_real=bl.assess(uniform(F(1,1100),F(0)),zero_target,'g')
    need(not both['arithmetic_budget_pass'] and only_real['arithmetic_budget_pass'], 'Complex g requires both real and imaginary uncertainty')
    residual_target=[[[0,0] for _ in range(12)] for _ in range(12)]
    residual_target[0][1]=['1/1000','0']
    residual_target[1][0]=['1/1000','0']
    rr=bl.assess(uniform(F(0),F(0)),residual_target,'Q')
    need(rr['scaled_operator_upper']>=F(169,100000), 'Actual off-template b0 target residual retained')
    records=[{'id':'x','intervals':uniform(F(1,4000),F(0)),'multiplicity':2}]
    ranked=bl.rank_regions(records)
    need(ranked[0]['scaled_uncertainty_operator_upper']==2*good['scaled_operator_upper'], 'Physical symmetry reuse counted twice without independence assumption')
    rejects('Fractional physical multiplicity rejected',lambda:bl.rank_regions([dict(records[0],multiplicity=F(1,2))]))
    rejects('Float endpoints rejected',lambda:bl.assess(uniform(0.1,0),zero_target,'Q'))
    for kind in ('g','Q'):
        allocation=bl.recommend_allocations(F(1,10000),kind,[str(F(k+1,5)) for k in range(17)])
        radii=allocation['global_moment_modulus_radii']
        actual=max(sum(row[k]*radii[k] for k in range(17)) for row in mb.multiplicities()['rows'])
        need(mb.KSTAR*(F(1,10000)+actual)<=mb.THRESHOLDS[kind],kind+'nonuniform all-row constraints close exactly')
    rejects('Nonpositive allocation cost rejected',lambda:bl.recommend_allocations(F(0),'Q',[0]*17))
    rejects('Full reserve rejected',lambda:bl.recommend_allocations(F(0),'Q',reserve=1))
    need(bl.recommend_allocations(F(1),'Q')['allocation_exists'] is False,'Insufficient center reserve not reported as feasible allocation')
    regions=[{'id':'c0','area':'1','multiplicity':1,'relative_cost':'1'},
             {'id':'c1','area':'1','multiplicity':2,'relative_cost':'3'}]
    for kind in ('g','Q'):
        allocated=bl.allocate_regions(regions,[F(1,10000)]*17,kind)
        need(all(sum(x['multiplicity']*x['moment_modulus_radii'][k] for x in allocated)==F(1,10000) for k in range(17)),kind+'weighted region allocation counts all3physicalcharts')
        need(allocated[1]['component_density_tolerances'][0]==3*allocated[0]['component_density_tolerances'][0],kind+'costweighted region density ratio retained')
    gx=bl.allocate_regions(regions,[F(1,10000)]*17,'g')[0]
    need(2*gx['component_density_tolerances'][0]**2<=gx['moment_modulus_radii'][0]**2,'Green region component allocation retains both complex components')
    rejects('Duplicate allocation region rejected',lambda:bl.allocate_regions([regions[0],regions[0]],[F(1,10)]*17,'Q'))
    return {'status':'passed','check_count':len(checks),'checks':checks,
            'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in ('budget_ledger.py','matrix_budget.py','check_budget_ledger.py')}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=run();a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
