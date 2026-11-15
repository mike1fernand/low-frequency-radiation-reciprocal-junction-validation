"""Independent R4 sum and direct operator-bound calculation.

The completed tile rectangles are evidence inputs. This script re-adds their
canonical endpoints, reconstructs every matrix entry without the budget
module's assembler, and proves the sufficient direct norm budgets using
integer/rational arithmetic. It does not rerun each interval integral.
"""
import hashlib
import json
import math
import subprocess
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def dyadic(end):return F(int(end[0]))*F(2)**int(end[1])


def upper_sqrt(q):
    assert q>=0
    numerator=q.numerator<<360
    k=math.isqrt(numerator//q.denominator)
    if k*k*q.denominator<numerator:k+=1
    return F(k,1<<180)


def main():
    ledger_path=ROOT/'budget/FINAL_CAMPAIGN_LEDGER.json'
    ledger=json.loads(ledger_path.read_text())
    campaign_path=ROOT/'CAMPAIGN.json'
    campaign=json.loads(campaign_path.read_text())
    assert sha(campaign_path)==ledger['input_sha256']['campaign']
    prior_partition=json.loads(Path(__file__).with_name('CAMPAIGN_CHECKS.json').read_text())
    assert prior_partition['campaign_sha256']==sha(campaign_path)
    assert len(campaign['jobs'])==222 and len(ledger['job_audits'])==222
    assert all(x['status']=='rigorous_complete_tile_enclosure' for x in ledger['job_audits'])
    baseline=(ROOT/campaign['baseline']).resolve()
    inherited=subprocess.run([sys.executable,'-B',str(ROOT/'budget/inherited_reader.py'),
                              '--root',str(baseline)],capture_output=True,text=True,check=True)
    inherited=json.loads(inherited.stdout)['records']
    sums={k:[[[F(0),F(0)],[F(0),F(0)]] for _ in range(17)] for k in ('g','Q')}
    def add(kind,values,multiplicity=1):
        assert len(values)==17
        for m,row in enumerate(values):
            for component,pair in enumerate(row):
                for endpoint,value in enumerate(pair):
                    sums[kind][m][component][endpoint]+=multiplicity*F(value)
    for record in inherited:add(record['evidence']['kind'],record['values'])
    saved_audits={x['id']:x for x in ledger['job_audits']}
    accepted=0
    stages=0
    for job in campaign['jobs']:
        path=ROOT/job['state']
        assert sha(path)==saved_audits[job['id']]['state_sha256']
        state=json.loads(path.read_text())
        assert state['status']=='rigorous_complete_tile_enclosure' and not state['pending_ids']
        assert state['leg_count']==(3 if job['kind']=='Q' else 4)
        accepted+=len(state['accepted_ids']);stages+=len(state['stages'])
        values=[]
        for row in state['aggregate_endpoints']:
            a,b,c,d=map(dyadic,row)
            values.append([[c,d],[F(0),F(0)]] if job['kind']=='Q' else [[a,b],[c,d]])
        add(job['kind'],values,len(job['physical_entries']))
    targets_path=ROOT/'budget/TARGETS.json'
    assert sha(targets_path)==ledger['input_sha256']['targets']
    targets=json.loads(targets_path.read_text())['targets']
    answers={}
    for kind in ('g','Q'):
        expected=[[[F(x) for x in pair] for pair in row] for row in ledger['complete_moments'][kind]]
        assert sums[kind]==expected
        result=ledger['results'][kind]
        assert F(result['fallback_parameter_area'])==0 and F(result['refined_parameter_area'])==3
        assert all(F(c['refined_area'])==1 and F(c['fallback_area'])==0 for c in result['coverage'].values())
        envelope=[]
        for i in range(12):
            row=[]
            for j in range(12):
                segment_i,channel_i=divmod(i,2)
                segment_j,channel_j=divmod(j,2)
                delta=segment_i-segment_j
                distance=abs(delta)
                if channel_i==channel_j:
                    slot=distance+(6 if channel_i else 0);sign=1
                    value=sums[kind][slot]
                elif not delta:
                    value=[[F(0),F(0)],[F(0),F(0)]];sign=1
                else:
                    slot=11+distance
                    sign=(1 if delta>0 else -1)*(1 if channel_i==0 else -1)
                    value=sums[kind][slot]
                target=list(map(F,targets[kind][i][j]))
                squared=F(0)
                for part in range(2):
                    lo,hi=value[part]
                    if sign<0:lo,hi=-hi,-lo
                    bound=max(abs(lo-target[part]),abs(hi-target[part]))
                    squared+=bound*bound
                row.append(upper_sqrt(squared))
            envelope.append(row)
        rowbound=max(sum(row) for row in envelope)
        colbound=max(sum(envelope[i][j] for i in range(12)) for j in range(12))
        frobenius=upper_sqrt(sum(x*x for row in envelope for x in row))
        direct=F(169,100)*min(upper_sqrt(rowbound*colbound),frobenius)
        limit=F(1,50) if kind=='g' else F(1,200)
        assert direct<=limit
        assert direct==F(result['operator_budget']['direct_envelope']['scaled_operator_upper'])
        answers[kind]={'independent_scaled_direct_upper':str(direct),
                       'diagnostic_decimal':float(direct),'threshold':str(limit),'passes':True}
    output={'status':'passed','ledger_sha256':sha(ledger_path),'jobs':222,
            'accepted_leaves':accepted,'completed_stages':stages,'results':answers,
            'claim':'Independent rational summation and144-entry assembly reproduce all34 complete moments and prove both direct operator budgets. The new center-LDL sharpening is not needed for acceptance.',
            'limits':'Inherited evidence parsed by unchanged R3 auditor; interval quadrature itself not rerun. Exact-Q positivity and bridge remain separate reviewed premises.'}
    Path(__file__).with_name('FINAL_GATE_CHECKS.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='results'},indent=2))
    print(json.dumps({k:{x:y for x,y in v.items() if x!='independent_scaled_direct_upper'} for k,v in answers.items()}))


if __name__=='__main__':main()
