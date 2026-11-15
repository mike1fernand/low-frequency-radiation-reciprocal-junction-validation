"""Audit first eight complete R4 production outputs against canonical sources."""
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'integrator'))
import compact_controller as cc
from flint import ctx


def main():
    jobs={j['id']:j for j in json.loads((ROOT/'CAMPAIGN.json').read_text())['jobs']}
    selected=[]
    for line in (ROOT/'campaign/RUN_LOG.jsonl').read_text().splitlines():
        item=json.loads(line)
        if item.get('stage_summary',{}).get('status')=='rigorous_complete_tile_enclosure' and item['id'] not in selected:
            selected.append(item['id'])
        if len(selected)==8:break
    assert len(selected)==8
    current=cc.identity()
    records=[]
    for name in selected:
        job=jobs[name]
        path=ROOT/job['state']
        state=json.loads(path.read_text())
        assert state['identity']==current
        assert state['settings']==cc.normalized_settings(json.loads((ROOT/job['settings']).read_text()))
        ctx.prec=state['settings']['bits']
        t,certificate=cc.endpoint_certificate(state['settings']['bounds'],state['settings']['chart'],state['settings']['endpoint_candidates'])
        assert certificate==state['endpoint_certificate']
        cc.audit_state(state)
        assert not state['pending_ids']
        expected=[[F(0)]*4 for _ in range(17)]
        for stage in state['stages']:
            for m,row in enumerate(stage['increment_endpoints']):
                for c,(mantissa,exponent) in enumerate(row):
                    expected[m][c]+=F(int(mantissa))*F(2)**int(exponent)
        actual=[[F(int(m))*F(2)**int(e) for m,e in row] for row in state['aggregate_endpoints']]
        assert actual==expected
        records.append({'id':name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                        'accepted_leaves':len(state['accepted_ids']),
                        'legs':state['leg_count'],
                        'normalized_maximum_radius_diagnostic':state['normalized_maximum_radius_diagnostic']})
    # Integrand and Cauchy proof code are exact retained R3 bytes.
    baseline=ROOT.parents[1]/'output_r3/paper1-complete-domain-validation-r3/quadrature/preconditioned'
    unchanged=[]
    for name in ('arbitrary_tile_sharp.py','certified_tensor_sharp.py','source_moments.py'):
        assert (ROOT/'integrator'/name).read_bytes()==(baseline/name).read_bytes()
        unchanged.append(name)
    out={'status':'passed','production_states':records,'unchanged_R3_primitives':unchanged,
         'checks':'Settings, runtime/source identities, tile endpoint certificates, compact split-tree coverage, public moments/radii, and independent exact stage sums checked.',
         'interpretation':'These are complete local integral enclosures for the unchanged reference model. Global matrix budgets and the bridge are not tested by this local audit.'}
    Path(__file__).with_name('FIRST_PRODUCTION_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))


if __name__=='__main__':main()
