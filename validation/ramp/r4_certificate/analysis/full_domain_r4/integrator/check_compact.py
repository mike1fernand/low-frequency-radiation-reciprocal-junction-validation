#!/usr/bin/env python3
"""Meaningful compact-controller checks; no large global campaign required."""
import argparse
import copy
import hashlib
import json
import tempfile
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx
import compact_controller as c
from arbitrary_tile_sharp import decode_ball
from source_moments import ball_serialization


def run(reference, output):
    checks = []
    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)
    def seal(state):
        state.pop('state_sha256',None)
        state['state_sha256']=c.digest(state)
    def reject(name, function):
        try:
            function()
        except (ValueError,KeyError):
            checks.append(name)
        else:
            raise AssertionError('unexpected acceptance: '+name)

    ctx.prec=128
    settings={'bounds':['31/32','1','31/32','1'],'chart':1,'mode':'qt','height':'1/4',
              'orders':[5,5,5],'density_tolerance':'1/25000','bits':128}
    with tempfile.TemporaryDirectory(prefix='paper1-compact-check-') as td:
        root=Path(td)
        partial=c.integrate_stage(settings,root/'resumed.json',0.05)
        check('deliberately interrupted stage remains partial', bool(partial['pending_ids']))
        reject('partial state cannot supply completed moments', lambda:c.completed_moments(partial))
        resumed=c.integrate_stage(settings,root/'resumed.json',20)
        if resumed['pending_ids']:
            resumed=c.integrate_stage(settings,root/'resumed.json',20)
        check('resumed complete body cell',not resumed['pending_ids'])
        direct=c.integrate_stage(settings,root/'direct.json',20)
        if direct['pending_ids']:
            direct=c.integrate_stage(settings,root/'direct.json',20)
        check('uninterrupted complete body cell',not direct['pending_ids'])
        check('restart preserves exact endpoint sums',resumed['aggregate_endpoints']==direct['aggregate_endpoints'])
        check('restart preserves exact split tree',resumed['splits']==direct['splits'])
        check('restart preserves accepted leaves',resumed['accepted_ids']==direct['accepted_ids'])
        again=c.integrate_stage(settings,root/'resumed.json',0.01)
        check('complete replay does not recompute',again==resumed)

        original=json.loads(Path(reference).read_text())
        check('reference is complete same physical cell',original['status']=='rigorous_complete_tile_enclosure'
              and original['settings']==settings)
        old_sum=c.zero_sum()
        for leaf in original['accepted']:
            old_sum=c.add_sums(old_sum,c.ball_endpoints([decode_ball(v) for v in leaf['moments']]))
        # Decode endpoints directly also: reading a ball can widen it slightly.
        old_direct=c.zero_sum()
        for leaf in original['accepted']:
            exact_rows=[]
            for moment in leaf['moments']:
                exact_rows.append([c.normalize_endpoint(moment[part][end]['mantissa'],moment[part][end]['exponent2'])
                                  for part,end in [('real','lower'),('real','upper'),('imag','lower'),('imag','upper')]])
            old_direct=c.add_sums(old_direct,exact_rows)
        raw_fraction=[[sum((c.endpoint_fraction([leaf['moments'][m][part][end]['mantissa'],leaf['moments'][m][part][end]['exponent2']]) for leaf in original['accepted']),Fraction())
                       for part,end in [('real','lower'),('real','upper'),('imag','lower'),('imag','upper')]] for m in range(17)]
        check('compact exact sums equal independent Fraction sum of all R3 leaves',
              [[c.endpoint_fraction(e) for e in row] for row in old_direct]==raw_fraction)
        old_balls=c.aggregate_balls(old_direct)
        new_balls=c.aggregate_balls(resumed['aggregate_endpoints'])
        check('all 17 fresh whole-cell qT intervals overlap completed R3 intervals',
              all(a.imag.overlaps(b.imag) for a,b in zip(old_balls,new_balls)))
        check('fresh whole-cell actual density radii meet requested scalar target',
              max(resumed['normalized_radii_diagnostic'])<float(Fraction(settings['density_tolerance'])))

        for name, mutation in [
            ('missing pending leaf', lambda x:x['pending_ids'].pop()),
            ('missing contour leg', lambda x:x.update(leg_count=2)),
            ('duplicate accepted leaf',lambda x:x['accepted_ids'].append(x['accepted_ids'][0])),
            ('corrupt exact aggregate',lambda x:x['aggregate_endpoints'][0][0].__setitem__(0,'123')),
            ('corrupt stage increment',lambda x:x['stages'][0]['increment_endpoints'][0][0].__setitem__(0,'123')),
        ]:
            test=copy.deepcopy(partial if name=='missing pending leaf' else resumed)
            mutation(test);seal(test)
            reject(name,lambda test=test:c.audit_state(test))
        test=copy.deepcopy(resumed);test['moments'][0]['lower']['mantissa']='1';seal(test)
        reject('public moments inconsistent with aggregate',lambda:c.audit_state(test))
        test=copy.deepcopy(resumed);test['maximum_radius_diagnostic']=0;seal(test)
        reject('false zero reported radius',lambda:c.audit_state(test))

        for name, identity_key in [('changed source rejected','sources'),('changed runtime rejected','python')]:
            test=copy.deepcopy(partial)
            test['identity'][identity_key]={} if identity_key=='sources' else 'different runtime'
            seal(test);c.atomic_json(root/'mismatch.json',test)
            reject(name,lambda:c.integrate_stage(settings,root/'mismatch.json',0.01))
        changed=dict(settings);changed['height']='1/8'
        reject('changed settings rejected',lambda:c.integrate_stage(changed,root/'resumed.json',0.01))
        reject('zero tolerance rejected',lambda:c.normalized_settings(dict(settings,density_tolerance='0')))
        bad=dict(settings);bad.pop('density_tolerance');bad['density_tolerances']=['1']*16
        reject('wrong vector length rejected',lambda:c.normalized_settings(bad))
        green_settings=dict(settings,mode='green',bounds=['0','1/1024','0','1/1024'],density_tolerance='1/100')
        green=c.integrate_stage(green_settings,root/'green.json',0.01)
        check('Green initializes all four contour legs including exterior',green['leg_count']==4
              and any(n.startswith('3:') for n in green['pending_ids']+green['accepted_ids']))
        broken=copy.deepcopy(green)
        broken['pending_ids']=[x for x in broken['pending_ids'] if not x.startswith('3:')]
        broken['accepted_ids']=[x for x in broken['accepted_ids'] if not x.startswith('3:')]
        seal(broken)
        reject('omitted Green exterior rejected',lambda:c.audit_state(broken))

        # Heterogeneous budgets exercised on a complete cell, not merely parsed.
        vector=dict(settings,bounds=['31/32','1','31/32','1'])
        vector.pop('density_tolerance')
        vector['density_tolerances']=['1/25000' if m%2==0 else '1/12500' for m in range(17)]
        vs=c.integrate_stage(vector,root/'vector.json',20)
        if vs['pending_ids']:
            vs=c.integrate_stage(vector,root/'vector.json',20)
        check('heterogeneous vector budgets complete all contour legs',not vs['pending_ids'])
        check('all17 vector actual radii meet respective requested targets',
              all(r<float(Fraction(t)) for r,t in zip(vs['normalized_radii_diagnostic'],vector['density_tolerances'])))
        check('vector and scalar whole-cell enclosures overlap',
              all(a.overlaps(b) for a,b in zip(c.completed_moments(vs),c.completed_moments(resumed))))
        result={'status':'passed','checks':checks,'check_count':len(checks),
                'reference_sha256':hashlib.sha256(Path(reference).read_bytes()).hexdigest(),
                'source_identity':c.identity(),
                'fresh_scalar_summary':c.summary(resumed),'fresh_vector_summary':c.summary(vs),
                'old_reference_bytes':Path(reference).stat().st_size,
                'fresh_compact_bytes':(root/'resumed.json').stat().st_size,
                'roundoff_and_error_note':'Acceptance diagnostics use complete resulting outward interval widths, not only analytic remainder targets.'}
        Path(output).write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'status':result['status'],'check_count':len(checks),
                          'fresh_compact_bytes':result['fresh_compact_bytes'],'old_reference_bytes':result['old_reference_bytes']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.reference,a.output)
