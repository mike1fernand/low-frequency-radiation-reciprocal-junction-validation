"""Audit complete R4 checkpoints, retain R3 evidence, and account for every gap.

Canonical compact endpoint sums are used; redundant display moments never
substitute for the saved exact stage sums. Partial tiles receive full-region
analytic fallback until completed. Final norm gates use all actual targets.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys

import matrix_budget as mb
import budget_ledger as bl


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local_path(root, relative):
    value=Path(relative)
    if value.is_absolute() or '..' in value.parts:
        raise ValueError('Campaign inputs require paths within campaign root')
    result=(root/value).resolve()
    if root not in result.parents or not result.is_file():
        raise ValueError('Campaign input missing or outside root')
    return result


def derive_physical(settings, entry):
    chart=settings['chart'];bounds=list(map(str,map(F,settings['bounds'])))
    transform=entry.get('transform')
    if transform=='swap_uv_chart0':
        if chart!=0:raise ValueError('Swap u/v only valid in chart0')
        bounds=bounds[2:]+bounds[:2]
    elif transform=='chart1_to_chart2':
        if chart!=1:raise ValueError('Chart2 reuse requires chart1 source')
        chart=2
    elif transform is not None:raise ValueError('Unknown symmetry transform')
    if type(entry['chart']) is not int or entry['chart']!=chart or list(map(F,entry['bounds']))!=list(map(F,bounds)):
        raise ValueError('Physical entry does not match source geometry/symmetry')
    return {'chart':chart,'bounds':bounds}


def exact_interval(endpoint):
    if set(endpoint)!={'mantissa','exponent2'} or type(endpoint['exponent2']) is not int:
        raise ValueError('Invalid fallback exact endpoint')
    return F(int(endpoint['mantissa']))*F(2)**endpoint['exponent2']


def read_compact(path, settings, cc):
    from flint import ctx
    raw=path.read_bytes()
    obj=json.loads(raw)
    snapshot_sha256=hashlib.sha256(raw).hexdigest()
    if obj['settings']!=settings:
        raise ValueError('Saved checkpoint settings differ from campaign settings')
    current_identity=cc.identity()
    recorded_identity=obj['identity']
    if (recorded_identity['sources']!=current_identity['sources'] or
            recorded_identity['flint_version']!=current_identity['flint_version']):
        raise ValueError('Saved checkpoint mathematical source or reviewed Flint version mismatch')
    if not isinstance(recorded_identity['python'],str) or not recorded_identity['python'].strip():
        raise ValueError('Missing recorded producer Python identity')
    old_precision=ctx.prec
    try:
        ctx.prec=settings['bits']
        _,cert=cc.endpoint_certificate(settings['bounds'],settings['chart'],settings['endpoint_candidates'])
        if obj['endpoint_certificate']!=cert:
            raise ValueError('Saved endpoint certificate differs from recomputation')
        cc.audit_state(obj)
    finally:
        ctx.prec=old_precision
    if obj['status']!='rigorous_complete_tile_enclosure':
        return None,obj,snapshot_sha256
    values=[]
    for row in obj['aggregate_endpoints']:
        lo,hi,il,ih=map(cc.endpoint_fraction,row)
        values.append([[il,ih],[F(0),F(0)]] if settings['mode']=='qt' else [[lo,hi],[il,ih]])
    mb.parse_intervals(values)
    return values,obj,snapshot_sha256


def assemble(root, manifest_path):
    root=root.resolve();manifest_path=manifest_path.resolve()
    manifest=json.loads(manifest_path.read_text())
    if manifest['schema']!='paper1-r4-campaign-v1':raise ValueError('Unknown campaign schema')
    baseline=(root/manifest['baseline']).resolve()
    sealed_manifest_path=baseline/'MANIFEST.json'
    sealed=json.loads(sealed_manifest_path.read_text())
    for relative in ('complete_ledger.py','coverage.py','matrix/matrix_budget.py','matrix/TARGETS.json',
                     'math/COMPLETE_FALLBACK_CHECKS.json','math/COMPLETE_FALLBACK_PROOF.md'):
        if sha(baseline/relative)!=sealed['files'][relative]:
            raise ValueError('Inherited auditor/target/fallback file differs from sealed R3 manifest')
    baseline_manifest=baseline/'EVIDENCE_MANIFEST.json'
    original=json.loads(baseline_manifest.read_text())
    if sha(baseline_manifest)!=manifest['baseline_manifest_sha256']:
        raise ValueError('Inherited evidence manifest changed')
    if manifest['inherited_evidence']!=original['evidence']:
        raise ValueError('Inherited records differ from selected R3 manifest')
    response=subprocess.run([sys.executable,'-B',str(Path(__file__).with_name('inherited_reader.py')),
                             '--root',str(baseline)],capture_output=True,text=True,check=True)
    inherited=json.loads(response.stdout)
    if inherited['manifest_sha256']!=manifest['baseline_manifest_sha256']:
        raise ValueError('Inherited audit manifest identity mismatch')
    for rec in inherited['records']:
        rec['values']=mb.parse_intervals(rec['values'])
    coverage_spec=importlib.util.spec_from_file_location('r4_exact_coverage',baseline/'coverage.py')
    coverage=importlib.util.module_from_spec(coverage_spec);coverage_spec.loader.exec_module(coverage)
    sys.path.insert(0,str(root/'integrator'))
    import compact_controller as cc
    jobs=manifest['jobs']
    if len({x['id'] for x in jobs})!=len(jobs) or len({x['state'] for x in jobs})!=len(jobs):
        raise ValueError('Duplicate job or checkpoint identity')
    completed=[];planned=[];job_audits=[];priorities={'g':[],'Q':[]};bindings={}
    inherited_groups={}
    for record in inherited['records']:
        evidence=record['evidence'];key=(evidence['kind'],evidence['path'],evidence['sha256'])
        if key not in inherited_groups:
            inherited_groups[key]={'id':'R3:'+evidence['path'],'intervals':record['values'],'multiplicity':0}
        elif inherited_groups[key]['intervals']!=record['values']:
            raise ValueError('Inherited symmetry copies have inconsistent enclosures')
        inherited_groups[key]['multiplicity']+=1
    for (kind,_,_),record in inherited_groups.items():
        priorities[kind].append(record)
    for job in jobs:
        kind=job['kind']
        if kind not in ('g','Q'):raise ValueError('Invalid job observable')
        settings_path=local_path(root,job['settings'])
        settings=cc.normalized_settings(json.loads(settings_path.read_text()))
        if settings['mode']!={'g':'green','Q':'qt'}[kind]:raise ValueError('Job mode/observable mismatch')
        if type(settings['chart']) is not int or settings['chart'] not in (0,1):
            raise ValueError('Independent jobs must be on chart0 or chart1')
        physical=[derive_physical(settings,e) for e in job['physical_entries']]
        if type(job['multiplicity']) is not int or job['multiplicity']!=len(physical) or not physical:
            raise ValueError('Physical multiplicity mismatch')
        if F(job['area'])!=coverage.area(settings['bounds']):raise ValueError('Job parameter area mismatch')
        for rec in physical:
            coverage.area(rec['bounds']);planned.append(dict(rec,kind=kind))
        state_rel=Path(job['state'])
        if state_rel.is_absolute() or '..' in state_rel.parts:raise ValueError('State outside campaign root')
        state_path=(root/state_rel).resolve()
        if root not in state_path.parents:raise ValueError('State outside campaign root')
        audit={'id':job['id'],'kind':kind,'multiplicity':job['multiplicity'],
               'physical_parameter_area':str(job['multiplicity']*F(job['area'])),
               'settings_sha256':sha(settings_path),'state_sha256':None,'status':'not_started'}
        if state_path.is_file():
            values,state,snapshot_sha256=read_compact(state_path,settings,cc)
            audit['state_sha256']=snapshot_sha256;audit['status']=state['status']
            audit['completed_stages']=len(state['stages'])
            if values is not None:
                for rec in physical:completed.append(dict(rec,values=values,evidence={'kind':kind,'id':job['id']}))
                priorities[kind].append({'id':job['id'],'intervals':values,'multiplicity':job['multiplicity']})
        job_audits.append(audit)
    targets_path=Path(__file__).with_name('TARGETS.json')
    if sha(targets_path)!=sha(baseline/'matrix/TARGETS.json'):
        raise ValueError('Exact selected targets changed')
    targets=json.loads(targets_path.read_text())
    fallback_path=baseline/'math/COMPLETE_FALLBACK_CHECKS.json'
    fallback=json.loads(fallback_path.read_text())
    if fallback['status']!='proved_complete_chart_fallback_bounds':raise ValueError('Fallback proof status mismatch')
    if [x.replace('_','') for x in fallback['moment_order']]!=list(mb.NAMES):raise ValueError('Fallback moment ordering differs')
    results={};complete_moments={}
    records=inherited['records']+completed
    for kind in ('g','Q'):
        inherited_kind=[x for x in inherited['records'] if x['evidence']['kind']==kind]
        plan_cover=coverage.full_cover(inherited_kind+[x for x in planned if x['kind']==kind])
        if any(F(c['fallback_area']) for c in plan_cover.values()):raise ValueError('Campaign plan omits part of complete domain')
        selected=[x for x in records if x['evidence']['kind']==kind]
        cover=coverage.full_cover(selected)
        remainder=sum((F(x['fallback_area']) for x in cover.values()),F(0))
        total=[[[F(0),F(0)],[F(0),F(0)]] for _ in mb.NAMES]
        for record in selected:
            for k,value in enumerate(record['values']):
                for part in range(2):
                    total[k][part]=[total[k][part][j]+value[part][j] for j in range(2)]
        rates=fallback['per_unit_angular_area']
        for k in range(17):
            for part,key in enumerate(['g_real','g_imag'] if kind=='g' else ['qt',None]):
                if key is None:continue
                v=remainder*exact_interval(rates[key][k]['upper'])
                nonnegative=k in (0,6) and (kind=='Q' or part==1)
                total[k][part][0]+=F(0) if nonnegative else -v
                total[k][part][1]+=v
        actual=bl.assess(total,targets['targets'][kind],kind)
        results[kind]={'coverage':cover,'refined_parameter_area':3-remainder,'fallback_parameter_area':remainder,
                       'operator_budget':actual,'refinement_priorities':bl.rank_regions(priorities[kind]),
                       'all_new_jobs_complete':all(x['status']=='rigorous_complete_tile_enclosure' for x in job_audits if x['kind']==kind)}
        complete_moments[kind]=total
    return {'schema':'paper1-complete-domain-budget-r4.v1','results':results,'complete_moments':complete_moments,
            'both_arithmetic_budgets_pass':all(x['operator_budget']['arithmetic_budget_pass'] for x in results.values()),
            'exact_lattice_floor_certified':False,'job_audits':job_audits,
            'input_sha256':{'campaign':sha(manifest_path),'baseline_evidence_manifest':sha(baseline_manifest),
                            'baseline_sealed_manifest':sha(sealed_manifest_path),
                            'targets':sha(targets_path),'fallback':sha(fallback_path),
                            'adapter':sha(Path(__file__))},
            'claim':'Complete-domain enclosures include analytic fallback for all incomplete regions. Floor needs both budgets, exact-Q PSD premise, and bridge replay.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--manifest',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=assemble(a.root,a.manifest or a.root/'CAMPAIGN.json')
    a.output.write_text(json.dumps(mb.json_exact(result),indent=2)+'\n')
    print(json.dumps({'both_budgets_pass':result['both_arithmetic_budgets_pass'],
                     'results':{k:{'remaining_area':str(v['fallback_parameter_area']),
                                   'scaled_upper':float(v['operator_budget']['scaled_operator_upper'])}
                                for k,v in result['results'].items()}}))


if __name__=='__main__':main()
