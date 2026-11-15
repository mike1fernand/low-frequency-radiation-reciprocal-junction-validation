#!/usr/bin/env python3
"""Negative controls for optional bridge gating and process/result reporting.

These controls use manufactured ledgers and deliberately unsuccessful helper
programs. They do not run or claim a successful scientific bridge.
"""
import argparse
import copy
import json
from pathlib import Path
import tempfile
from verify_saved import replay_optional_bridge, require_bridge_gate


def run(output):
    checks=[]
    good={'both_arithmetic_budgets_pass':True,'results':{
        'g':{'fallback_parameter_area':'0','operator_budget':{'arithmetic_budget_pass':True,'threshold':'1/50','scaled_operator_upper':'1/100'}},
        'Q':{'fallback_parameter_area':'0','operator_budget':{'arithmetic_budget_pass':True,'threshold':'1/200','scaled_operator_upper':'1/400'}}}}
    require_bridge_gate(good)
    checks.append('manufactured zero-fallback exact budget prerequisites accepted')
    with tempfile.TemporaryDirectory(prefix='paper1-bridge-runner-') as temporary:
        base=Path(temporary)
        root=base/'analysis/full_domain_r4';(root/'budget').mkdir(parents=True)
        source=root/'budget/replay_bridge.py'
        def rejection(name, ledger, script):
            destination=base/name;destination.mkdir()
            source.write_text(script)
            report={'status':'in_progress','stages':[],'exact_lattice_floor_certified':False}
            def save():
                (destination/'VERIFICATION.json').write_text(json.dumps(report))
            try:replay_optional_bridge(root,destination,ledger,report,save,5)
            except (ValueError,KeyError,FileNotFoundError):pass
            else:raise AssertionError(name+' unexpectedly passed')
            saved=json.loads((destination/'VERIFICATION.json').read_text())
            assert saved['status']=='failed' and saved['exact_lattice_floor_certified'] is False
            checks.append(name+' leaves terminal failed status without floor promotion')
            return saved
        unavailable=copy.deepcopy(good);unavailable['both_arithmetic_budgets_pass']=False
        state=rejection('unavailable_budget',unavailable,'raise AssertionError("must not run")\n')
        assert state['stages']==[]
        checks.append('unavailable budget rejected before process launch')
        fallback=copy.deepcopy(good);fallback['results']['g']['fallback_parameter_area']='1/1000'
        rejection('nonzero_fallback',fallback,'raise AssertionError("must not run")\n')
        false_boolean=copy.deepcopy(good)
        false_boolean['results']['Q']['operator_budget']['scaled_operator_upper']='1/100'
        rejection('false_budget_boolean',false_boolean,'raise AssertionError("must not run")\n')
        failed=rejection('failed_process',good,'import sys\nsys.exit(7)\n')
        assert failed['stages'][-1]['status']=='failed' and failed['stages'][-1]['returncode']==7
        checks.append('actual failing subprocess exit code 7 retained')
        malformed='''import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root');p.add_argument('--output-dir');a=p.parse_args()
d=Path(a.output_dir);d.mkdir(parents=True)
(d/'EXACT_LATTICE_BRIDGE.json').write_text(json.dumps({'status':'pretend_success','exact_lattice_floor_certified':True,'rational_floor':'1037/10000'}))
'''
        bad=rejection('malformed_success_result',good,malformed)
        assert bad['stages'][-1]['status']=='passed' and bad['status']=='failed'
        checks.append('zero process exit does not override malformed bridge contract')
    result={'status':'passed','checks':checks,'check_count':len(checks),
            'successful_scientific_bridge_executed':False,
            'scope':'Manufactured prerequisites and deliberately unsuccessful subprocess/results; no positive lattice-floor claim.'}
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','check_count':len(checks)}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.output)
