#!/usr/bin/env python3
"""Bounded audit of shipped evidence; no fresh global quadrature by default.

An unresolved but correctly audited budget is a successful reproduction of an
unresolved result. It must never be reported as an exact-lattice certificate.
"""
import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from bootstrap import PACKAGE_ROOT, bootstrap, sha, verify_manifest


def runtime_audit_record(producer_runtimes, audit_runtime):
    """Audit saved evidence without requiring the producer's Python build.

    Source identity, precision, certificates, exact sums and coverage remain
    checked by campaign_ledger. Integration resume keeps the original stricter
    runtime comparison in the unchanged compact controller.
    """
    if audit_runtime['python_flint']!='0.8.0' or audit_runtime['flint']!='3.3.1':
        raise ValueError('Evidence audit requires reviewed python-flint 0.8.0 / FLINT 3.3.1')
    if any(x['python_flint']!=audit_runtime['python_flint'] for x in producer_runtimes):
        raise ValueError('Saved checkpoint python-flint version differs from reviewed audit runtime')
    return {'audit_runtime':audit_runtime,'checkpoint_producer_runtimes':producer_runtimes,
            'python_build_matches_all_producers':all(x['python']==audit_runtime['python'] for x in producer_runtimes),
            'runtime_policy':'Read-only saved-evidence audit permits a different Python build; unchanged source/Flint and all certificate, exact-sum and coverage checks remain required. Checkpoint resume retains exact producer-runtime identity.'}


def require_equal_ledgers(current, saved):
    """Compare every scientific-ledger field; runtime reports live separately."""
    if current!=saved:
        raise ValueError('fresh audited ledger differs from selected saved result')


def require_bridge_gate(ledger):
    if ledger.get('both_arithmetic_budgets_pass') is not True:
        raise ValueError('Bridge requires both actual global arithmetic budgets to pass')
    for kind,limit in [('g',Fraction(1,50)),('Q',Fraction(1,200))]:
        result=ledger['results'][kind]
        gate=result['operator_budget']
        if (gate.get('arithmetic_budget_pass') is not True or
                Fraction(gate['threshold'])!=limit or
                Fraction(gate['scaled_operator_upper'])>limit):
            raise ValueError('Bridge actual exact operator budget does not pass: '+kind)
        if Fraction(result['fallback_parameter_area'])!=0:
            raise ValueError('Bridge requires zero fallback parameter area: '+kind)


def execute_stage(name, command, destination, report, save, timeout):
    """Record success, failed process, and timeout before returning/raising."""
    started=time.monotonic()
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    try:
        result=subprocess.run([sys.executable,'-B',*command],capture_output=True,text=True,
                              env=env,timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        report['stages'].append({'name':name,'returncode':None,
                                'wall_seconds':time.monotonic()-started,'status':'timeout'})
        report['status']='failed'
        report['error']='replay stage timed out: '+name
        (destination/(name+'.stderr.txt')).write_text(str(exc))
        save()
        raise
    (destination/(name+'.stdout.txt')).write_text(result.stdout)
    (destination/(name+'.stderr.txt')).write_text(result.stderr)
    record={'name':name,'returncode':result.returncode,
            'wall_seconds':time.monotonic()-started,'status':'passed' if result.returncode==0 else 'failed'}
    report['stages'].append(record)
    if result.returncode:
        report['status']='failed'
        report['error']='replay stage failed: '+name
    save()
    if result.returncode:
        raise ValueError('replay stage failed: '+name)


def validate_bridge_result(root, directory, ledger):
    """Bind the successful original bridge to this exact audited global gate."""
    require_bridge_gate(ledger)
    summary_path=directory/'EXACT_LATTICE_BRIDGE.json'
    summary=json.loads(summary_path.read_text())
    if (summary.get('status')!='exact_reference_lattice_floor_certified_from_complete_enclosures_and_original_bridge' or
            summary.get('exact_lattice_floor_certified') is not True or
            summary.get('rational_floor')!='1037/10000' or
            summary.get('fresh_original_process_exit')!=0):
        raise ValueError('Malformed or unsuccessful exact-lattice bridge summary')
    gate_path=directory/'GLOBAL_GATE.json'
    require_equal_ledgers(json.loads(gate_path.read_text()),ledger)
    proof_path=directory/'finite_proof/finite_kernel_certificate.json'
    proof=json.loads(proof_path.read_text())
    bridge=proof['hypothetical_kernel_error_conclusion']
    if (proof['status']!='finite_cached_kernel_certificate_passed' or
            proof['certificate_revision']!='basis_invariant_n192_sobol23_m256_v1' or
            proof['finite_rational_floor']!='27/250' or len(proof['vertices'])!=64 or
            bridge['eG']!='1/50' or bridge['eQ']!='1/200' or
            bridge['rational_floor']!='1037/10000' or
            summary['original_conditional_floor_enclosure']!=bridge['floor_enclosure']):
        raise ValueError('Fresh original proof violates exact bridge contract')
    from flint import arb
    if not arb(bridge['floor_enclosure'])>arb(1037)/10000:
        raise ValueError('Original rigorous bridge enclosure does not prove strict stated floor')
    package=root/'budget'
    inputs=package/'bridge_inputs'
    readiness=json.loads((package/'BRIDGE_READINESS.json').read_text())
    expected_sources={key[len('staged_repo/'):]:value for key,value in readiness['sha256'].items()
                      if key.startswith('staged_repo/')}
    if proof['source_sha256']!=expected_sources:
        raise ValueError('Original bridge selected-cache identities differ')
    baseline=(root/json.loads((root/'CAMPAIGN.json').read_text())['baseline']).resolve()
    expected={'global_gate_sha256':sha(gate_path),'fresh_original_proof_sha256':sha(proof_path),
              'original_producer_sha256':sha(inputs/'finite_kernel_certificate.py'),
              'exact_Q_psd_proof_sha256':sha(baseline/'math/COMPLETE_FALLBACK_PROOF.md'),
              'targets_sha256':sha(package/'TARGETS.json')}
    if summary.get('bindings')!=expected:
        raise ValueError('Bridge gate/proof/source/PSD/target binding mismatch')
    return summary


def replay_optional_bridge(root, destination, ledger, report, save, timeout):
    """Keep every prerequisite, process and output failure terminally failed."""
    try:
        require_bridge_gate(ledger)
        execute_stage('exact_lattice_bridge',[str(root/'budget/replay_bridge.py'),
                      '--root',str(root),'--output-dir',str(destination/'bridge')],
                      destination,report,save,min(timeout,55))
        summary=validate_bridge_result(root,destination/'bridge',ledger)
        report['bridge_summary_sha256']=sha(destination/'bridge/EXACT_LATTICE_BRIDGE.json')
        report['bridge']=summary
        report['exact_lattice_floor_certified']=True
        report['rational_floor']=summary['rational_floor']
        save()
    except Exception as exc:
        report['status']='failed'
        report['exact_lattice_floor_certified']=False
        report['error']=type(exc).__name__+': '+str(exc)
        save()
        raise


def run(args):
    package=args.package_root.resolve()
    root=package/'analysis/full_domain_r4'
    destination=args.output.absolute()
    if destination==package or package in destination.parents:
        raise ValueError('verification output must be outside the sealed package')
    if destination.exists() and any(destination.iterdir()):
        raise ValueError('verification output directory must be new or empty')
    destination.mkdir(parents=True,exist_ok=True)
    report={'schema':'paper1-saved-evidence-replay-r4.v1','status':'in_progress',
            'fresh_global_quadrature_performed':False,'stages':[]}
    def save():
        (destination/'VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    try:
        report['outer_manifest']=verify_manifest(package,args.manifest or package/'MANIFEST.json')
        report['inherited']=bootstrap(package)
        import flint
        audit_runtime={'python':sys.version,'python_flint':flint.__version__,
                       'flint':flint.__FLINT_VERSION__}
        campaign=json.loads((root/'CAMPAIGN.json').read_text())
        identities=[]
        for job in campaign['jobs']:
            state=root/job['state']
            if state.is_file():
                identity=json.loads(state.read_text())['identity']
                expected={'python':identity['python'],'python_flint':identity['flint_version']}
                if expected not in identities:
                    identities.append(expected)
        report.update(runtime_audit_record(identities,audit_runtime))
        save()
        stages=[('budget_controls',[str(root/'budget/check_budget_ledger.py'),
                                    '--output',str(destination/'budget_checks.json')]),
                ('audited_campaign_ledger',[str(root/'budget/campaign_ledger.py'),
                                             '--root',str(root),'--output',str(destination/'ledger.json')])]
        if args.physical_check:
            stages.append(('compact_physical_controls',[str(root/'integrator/check_compact.py'),
                           '--reference',str(package/'output_r3/paper1-complete-domain-validation-r3/quadrature/qt_chart1_body_preconditioned_gl5.json'),
                           '--output',str(destination/'compact_checks.json')]))
        for name,command in stages:
            execute_stage(name,command,destination,report,save,args.stage_timeout)
        ledger_path=destination/'ledger.json'
        ledger=json.loads(ledger_path.read_text())
        report['ledger_sha256']=sha(ledger_path)
        report['both_arithmetic_budgets_pass']=ledger['both_arithmetic_budgets_pass']
        report['exact_lattice_floor_certified']=False
        report['physical_replay_performed']=args.physical_check
        if args.expected_ledger:
            expected=args.expected_ledger
            if not expected.is_absolute():
                expected=package/expected
            require_equal_ledgers(ledger,json.loads(expected.read_text()))
            report['saved_ledger_reproduced_exactly']=True
            report['selected_saved_ledger_sha256']=sha(expected)
        else:
            report['saved_ledger_reproduced_exactly']=None
            report['comparison_note']='No selected saved-ledger comparison requested; current shipped states were audited and aggregated.'
        if args.bridge:
            replay_optional_bridge(root,destination,ledger,report,save,args.stage_timeout)
        report['status']='passed'
        report['scientific_result']=('exact_reference_lattice_floor_certified' if report['exact_lattice_floor_certified'] else
                                     'arithmetic_budgets_pass_floor_requires_separate_bridge_and_PSD' if ledger['both_arithmetic_budgets_pass'] else
                                     'arithmetic_budgets_unresolved')
        report['limits']='Saved evidence audit and exact budget recomputation; optional bridge only if all prerequisites and output bindings pass. No global quadrature rerun, second-platform validation, packet-average or frequency-band guarantee.'
        save()
        print(json.dumps({'status':'passed','scientific_result':report['scientific_result'],
                          'saved_ledger_reproduced_exactly':report['saved_ledger_reproduced_exactly']}))
        return 0
    except Exception as exc:
        report['status']='failed'
        report['error']=type(exc).__name__+': '+str(exc)
        save()
        print(report['error'],file=sys.stderr)
        return 2


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package-root',type=Path,default=PACKAGE_ROOT)
    p.add_argument('--manifest',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--expected-ledger',type=Path)
    p.add_argument('--physical-check',action='store_true')
    p.add_argument('--bridge',action='store_true',help='after zero-fallback budget success, rerun and validate the original exact-lattice bridge')
    p.add_argument('--stage-timeout',type=float,default=55)
    a=p.parse_args()
    if not 0<a.stage_timeout<=60:
        p.error('stage timeout must be positive and <=60 seconds')
    raise SystemExit(run(a))
