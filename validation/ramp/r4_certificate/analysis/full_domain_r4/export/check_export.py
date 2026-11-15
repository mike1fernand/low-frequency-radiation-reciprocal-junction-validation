#!/usr/bin/env python3
"""Test portable archive bootstrap and manifest rejection controls."""
import argparse
import copy
import json
from pathlib import Path
import tempfile
import time
from bootstrap import bootstrap, safe_relative, sha, verify_manifest
from verify_saved import runtime_audit_record, require_equal_ledgers


def run(archive, output):
    checks=[]
    def rejected(name, fn):
        try:fn()
        except (ValueError,FileNotFoundError):checks.append(name)
        else:raise AssertionError(name)
    for bad in ('../x','/absolute','x/../../y','x\\y',''):
        rejected('reject unsafe path '+repr(bad),lambda bad=bad:safe_relative(bad))
    producer=[{'python':'3.12.14 (producer compiler)','python_flint':'0.8.0'}]
    auditor={'python':'3.12.14 (different compiler)','python_flint':'0.8.0','flint':'3.3.1'}
    before=copy.deepcopy(producer)
    record=runtime_audit_record(producer,auditor)
    assert not record['python_build_matches_all_producers'] and record['checkpoint_producer_runtimes']==before
    assert producer==before
    checks.append('read-only audit records differing Python build without rewriting producer identity')
    rejected('unreviewed python-flint audit version rejected',lambda:runtime_audit_record(producer,dict(auditor,python_flint='0.9.0')))
    rejected('unreviewed FLINT library audit version rejected',lambda:runtime_audit_record(producer,dict(auditor,flint='3.4.0')))
    rejected('producer python-flint mismatch rejected',lambda:runtime_audit_record([dict(producer[0],python_flint='0.7.0')],auditor))
    ledger={'results':{'Q':{'upper':'1/200'}},'job_audits':[{'state_sha256':'abc'}],
            'input_sha256':{'adapter':'source','targets':'target'}}
    original_report={'ledger':copy.deepcopy(ledger),'audit_runtime':dict(auditor,python=producer[0]['python'])}
    foreign_report={'ledger':copy.deepcopy(ledger),'audit_runtime':auditor}
    require_equal_ledgers(foreign_report['ledger'],original_report['ledger'])
    checks.append('separate audit-runtime change leaves every scientific ledger field comparable')
    for name,change in [
        ('changed mathematical bound rejected',lambda x:x['results']['Q'].update(upper='1/100')),
        ('changed state identity rejected',lambda x:x['job_audits'][0].update(state_sha256='def')),
        ('changed source identity rejected',lambda x:x['input_sha256'].update(adapter='new source')),
        ('changed selected target identity rejected',lambda x:x['input_sha256'].update(targets='new target'))]:
        changed=copy.deepcopy(ledger);change(changed)
        rejected(name,lambda changed=changed:require_equal_ledgers(changed,ledger))
    with tempfile.TemporaryDirectory(prefix='paper1-export-check-') as temporary:
        root=Path(temporary)
        payload=root/'payload.txt';payload.write_text('bound source\n')
        manifest=root/'MANIFEST.json'
        manifest.write_text(json.dumps({'files':{'payload.txt':sha(payload)}}))
        assert verify_manifest(root,manifest)['verified_payloads']==1
        checks.append('correct manifest payload accepted')
        payload.write_text('changed source\n')
        rejected('changed manifest payload rejected',lambda:verify_manifest(root,manifest))
        manifest.write_text(json.dumps({'files':{'../outside':'0'}}))
        rejected('traversal manifest key rejected',lambda:verify_manifest(root,manifest))
        start=time.monotonic()
        extracted=bootstrap(root,archive)
        checks.append('real immutable archive extracted and all inherited manifest payloads verified')
        reused=bootstrap(root,archive)
        assert reused==extracted
        checks.append('already materialized immutable baseline verified without re-extraction')
        baseline=Path(extracted['baseline'])
        target=baseline/'README.md'
        if not target.exists():raise AssertionError('expected R3 README missing')
        target.write_text(target.read_text()+'\nchanged\n')
        rejected('changed extracted baseline payload rejected',lambda:bootstrap(root,archive))
        report={'status':'passed','checks':checks,'check_count':len(checks),
                'archive_sha256':sha(archive),'verified_inherited_payloads':extracted['verified_payloads'],
                'bootstrap_and_recheck_seconds':time.monotonic()-start,
                'runtime_policy_test_scope':'Synthetic producer/auditor build labels exercise policy only; no second Python build or platform was executed.',
                'limits':'Tests bootstrap and manifest guards. A complete outer-package cold replay and current campaign ledger comparison are separate final export checks.'}
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','checks':len(checks)}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.archive,a.output)
