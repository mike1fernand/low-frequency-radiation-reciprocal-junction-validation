"""Adapter-specific negative controls, using the completed body-diagonal tile."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import sys
import tempfile

import campaign_ledger as cl


def run(root):
    sys.path.insert(0,str(root/'integrator'))
    import compact_controller as cc
    sample=root/'integrator/BODY_SAMPLE.json'
    source=json.loads(sample.read_text());settings=source['settings']
    checks=[]
    def need(ok,name):
        if not ok:raise ArithmeticError(name)
        checks.append(name)
    def reject(name,call):
        try:call()
        except (ValueError,KeyError,ArithmeticError):checks.append(name)
        else:raise ArithmeticError(name+' was accepted')
    values,state,binding=cl.read_compact(sample,settings,cc)
    need(values is not None and len(values)==17,'Complete compact sample accepted')
    need(binding==hashlib.sha256(sample.read_bytes()).hexdigest(),'Hash identifies same byte snapshot used for audit')
    for k,row in enumerate(state['aggregate_endpoints']):
        need(values[k][0]==list(map(cc.endpoint_fraction,row[2:])),'Q moment'+str(k)+' uses canonical imaginary endpoints')
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'changed.json'
        def edited(mutator):
            obj=copy.deepcopy(source);mutator(obj)
            obj['state_sha256']=cc.digest({k:v for k,v in obj.items() if k!='state_sha256'})
            path.write_text(json.dumps(obj))
            return lambda:cl.read_compact(path,settings,cc)
        reject('Producer change rejected despite recomputed checkpoint hash',edited(lambda x:x['identity']['sources'].update({'source_moments.py':'0'*64})))
        changed_python=edited(lambda x:x['identity'].__setitem__('python',x['identity']['python']+' [different audited build]'))
        values,changed,changed_binding=changed_python()
        need(values is not None and changed['identity']['python']==source['identity']['python']+' [different audited build]',
             'Different producer Python build accepted for audit with original identity preserved')
        reject('Strict computation resume rejects changed producer Python build',lambda:cc.integrate_stage(settings,path,max_seconds=0.01))
        reject('Changed Flint version rejected despite recomputed checkpoint hash',edited(lambda x:x['identity'].__setitem__('flint_version','unreviewed-version')))
        reject('Endpoint certificate change rejected despite recomputed checkpoint hash',edited(lambda x:x['endpoint_certificate'].update({'fake':True})))
        reject('Canonical endpoint change rejected despite recomputed checkpoint hash',edited(lambda x:x['aggregate_endpoints'][0][0].__setitem__(0,'12345')))
        reject('Settings change rejected despite recomputed checkpoint hash',edited(lambda x:x['settings'].__setitem__('height','1/8')))
    s={'chart':0,'bounds':['0','1/4','1/2','3/4']}
    entry={'chart':0,'bounds':['1/2','3/4','0','1/4'],'transform':'swap_uv_chart0'}
    need(cl.derive_physical(s,entry)['bounds']==entry['bounds'],'Permitted chart0 reflected rectangle accepted')
    reject('Wrong reflected rectangle rejected',lambda:cl.derive_physical(s,dict(entry,bounds=s['bounds'])))
    reject('Chart1 reflected u/v rejected',lambda:cl.derive_physical(dict(s,chart=1),entry))
    entry2={'chart':2,'bounds':s['bounds'],'transform':'chart1_to_chart2'}
    need(cl.derive_physical(dict(s,chart=1),entry2)['chart']==2,'Permitted chart1to2 identity accepted')
    reject('Chart0to2 unsupported physical exchange rejected',lambda:cl.derive_physical(s,entry2))
    return {'status':'passed','check_count':len(checks),'checks':checks,
            'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in ('campaign_ledger.py','inherited_reader.py','check_campaign_adapter.py')},
            'sample_sha256':binding}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=run(a.root.resolve());a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
