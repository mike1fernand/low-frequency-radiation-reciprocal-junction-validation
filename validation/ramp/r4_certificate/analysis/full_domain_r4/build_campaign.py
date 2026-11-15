"""Exact complement of inherited complete tiles, partitioned for resumable work."""
import json,sys,hashlib
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
R3=ROOT.parents[1]/'output_r3/paper1-complete-domain-validation-r3'
sys.path.insert(0,str(R3))
from coverage import rect,area,intersection_area,full_cover

def build():
    inherited=json.loads((R3/'EVIDENCE_MANIFEST.json').read_text())['evidence']
    jobs=[]
    for kind in ['Q','g']:
        for chart in [0,1]:
            old=[rect(e['bounds']) for e in inherited if e['kind']==kind and e['chart']==chart]
            new=[]
            for i in range(8):
                for j in range(8):
                    base=(F(i,8),F(i+1,8),F(j,8),F(j+1,8))
                    cuts=[r for r in old if intersection_area(r,base)>0]
                    xs=sorted({base[0],base[1],*(v for r in cuts for v in r[:2] if base[0]<v<base[1])})
                    ys=sorted({base[2],base[3],*(v for r in cuts for v in r[2:] if base[2]<v<base[3])})
                    for a,b in zip(xs,xs[1:]):
                        for c,d in zip(ys,ys[1:]):
                            r=(a,b,c,d)
                            if not any(intersection_area(r,o)>0 for o in old):new.append(r)
            for index,r in enumerate(sorted(new)):
                mirror=(r[2],r[3],r[0],r[1])
                if chart==0:
                    assert mirror in new,'chart0 partition is not symmetric'
                    if r>mirror:continue
                name=f'{kind}_c{chart}_{index:03d}'
                bounds=list(map(str,r))
                physical=[{'chart':chart,'bounds':bounds}]
                if chart==1:physical.append({'chart':2,'bounds':bounds,'transform':'chart1_to_chart2'})
                elif r!=mirror:physical.append({'chart':0,'bounds':list(map(str,mirror)),'transform':'swap_uv_chart0'})
                settings={'bounds':bounds,'chart':chart,'mode':'qt' if kind=='Q' else 'green','height':'1/4','orders':[5,5,5],'bits':128,'density_tolerance':'1/6000' if kind=='Q' else '1/1500'}
                jobs.append({'id':name,'kind':kind,'settings':f'campaign/settings/{name}.json','state':f'campaign/states/{name}.json','physical_entries':physical,'area':str(area(r)),'multiplicity':len(physical)})
                dest=ROOT/'campaign/settings'/f'{name}.json';dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_text(json.dumps(settings,indent=2)+'\n')
    for kind in ['Q','g']:
        records=[e for e in inherited if e['kind']==kind]
        records += [p for job in jobs if job['kind']==kind for p in job['physical_entries']]
        cover=full_cover(records)
        assert all(v['fallback_area']=='0' for v in cover.values())
    payload={'schema':'paper1-r4-campaign-v1','baseline':'../../output_r3/paper1-complete-domain-validation-r3','baseline_manifest_sha256':hashlib.sha256((R3/'EVIDENCE_MANIFEST.json').read_bytes()).hexdigest(),'inherited_evidence':inherited,'jobs':jobs,'initial_tolerances_are_exploratory':True}
    (ROOT/'CAMPAIGN.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'jobs':len(jobs),'groups':{f'{k}_c{c}':sum(j['kind']==k and j['id'].startswith(f'{k}_c{c}') for j in jobs) for k in ['Q','g'] for c in [0,1]},'physical_new_areas':{k:str(sum(F(j['area'])*j['multiplicity'] for j in jobs if j['kind']==k)) for k in ['Q','g']},'coverage':'exact complete partition; pending jobs not yet enclosures'}))
if __name__=='__main__':build()
