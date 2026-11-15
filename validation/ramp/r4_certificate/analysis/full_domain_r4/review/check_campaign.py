"""Independent exact R4 angular partition and physical-symmetry audit."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def rectangle(values):
    a,b,c,d=map(F,values)
    assert 0<=a<b<=1 and 0<=c<d<=1
    return a,b,c,d


def area(values):
    a,b,c,d=rectangle(values)
    return (b-a)*(d-c)


def main():
    data=json.loads((ROOT/'CAMPAIGN.json').read_text())
    jobs=data['jobs']
    assert len({j['id'] for j in jobs})==len(jobs)
    assert len({j['state'] for j in jobs})==len(jobs)
    assert len({j['settings'] for j in jobs})==len(jobs)
    checked=[]
    for job in jobs:
        settings=json.loads((ROOT/job['settings']).read_text())
        r=rectangle(settings['bounds'])
        assert settings['mode']=={'Q':'qt','g':'green'}[job['kind']]
        assert job['area']==str(area(r))
        assert job['multiplicity']==len(job['physical_entries'])
        original=job['physical_entries'][0]
        assert original['chart']==settings['chart']
        assert rectangle(original['bounds'])==r
        assert original.get('transform') is None
        assert settings['chart'] in (0,1)
        expected=1 if settings['chart']==0 and r==(r[2],r[3],r[0],r[1]) else 2
        assert job['multiplicity']==expected
        if expected==2:
            copied=job['physical_entries'][1]
            if settings['chart']==0:
                assert copied['chart']==0 and copied['transform']=='swap_uv_chart0'
                assert rectangle(copied['bounds'])==(r[2],r[3],r[0],r[1])
            else:
                assert copied['chart']==2 and copied['transform']=='chart1_to_chart2'
                assert rectangle(copied['bounds'])==r
    checked.append(f'{len(jobs)} job settings and physical copies independently matched')
    records=data['inherited_evidence']+[
        dict(p,kind=j['kind']) for j in jobs for p in j['physical_entries']]
    areas={}
    for kind in ('Q','g'):
        for chart in (0,1,2):
            boxes=[rectangle(x['bounds']) for x in records if x['kind']==kind and x['chart']==chart]
            for index,(a,b,c,d) in enumerate(boxes):
                for e,f,g,h in boxes[:index]:
                    assert min(b,f)<=max(a,e) or min(d,h)<=max(c,g)
            total=sum((area(box) for box in boxes),F(0))
            assert total==1
            areas[f'{kind}_chart{chart}']={'rectangles':len(boxes),'exact_area':str(total)}
    checked.append('all six observable/chart covers have unit area and pairwise disjoint interiors')
    out={'status':'passed','checks':checked,'coverage':areas,
         'campaign_sha256':hashlib.sha256((ROOT/'CAMPAIGN.json').read_bytes()).hexdigest(),
         'limits':'Partition only; scheduled jobs remain unperformed until separately validated complete checkpoints exist.'}
    Path(__file__).with_name('CAMPAIGN_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))


if __name__=='__main__':main()
