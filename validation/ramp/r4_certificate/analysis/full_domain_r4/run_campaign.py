"""Bounded process scheduler; each tile is independently restartable.
Run repeatedly; existing complete tiles are inspected and skipped.
"""
import argparse,json,sys,time,subprocess,concurrent.futures,collections,os
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent

def load_jobs(manifest):
    groups=collections.defaultdict(list)
    for job in manifest['jobs']:groups[(job['kind'],job['physical_entries'][0]['chart'])].append(job)
    for g in groups.values():g.sort(key=lambda j:(float(F(j['physical_entries'][0]['bounds'][1])+F(j['physical_entries'][0]['bounds'][3])),j['id']))
    jobs=[]
    while any(groups.values()):
        for key in sorted(groups):
            if groups[key]:jobs.append(groups[key].pop(0))
    return jobs

def status(job):
    p=ROOT/job['state']
    if not p.exists():return None
    return json.loads(p.read_text())['status']

def stage(job,seconds):
    command=[sys.executable,str(ROOT/'integrator/run_stage.py'),'--settings',str(ROOT/job['settings']),'--output',str(ROOT/job['state']),'--max-seconds',str(seconds)]
    started=time.monotonic()
    try:
        r=subprocess.run(command,capture_output=True,text=True,timeout=45)
        result={'id':job['id'],'wall_seconds':time.monotonic()-started,'returncode':r.returncode}
        if r.returncode==0:result['stage_summary']=json.loads(r.stdout)
        else:result['error']=(r.stderr or r.stdout)[-3000:]
    except subprocess.TimeoutExpired:
        result={'id':job['id'],'wall_seconds':time.monotonic()-started,'returncode':124,'error':'45-second process limit; last atomic checkpoint retained'}
    return result

def progress(jobs,started,stages):
    states={j['id']:status(j) for j in jobs}
    return {'wall_seconds':round(time.monotonic()-started,2),'stages_this_run':stages,'jobs':len(jobs),'complete':sum(v=='rigorous_complete_tile_enclosure' for v in states.values()),'partial':sum(v is not None and v!='rigorous_complete_tile_enclosure' for v in states.values()),'new_complete_physical_area':{k:str(sum(F(j['area'])*j['multiplicity'] for j in jobs if j['kind']==k and states[j['id']]=='rigorous_complete_tile_enclosure')) for k in ['Q','g']}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--max-seconds',type=float,default=600);p.add_argument('--stage-seconds',type=float,default=20);p.add_argument('--limit-jobs',type=int);a=p.parse_args()
    if not 1<=a.workers<=6:raise ValueError('workers must be 1..6')
    manifest=json.loads((ROOT/'CAMPAIGN.json').read_text());jobs=load_jobs(manifest)
    if a.limit_jobs:jobs=jobs[:a.limit_jobs]
    todo=collections.deque(j for j in jobs if status(j)!='rigorous_complete_tile_enclosure')
    started=time.monotonic();stages=0;lastprint=started;failed=False
    logpath=ROOT/'campaign/RUN_LOG.jsonl';logpath.parent.mkdir(parents=True,exist_ok=True)
    with logpath.open('a') as log,concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        active={}
        while todo or active:
            while todo and len(active)<a.workers and time.monotonic()-started<a.max_seconds and not failed:
                j=todo.popleft();active[pool.submit(stage,j,a.stage_seconds)]=j
            if not active:break
            done,_=concurrent.futures.wait(active,timeout=2,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in done:
                j=active.pop(f);r=f.result();stages+=1
                log.write(json.dumps(r,separators=(',',':'))+'\n');log.flush()
                if r['returncode']!=0 or r.get('stage_summary',{}).get('status')=='blocked_checkpoint_no_tile_enclosure':
                    failed=True;print(json.dumps({'error':r}),flush=True)
                elif r['stage_summary']['status']!='rigorous_complete_tile_enclosure':
                    # Finish a started tile before launching further ones: preserves useful complete enclosures.
                    todo.appendleft(j)
            if time.monotonic()-lastprint>=20:
                print(json.dumps(progress(jobs,started,stages)),flush=True);lastprint=time.monotonic()
    final=progress(jobs,started,stages);final['failed']=failed;final['all_selected_complete']=final['complete']==len(jobs)
    (ROOT/'campaign/LAST_RUN.json').write_text(json.dumps(final,indent=2)+'\n');print(json.dumps(final),flush=True)
    if failed:sys.exit(1)
if __name__=='__main__':main()
