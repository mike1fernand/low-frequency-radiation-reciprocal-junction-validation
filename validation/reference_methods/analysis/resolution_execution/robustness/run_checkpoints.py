#!/usr/bin/env python3
"""Sequential checkpoint runner; completed matching records may be resumed."""
from pathlib import Path
import hashlib,json,subprocess,sys,time,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[2];WORKER=HERE/'worker.py'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
expected={str(p.relative_to(BASE)):digest(p) for p in [WORKER,*sorted((HERE/'src').glob('*.py'))]}
jobs=[]
for i in range(49):
 w=.6+.025*i;n=10 if i==24 else 6;tag=f'omega{w:.6f}_N80_nseg{n}'
 jobs.append((HERE/'kernels'/tag,['--family','shell','--omega',format(w,'.6f'),'--grid','80','--nseg',str(n)]))
for kd in (.75,.9,1.,1.1,1.25):
 tag=f'kdiag{kd:g}_N80';jobs.append((HERE/'stiffness'/tag,['--family','stiffness','--kdiag',str(kd),'--grid','80']))
records=[];start=time.perf_counter()
for stem,args in jobs:
 p=stem.with_suffix('.json') if '.' not in stem.name else Path(str(stem)+'.json')
 # str addition preserves decimal points in the scientifically meaningful stem.
 p=Path(str(stem)+'.json')
 resumed=False
 if p.exists():
  record=json.loads(p.read_text());need(record.get('status')=='completed_for_declared_finite_sample_scope','Existing checkpoint not complete')
  need(record['source_sha256']==expected,'Existing checkpoint belongs to different producer/source')
  matrix=p.parent/record['matrix_file'];need(digest(matrix)==record['matrix_sha256'],'Existing matrix hash mismatch');resumed=True
 else:
  log=Path(str(stem)+'.log')
  with log.open('w') as f:
   process=subprocess.run([sys.executable,'-B',str(WORKER),*args],stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1'})
  need(process.returncode==0,f'Checkpoint failed: {log}')
  record=json.loads(p.read_text());need(record['source_sha256']==expected,'Producer changed during run')
 records.append(dict(record=str(p.relative_to(BASE)),sha256=digest(p),matrix_sha256=record['matrix_sha256'],seconds=record['seconds'],resumed=resumed))
 result=dict(status='in_progress',completed=len(records),total=len(jobs),elapsed_seconds=time.perf_counter()-start,records=records,source_sha256=expected)
 (HERE/'CHECKPOINT_RUN.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(dict(completed=len(records),total=len(jobs),record=p.name,resumed=resumed,seconds=record['seconds'])),flush=True)
result['status']='completed_all_54_kernel_jobs';(HERE/'CHECKPOINT_RUN.json').write_text(json.dumps(result,indent=2)+'\n')
