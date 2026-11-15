#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(cmd,cwd):
 p=subprocess.run(cmd,cwd=cwd,text=True,capture_output=True,timeout=180,env={**os.environ,'OMP_NUM_THREADS':'2'})
 if p.returncode:raise RuntimeError(str(cmd)+'\n'+p.stdout+'\n'+p.stderr)
 return p
with tempfile.TemporaryDirectory(prefix='p12_checkpoint_') as td:
 d=Path(td);orig=d/'orig';chk=d/'chk'
 subprocess.run(['g++','-O2','-std=c++17','-fopenmp',str(ROOT/'run/p07_trajectory.cpp'),'-o',str(orig)],check=True)
 subprocess.run(['g++','-O2','-std=c++17','-fopenmp',str(ROOT/'run/p07_trajectory_checkpoint.cpp'),'-o',str(chk)],check=True)
 run([str(orig),'2',str(d/'mono.csv'),str(d/'mono.bin'),'80','12','12','0.04','2','1.2','0.4'],d)
 run([str(chk),'2','-',str(d/'c1.bin'),'-','-','80','12','12','0.04','1','1.2','0.4'],d)
 run([str(chk),'2',str(d/'c1.bin'),'-',str(d/'seg.csv'),str(d/'seg.bin'),'80','12','12','0.04','2','1.2','0.4'],d)
 checks={'csv':(d/'mono.csv').read_bytes()==(d/'seg.csv').read_bytes(),'snapshots':(d/'mono.bin').read_bytes()==(d/'seg.bin').read_bytes()}
 out={'status':'passed' if all(checks.values()) else 'failed','checks':checks,'hashes':{'csv':sha(d/'mono.csv'),'snapshots':sha(d/'mono.bin')}}
 print(json.dumps(out,indent=2,sort_keys=True))
 if out['status']!='passed':raise SystemExit(1)
