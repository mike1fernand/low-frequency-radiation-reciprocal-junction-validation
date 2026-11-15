"""Recreate the original 120-site free-boundary junction in a temporary directory."""
import json
import os
from math import isclose
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT=Path(__file__).resolve().parent
def check(expected, got, path):
    if isinstance(expected,dict):
        if set(expected)!=set(got):raise ValueError(path)
        for key in expected:check(expected[key],got[key],path+'/'+key)
    elif isinstance(expected,list):
        if len(expected)!=len(got):raise ValueError(path)
        for a,b in zip(expected,got):check(a,b,path)
    elif isinstance(expected,float):
        if not isclose(expected,got,rel_tol=1e-8,abs_tol=2e-13):raise ValueError((path,expected,got))
    elif expected!=got:raise ValueError((path,expected,got))
def main():
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4')
    with TemporaryDirectory() as td:
        path=Path(td); solver=path/'evolve'
        subprocess.run(['g++','-std=c++17','-O3','-fopenmp',str(ROOT/'evolve.cpp'),'-o',str(solver)],check=True)
        data=subprocess.check_output([sys.executable,str(ROOT/'prepare.py')],env=env)
        subprocess.run([str(solver),td],input=data,env=env,check=True)
        results=json.loads(subprocess.check_output([sys.executable,str(ROOT/'analyze.py'),td],env=env))
    check(json.loads((ROOT/'results.json').read_text()),results,'root')
    print(json.dumps({k:results[k] for k in ('longest','mean_low_high_total','percent_current_balance_residual',
                                           'open_loaded_vs_free_periodic_velocity_relative_RMS')},indent=2))
if __name__=='__main__':main()
