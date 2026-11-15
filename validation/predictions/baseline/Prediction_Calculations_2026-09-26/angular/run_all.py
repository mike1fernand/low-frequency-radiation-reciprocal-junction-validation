"""Replay from the copied inputs, with the currently selected Python environment."""
import os,subprocess,sys
from pathlib import Path

P=Path(__file__).resolve().parent
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
jobs=[
 ['integrate.py','--nmu','128','--nphi','64','--scan','513','--angles','.5','1','2','5','10','20','30','--output','integral_medium.json'],
 ['integrate.py','--nmu','256','--nphi','128','--scan','1025','--angles','.5','1','2','5','10','20','30','--output','integral_fine.json'],
 ['check_selection.py'],['sensitivity.py'],['summarize.py']]
for job in jobs:
 print('Running',job[0],flush=True)
 subprocess.run([sys.executable,str(P/job[0]),*job[1:]],env=env,check=True)
