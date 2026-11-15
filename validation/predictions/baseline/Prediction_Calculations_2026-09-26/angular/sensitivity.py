"""Change only the saved imaginary kernel in the dressed source, N160 vs N192."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np
import core,integrate

P=Path(__file__).resolve().parent
A=P/'inputs/omega1.2_N160.npz'
if not A.exists():
 A=P.parents[1]/'review_r8/input/paper1-validation-integrated-r8/legacy_rc4/validation/current_numerics/analysis/resolution_execution/surface/results/omega1.2_N160.npz'
record=json.loads(A.with_suffix('.json').read_text());h=hashlib.sha256(A.read_bytes()).hexdigest()
assert h==record['matrix_sha256'] and record['omega']==1.2 and record['grid']==160
if A.resolve()!=(P/'inputs/omega1.2_N160.npz').resolve():
 shutil.copyfile(A,P/'inputs/omega1.2_N160.npz');shutil.copyfile(A.with_suffix('.json'),P/'inputs/omega1.2_N160.json')
z=np.load(A);selected=np.load(P/'inputs/REFERENCE_KERNEL.npz');pv=selected['PV'][:12,:12]
outpath=P/'inputs/REFERENCE_KERNEL_N160_SAME_PV.npz'
np.savez_compressed(outpath,G0=pv+1j*z['Im_total'],**{f'Im_branch_{j}':z[f'Im{j}'] for j in range(3)})
old=core.states
integrate.states=lambda:old(outpath)
result=integrate.run(128,64,513,[10.,20.,30.])
result['input']={'source':str(A),'source_sha256':h,'derived_kernel_sha256':hashlib.sha256(outpath.read_bytes()).hexdigest(),
                 'change':'Only imaginary source-state Green kernel replaced with saved N160 shell; selected M256 PV exactly unchanged.'}
(P/'integral_source_N160.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result['total_relative_difference_from_selected_kernel'],indent=2))
