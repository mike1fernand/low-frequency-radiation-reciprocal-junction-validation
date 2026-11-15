"""Combine retained samples and fresh held-out kernels; check independent tetra."""
from pathlib import Path
import sys,json,itertools,csv
from dataclasses import asdict
import numpy as np
from scipy.stats import qmc
from screen import response
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'tetra_src'))
from combined_tetra_diagnostic import Parameters,TetrahedronResult,independent_scatter
base=np.array([.1,.2,.2,.2,.2,.1])
factors=[('nominal',np.ones(6))]+[('corner_'+str(i),1+.05*np.array(s)) for i,s in enumerate(itertools.product([-1,1],repeat=6))]
factors += [('interior_'+str(i),.95+.10*u) for i,u in enumerate(qmc.Sobol(6,scramble=True,seed=20260926).random_base2(7))]
prior=json.loads((HERE/'screen_expanded.json').read_text())
rows=[dict(shell_grid=80,PV_order=256,**r) for r in prior['rows'] if abs(r['omega']-1.3)>1e-10]
tetra_rows=[];comparisons=[]
for j,w in enumerate([1.1125,1.2375,1.3]):
    stem=f'omega{w:.6f}_N128'.replace('.','p')
    with np.load(HERE/'fresh'/f'{stem}.npz',allow_pickle=False) as z:
        pv=z['PV_256'];bs=[z[f'Im_branch_{i}'] for i in range(3)]
    current=[dict(shell_grid=128,PV_order=256,setting=label,**response(w,base*f,pv,bs)) for label,f in factors]
    rows+=current
    with np.load(HERE/'tetra_N64'/f'omega_{j:03d}.npz',allow_pickle=False) as z:
        t=TetrahedronResult(w,1e-10,64,6,z['G0'],[z[f'Im{i}'] for i in range(3)],z['ImL'],[z[f'ImBL{i}'] for i in range(3)],0.,float(np.linalg.eigvalsh(z['G0'].imag).min()),0.,{})
    for label,f in factors:
        result=asdict(independent_scatter(w,base*f,Parameters(),t,eta_lead=0))
        assert result['closure_amplitude']<1e-10 and result['component_closure_modal']<1e-10
        tetra_rows.append(dict(setting=label,g=(base*f).tolist(),**result))
    surf=current[1];ind=tetra_rows[-len(factors)+1]
    comparisons.append({'omega':w,'setting':'all factors 0.95','surface_N128':surf['TqT'],'tetra_N64':ind['TT_modal'],'absolute_method_difference':abs(surf['TqT']-ind['TT_modal'])})
rows.sort(key=lambda r:(r['omega'],r['setting']))
summary={'scope':'Finite sampled evidence; no frequency-continuous or whole-box certificate. Observed method differences are not rigorous error bounds.',
 'samples':len(rows),'frequencies':sorted(set(r['omega'] for r in rows)),'couplings_per_frequency':len(factors),
 'minimum_surface':min(rows,key=lambda r:r['TqT']),'independent_samples':len(tetra_rows),
 'minimum_tetra':min(tetra_rows,key=lambda r:r['TT_modal']),
 'max_surface_closure':max(r['closure'] for r in rows),'max_independent_closure':max(r['closure_amplitude'] for r in tetra_rows),
 'comparisons':comparisons}
(HERE/'FINAL_SCREEN.json').write_text(json.dumps({'summary':summary,'surface_rows':rows,'independent_rows':tetra_rows},indent=2)+'\n')
(HERE/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
with (HERE/'FINAL_SCREEN.csv').open('w',newline='') as out:
    fields=['omega','setting','shell_grid','PV_order','TqT','TqL','Rq','Tline','closure'];writer=csv.DictWriter(out,fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
print(json.dumps(summary,indent=2))
