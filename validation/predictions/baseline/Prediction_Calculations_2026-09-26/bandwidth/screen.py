from pathlib import Path
import os,sys,json,itertools,time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
sys.dont_write_bytecode=True
import numpy as np
from scipy.stats import qmc
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'src'))
from model_api_v1_1 import ModelParameters
from reduced_core import scatter_from_reduced_green
from spectral_histogram import SpectralHistogram

def histories():
    out=[]
    for p in sorted((HERE/'inputs').glob('histogram_*.npz')):
        with np.load(p,allow_pickle=False) as z:
            out.append(SpectralHistogram(z['edges'],z['mass'],int(z['sample_count']),'bound R8 Sobol',{}))
    return out

def response(w,g,pv,bs):
    r,x=scatter_from_reduced_green(w,g,ModelParameters(),pv+1j*sum(bs),eta_lead=0,return_internal=True)
    powers=[float(2*x['gamma_in']*np.vdot(x['t'],b@x['t']).real) for b in bs]
    row={'omega':w,'g':g.tolist(),'Rq':r.Rq,'Tline':r.Tq,'TqT':sum(powers[:2]),'TqL':powers[2],
         'closure':r.closure_error,'modal_residual':abs(sum(powers)-r.Tbulk)}
    assert row['closure']<1e-10 and row['modal_residual']<1e-11,row
    assert min(row[k] for k in ['Rq','Tline','TqT','TqL'])>=-1e-12,row
    return row

def run(interior=False):
    hs=histories();base=np.array([.1,.2,.2,.2,.2,.1]);rows=[]
    factors=[('nominal',np.ones(6))]+[('corner_'+str(i),1+.05*np.array(signs)) for i,signs in enumerate(itertools.product([-1,1],repeat=6))]
    if interior:
        factors += [('interior_'+str(i),.95+.10*u) for i,u in enumerate(qmc.Sobol(6,scramble=True,seed=20260926).random_base2(7))]
    for path in sorted((HERE/'inputs').glob('shell_*_N80.npz')):
        w=float(path.name.split('_')[1])
        if not interior and all(abs(w-v)>1e-10 for v in [1.10,1.20,1.30]):
            continue
        pv=np.mean([h.principal_value_jackson(w*w,256) for h in hs],axis=0)
        with np.load(path,allow_pickle=False) as z:bs=[z[f'Im_branch_{i}'] for i in range(3)]
        for label,f in factors:rows.append(dict(setting=label,**response(w,base*f,pv,bs)))
    worst=min(rows,key=lambda x:x['TqT']);summary={'rows':len(rows),'minimum_TqT':worst,'all_sampled_TqT_at_least_0_10':worst['TqT']>=.1,'max_closure':max(r['closure'] for r in rows),'scope':'Finite sampled numerical screen, not a domain certificate.'}
    result={'summary':summary,'rows':rows}
    output=HERE/('screen_expanded.json' if interior else 'screen_pilot.json');output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':run('--interior' in sys.argv)
