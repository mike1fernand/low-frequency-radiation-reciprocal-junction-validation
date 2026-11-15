"""Compare fixed and refitted PV and directional spectral matrices of the omega=2.6 control (24 complete motifs)."""
import sys, ast, json
from pathlib import Path
import numpy as np
base=Path(__file__).resolve().parent
sys.path.insert(0,str(base/'replay_sources'))
from p1solver import assemble,assemble_im,solve_junction,reference_profile
nc=24;w=2.6
p0=base/'prior_staged'
def fit(eps,a,degree):
    return np.polynomial.polynomial.polyfit(eps,a.reshape(len(eps),-1),degree)[0].reshape(a.shape[1:])
def calc(real,shell):
    qm,qp=shell
    q=qm+qp
    g=assemble(real,nc) + 1j*q
    s=solve_junction(reference_profile(nc),w,g,q,0*q,Kq=2.25)
    tm=float(2*(s['eta'].conj()@qm@s['eta']).real);tb=float(s['Tbulk'])
    return dict(Tbulk=tb,Tminus=tm,Bminus=tm/tb,closure=float(s['closure']))
a={ast.literal_eval(k):v for k,v in np.load(p0/'replay_sources/r3_eps_scan.npy',allow_pickle=True).item().items()}
eps=[.006,.008,.012,.016,.024]
orig_re=fit(1j*np.array(eps),np.array([a[w,e] for e in eps]),2).real
def blocks(I):
    return [assemble_im(I.sum(axis=0)[direction],nc,nc) for direction in (1,0)]
with np.load(p0/'dressed_26/shell_N2048.npz', allow_pickle=False) as z:
    orig_q=blocks(z['Im'])
report={'historical_control':calc(orig_re,orig_q),'same_pv_cartesian_shell':{},'cartesian_pv_same_shell':{},'cartesian_all':{},'combined_finite_powers':{},'realaxis_all':{}}
for f in sorted((base/'cartesian_quadrant').glob('*/combined_inputs.npz')):
    with np.load(f,allow_pickle=False) as d:
        for deg in (2,3,4):
            k=f'{int(d["N"])}_{deg}'; re=fit(d['epsilon'],d['G'].real,deg)[nc-1:]
            q=[assemble_im(fit(d['epsilon'],d[v],deg),nc,nc) for v in ['Qm','Qp']]
            report['same_pv_cartesian_shell'][k]=calc(orig_re,q)
            report['cartesian_pv_same_shell'][k]=calc(re,orig_q)
            report['cartesian_all'][k]=calc(re,q)
            rows=json.loads(f.with_name('results.json').read_text())['rows'];t=np.array([[s['T_minus_broadened_diagnostic'],s['T_plus_broadened_diagnostic']] for s in rows]);ends=fit(d['epsilon'],t,deg)
            report['combined_finite_powers'][k]={'Tbulk':ends.sum(),'Tminus':ends[0],'Bminus':ends[0]/ends.sum()}
for f in sorted((base/'realaxis').glob('*_result.json')):
    v=json.loads(f.read_text());report['realaxis_all'][f"{v['omega']:.2f}_N{v['Nperp']}"]={k:v[k] for k in ['Tbulk','T_minus','B_minus','closure']}
(base/'method-comparison.json').write_text(json.dumps(report,indent=2)+'\n')
for k in report:print(k,report[k].get('512_3',report[k]))
