from pathlib import Path
import os,sys,json,time,hashlib,argparse
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1');sys.dont_write_bytecode=True
import numpy as np
from scipy.ndimage import map_coordinates
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'src'))
import on_shell_basis as shell
from model_api_v1_1 import ModelParameters
from screen import histories,response
geometry=[];original=shell.marching_cubes
def strict(volume,*args,**kwargs):
    assert np.isfinite(volume).all() and volume.min()<0<volume.max()
    ans=original(volume,*args,**kwargs);v,f=ans[:2];cent=v[f].mean(axis=1)
    raw=volume[:-1,:-1,:-1];dk=2*np.pi/raw.shape[0]
    grad=[(np.roll(raw,-1,axis=i)-np.roll(raw,1,axis=i))/(2*dk) for i in range(3)]
    gn=np.sqrt(sum(map_coordinates(shell._extend(g),cent.T,order=1,mode='wrap')**2 for g in grad))
    area=.5*np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)*dk**2
    assert np.isfinite(gn).all() and np.all(gn>1e-11)
    assert np.isfinite(area).all() and np.all(area>0)
    geometry.append({'triangles':len(f),'min_gradient':float(gn.min()),'min_area':float(area.min()),'omitted_cells':0})
    return ans
shell.marching_cubes=strict
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--omega',type=float,required=True);ap.add_argument('--grid',type=int,required=True);a=ap.parse_args()
    start=time.monotonic();tag=f'omega{a.omega:.6f}_N{a.grid}'.replace('.','p');p=HERE/'fresh'/tag;p.parent.mkdir(exist_ok=True)
    if p.with_suffix('.json').exists():raise RuntimeError('Checkpoint exists; inspect/reuse instead of overwriting')
    on=shell.integrate_on_shell(a.omega,ModelParameters(),nseg=6,grid=a.grid)
    assert len(geometry)==3
    for m in [on.Im_total,*on.Im_branches]:
        assert np.isfinite(m).all() and np.array_equal(m,m.T)
        assert np.linalg.eigvalsh(m).min()>-1e-13
    hs=histories();pvs={M:np.mean([h.principal_value_jackson(a.omega*a.omega,M) for h in hs],axis=0) for M in [192,256,384]}
    base=np.array([.1,.2,.2,.2,.2,.1]);rows=[]
    for M,pv in pvs.items():
        for scale in [.95,1,1.05]:rows.append({'PV_order':M,'common_scale':scale,**response(a.omega,base*scale,pv,on.Im_branches)})
    payload={'Im_total':on.Im_total,**{f'Im_branch_{i}':m for i,m in enumerate(on.Im_branches)},**{f'Im_branch_long_{i}':m for i,m in enumerate(on.Im_branch_longitudinal)},**{f'PV_{M}':pv for M,pv in pvs.items()}}
    np.savez_compressed(p.with_suffix('.npz'),**payload)
    record={'omega':a.omega,'grid':a.grid,'seconds':time.monotonic()-start,'geometry':geometry,'rows':rows,'raw_reciprocity_residual':on.symmetry_residual_before_projection,'matrix_sha256':hashlib.sha256(p.with_suffix('.npz').read_bytes()).hexdigest(),'scope':'Fresh finite-resolution numerical sample; observed convergence, not a rigorous integration enclosure.'}
    p.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'omega':a.omega,'grid':a.grid,'seconds':record['seconds'],'M256_rows':[x for x in rows if x['PV_order']==256]},indent=2))
if __name__=='__main__':main()
