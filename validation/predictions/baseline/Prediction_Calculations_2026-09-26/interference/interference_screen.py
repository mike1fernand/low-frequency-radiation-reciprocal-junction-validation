"""Fresh finite solves on read-only archived kernels; not fresh bulk quadrature.

Run: OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
     PYTHONDONTWRITEBYTECODE=1 python interference_screen.py [--extended]
No archive implementation imports. Scalar line mass, bulk mass and spacing = 1.
"""
from pathlib import Path
import argparse, csv, hashlib, json, os, platform
from itertools import product
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ARCHIVE=(HERE/'inputs') if (HERE/'inputs').is_dir() else ROOT/'review_r8/input/paper1-validation-integrated-r8/legacy_rc4/validation/current_numerics/analysis/resolution_execution'
MAIN=ROOT/'paper_update_20260926/build/paper1_main.tex'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def path_id(p):
    p=Path(p)
    return str(p.relative_to(HERE)) if p.is_relative_to(HERE) else str(p.relative_to(ROOT))
def resolve_id(rel):
    return (HERE/rel) if (HERE/rel).exists() else ROOT/rel
def main_binding():
    binding=HERE/'MAIN_SOURCE_BINDING.json'
    return json.loads(binding.read_text()) if binding.exists() else {'path':path_id(MAIN),'sha256':sha(MAIN)}

def assemble(profile):
    """Assemble all eight central springs per segment in Cartesian coordinates."""
    ns=len(profile); nq=ns+1; nb=12*ns
    B=np.zeros((nb,2*ns)); HX=np.zeros((nq+nb,nq+nb)); HR=HX.copy(); HP=HX.copy()
    for s,g in enumerate(profile):
        for c,(sy,sz) in enumerate(product([-1,1],repeat=2)):
            rb=12*s+3*c; r=nq+rb
            B[rb,2*s]=.5
            B[rb+1,2*s+1]=sy/(2*np.sqrt(2))
            B[rb+2,2*s+1]=sz/(2*np.sqrt(2))
            x=np.zeros(nq+nb);x[r]=1;x[s]=x[s+1]=-.5
            v=np.zeros(nq+nb);v[r+1]=sy;v[r+2]=sz;v[s]=-.5;v[s+1]=.5
            HX+=2*g/3*np.outer(x,x);HR+=2*g/3*np.outer(v,v)
            for ep,sign in [(s,1),(s+1,-1)]:
                n=np.array([sign,sy,sz])/np.sqrt(3)
                b=np.zeros(nq+nb);b[ep]=-n[0];b[r:r+3]=n
                HP+=g*np.outer(b,b)
    def blocks(h): return B.T@h[nq:,nq:]@B,B.T@h[nq:,:nq],h[:nq,:nq]
    return blocks(HX),blocks(HR),float(abs(HP-HX-HR).max())

def load_kernel(path,ns=6):
    with np.load(path,allow_pickle=False) as z:
        G=z['G0'][:2*ns,:2*ns].copy()
        Q=[z[f'Im_branch_{i}' if f'Im_branch_{i}' in z.files else f'Im{i}'][:2*ns,:2*ns].copy() for i in range(3)]
    meta=json.loads(path.with_suffix('.json').read_text())
    return G,Q,meta

def solve(G,Q,omega,cq,profile,ax=1.,ar=1.,state=False):
    (Ux,Vx,Kx),(Ur,Vr,Kr),herr=assemble(profile)
    U=ax*Ux+ar*Ur;V=ax*Vx+ar*Vr;Ks=ax*Kx+ar*Kr
    nc=len(profile);nq=nc+1;Kq=cq**2
    kap=2*np.arcsin(omega/(2*cq));gamma=2*Kq*np.sin(kap)
    KD=2*Kq*np.eye(nq)-Kq*(np.eye(nq,k=1)+np.eye(nq,k=-1))
    Sl=np.zeros((nq,nq),complex);Sl[0,0]=Sl[-1,-1]=-Kq*np.exp(1j*kap)
    gc=np.linalg.solve(np.eye(2*nc)+G@U,G)
    inverse=omega**2*np.eye(nq)-KD-Ks-Sl+V.T@gc@V
    GD=np.linalg.inv(inverse);chi=GD[:,0]*np.sqrt(gamma)
    A=np.linalg.inv(np.eye(2*nc)+U@G);eta=A@V@chi
    ts=np.array([2*np.real(eta.conj()@q@eta) for q in Q])
    tq=gamma**2*abs(GD[-1,0])**2;rq=abs(1-1j*gamma*GD[0,0])**2
    # Independent simultaneous line/core equations test the elimination arithmetic.
    aug=np.block([[omega**2*np.eye(nq)-KD-Ks-Sl,-V.T],[G@V,np.eye(2*nc)+G@U]])
    rhs=np.zeros(nq+2*nc,complex);rhs[0]=np.sqrt(gamma)
    direct=np.linalg.solve(aug,rhs)
    eta_direct=V@direct[:nq]+U@direct[nq:]
    result={'Rq':float(rq),'Tline':float(tq),'Tbulk':float(ts.sum()),'TqT':float(ts[:2].sum()),'TqL':float(ts[2]),
        'closure':float(abs(rq+tq+ts.sum()-1)), 'hessian_pair_error':herr,
        'augmented_eta_error':float(np.max(abs(eta_direct-eta))),
        'min_core_singular':float(np.linalg.svd(np.eye(2*nc)+G@U,compute_uv=False)[-1]),
        'min_device_singular':float(np.linalg.svd(inverse,compute_uv=False)[-1])}
    if state:
        ex=A@(ax*Vx)@chi;er=A@(ar*Vr)@chi
        def terms(q):
            xx=float(2*np.real(ex.conj()@q@ex));rr=float(2*np.real(er.conj()@q@er));cross=float(4*np.real(ex.conj()@q@er))
            return {'XX':xx,'RR':rr,'XR':cross,'sum':xx+rr+cross,'normalized_cross':cross/(xx+rr+cross),
                    'overlap_cosine':cross/(2*np.sqrt(xx*rr)) if xx*rr>0 else None}
        result['same_state_qT']=terms(Q[0]+Q[1]);result['same_state_qL']=terms(Q[2])
        result['decomposition_residual']=abs(result['same_state_qT']['sum']-result['TqT'])
    return result

def row(path,cq,profile):
    G,Q,m=load_kernel(path,len(profile));om=float(m['omega'])
    if om>=2*cq: return None
    full=solve(G,Q,om,cq,profile,state=True)
    controls={str(a):solve(G,Q,om,cq,profile,ar=a) for a in [0.,.5,1.]}
    controls['R_only']=solve(G,Q,om,cq,profile,ax=0.)
    return {'kernel':path_id(path),'kernel_sha256':sha(path),'metadata_sha256':sha(path.with_suffix('.json')),
            'omega':om,'cq':cq,'profile':list(profile),'grid':m.get('grid',m.get('surface_grid')),
            'kernel_spectral_identity_residual':float(np.max(abs(sum(Q)-(G-G.conj().T)/(2j)))),
            'qT_min_eigenvalue':float(np.linalg.eigvalsh(Q[0]+Q[1])[0]),
            'full':full,'controls':controls,
            'redressed_full_minus_X':full['TqT']-controls['0.0']['TqT'],
            'alpha_half_minus_zero':controls['0.5']['TqT']-controls['0.0']['TqT'],
            'alpha_one_minus_half':full['TqT']-controls['0.5']['TqT']}

def surface(om,n=128):
    if om==1.2: return ARCHIVE/'selected/REFERENCE_KERNEL.npz'
    d=ARCHIVE/f'surface/selected_nodes/N{n}'
    candidates=[p for p in d.glob('omega*.npz') if abs(float(json.loads(p.with_suffix('.json').read_text())['omega'])-om)<1e-12]
    if not candidates and om==1.:
        return ARCHIVE/'surface/selected_nodes/N192/omega1.npz'
    return candidates[0]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--extended',action='store_true');args=ap.parse_args()
    profile=[.1,.2,.2,.2,.2,.1]
    oms=[.8,1.2,1.6] if not args.extended else [.8,1.,1.2,1.4,1.6,1.7]
    speeds=[.8,1.3,2.] if not args.extended else [.65,.8,1.,1.3,1.7,2.]
    rows=[]
    for om,cq in product(oms,speeds):
        r=row(surface(om),cq,profile)
        if r: rows.append(r)
    name='extended' if args.extended else 'pilot'
    out={'status':'completed_saved_kernel_screen','scope':'Fresh Cartesian spring assembly and finite solves on archived bulk kernels; not a new bulk numerical certificate.',
         'script_sha256':sha(__file__),'main_source':main_binding(),'environment':{'python':platform.python_version(),'numpy':np.__version__,'threads':{k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}},'rows':rows}
    (HERE/f'{name}.json').write_text(json.dumps(out,indent=2))
    flat=[{'omega':r['omega'],'cq':r['cq'],'grid':r['grid'],'TqT':r['full']['TqT'],'XR':r['full']['same_state_qT']['XR'],
           'XR_fraction':r['full']['same_state_qT']['normalized_cross'],'full_minus_X':r['redressed_full_minus_X'],
           'half_minus_zero':r['alpha_half_minus_zero'],'one_minus_half':r['alpha_one_minus_half'],'kernel_sha256':r['kernel_sha256']} for r in rows]
    with (HERE/f'{name}.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    print(json.dumps({'run':name,'rows':flat,'max_closure':max(r['full']['closure'] for r in rows)},indent=2))

if __name__=='__main__':main()
