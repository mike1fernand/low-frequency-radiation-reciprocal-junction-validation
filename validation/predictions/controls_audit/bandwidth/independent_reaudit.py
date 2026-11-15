"""Bounded independent audit of archived bandwidth rows; baseline is read-only.

Reassembles paired Cartesian spring Hessians and solves the joint line/core
system. The archived bulk quadratures are treated as inputs, not regenerated.
"""
from pathlib import Path
import os, sys, json, hashlib, itertools
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.dont_write_bytecode=True
import numpy as np
from scipy.stats import qmc

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'baseline/Prediction_Calculations_2026-09-26/bandwidth'
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'src'))
from spectral_histogram import SpectralHistogram

def cartesian_matrices(g):
    n=len(g); B=np.zeros((12*n,2*n)); H=np.zeros((12*n+n+1,12*n+n+1))
    for s,gs in enumerate(g):
        for c,(sy,sz) in enumerate(itertools.product([-1,1],repeat=2)):
            a=12*s+3*c
            B[a,2*s]=.5
            B[a+1,2*s+1]=sy/(2*np.sqrt(2))
            B[a+2,2*s+1]=sz/(2*np.sqrt(2))
            for end,sx in [(s,1),(s+1,-1)]:
                direction=np.array([sx,sy,sz])/np.sqrt(3)
                extension=np.zeros(len(H))
                extension[a:a+3]=direction
                extension[12*n+end]=-direction[0]
                H+=gs*np.outer(extension,extension)
    return B.T@H[:12*n,:12*n]@B,B.T@H[:12*n,12*n:],H[12*n:,12*n:]

def physical_solve(w,g,G0,branch_matrices):
    n=len(g); U,V,Kself=cartesian_matrices(g)
    kq=2*np.arcsin(w/(2*np.sqrt(1.69)))
    sigma=-1.69*np.exp(1j*kq)
    gamma=2*1.69*np.sin(kq)
    Kline=np.diag(np.full(n+1,2*1.69))+np.diag(np.full(n,-1.69),1)+np.diag(np.full(n,-1.69),-1)
    D=w*w*np.eye(n+1)-Kline-Kself
    D=D.astype(complex); D[0,0]-=sigma; D[-1,-1]-=sigma
    # eta = V q - U G0 eta; line equation D q + V^T G0 eta = drive.
    joint=np.block([[D,V.T@G0],[-V,np.eye(2*n)+U@G0]])
    drive=np.zeros(3*n+1,complex);drive[0]=1j*gamma
    answer=np.linalg.solve(joint,drive);q=answer[:n+1];force=-answer[n+1:]
    incoming=w*gamma/4
    powers=np.array([w/2*np.vdot(force,b@force).real/incoming for b in branch_matrices])
    reflection=abs(1-q[0])**2
    transmission=abs(q[-1])**2
    return dict(TqT=float(sum(powers[:2])),TqL=float(powers[2]),Rq=float(reflection),Tline=float(transmission),closure=float(abs(reflection+transmission+sum(powers)-1)))

def main():
    record=json.loads((BASE/'FINAL_SCREEN.json').read_text())
    hs=[]
    for p in sorted((BASE/'inputs').glob('histogram_*.npz')):
        with np.load(p) as z:hs.append(SpectralHistogram(z['edges'],z['mass'],int(z['sample_count']),'audit',{}))
    base=np.array([.1,.2,.2,.2,.2,.1])
    settings={'nominal':base}|{f'corner_{i}':base*(1+.05*np.array(signs)) for i,signs in enumerate(itertools.product([-1,1],repeat=6))}
    settings.update({f'interior_{i}':base*(.95+.1*u) for i,u in enumerate(qmc.Sobol(6,scramble=True,seed=20260926).random_base2(7))})
    kernels={}
    for w in sorted(set(r['omega'] for r in record['surface_rows'])):
        if w in [1.1125,1.2375,1.3]:
            p=BASE/'fresh'/(f'omega{w:.6f}_N128'.replace('.','p')+'.npz')
            with np.load(p) as z:pv=z['PV_256'];bs=[z[f'Im_branch_{i}'] for i in range(3)]
        else:
            with np.load(BASE/'inputs'/f'shell_{w:.6f}_N80.npz') as z:bs=[z[f'Im_branch_{i}'] for i in range(3)]
            pv=np.mean([h.principal_value_jackson(w*w,256) for h in hs],axis=0)
        kernels[w]=(pv+1j*sum(bs),bs)
    surface=[];maxdiff=0;maxclosure=0
    for row in record['surface_rows']:
        assert np.array_equal(np.array(row['g']),settings[row['setting']])
        fresh=physical_solve(row['omega'],row['g'],*kernels[row['omega']])
        maxdiff=max(maxdiff,max(abs(fresh[k]-row[k]) for k in ['TqT','TqL','Rq','Tline']))
        maxclosure=max(maxclosure,fresh['closure'])
        surface.append(dict(omega=row['omega'],setting=row['setting'],**fresh))
    tetra_kernels={}
    for i,w in enumerate([1.1125,1.2375,1.3]):
        with np.load(BASE/'tetra_N64'/f'omega_{i:03d}.npz') as z:tetra_kernels[w]=(z['G0'],[z[f'Im{i}'] for i in range(3)])
    tetra=[];max_tetra_diff=0
    for row in record['independent_rows']:
        fresh=physical_solve(row['omega'],row['g'],*tetra_kernels[row['omega']])
        max_tetra_diff=max(max_tetra_diff,max(abs(fresh[a]-row[b]) for a,b in [('TqT','TT_modal'),('TqL','TL_modal'),('Rq','Rq_amplitude'),('Tline','Tq')]))
        maxclosure=max(maxclosure,fresh['closure'])
        tetra.append(dict(omega=row['omega'],setting=row['setting'],**fresh))
    grouped={str(w):len({r['setting'] for r in surface if r['omega']==w}) for w in sorted(kernels)}
    spec=json.loads((BASE/'INPUTS.json').read_text())
    bindings={}
    for category,folder in [('copied_sources','src'),('derived_inputs','inputs')]:
        for name,info in spec[category].items():bindings[str(Path(folder)/name)]=hashlib.sha256((BASE/folder/name).read_bytes()).hexdigest()==info['sha256']
    long_lower=[sum(r['branch_longitudinal_powers'][:2])/r['TT_modal'] for r in record['independent_rows']]
    high_long=[r['branch_longitudinal_fractions'][2] for r in record['independent_rows']]
    out={
      'scope':'Independent finite scattering solve from Cartesian bond Hessian and analytical propagating lead; archived bulk quadratures reused. No new quadrature or certified error bound.',
      'surface_count':len(surface),'tetra_count':len(tetra),'unique_settings_per_frequency':grouped,
      'surface_max_abs_coefficient_difference':maxdiff,'tetra_max_abs_coefficient_difference':max_tetra_diff,'independent_max_closure':maxclosure,
      'minimum_surface':min(surface,key=lambda r:r['TqT']),'minimum_tetra':min(tetra,key=lambda r:r['TqT']),
      'surface_all_at_least_0p10':all(r['TqT']>=.1 for r in surface),'tetra_all_at_least_0p10':all(r['TqT']>=.1 for r in tetra),
      'all_declared_input_hashes_match':all(bindings.values()),'input_hash_checks':bindings,
      'saved_tetra_lower_cluster_longitudinal_power_fraction_range':[min(long_lower),max(long_lower)],
      'saved_tetra_upper_branch_longitudinal_power_fraction_range':[min(high_long),max(high_long)],
      'max_modal_vs_wave_normal_transverse_fraction_difference':max(abs(r['TT_modal']-r['TT_projector']) for r in record['independent_rows']),
      'minimum_wave_normal_transverse_fraction':min(r['TT_projector'] for r in record['independent_rows'])}
    (OUT/'INDEPENDENT_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
