"""Independent finite-matrix and analytic-sensitivity audit on retained kernels.

Run from any directory: python interference_reaudit.py
Requires NumPy only; never writes into the immutable baseline package.
No fresh bulk quadrature, interval certificate or physical alpha actuator is implied.
"""
from pathlib import Path
import hashlib, json, csv
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'baseline'/'Prediction_Calculations_2026-09-26'
INT=BASE/'interference'
PROFILE=np.array([.1,.2,.2,.2,.2,.1])

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def blocks(profile):
    """Closed-form motif Hessians from paired energy, independent of Cartesian producer."""
    n=len(profile); nq=n+1
    Ux=np.zeros((2*n,2*n)); Ur=Ux.copy()
    Vx=np.zeros((2*n,nq)); Vr=Vx.copy()
    Kx=np.zeros((nq,nq)); Kr=Kx.copy()
    for i,g in enumerate(profile):
        a=2*g/3
        Ux[2*i,2*i]=a; Ur[2*i+1,2*i+1]=2*a
        Vx[2*i,i:i+2]=[-a,-a]
        Vr[2*i+1,i:i+2]=[-np.sqrt(2)*a,np.sqrt(2)*a]
        Kx[i:i+2,i:i+2]+=a*np.ones((2,2))
        Kr[i:i+2,i:i+2]+=a*np.array([[1.,-1.],[-1.,1.]])
    return (Ux,Vx,Kx),(Ur,Vr,Kr)

def load(path,n):
    with np.load(path,allow_pickle=False) as z:
        G=z['G0'][:2*n,:2*n].copy()
        Q=[z[f'Im_branch_{j}' if f'Im_branch_{j}' in z.files else f'Im{j}'][:2*n,:2*n].copy() for j in range(3)]
    return G,Q

def calc(G,Q,om,cq,profile,alpha=1.,scale=1.,deriv='alpha',ax=1.):
    (Ux,Vx,Kx),(Ur,Vr,Kr)=blocks(profile)
    U=scale*(ax*Ux+alpha*Ur); V=scale*(ax*Vx+alpha*Vr); Ks=scale*(ax*Kx+alpha*Kr)
    if deriv=='alpha': dU=scale*Ur; dV=scale*Vr; dK=scale*Kr
    elif deriv=='scale': dU=ax*Ux+alpha*Ur; dV=ax*Vx+alpha*Vr; dK=ax*Kx+alpha*Kr
    else: raise ValueError(deriv)
    nc=U.shape[0]; nq=V.shape[1]
    kap=2*np.arcsin(om/(2*cq)); gamma=2*cq*cq*np.sin(kap)
    KD=cq*cq*(2*np.eye(nq)-np.eye(nq,k=1)-np.eye(nq,k=-1))
    Sigma=np.zeros((nq,nq),complex); Sigma[0,0]=Sigma[-1,-1]=-cq*cq*np.exp(1j*kap)
    L=om*om*np.eye(nq)-KD-Ks-Sigma
    M=np.block([[L,-V.T],[G@V,np.eye(nc)+G@U]])
    dM=np.block([[-dK,-dV.T],[G@dV,G@dU]])
    b=np.zeros(nq+nc,complex); b[0]=np.sqrt(gamma)
    x=np.linalg.solve(M,b); dx=np.linalg.solve(M,-dM@x)
    chi=x[:nq]; u=x[nq:]; dchi=dx[:nq]; du=dx[nq:]
    eta=V@chi+U@u
    deta=dV@chi+V@dchi+dU@u+U@du
    A=np.linalg.solve(np.eye(nc)+U@G,np.eye(nc))
    drive=A@dV@chi
    feedback=A@(V@dchi-dU@G@eta)
    ex=A@(scale*ax*Vx)@chi; er=A@(scale*alpha*Vr)@chi
    q=Q[0]+Q[1]
    quad=lambda e,q:float(2*np.real(e.conj()@q@e))
    T=quad(eta,q); XX=quad(ex,q); RR=quad(er,q); XR=float(4*np.real(ex.conj()@q@er))
    ds=float(4*np.real(eta.conj()@q@deta))
    dd=float(4*np.real(eta.conj()@q@drive))
    df=float(4*np.real(eta.conj()@q@feedback))
    # Schur derivative is a separately formed check on the simultaneous solve.
    gc=np.linalg.solve(np.eye(nc)+G@U,G)
    D=L+V.T@gc@V
    Dp=-dK+dV.T@gc@V+V.T@gc@dV-V.T@gc@dU@gc@V
    chi_s=np.linalg.solve(D,b[:nq]); dchi_s=np.linalg.solve(D,-Dp@chi_s)
    deta_s=A@(dV@chi_s+V@dchi_s-dU@G@(A@V@chi_s))
    Tbulk=sum(quad(eta,qq) for qq in Q)
    r=1-1j*np.sqrt(gamma)*chi[0]; t=-1j*np.sqrt(gamma)*chi[-1]
    R=float(abs(r)**2)
    Tl=float(abs(t)**2)
    seven=float(abs(r+t)**2); sodd=float(abs(r-t)**2)
    return {'TqT':T,'Tbulk':Tbulk,'Rq':R,'Tline':Tl,'XX':XX,'RR':RR,'XR':XR,
      'scattering_even_abs2':seven,'scattering_odd_abs2':sodd,
      'derivative':ds,'direct_source_derivative':dd,'redressing_derivative':df,
      'decomposition_residual':abs(T-XX-RR-XR),'flux_residual':abs(R+Tl+Tbulk-1),
      'source_identity_residual':float(np.max(abs(eta-A@V@chi))),
      'derivative_identity_residual':abs(ds-dd-df),
      'schur_derivative_source_residual':float(np.max(abs(deta-deta_s)))}

def audit():
    results={'scope':'Independent finite-matrix replay and exact analytic sensitivity; retained bulk kernels, no fresh quadrature/certificate.',
      'source_script_sha256':sha(__file__),'numpy':np.__version__}
    ext=json.loads((INT/'extended.json').read_text())
    replay=[]
    for row in ext['rows']:
        path=INT/row['kernel']; p=row['profile']; G,Q=load(path,len(p))
        full=calc(G,Q,row['omega'],row['cq'],p)
        xonly=calc(G,Q,row['omega'],row['cq'],p,alpha=0.)
        replay.append({'omega':row['omega'],'cq':row['cq'],'kernel':row['kernel'],
          'TqT_error':abs(full['TqT']-row['full']['TqT']),
          'XR_error':abs(full['XR']-row['full']['same_state_qT']['XR']),
          'delta_error':abs(full['TqT']-xonly['TqT']-row['redressed_full_minus_X']),
          'flux_residual':full['flux_residual']})
    results['replayed_31_rows']=replay
    selected=[('surface_0p8','inputs/surface/selected_nodes/N160/omega0.8.npz',.8),
      ('tetra_0p8','inputs/tetra/packet_N128/omega_000.npz',.8),
      ('surface_reference','inputs/selected/REFERENCE_KERNEL.npz',1.2),
      ('tetra_reference','inputs/tetra/reference_N192/omega_000.npz',1.2)]
    rows=[]
    for label,rel,om in selected:
        path=INT/rel; G,Q=load(path,6)
        for cq in ([.8,1.3] if om==.8 else [1.3]):
            full=calc(G,Q,om,cq,PROFILE); xo=calc(G,Q,om,cq,PROFILE,alpha=0.)
            fd=[]
            for h in [1e-4,5e-5]:
                num=(calc(G,Q,om,cq,PROFILE,alpha=1+h)['TqT']-calc(G,Q,om,cq,PROFILE,alpha=1-h)['TqT'])/(2*h)
                fd.append({'step':h,'derivative':num,'abs_error':abs(num-full['derivative'])})
            physical=calc(G,Q,om,cq,PROFILE,deriv='scale')
            h=5e-5
            num=(calc(G,Q,om,cq,PROFILE,scale=1+h)['TqT']-calc(G,Q,om,cq,PROFILE,scale=1-h)['TqT'])/(2*h)
            rows.append({'case':label,'kernel':rel,'kernel_sha256':sha(path),'omega':om,'cq':cq,
              'full':full,'X_only':xo,'full_minus_X':full['TqT']-xo['TqT'],
              'X_self_loading_change':full['XX']-xo['TqT'],
              'radial_and_cross_gain':full['RR']+full['XR'],
              'finite_difference_alpha':fd,'physical_common_scale_derivative':physical['derivative'],
              'physical_scale_fd_abs_error':abs(num-physical['derivative']),
              'physical_common_scale_finite':{str(s):calc(G,Q,om,cq,PROFILE,scale=s)['TqT'] for s in [.99,1.,1.01]}})
    results['selected_rows']=rows
    # The one-motif parity null can be strengthened to exact redressed additivity.
    # Orthogonal parity + left/right symmetric 2-site line diagonalizes even/odd channels.
    motif=[]
    saved=json.loads((INT/'confirmation.json').read_text())
    for row in saved['single_motif_rows']:
        path=INT/row['kernel']; G,Q=load(path,1); p=row['profile']; om=row['omega']; cq=row['cq']
        full=calc(G,Q,om,cq,p); xo=calc(G,Q,om,cq,p,alpha=0.); ro=calc(G,Q,om,cq,p,ax=0.)
        motif.append({'kernel':row['kernel'],'omega':om,'cq':cq,'g':p[0],
          'full_X_self_minus_Xonly':full['XX']-xo['TqT'],'XR':full['XR'],
          'full_minus_X':full['TqT']-xo['TqT'],'RR':full['RR'],
          'Ronly':ro['TqT'],'Xonly':xo['TqT'],
          'Xonly_bulk':xo['Tbulk'],'Ronly_bulk':ro['Tbulk'],
          'Xonly_uncoupled_parity_unitarity_residual':abs(xo['scattering_odd_abs2']-1),
          'Ronly_uncoupled_parity_unitarity_residual':abs(ro['scattering_even_abs2']-1),
          'Xonly_half_channel_flux_residual':abs(xo['Tbulk']-(1-xo['scattering_even_abs2'])/2),
          'Ronly_half_channel_flux_residual':abs(ro['Tbulk']-(1-ro['scattering_odd_abs2'])/2),
          'full_R_self_minus_Ronly':full['RR']-ro['TqT'],
          'qT_redressed_additivity_residual':abs(full['TqT']-xo['TqT']-ro['TqT']),
          'all_bulk_redressed_additivity_residual':abs(full['Tbulk']-xo['Tbulk']-ro['Tbulk'])})
    results['single_motif_checks']=motif
    checks={
      'max_replay_TqT_error':max(x['TqT_error'] for x in replay),
      'max_replay_XR_error':max(x['XR_error'] for x in replay),
      'max_replay_delta_error':max(x['delta_error'] for x in replay),
      'max_replay_flux_residual':max(x['flux_residual'] for x in replay),
      'max_analytic_derivative_identity_residual':max(x['full']['derivative_identity_residual'] for x in rows),
      'max_schur_derivative_source_residual':max(x['full']['schur_derivative_source_residual'] for x in rows),
      'max_alpha_finite_difference_error':max(f['abs_error'] for x in rows for f in x['finite_difference_alpha']),
      'max_physical_scale_finite_difference_error':max(x['physical_scale_fd_abs_error'] for x in rows),
      'max_motif_X_self_minus_Xonly':max(abs(x['full_X_self_minus_Xonly']) for x in motif),
      'max_motif_delta_minus_RR':max(abs(x['full_minus_X']-x['RR']) for x in motif),
      'max_motif_R_self_minus_Ronly':max(abs(x['full_R_self_minus_Ronly']) for x in motif),
      'max_motif_qT_redressed_additivity_residual':max(x['qT_redressed_additivity_residual'] for x in motif),
      'max_motif_all_bulk_redressed_additivity_residual':max(x['all_bulk_redressed_additivity_residual'] for x in motif),
      'max_motif_Xonly_bulk':max(x['Xonly_bulk'] for x in motif),
      'max_motif_Ronly_bulk':max(x['Ronly_bulk'] for x in motif),
      'max_motif_Xonly_qT':max(x['Xonly'] for x in motif),
      'max_motif_Ronly_qT':max(x['Ronly'] for x in motif),
      'max_motif_uncoupled_parity_unitarity_residual':max(x[k] for x in motif for k in ['Xonly_uncoupled_parity_unitarity_residual','Ronly_uncoupled_parity_unitarity_residual']),
      'max_motif_half_channel_flux_residual':max(x[k] for x in motif for k in ['Xonly_half_channel_flux_residual','Ronly_half_channel_flux_residual']),
    }
    results['checks']=checks
    (HERE/'interference_reaudit.json').write_text(json.dumps(results,indent=2))
    with (HERE/'alpha_sensitivity.csv').open('w',newline='') as f:
        fields=['case','omega','cq','TqT','XR','RR','full_minus_X','X_self_loading_change','radial_and_cross_gain','alpha_derivative','direct_source_derivative','redressing_derivative','physical_common_scale_derivative']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for x in rows:
            w.writerow(dict(case=x['case'],omega=x['omega'],cq=x['cq'],TqT=x['full']['TqT'],XR=x['full']['XR'],RR=x['full']['RR'],full_minus_X=x['full_minus_X'],X_self_loading_change=x['X_self_loading_change'],radial_and_cross_gain=x['radial_and_cross_gain'],alpha_derivative=x['full']['derivative'],direct_source_derivative=x['full']['direct_source_derivative'],redressing_derivative=x['full']['redressing_derivative'],physical_common_scale_derivative=x['physical_common_scale_derivative']))
    print(json.dumps(checks,indent=2))
    for x in rows:
        print(x['case'],x['omega'],x['cq'],'alpha slope',x['full']['derivative'],'direct',x['full']['direct_source_derivative'],'feedback',x['full']['redressing_derivative'],'physical scale',x['physical_common_scale_derivative'])

if __name__=='__main__': audit()
