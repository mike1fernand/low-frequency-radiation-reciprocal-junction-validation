"""Independent radial-shell totals and ideal wave/group-direction cone powers.

No Monte Carlo, broadened resolvent, paraxial assumption or fitted theta law.
The reference radial census proves one root of each ordered branch on every ray.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from core import roots,states,density

P=Path(__file__).resolve().parent
def dirs(theta,phi,sign=1):
    theta,phi=np.broadcast_arrays(theta,phi)
    return np.stack([sign*np.cos(theta),np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi)],axis=-1)

def weights(rr,eta):
    vals={name:density(rr,t) for name,t in eta.items()}
    vals['same_cross']=vals['full']-vals['same_X']-vals['same_R']
    return vals

def run(nmu,nphi,scan_count,angles):
    start=time.monotonic();eta,ref=states();phi,wp=np.polynomial.legendre.leggauss(nphi)
    phi=(phi+1)*np.pi/8;wp=wp*np.pi/8*8 # full transverse circle by exact C4v symmetry
    mu,wm=np.polynomial.legendre.leggauss(nmu)
    mm,pp=np.meshgrid(mu,phi,indexing='ij');n=dirs(np.arccos(mm),pp)
    full={name:[] for name in list(eta)+['same_cross']};globaldiag=[]
    for branch in range(3):
        rr=roots(n,branch);wt=wm[:,None]*wp[None,:]
        for name,v in weights(rr,eta).items():full[name].append(float(np.sum(wt*v)))
        globaldiag.append({'branch':branch,'max_shell_residual':float(abs(rr['values'][...,branch]-1.44).max()),
                           'min_radial_derivative':float(rr['radial'].min())})
    # Globally scan polar rays in each symmetry-sector azimuth to discover all
    # finite-cone indicator components. Never infer support from a local node.
    scan=np.linspace(0,np.pi/2,scan_count)
    tt,pp=np.meshgrid(scan,phi,indexing='ij')
    cone={};scan_summary=[]
    localn=max(24,nmu//2)
    u,wu=np.polynomial.legendre.leggauss(localn)
    for branch in range(3):
        rr=roots(dirs(tt,pp),branch)
        cosg=rr['grad'][...,0]/np.linalg.norm(rr['grad'],axis=-1)
        for angle in angles:
            threshold=np.cos(np.deg2rad(angle));inside=cosg>=threshold
            transitions=np.diff(inside.astype(int),axis=0)
            exits=(transitions==-1).sum(axis=0);entries=(transitions==1).sum(axis=0)
            record={'branch':branch,'angle_degrees':angle,'exits_min':int(exits.min()),'exits_max':int(exits.max()),
                    'entries_max':int(entries.max()),'all_start_inside':bool(inside[0].all()),'all_equator_outside':bool((~inside[-1]).all()),
                    'opposite_hemisphere_cone_nodes':int((cosg<=-threshold).sum()),
                    'min_margin_beyond_2theta':float((threshold-cosg[scan>=min(np.pi/2,2*np.deg2rad(angle))]).min())}
            scan_summary.append(record)
            # This run may only use a simple boundary after explicitly observing
            # that structure over the whole polar range at every quadrature azimuth.
            if not (np.all(exits==1) and np.all(entries==0) and inside[0].all() and (~inside[-1]).all() and not np.any(cosg<=-threshold)):
                raise RuntimeError(f'Multiple cone preimages require general interval integration: {record}')
            index=np.argmax(transitions==-1,axis=0);lo=scan[index];hi=scan[index+1]
            for _ in range(36):
                mid=(lo+hi)/2;z=roots(dirs(mid,phi),branch,iterations=38)
                cg=z['grad'][...,0]/np.linalg.norm(z['grad'],axis=-1)
                above=cg>=threshold;lo=np.where(above,mid,lo);hi=np.where(above,hi,mid)
            boundary=(lo+hi)/2
            for kind,bd in [('group',boundary),('wave',np.full_like(phi,np.deg2rad(angle)))]:
                # integrate dmu, so the sphere's sin(theta) Jacobian is exact
                lower=np.cos(bd);mus=lower[None,:]+(u[:,None]+1)*(1-lower[None,:])/2
                wt=wu[:,None]*(1-lower[None,:])/2*wp[None,:]
                for side,sign in [('forward',1),('backward',-1)]:
                    r=roots(dirs(np.arccos(mus),phi[None,:],sign),branch)
                    key=f'{kind}_{side}_{angle:g}'
                    if key not in cone:cone[key]={name:[] for name in full}
                    for name,v in weights(r,eta).items():cone[key][name].append(float(np.sum(wt*v)))
            record['boundary_theta_range_degrees']=np.rad2deg([boundary.min(),boundary.max()]).tolist()
    error={name:{'qT_relative':(sum(full[name][:2])-ref[name]['qT'])/ref[name]['qT'],
                 'qL_relative':(full[name][2]-ref[name]['qL'])/ref[name]['qL']} for name in ref}
    return {'quadrature':{'nmu':nmu,'nphi_sector':nphi,'polar_scan_count':scan_count,'cone_local_nmu':localn},
            'total_branch_powers':full,'selected_kernel_reference':ref,'total_relative_difference_from_selected_kernel':error,
            'cone_branch_powers':cone,'global_root_diagnostics':globaldiag,'whole_polar_scan':scan_summary,
            'scope':'Finite-accuracy integration of full reference shells; fixed selected-kernel dressed sources; ideal angular power collection, no physical receiver loading. Scan is numerical, not a certified exclusion of all infinitesimal ray preimages.',
            'elapsed_seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--nmu',type=int,default=64);p.add_argument('--nphi',type=int,default=40)
    p.add_argument('--scan',type=int,default=257);p.add_argument('--angles',type=float,nargs='+',default=[2,5,10,20,30]);p.add_argument('--output',required=True)
    a=p.parse_args();out=run(a.nmu,a.nphi,a.scan,a.angles)
    (P/a.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ['quadrature','total_branch_powers','total_relative_difference_from_selected_kernel','elapsed_seconds']},indent=2))
