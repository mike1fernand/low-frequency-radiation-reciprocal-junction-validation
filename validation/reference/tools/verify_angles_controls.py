"""Recompute Section VI.B/C observables and compare recorded shell/sector outputs.

Supply --evidence PATH for a directory containing the original records in SOURCES.json.
The continuum integration is independent. Finite-frequency cone checks reuse stored
branch integrals; none of these operations independently recomputes the lattice PV.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

STRENGTHS = np.array([.1, .2, .2, .2, .2, .1])
OMEGA, C_LINE = 1.2, 1.3
KAPPA = 2*np.arcsin(OMEGA/(2*C_LINE))


def read(evidence, filename):
    return json.loads((evidence/filename).read_text())


def sphere_directions(ntheta):
    mu, weights = leggauss(ntheta)
    nphi = 2*ntheta
    phi = np.arange(nphi)*2*np.pi/nphi
    radius = np.sqrt(1-mu**2)[:, None]
    directions = np.stack(np.broadcast_arrays(mu[:, None], radius*np.cos(phi),
                                                 radius*np.sin(phi)), axis=-1)
    return directions.reshape(-1, 3), np.repeat(weights, nphi)*2*np.pi/nphi


def continuum_shares(ntheta):
    """Project both Cartesian force patterns on their own constant-frequency spheres."""
    corners, f_x, f_r = [], [], []
    q = np.exp(1j*KAPPA*np.arange(7))
    for i,g in enumerate(STRENGTHS):
        for y,z in itertools.product((-1,1),repeat=2):
            force = np.zeros(3,dtype=complex)
            for j,dx in enumerate((1,-1)):
                unit_bond = np.array([dx,y,z])/np.sqrt(3)
                force += g*q[i+j]*unit_bond[0]*unit_bond
            corners.append([i+.5,y/2,z/2])
            f_x.append([force[0],0,0])
            f_r.append([0,force[1],force[2]])
    points=np.array(corners)
    directions, weights=sphere_directions(ntheta)
    channel_power=[]
    for mode,speed in [('S',1.),('P',np.sqrt(3))]:
        phase=np.exp(-1j*(OMEGA/speed)*directions@points.T)
        Fx,Fr=phase@f_x,phase@f_r
        powers=[]
        for force in (Fx,Fr,Fx+Fr):
            parallel=(force*directions).sum(axis=1)
            projected=parallel[:,None]*directions if mode=='P' else force-parallel[:,None]*directions
            powers.append(weights@np.sum(abs(projected)**2,axis=1)/speed**3)
        channel_power.append(powers)
    s,p=np.array(channel_power)
    np.testing.assert_allclose([*s[:2],s[2]-s[:2].sum(),*p[:2],p[2]-p[:2].sum()],
                               [9.81946085,.35104365,3.29640150,1.40339500,.03147380,-.25019629],
                               rtol=0,atol=6e-8)
    return np.r_[s/(s+p),s[:2].sum()/(s[:2]+p[:2]).sum()]


def force_terms(evidence):
    """Reassemble both square energies, solve the full state, and split its source."""
    vectors=np.zeros((12,7))
    for i in range(6):
        vectors[2*i,i:i+2] = [-1,-1]
        vectors[2*i+1,i:i+2] = [-1,1]
    strength=np.repeat(2*STRENGTHS/3,2)
    scale=np.tile([1,np.sqrt(2)],6)
    V=(strength*scale)[:,None]*vectors
    U=np.diag(strength*scale**2)
    Kself=vectors.T@(strength[:,None]*vectors)
    with np.load(evidence/'REFERENCE_KERNEL.npz') as kernel:
        G=kernel['G0'][:12,:12]
        shell=kernel['Im_qT'][:12,:12]
        total_shell=kernel['Im_total'][:12,:12]
    endpoint=np.eye(7)[0]
    lead=np.diag([1,0,0,0,0,0,1])
    kline=C_LINE**2
    bare=OMEGA**2*np.eye(7)-kline*(2*np.eye(7)-np.eye(7,k=1)-np.eye(7,k=-1))+kline*np.exp(1j*KAPPA)*lead
    amp=np.sqrt(2*kline*np.sin(KAPPA))*endpoint
    corner_drive=np.linalg.solve(np.eye(12)+U@G,V)
    chi=np.linalg.solve(bare-Kself+V.T@G@corner_drive,amp)
    fields=[]
    for parity in (0,1):
        piece=np.zeros_like(V)
        piece[parity::2]=V[parity::2]
        fields.append(np.linalg.solve(np.eye(12)+U@G,piece@chi))
    x,r=fields
    weights=np.array([2*np.vdot(x,shell@x).real,2*np.vdot(r,shell@r).real,4*np.vdot(x,shell@r).real])
    np.testing.assert_allclose(weights,[.10147813,.00179177,.02238935],rtol=0,atol=6e-9)
    np.testing.assert_allclose(weights.sum(),2*np.vdot(x+r,shell@(x+r)).real,rtol=0,atol=2e-13)
    weak=V@np.linalg.solve(bare,amp)
    weak_share=np.vdot(weak,shell@weak).real/np.vdot(weak,total_shell@weak).real
    assert round(weak_share,4)==.9271
    return weights, weak_share


def cone_results(evidence):
    fine=read(evidence,'integral_fine.json')['cone_branch_powers']
    other={row['half_angle_deg']:row['forward'] for row in read(evidence,'INDEPENDENT_ANGULAR_CHECK.json')['cones']}
    finite,low=[],[]
    for angle,qT,qL,pct in [(10,9.17e-6,7.61e-4,'1.19'),(20,1.99e-4,2.69e-3,'6.89'),(30,2.87e-3,5.09e-3,'36.1')]:
        f=fine[f'group_forward_{angle}']['full']
        part=sum(f[:2]),f[2]
        np.testing.assert_allclose(part,[other[angle]['qT'],other[angle]['qL']],rtol=3e-10,atol=0)
        assert [float(format(v,'.3g')) for v in part]==[qT,qL]
        percent=100*part[0]/sum(part)
        assert format(percent,'.3g')==pct
        finite.append(percent)
        c=np.cos(np.deg2rad(angle))
        T=2*np.sqrt(3)**3*(2-3*c+c**3)/4
        L=(1-c**3)/2
        low.append(100*T/(T+L))
    assert [format(v,'.3g') for v in low]==['7.39','24.6','43.3']
    radial=fine['group_forward_10']['R_only']
    assert format(sum(radial[:2]),'.3g')=='8.84e-07' and format(100*sum(radial[:2])/sum(radial),'.3g')=='99.2'
    return finite,low


def recorded_controls(evidence,weights):
    sweep=read(evidence,'R_only_frequency_N48.json')['rows']
    assert [r['omega'] for r in sweep]==[.8,1.,1.2,1.4,1.7]
    assert all(r['TT_modal_R_only']>0 and r['TT_modal_R_only']<=r['Tbulk_R_only'] for r in sweep)
    data=read(evidence,'DOWNSTREAM_TETRA_N192.json')
    for name,expect in [('full',(.1366,.1256,.9194)),('X_only',(.1152,.1022,.8873)),('R_only',(.002430,.002231,.9179))]:
        row=data['controls'][name]
        actual=(row['Tbulk'],row['TT_modal'],row['TT_modal']/row['Tbulk'])
        assert [float(format(x,'.4g')) for x in actual]==list(expect)
    steps=data['alpha']['rows'];alpha=[r['alpha'] for r in steps];values=[r['TT_modal'] for r in steps]
    np.testing.assert_allclose(alpha,np.linspace(0,1,21),rtol=0,atol=1e-15)
    assert min(np.diff(values))>0 and round(100*(values[-1]-values[0]),2)==2.34
    other_cross=data['same_state']['modal_qT']['XR_cross']
    assert round(100*abs(other_cross-weights[2])/weights[2],3)==.026
    detune=read(evidence,'interference_reaudit.json')
    results={}
    for r in detune['selected_rows']:
        if (r['omega'],r['cq'])!=(.8,.8):continue
        f,x=r['full'],r['X_only']
        delta=f['TqT']-x['TqT']
        assert delta<0<f['XR'] and round(100*delta,2)==-1.18
        np.testing.assert_allclose(delta,f['XX']-x['TqT']+f['RR']+f['XR'],rtol=0,atol=1e-12)
        results[r['case']]=[f['TqT'],x['TqT'],f['XR']]
    assert len(results)==2
    single=detune['single_motif_checks']
    assert len(single)==45
    assert max(abs(r['qT_redressed_additivity_residual']) for r in single)<2e-17
    assert max(max(r['Xonly_bulk'],r['Ronly_bulk']) for r in single)<.5
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path,default=Path(__file__).resolve().parents[1]/'evidence/controls')
    e=parser.parse_args().evidence
    for name,source in read(e,'SOURCES.json').items():
        assert hashlib.sha256((e/name).read_bytes()).hexdigest()==source['sha256'],name
    a,b=continuum_shares(16),continuum_shares(24)
    np.testing.assert_allclose(a,b,rtol=0,atol=1e-12)
    finite,low=cone_results(e)
    weights,weak=force_terms(e)
    detuned=recorded_controls(e,weights)
    print(json.dumps({'status':'PASS','finite_forward_qT_shares_percent':finite,
                      'low_frequency_forward_qT_shares_percent':low,'selected_kernel_full_state_XX_RR_XR':weights.tolist(),
                      'continuum_shares_X_R_coherent_incoherent':b.tolist(),'discrete_weak_coupling_emitted_qT_share':weak,
                      'detuned_full_Xonly_cross':detuned,
                      'scope':'Archived finite-cone branch integrals and independent angular output; direct prescribed-force continuum quadrature; new solve with saved bulk response; 45 one-motif receipts; no new exact-lattice enclosure.'},indent=2))
