"""Rebuild the initial packet and solve both reciprocal contacts at its retained frequencies."""
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check(actual, expected, tolerance=1e-11):
    np.testing.assert_allclose(actual, expected, atol=tolerance, rtol=0)


def initial_packet(cfg):
    n, a, mass, stiffness = (cfg[k] for k in ("n", "a", "m_q", "K_q"))
    sites = a*np.arange(n) - cfg["center"]
    k = 2*np.pi*np.fft.fftfreq(n, d=a)
    max_frequency = 2*np.sqrt(stiffness/mass)
    carrier = 2*np.arcsin(cfg["omega0"]/max_frequency)/a
    omega = max_frequency*np.abs(np.sin(k*a/2))
    q = cfg["amplitude"] * np.exp(-0.5*(sites/cfg["sigma"])**2)*np.cos(carrier*sites)
    qhat = np.fft.fft(q)
    vhat = -1j*np.sign(k)*omega*qhat
    vhat[[0, n//2]] = 0
    v = np.fft.ifft(vhat)
    require(np.max(np.abs(v.imag)) < 1e-15, "Initial velocity must be real")
    weights = mass*(abs(vhat)**2 + omega**2*abs(qhat)**2)/(2*n)
    open_energy = mass*np.dot(v.real, v.real)/2 + stiffness*np.dot(np.diff(q), np.diff(q))/2
    check(weights.sum(), open_energy + stiffness*(q[0]-q[-1])**2/2, 1e-18)
    return q, omega, weights, open_energy


def contact_hessian(strengths):
    """Add the eight Cartesian central-spring extensions of each complete motif."""
    nline = len(strengths)+1
    size = nline + 2*len(strengths)
    hessian = np.zeros((size,size))
    for j, g in enumerate(strengths):
        for sign, line_site in ((1, j), (-1, j+1)):
            # Each of the four corners has the same projected extension in the centered basis.
            extension = np.zeros(size)
            extension[line_site] = -sign/np.sqrt(3)
            extension[nline+2*j] = sign/(2*np.sqrt(3))
            extension[nline+2*j+1] = 1/np.sqrt(6)
            hessian += 4*g*np.outer(extension,extension)
    return hessian


def powers(omega, hessian, G, shell_T, shell_L, kq, mass):
    """Simultaneous seven-line/twelve-corner solve with the physical outgoing lead limit."""
    n = 7
    line_k = 2*np.arcsin(omega/(2*np.sqrt(kq/mass)))
    gamma = 2*kq*np.sin(line_k)
    line = kq*(2*np.eye(n)-np.eye(n,k=1)-np.eye(n,k=-1))
    device = (mass*omega**2*np.eye(n) - line - hessian[:n,:n]).astype(complex)
    device[[0,-1],[0,-1]] += kq*np.exp(1j*line_k)
    coupling, core = hessian[n:,:n], hessian[n:,n:]
    response = np.block([[device, -coupling.T], [G@coupling, np.eye(len(G))+G@core]])
    incoming = np.zeros(len(hessian)); incoming[0] = np.sqrt(gamma)
    motion = np.linalg.solve(response,incoming)
    force = -coupling@motion[:n] - core@motion[n:]
    check(motion[n:], G@force, 1e-12)
    result = np.array([2*np.vdot(force,shell@force).real for shell in (shell_T,shell_L)] +
                      [gamma*abs(motion[n-1])**2, abs(1-1j*np.sqrt(gamma)*motion[0])**2])
    check(result.sum(), 1, 5e-12)
    require(np.all(result >= -1e-12), "Negative outgoing power")
    return result


def verify(root):
    binding = json.loads((root/'PROVENANCE.json').read_bytes())['sha256_by_original_member_suffix']
    def authenticated(data, member):
        require(sha256(data).hexdigest()==binding[member], f"Altered original input {member}")
        return data
    source = json.loads(authenticated((root/'PACKET_SPECTRUM.json').read_bytes(),'PACKET_SPECTRUM.json'))
    contract = json.loads(authenticated((root/'FREQUENCY_CONTRACT.json').read_bytes(),'inputs/FREQUENCY_CONTRACT.json'))
    cfg = contract['preparation']
    q,omega,w,eopen = initial_packet(cfg)
    poly = np.argsort(omega)
    percentiles = omega[poly[np.searchsorted(np.cumsum(w[poly])/w.sum(),[.01,.99])]]
    check(percentiles,[.93175667,1.47265622],6e-9)
    positive_bins = np.arange(1,cfg['n']//2)
    lo,hi = cfg['band']
    bins = positive_bins[(omega[positive_bins]>=lo)&(omega[positive_bins]<=hi)]
    pairs = (w[bins]+w[-bins])/w.sum()
    coverage = float(pairs.sum())
    check(coverage,contract['inside_coverage'],2e-15)
    check(w.sum()/eopen,source['periodic_to_initial_open_total_energy_factor'],2e-14)
    check(w[cfg['n']//2]/w.sum(),2.5768734723116623e-7,1e-17)
    require(w[0]==0, "Unexpected zero-mode weight")
    check(bins,[r['positive_k_index'] for r in contract['central_nodes']],0)
    check(pairs,[r['weight_over_full_periodic_energy'] for r in contract['central_nodes']],1e-15)
    profiles = {'reference': [.1,.2,.2,.2,.2,.1], 'ramp': np.arange(1,7)/7}
    blocks = {name:contact_hessian(g) for name,g in profiles.items()}
    family_results = {}
    with ZipFile(root/'kernels.zip') as archive:
        require(len(archive.namelist())==len(bins)*len(source['families']), "Missing or extra kernels")
        for family in source['families']:
            sums = {name:np.zeros(4) for name in profiles}
            method = family['method']
            for b,weight,node in zip(bins,pairs,family['rows']):
                check(omega[b],node['omega'],1e-13)
                member = f"{family['label']}/dft_{b:03d}.npz"
                raw = authenticated(archive.read(member),'inputs/'+member)
                with np.load(BytesIO(raw),allow_pickle=False) as data:
                    G = data['G0']
                    keys = [f'Im_branch_{i}' if method=='surface' else f'Im{i}' for i in range(3)]
                    qT, qL = data[keys[0]]+data[keys[1]],data[keys[2]]
                    check((G-G.conj().T)/(2j),qT+qL,2e-12)
                    for name,block in blocks.items():
                        calculated = powers(omega[b],block,G,qT,qL,cfg['K_q'],cfg['m_q'])
                        state = node['states'][name]
                        check(calculated,[state[k] for k in ['qT','qL','Tline','Rq']],1e-10)
                        sums[name] += weight*calculated
            for result in sums.values():
                check(result.sum(),coverage,3e-12)
            central = {name:float(total[0]) for name,total in sums.items()}
            check(central['reference'],family['original_reference_qT_in_band_contribution'])
            check(central['ramp'],family['new_ramp_full_energy_in_band_contribution']['qT'])
            family_results[family['label']] = {
                'reference_qT_per_periodic_packet_energy':central['reference'],
                'ramp_qT_per_open_line_energy':central['ramp']*w.sum()/eopen,
            }
    selected = family_results['surface_N160']
    outside = 1-coverage
    return {'status':'PASS','direct_stationary_contact_solves':len(profiles)*len(bins)*len(family_results),
            'positive_bin_range':[int(bins.min()),int(bins.max())],
            'distinct_evaluated_frequencies':len(bins),'signed_grid_size':cfg['n'],
            'energy_quantiles_01_99':percentiles.tolist(),'energy_coverage':coverage,
            'nyquist_fraction':float(w[cfg['n']//2]/w.sum()),'periodic_to_open_factor':float(w.sum()/eopen),
            'omitted_weight_at_most':outside,'reference_qT_surface_tail_only_interval':
            [selected['reference_qT_per_periodic_packet_energy'],selected['reference_qT_per_periodic_packet_energy']+outside],
            'stationary_in_band_contributions':family_results,
            'scope':'Independent coupled line-and-corner solves from archived on-frequency Green and branch matrices; no fresh bulk quadrature. Omitted-weight intervals take computed in-band values as given.'}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    base=Path(__file__).resolve().parents[1]/'evidence/packet_spectrum'
    parser.add_argument('--evidence',type=Path,default=base)
    print(json.dumps(verify(parser.parse_args().evidence),indent=2))
