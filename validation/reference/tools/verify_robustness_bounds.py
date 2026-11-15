"""Replay Section VIII's positive-overlap example and its fixed-frequency tolerance bounds.

The exact rational ramp calculation starts from the selected bulk matrices. The
reference calculation starts from its archived outward-rounded vertex bounds.
Both use the saved contour ledger and its physical-validation record. With
--campaign, also check the 222 original contour state and settings files.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def witness():
    bonds = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0),
                      (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)])
    directions = bonds / np.linalg.norm(bonds, axis=1)[:, None]
    seed = np.array([0.6318952450, -0.6028746482, -0.9641061600])  # Rounded SM3 locator.

    def modes(k):
        matrix = np.einsum('b,bi,bj->ij', 2 - 2*np.cos(bonds @ k), directions, directions)
        return np.linalg.eigh(matrix)

    lo, hi = 1 - 1e-6, 1 + 1e-6
    require(modes(lo*seed)[0][0] < 1.44 < modes(hi*seed)[0][0], 'Unbracketed witness')
    for _ in range(40):
        mid = (lo + hi)/2
        if modes(mid*seed)[0][0] < 1.44:
            lo = mid
        else:
            hi = mid
    k = seed * ((lo+hi)/2)
    squared, vectors = modes(k)
    polarization = vectors[:, 0]
    gradient = np.einsum('b,bi->i', 2*np.sin(bonds@k)*(directions@polarization)**2, bonds)
    # Caroli-normalized free right-going input. A common phase has no effect.
    k_line = 2*np.arcsin(1.2 / (2*1.3))
    line = np.exp(1j*k_line*np.arange(7)) / np.sqrt(2*1.69*np.sin(k_line))
    overlap = 0j
    for segment, strength in enumerate((0.1, 0.2, 0.2, 0.2, 0.2, 0.1)):
        for y, z in product((-1, 1), repeat=2):
            point = np.array([segment+0.5, y/2, z/2])
            corner_force = sum(strength * line[segment+j] * dx*np.array([dx,y,z])/3
                               for j, dx in enumerate((1, -1)))
            overlap += np.exp(-1j*point@k) * (polarization @ corner_force)
    result = [np.linalg.norm(gradient), *np.diff(squared), abs(overlap)**2]
    np.testing.assert_allclose(result, [1.8578005, 0.1007240, 2.6995796, 0.7492079],
                               atol=6e-8, rtol=0)
    require(abs(squared[0] - 1.44) < 1e-12, 'Shell residual')
    return dict(zip(('gradient_norm', 'lower_split', 'upper_gap', 'squared_source_overlap'), result))


def rational_ball(text):
    match = re.fullmatch(r'\[([^ ]+) \+/- ([^\]]+)\]', text)
    require(match is not None, 'Malformed printed Arb interval')
    center, radius = map(F, match.groups())
    return center-radius, center+radius


def positive_definite(matrix):
    """Exact symmetric elimination; strict positive pivots prove positive definiteness."""
    a = [row[:] for row in matrix]
    for k in range(len(a)):
        require(a[k][k] > 0, 'Non-positive exact pivot')
        for i in range(k+1, len(a)):
            multiplier = a[i][k] / a[k][k]
            for j in range(k+1, len(a)):
                a[i][j] -= multiplier * a[k][j]


def reference_floor(ref, metric):
    # Inflate the archived source-norm upper bound and independently certify
    # ||F|| < 1.61 from the exact dyadic lower-pair matrix (F*F = 2 Kq Q).
    require(len(ref['vertices']) == 64, 'Incorrect reference vertex count')
    require({tuple(v['signs']) for v in ref['vertices']} == set(product((-1,1), repeat=6)),
            'Incomplete reference box')
    caps = {k: max(F(v['bounds'][k]) for v in ref['vertices'])
            for k in ref['uniform_norm_upper_rationals']}
    require(caps == {k:F(v) for k,v in ref['uniform_norm_upper_rationals'].items()},
            'Reference vertex maxima')
    norm_bound = F(161,100)
    positive_definite([[int(i==j)*norm_bound**2 - 2*F(169,100)*x for j,x in enumerate(row)]
                       for i,row in enumerate(metric)])
    t_lower, _ = rational_ball(ref['T0_enclosure'])
    unit = 10**24
    amp = F(isqrt((t_lower.numerator * unit**2) // t_lower.denominator), unit)
    r, h = caps['r'], caps['h']
    amp -= caps['d'] + caps['w'] * caps['v'] / (1-r)
    eg, eq = F(1,50), F(1,200)
    rho = r+h*eg
    require(rho < 1 and caps['beta']+caps['hcore']*eg < 1, 'Reference inverse bound')
    _, z = rational_ball(ref['hypothetical_kernel_error_conclusion']['source_norm_upper'])
    state = h*eg*z/(1-rho)
    amplitude = amp - norm_bound*state
    require(amplitude > 0, 'Reference lower amplitude')
    floor = amplitude**2 - 2*eq*(z+state)**2
    require(floor > F(1037,10000), 'Reference whole-box floor')
    return float(floor)


def check_campaign(folder, aud):
    campaign = read(folder/'CAMPAIGN.json')
    jobs = {j['id']:j for j in campaign['jobs']}
    require(len(jobs) == len(aud) == 222, 'Unexpected campaign size')
    for a in aud:
        job = jobs[a['id']]
        for data in ('settings', 'state'):
            require(digest(folder/job[data]) == a[data+'_sha256'], 'Raw campaign binding')
        require(a['status']=='rigorous_complete_tile_enclosure', 'Incomplete tile')
        require(job['kind']==a['kind'] and F(job['area'])*job['multiplicity']==F(a['physical_parameter_area']),
                'Inconsistent campaign tile area')


def verify(root, campaign):
    e=root/'evidence/robustness_bounds'
    chain=read(root/'evidence/vi_results/RAMP_CERTIFICATE_CHAIN.json')
    paths=e/'ramp_inputs'
    ledger=read(e/'FINAL_CAMPAIGN_LEDGER.json')
    verification=read(e/'VERIFICATION.json')
    ref=read(e/'REFERENCE_VERTEX_BOUNDS.json')
    provenance=verification['bridge']['bindings']
    identities={
        e/'FINAL_CAMPAIGN_LEDGER.json':provenance['global_gate_sha256'],
        e/'TARGETS.json':provenance['targets_sha256'],
        e/'REFERENCE_VERTEX_BOUNDS.json':provenance['fresh_original_proof_sha256'],
        e/'VERIFICATION.json':chain['files']['R4_fresh_saved_replay']['sha256'],
        paths/'inputs/REFERENCE_KERNEL.npz':chain['files']['bundled_cached_kernel']['sha256'],
        paths/'verify_point_rational.py':chain['files']['point_verifier']['sha256'],
        paths/'verify_box_rational.py':chain['files']['box_verifier']['sha256'],
        paths/'reproduce_response_certificates.py':chain['files']['portable_response_replayer']['sha256'],
    }
    require(all(digest(f)==v for f,v in identities.items()), 'Evidence hash binding')
    require(verification['exact_lattice_floor_certified'] and verification['physical_replay_performed']
            and verification['saved_ledger_reproduced_exactly'] and
            all(s['status']=='passed' for s in verification['stages']), 'Unverified archived domain')
    require(len(ledger['job_audits']) == 222, 'Unexpected contour tile count')
    if campaign is not None:
        check_campaign(campaign,ledger['job_audits'])
    targets=read(e/'TARGETS.json')['targets']
    with np.load(paths/'inputs/REFERENCE_KERNEL.npz') as data:
        metric = [[F(float(x)) for x in row] for row in data['Im_qT'][:12,:12]]
        for name, limit in [('g',F(1,250)),('Q',F(1,500))]:
            result=ledger['results'][name]
            require(result['all_new_jobs_complete'] and all(F(x['fallback_area'])==0 and
                    F(x['refined_area'])==1 for x in result['coverage'].values()), 'Chart coverage')
            require(F(result['operator_budget']['scaled_operator_upper']) < limit, 'Kernel norm budget')
            actual=data['G0' if name=='g' else 'Im_qT'][:12,:12]
            for i,j in product(range(12),repeat=2):
                z=actual[i,j]
                require([F(float(z.real)),F(float(z.imag))] == list(map(F,targets[name][i][j])),
                        'Incorrect kernel center')
    floor=reference_floor(ref,metric)
    with tempfile.TemporaryDirectory() as tmp:
        output=Path(tmp)/'response'
        subprocess.run([sys.executable,str(paths/'reproduce_response_certificates.py'),
                        '--output-dir',str(output)],stdout=subprocess.PIPE,text=True,check=True)
        point=read(output/'POINT_RATIONAL_CERTIFICATE.json')['tighter_conditional_interval']
        box=read(output/'BOX_RATIONAL_CERTIFICATE.json')
        lower,upper=F(point['lower_rational']),F(point['upper_rational'])
        all_box=F(box['tight_kernel_box_lower']['lower_rational'])
        require(lower==F(chain['nominal_infinite_lattice_interval']['lower_rational'])
                and upper==F(chain['nominal_infinite_lattice_interval']['upper_rational'])
                and all_box==F(chain['whole_six_motif_plus_minus5_percent_infinite_lattice_floor']['lower_rational']),
                'Ramp replay inconsistent with saved proof chain')
        require(F('0.732166') < lower < upper < F('0.739652') and all_box > F('.69050'),
                'Manuscript quoted ramp bounds')
    result={'status':'PASS', 'weak_reference_witness':witness(), 'saved_contour_jobs':len(ledger['job_audits']),
            'reference_auxiliary_vertices':64, 'independent_reference_lower_estimate':floor,
            'reference_reported_strict_floor':'1037/10000',
            'ramp_nominal_rational_endpoints':[str(lower),str(upper)],'ramp_replayed_box_floor':str(all_box),
            'ramp_auxiliary_vertices':box['vertex_count'], 'fresh_contour_integration':False,
            'scope':'The contour ledger and reference vertex enclosures are retained inputs. The ramp point and box rational response and the reference final lower bound are recalculated.'}
    if campaign is not None: result['original_settings_and_checkpoint_file_hashes_matched']=2*len(ledger['job_audits'])
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--campaign',type=Path,help='Optional unpacked full validation analysis/full_domain_r4 directory')
    opts=ap.parse_args();print(json.dumps(verify(opts.root,opts.campaign),indent=2))
