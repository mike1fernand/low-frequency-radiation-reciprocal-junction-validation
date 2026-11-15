"""Recompute the increasing-ramp packet selectors, modal averages and stationary comparison."""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

EVIDENCE = Path(__file__).resolve().parents[1] / 'evidence' / 'ramp_packet'
CASES = {
    'first_box_56': 'CANDIDATE_PACKET_RESULT.json',
    'primary_80': 'CANDIDATE_PACKET_BOX80_RESULT.json',
    'half_step_80': 'CONTROL_DT002_BOX80_RESULT.json',
    'larger_96': 'CONTROL_DT004_BOX96_RESULT.json',
}


def read(name):
    return json.loads((EVIDENCE / name).read_bytes())


def require(ok, message):
    if not ok:
        raise ValueError(message)


def verify_archive_inputs():
    originals = read('PROVENANCE.json')
    for name, record in originals['members'].items():
        require(hashlib.sha256((EVIDENCE / name).read_bytes()).hexdigest() == record['sha256'],
                'Archived input changed: ' + name)


def select(rows):
    accepted = []
    for row in rows:
        ok = (row['link'] <= 0.0025 and abs(row['J']) <= 0.0005 and row['shell'] <= 0.002
              and abs(row['periodic_difference']) <= 0.0025 and row['drift'] <= 1e-5)
        require(ok == row['passed'], 'Inconsistent selector at time ' + str(row['t']))
        accepted.append(ok)
    times = np.array([row['t'] for row in rows])
    np.testing.assert_allclose(np.diff(times), 0.2, atol=1e-12)
    groups = [(flag, list(group)) for flag, group in itertools.groupby(enumerate(accepted), key=lambda x: x[1])]
    runs = [[times[group[0][0]], times[group[-1][0]]] for passed, group in groups if passed]
    # Endpoints are cadence ticks. Rank lengths in integer ticks to handle ties exactly.
    valid = [run for run in runs if round((run[1]-run[0]) / 0.2) >= 10]
    best = max(valid, key=lambda run: (round((run[1]-run[0]) / 0.2), -run[0])) if valid else None
    return best


def compare_case(record, is_control):
    window = select(record['rows'])
    if record['status'] == 'NO_ACCEPTED_WINDOW':
        require(window is None, 'The initial pilot unexpectedly passed')
        return {'passing_window': None}
    if is_control:
        require(record['selection_mode'] == 'frozen_control', 'Missing preselected control interval')
        window = record['config']['frozen_window']
        require(all(r['passed'] for r in record['rows'] if window[0] <= r['t'] <= window[1]),
                'Control fails a selector condition on its frozen interval')
    require(window == [57.4, 64.2], 'Wrong ramp interval')
    samples = np.array(record['branch_samples'])
    np.testing.assert_allclose(samples[:, 0], np.arange(287, 322) * 0.2, rtol=0, atol=1e-12)
    means = np.trapezoid(samples[:, 1:], x=samples[:, 0], axis=0) / (window[1]-window[0])
    np.testing.assert_allclose(means, record['mean_low_high_total'], rtol=0, atol=1e-14)
    return {'passing_window': window, 'sample_count': len(samples),
            'qT_allocation_mean_endpoints': means[:2].tolist(), 'midpoint': float(np.mean(means[:2]))}


def stationary():
    record = read('PACKET_SPECTRUM.json')
    q, v = np.loadtxt(EVIDENCE / 'initial_qv.txt').reshape(2, 120)
    j = np.arange(q.size) - 18
    q_expected = 0.002 * np.exp(-0.5 * (j/7)**2) * np.cos(2 * np.arcsin(1.2/2.6) * j)
    np.testing.assert_allclose(q, q_expected, atol=1e-17, rtol=0)
    k = 2 * np.pi * np.fft.fftfreq(q.size)
    omega = 2.6 * abs(np.sin(k/2))
    sense = np.sign(k)
    sense[q.size//2] = 0
    qhat = np.fft.fft(q)
    np.testing.assert_allclose(v, np.fft.ifft(-1j * sense * omega * qhat).real, atol=1e-17, rtol=0)
    energy = (abs(np.fft.fft(v))**2 + (omega * abs(qhat))**2)/(2*q.size)
    open_line = (v @ v + 1.69 * np.diff(q) @ np.diff(q)) / 2
    frequencies = next(f for f in record['families'] if f['label'] == 'surface_N160')
    pairs, contributions = [], []
    for row in frequencies['rows']:
        i = row['positive_index']
        require(row['negative_index'] == q.size-i, 'Wrong conjugate bin')
        require(abs(omega[i] - row['omega']) < 1e-14, 'Wrong stationary frequency')
        pair = energy[i] + energy[-i]
        pairs.append(pair)
        contributions.append(pair * row['states']['ramp']['qT'])
    require(len(pairs) == 14, 'Unexpected number of evaluated frequencies')
    np.testing.assert_allclose(sum(pairs), sum(energy[(omega >= .8) & (omega <= 1.6)]),
                               rtol=1e-13, atol=0)
    estimate = float(sum(contributions) / open_line)
    np.testing.assert_allclose(estimate, frequencies['qT_in_band_contribution_per_initial_open_total_energy'],
                               rtol=0, atol=1e-13)
    return {'evaluated_distinct_frequencies': len(pairs),
            'periodic_energy_coverage': float(sum(pairs)/sum(energy)),
            'periodic_to_open_line_energy_factor': float(sum(energy)/open_line),
            'qT_in_band_per_open_line_energy': estimate,
            'unevaluated_full_weight_allowance_per_open_line_energy': float((sum(energy)-sum(pairs))/open_line)}


def main():
    verify_archive_inputs()
    cases = {name: compare_case(read(record), name in ('half_step_80','larger_96'))
             for name, record in CASES.items()}
    weighted = stationary()
    direct = cases['primary_80']['midpoint']
    diff = 100 * (direct - weighted['qT_in_band_per_open_line_energy'])
    require(abs(diff-0.0149932255418) < 1e-10, 'Incorrect direct minus stationary comparison')
    print(json.dumps({'status':'PASS','cases':cases,'stationary_from_signed_bins':weighted,
                      'primary_minus_in_band_percentage_points':diff,
                      'scope':'Five-predicate re-selection, archived modal allocations, and reweighted stationary samples. Rerun the primary evolution separately with replay_ramp_packet.py.'}, indent=2))


if __name__ == '__main__':
    main()
