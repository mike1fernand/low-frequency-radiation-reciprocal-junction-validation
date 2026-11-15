#!/usr/bin/env python3
"""Bounded engine checks of the derived old-material/H18 bulk correspondence.

No junction solve, quantum/clock identification, production trajectory, or
claim of physical ancestry is performed. The proof is in CONTROLLED_BULK_BRIDGE.md.
"""
from pathlib import Path
import hashlib
import json
import sys
import tempfile
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
INPUTS = HERE.parent/'recheck/purity_inputs'
if not INPUTS.is_dir():
    INPUTS = HERE/'package/paper1-significance/recheck/purity_inputs'


def old_symbol(theta, omega2):
    # S0335's centered divergence and three componentwise axial differences.
    s = np.sin(theta)
    Q = float(np.sum(4*np.sin(theta/2)**2))
    return omega2*(Q*np.eye(3)+2*np.outer(s, s))


def main():
    known = {x['path']: x for x in json.loads((INPUTS/'INPUT_PROVENANCE.json').read_text())['files']}
    payload = (INPUTS/'model_api_v1_1.py').read_bytes()
    assert hashlib.sha256(payload).hexdigest() == known['model_api_v1_1.py']['sha256']
    rng = np.random.default_rng(20261107)
    directions = np.vstack([np.eye(3), [1, 1, 0], [1, 1, 1], [1, 2, 3], rng.normal(size=(18, 3))])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    fixtures = [(1., 1., 1.), (.7, 1.3, .8)]  # K, m, a: nonunit fixture detects convention errors.
    records, flows, tails = [], [], []
    max_identity_error = 0.
    max_metric_residual = 0.
    sys.dont_write_bytecode = True
    with tempfile.TemporaryDirectory(prefix='p1-bulk-bridge-') as tmp:
        (Path(tmp)/'model_api_v1_1.py').write_bytes(payload)
        sys.path.insert(0, tmp)
        from model_api_v1_1 import ModelParameters, bulk_dynamical_matrix
        for K, mass, a in fixtures:
            par = ModelParameters(a=a, m_bulk=mass, k_axis=K, k_diag=K)
            omega = np.sqrt(K/mass)
            assert abs(par.c_bulk_T_longwave-a*omega) < 1e-14
            theta_set = [rho*d for rho in [.01, .05, .2, .5, 1., 1.4] for d in directions]
            theta_set += [np.array([np.pi, 0, 0]), np.array([np.pi]*3)]
            theta_set += list(rng.uniform(-np.pi, np.pi, (48, 3)))
            for theta in theta_set:
                physical_k = theta/a
                assert np.allclose(a*physical_k, theta, rtol=1e-14, atol=1e-14)
                D0 = old_symbol(theta, omega**2)
                D1 = bulk_dynamical_matrix(theta, par)
                q = 4*np.sin(theta/2)**2
                Q = float(np.sum(q))
                exact = omega**2*np.diag(q*(2*q-Q)/2)
                identity = np.linalg.norm(D1-D0-exact, 2)/omega**2
                max_identity_error = max(max_identity_error, float(identity))
                assert identity < 3e-14
                delta_norm = np.linalg.norm(D1-D0, 2)/omega**2
                assert delta_norm <= Q**2/2+3e-14
                eig, vec = np.linalg.eigh(D0)
                invS = (vec/np.sqrt(eig))@vec.T
                relative = np.linalg.eigvalsh(invS@D1@invS)
                assert relative.min() >= .25-1e-11
                assert relative.max() <= 7+1e-11
                rho = float(np.linalg.norm(theta))
                if rho < np.sqrt(2):
                    assert relative.min() >= 1-rho**2/2-1e-10
                    assert relative.max() <= 1+rho**2/2+1e-10
                records.append({'fixture': [K, mass, a], 'theta': theta.tolist(),
                    'rho': rho, 'Q': Q, 'scaled_identity_error': float(identity),
                    'relative_eigenvalues': relative.tolist(),
                    'scaled_symbol_error': float(delta_norm), 'Q_bound': Q**2/2})

            def flow_difference(theta, time):
                D0 = old_symbol(theta, omega**2)
                D1 = bulk_dynamical_matrix(theta, par)
                eig, vec = np.linalg.eigh(D0)
                S = (vec*np.sqrt(eig))@vec.T
                invS = (vec/np.sqrt(eig))@vec.T
                A0 = np.block([[np.zeros((3, 3)), S], [-S, np.zeros((3, 3))]])
                A1 = np.block([[np.zeros((3, 3)), S], [-D1@invS, np.zeros((3, 3))]])
                W1 = np.block([[invS@D1@invS, np.zeros((3, 3))], [np.zeros((3, 3)), np.eye(3)]])
                nonlocal max_metric_residual
                residual = np.linalg.norm(A1.T@W1+W1@A1, 2)/omega
                max_metric_residual = max(max_metric_residual, float(residual))
                assert residual < 1e-12
                F0, F1 = expm(time*A0), expm(time*A1)
                assert np.linalg.norm(F0.T@F0-np.eye(6), 2) < 1e-12
                assert np.linalg.norm(F1.T@W1@F1-W1, 2) < 1e-11
                _, sing, right = np.linalg.svd(F1-F0)
                x0 = right[0]
                u0, v0 = invS@x0[:3], x0[3:]
                return float(sing[0]), {'D_old': D0.tolist(), 'D_H18': D1.tolist(),
                    'initial_u': u0.tolist(), 'initial_p': (mass*v0).tolist(),
                    'energy_coordinate_initial': x0.tolist(), 'metric_residual': float(residual)}

            for rho in [.1, .2, .5]:
                delta = rho**2/2
                kappa = np.sqrt((1+delta)/(1-delta))
                epsilon = omega*rho**3/2
                for direction in directions[:6]:
                    theta = rho*direction
                    for scaled_time in [2., 10.]:
                        time = scaled_time/omega
                        error, details = flow_difference(theta, time)
                        bound = kappa*epsilon*time
                        assert error <= bound+1e-12
                        flows.append({'fixture': [K, mass, a], 'theta': theta.tolist(),
                            'structure_time': time, 'omega_structure_time': scaled_time,
                            'operator_flow_error': error, 'bound': float(bound), **details})

            rho = .2
            delta = rho**2/2
            kappa = np.sqrt((1+delta)/(1-delta))
            epsilon = omega*rho**3/2
            low = rho*np.array([1., 2., 3.])/np.sqrt(14)
            high = np.array([np.pi/2, 0, 0])
            for scaled_time in [2., 10.]:
                time = scaled_time/omega
                low_error, low_details = flow_difference(low, time)
                high_error, high_details = flow_difference(high, time)
                for fraction in [0., .0001, .01]:
                    error = np.sqrt((1-fraction)*low_error**2+fraction*high_error**2)
                    bound = np.sqrt((1-fraction)*(kappa*epsilon*time)**2+fraction*(np.sqrt(28)+1)**2)
                    assert error <= bound+1e-12
                    tails.append({'fixture': [K, mass, a], 'omega_structure_time': scaled_time,
                        'old_total_energy_tail_fraction': fraction, 'flow_error': float(error),
                        'bound': float(bound), 'low_theta': low.tolist(), 'high_theta': high.tolist(),
                        'low_unit_energy_coordinate_norm_initial_state': low_details,
                        'high_unit_energy_coordinate_norm_initial_state': high_details})

        # Discriminating controls: local-contact mixing does not preserve the band,
        # the exact lattice symbols are unequal, and a wrong material dictionary fails.
        theta = np.array([np.pi, 0, 0])
        par = ModelParameters()
        zone = bulk_dynamical_matrix(theta, par)-old_symbol(theta, 1.)
        assert np.allclose(np.diag(zone), [8., 0., 0.])
        small = np.array([.01, 0, 0])
        wrong_old = (np.sum(4*np.sin(small/2)**2)*np.eye(3)+np.outer(np.sin(small), np.sin(small)))
        wrong_error = np.linalg.norm(bulk_dynamical_matrix(small, par)-wrong_old, 2)
        assert wrong_error > 1000*np.linalg.norm(small)**4/2
        N = 64
        k = 2*np.pi*np.fft.fftfreq(N)
        low_mask = np.abs(k) <= .2
        force = np.zeros(N); force[0] = 1
        spectrum = np.fft.fft(force)
        assert np.any(abs(spectrum[~low_mask]) > .9)
        # The uniform displacement/velocity mode has identical free drift.
        assert np.array_equal(old_symbol(np.zeros(3), 1.), np.zeros((3, 3)))
        assert np.array_equal(bulk_dynamical_matrix(np.zeros(3), par), np.zeros((3, 3)))

    result = {'status': 'PASS', 'scope': 'Homogeneous harmonic bulk only; analytic proof is separate.',
        'dictionary': {'kappa_d': 'K', 'B': '2*K/a', 'm_d': 'm_b', 'v0': 'a**3'},
        'symbol_cases': len(records), 'flow_cases': len(flows), 'tail_cases': len(tails),
        'max_scaled_exact_identity_error': max_identity_error,
        'max_scaled_energy_metric_residual': max_metric_residual,
        'declared_controls': {'wrong_dictionary_rejected': True, 'exact_equality_rejected_at_zone_boundary': True,
            'localized_contact_can_excite_out_of_band': True, 'zero_mode_identical': True},
        'example': {'rho': .2, 'omega_structure_time': 10., 'relative_energy_form_bound': .02,
            'relative_energy_coordinate_flow_bound': float(np.sqrt(1.02/.98)*.004*10)},
        'engine_sha256': hashlib.sha256(payload).hexdigest(),
        'check_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'symbol_records': records, 'flow_records': flows, 'tail_records': tails,
        'limits': ['No exact original-DST-to-junction ancestry', 'No inferred tail for the existing Pulse',
            'No junction efficiency transfer', 'No nonlinear, quantum, observer-clock or fixed-c result',
            'Floating-point tests support formulas; do not certify them statistically']}
    (HERE/'CONTROLLED_BULK_CHECK.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if not k.endswith('_records')}, indent=2))


if __name__ == '__main__':
    main()
