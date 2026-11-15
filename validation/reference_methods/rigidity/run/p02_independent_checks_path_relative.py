#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'inputs'
OUT = ROOT / 'records/P02_INDEPENDENT_CHECKS_PASS16.json'

EXPECTED_HASHES = {
    'model_api_v1_1.py': 'aa4128940dc7d39682af53c87cf04ae2532de4af880d25393ef4a4fe0b067980',
    'phase12_implementation_audit_and_results.md': '35502fc0eedf64d1948b429fe40efbe2c490802e31c186a5466854ed964f4f3c',
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_model():
    path = SRC / 'model_api_v1_1.py'
    spec = importlib.util.spec_from_file_location('p1model_p02', path)
    if spec is None or spec.loader is None:
        raise RuntimeError('could not load model module')
    mod = importlib.util.module_from_spec(spec)
    sys.modules['p1model_p02'] = mod
    spec.loader.exec_module(mod)
    return mod


def svd_nullity(B: np.ndarray, tol: float = 1e-10) -> tuple[int, np.ndarray]:
    s = np.linalg.svd(B, compute_uv=False)
    rank = int(np.sum(s > tol))
    return int(B.shape[1] - rank), s


AX = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
FD = ((1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1))


def bulk_extension_matrix(nx: int, ny: int, nz: int, *, axial=True, diagonals=True):
    sites = list(product(range(nx), range(ny), range(nz)))
    index = {r: i for i, r in enumerate(sites)}
    rows = []
    labels = []
    directions = (AX if axial else ()) + (FD if diagonals else ())
    for r in sites:
        for d in directions:
            s = (r[0] + d[0], r[1] + d[1], r[2] + d[2])
            if s not in index:
                continue
            n = np.asarray(d, dtype=float)
            n /= np.linalg.norm(n)
            row = np.zeros(3 * len(sites), dtype=float)
            row[3 * index[r]:3 * index[r] + 3] -= n
            row[3 * index[s]:3 * index[s] + 3] += n
            rows.append(row)
            labels.append((r, s, d))
    return np.asarray(rows), sites, labels


def combined_extension_matrix(nx: int, ny: int = 2, nz: int = 2, *, coupling=0.2):
    if ny != 2 or nz != 2:
        raise ValueError('this fixture uses the four-corner transverse source cell')
    js = (-1, 0)
    ks = (-1, 0)
    bulk_sites = [(i, j, k) for i in range(nx) for j in js for k in ks]
    bindex = {r: i for i, r in enumerate(bulk_sites)}
    nline = nx + 1
    ndof = nline + 3 * len(bulk_sites)
    rows = []
    weights = []
    labels = []

    for i in range(nx):
        row = np.zeros(ndof)
        row[i] -= 1.0
        row[i + 1] += 1.0
        rows.append(row)
        weights.append(1.69)
        labels.append(('line', i))

    for r in bulk_sites:
        for d in AX + FD:
            s = (r[0] + d[0], r[1] + d[1], r[2] + d[2])
            if s not in bindex:
                continue
            n = np.asarray(d, dtype=float)
            n /= np.linalg.norm(n)
            row = np.zeros(ndof)
            off = nline
            row[off + 3 * bindex[r]:off + 3 * bindex[r] + 3] -= n
            row[off + 3 * bindex[s]:off + 3 * bindex[s] + 3] += n
            rows.append(row)
            weights.append(1.0)
            labels.append(('bulk', r, s, d))

    if coupling > 0:
        for i in range(nx):
            for sy, j in ((-1, -1), (1, 0)):
                for sz, k in ((-1, -1), (1, 0)):
                    r = (i, j, k)
                    for endpoint, nvec, side in (
                        (i, np.asarray((1, sy, sz), float) / math.sqrt(3), 'L'),
                        (i + 1, np.asarray((-1, sy, sz), float) / math.sqrt(3), 'R'),
                    ):
                        row = np.zeros(ndof)
                        row[endpoint] -= nvec[0]
                        off = nline + 3 * bindex[r]
                        row[off:off + 3] += nvec
                        rows.append(row)
                        weights.append(coupling)
                        labels.append(('source', i, sy, sz, side))

    return np.asarray(rows), np.asarray(weights), bulk_sites, nline, labels


def independent_bulk_D(k: np.ndarray, k_axis: float, k_diag: float, mass: float) -> np.ndarray:
    D = np.zeros((3, 3), float)
    for d, stiffness in [(d, k_axis) for d in AX] + [(d, k_diag) for d in FD]:
        v = np.asarray(d, float)
        n = v / np.linalg.norm(v)
        D += 2.0 * stiffness * (1.0 - math.cos(float(k @ v))) * np.outer(n, n)
    return D / mass


def main() -> None:
    checks = []

    def add(name: str, passed: bool, actual, expected=None, notes=''):
        checks.append({
            'name': name,
            'passed': bool(passed),
            'actual': actual,
            'expected': expected,
            'notes': notes,
        })

    # Exact-source bindings.
    for name, expected in EXPECTED_HASHES.items():
        actual = sha256(SRC / name)
        add(f'hash::{name}', actual == expected, actual, expected)

    model = load_model()
    p = model.ModelParameters()

    # Bond inventory.
    inventory = list(model.AXIAL_BONDS) + list(model.FACE_DIAGONAL_BONDS)
    canon = {tuple(int(x) for x in d) for d in inventory}
    signed_collision = any(tuple(-x for x in d) in canon for d in canon)
    add('bond_inventory_counts', len(model.AXIAL_BONDS) == 3 and len(model.FACE_DIAGONAL_BONDS) == 6,
        {'axial': len(model.AXIAL_BONDS), 'face_diagonal': len(model.FACE_DIAGONAL_BONDS)},
        {'axial': 3, 'face_diagonal': 6})
    add('bond_inventory_unique_positive_representatives', len(canon) == 9 and not signed_collision,
        {'unique': len(canon), 'signed_collision': signed_collision},
        {'unique': 9, 'signed_collision': False})

    # Symbolic source-link reduction.
    qi, qj, ux, uy, uz, g = sp.symbols('qi qj ux uy uz g', real=True)
    sy, sz = sp.symbols('sy sz', real=True)
    qbar = (qi + qj) / 2
    dq = qj - qi
    u = sp.Matrix((ux, uy, uz))
    ex = sp.Matrix((1, 0, 0))
    nL = sp.Matrix((1, sy, sz)) / sp.sqrt(3)
    nR = sp.Matrix((-1, sy, sz)) / sp.sqrt(3)
    eL = nL.dot(u - qi * ex)
    eR = nR.dot(u - qj * ex)
    H_links = sp.expand(g * (eL**2 + eR**2) / 2)
    H_reduced = sp.expand(g * ((ux - qbar)**2 + (sy * uy + sz * uz + dq / 2)**2) / 3)
    identity = sp.simplify(H_links - H_reduced)
    add('symbolic_source_link_rotation_identity', identity == 0, str(identity), '0')

    variables = (qi, qj, ux, uy, uz)
    Hlocal = sp.hessian(H_reduced, variables).subs({sy: 1, sz: 1, g: 1})
    local_eigs = sorted([float(v) for v in np.linalg.eigvalsh(np.asarray(Hlocal, float))])
    add('local_source_hessian_psd_rank2', min(local_eigs) > -1e-12 and sum(v > 1e-12 for v in local_eigs) == 2,
        local_eigs, 'three zeros and two positive eigenvalues')

    # Independent random force-gradient fixture.
    rng = np.random.default_rng(20260813)
    max_scaled_err = 0.0
    Hf = sp.lambdify((qi, qj, ux, uy, uz, g, sy, sz), H_reduced, 'numpy')
    grad = [sp.lambdify((qi, qj, ux, uy, uz, g, sy, sz), -sp.diff(H_reduced, v), 'numpy') for v in variables]
    for _ in range(100):
        x = rng.normal(size=5)
        gv = float(rng.uniform(0.05, 1.5))
        syv = int(rng.choice((-1, 1)))
        szv = int(rng.choice((-1, 1)))
        analytic = np.asarray([f(*x, gv, syv, szv) for f in grad], float)
        eps = 2e-6
        numeric = np.zeros(5)
        for j in range(5):
            xp = x.copy(); xm = x.copy()
            xp[j] += eps; xm[j] -= eps
            numeric[j] = -(Hf(*xp, gv, syv, szv) - Hf(*xm, gv, syv, szv)) / (2 * eps)
        scale = max(1.0, np.linalg.norm(analytic), np.linalg.norm(numeric))
        max_scaled_err = max(max_scaled_err, float(np.linalg.norm(analytic - numeric) / scale))
    add('independent_random_force_gradient', max_scaled_err < 2e-9, max_scaled_err, '<2e-9')

    # Production D versus independent bond-sum D.
    max_D_err = 0.0
    min_random_eig = float('inf')
    min_ratio = float('inf')
    for _ in range(5000):
        k = rng.uniform(-math.pi, math.pi, size=3)
        if np.linalg.norm(k) < 1e-8:
            continue
        D1 = model.bulk_dynamical_matrix(k, p)
        D2 = independent_bulk_D(k, p.k_axis, p.k_diag, p.m_bulk)
        max_D_err = max(max_D_err, float(np.max(np.abs(D1 - D2))))
        eig = np.linalg.eigvalsh(D1)
        min_random_eig = min(min_random_eig, float(eig[0]))
        # Periodic distance from Gamma, for a scale-aware positivity diagnostic.
        km = ((k + math.pi) % (2 * math.pi)) - math.pi
        min_ratio = min(min_ratio, float(eig[0] / max(np.dot(km, km), 1e-30)))
    add('production_bulk_D_matches_independent_bond_sum', max_D_err < 2e-14, max_D_err, '<2e-14')
    add('bulk_D_random_nonzero_k_positive', min_random_eig > 0, min_random_eig, '>0')
    add('bulk_D_random_scale_aware_floor_positive', min_ratio > 0, min_ratio, '>0')

    # Exact high-symmetry grid check, excluding Gamma.
    grid = np.linspace(-math.pi, math.pi, 17)
    min_grid_eig = float('inf')
    argmin = None
    for kx, ky, kz in product(grid, repeat=3):
        k = np.asarray((kx, ky, kz))
        if np.linalg.norm(k) < 1e-12:
            continue
        eig0 = float(np.linalg.eigvalsh(model.bulk_dynamical_matrix(k, p))[0])
        if eig0 < min_grid_eig:
            min_grid_eig = eig0; argmin = [float(kx), float(ky), float(kz)]
    add('bulk_D_grid_nonzero_k_positive', min_grid_eig > 1e-12,
        {'minimum_eigenvalue': min_grid_eig, 'k': argmin}, '>1e-12')

    # Mutations must fail in the predicted way.
    k_axis_test = np.asarray((1.1, 0.0, 0.0))
    D_axis_only = independent_bulk_D(k_axis_test, 1.0, 0.0, 1.0)
    axis_eigs = np.linalg.eigvalsh(D_axis_only)
    add('mutation_axis_only_has_two_transverse_floppy_modes', np.sum(np.abs(axis_eigs) < 1e-12) == 2,
        axis_eigs.tolist(), 'two zeros')
    D_negative = independent_bulk_D(np.asarray((0.8, 0.7, 0.3)), 1.0, -1.0, 1.0)
    neg_eigs = np.linalg.eigvalsh(D_negative)
    add('mutation_negative_diagonal_stiffness_breaks_nonnegativity', float(neg_eigs[0]) < -1e-6,
        neg_eigs.tolist(), 'at least one negative eigenvalue')

    # Fourier source formula check.
    max_source_rel = 0.0
    for _ in range(500):
        k = rng.uniform(-math.pi, math.pi, size=3)
        closed = model.coupling_vector(k, p).astype(complex)
        direct = model.fourier_transform_source_pattern(k, p)
        if np.linalg.norm(closed) < 1e-14 and np.linalg.norm(direct) < 1e-14:
            continue
        inner = np.vdot(direct, closed)
        phase = inner / abs(inner) if abs(inner) > 1e-14 else 1.0 + 0j
        rel = np.linalg.norm(direct * phase - closed) / max(np.linalg.norm(closed), 1e-14)
        max_source_rel = max(max_source_rel, float(rel))
    add('closed_form_source_matches_direct_DFT_up_to_global_phase', max_source_rel < 2e-14,
        max_source_rel, '<2e-14')

    # Incidence/Gram formulation and finite-cluster nullspaces.
    bulk_cluster_results = []
    for dims in ((2, 2, 2), (3, 2, 2), (3, 3, 3), (4, 3, 3)):
        B, sites, labels = bulk_extension_matrix(*dims)
        nullity, singular = svd_nullity(B)
        bulk_cluster_results.append({
            'shape': list(dims), 'dof': int(B.shape[1]), 'bonds': int(B.shape[0]),
            'nullity': nullity, 'smallest_nonzero_singular': float(singular[-7]) if len(singular) >= 7 else None,
        })
    add('representative_free_bulk_clusters_have_only_six_rigid_modes',
        all(r['nullity'] == 6 for r in bulk_cluster_results), bulk_cluster_results, 'nullity=6 each')

    combined_results = []
    for nx in (1, 2, 3, 4):
        B, W, sites, nline, labels = combined_extension_matrix(nx)
        Bw = np.sqrt(W)[:, None] * B
        nullity, singular = svd_nullity(Bw)
        combined_results.append({
            'active_segments': nx, 'dof': int(B.shape[1]), 'constraints': int(B.shape[0]),
            'nullity': nullity, 'smallest_nonzero_singular': float(singular[-3]) if len(singular) >= 3 else None,
        })
    add('representative_free_combined_clusters_have_two_exact_null_modes',
        all(r['nullity'] == 2 for r in combined_results), combined_results, 'nullity=2 each')

    # Direct residuals of the two predicted modes for the nx=3 fixture.
    B, W, bulk_sites, nline, labels = combined_extension_matrix(3)
    ndof = B.shape[1]
    axial = np.zeros(ndof); axial[:nline] = 1.0
    for bi in range(len(bulk_sites)):
        axial[nline + 3 * bi] = 1.0
    rotation = np.zeros(ndof)
    for bi, (i, j, k) in enumerate(bulk_sites):
        R = np.asarray((i + 0.5, j + 0.5, k + 0.5), float)
        disp = np.cross(np.asarray((1.0, 0.0, 0.0)), R)
        rotation[nline + 3 * bi:nline + 3 * bi + 3] = disp
    add('combined_axial_translation_residual', float(np.max(np.abs(B @ axial))) < 1e-14,
        float(np.max(np.abs(B @ axial))), '<1e-14')
    add('combined_bulk_rotation_about_line_residual', float(np.max(np.abs(B @ rotation))) < 1e-14,
        float(np.max(np.abs(B @ rotation))), '<1e-14')

    # Exact source-motif constraint theorem under a rigid bulk and connected line.
    tx, ty, tz, ox, oy, oz, c, a = sp.symbols('tx ty tz ox oy oz c a', real=True)
    t = sp.Matrix((tx, ty, tz)); Om = sp.Matrix((ox, oy, oz))
    equations = []
    xm = sp.Rational(1, 2) * a
    for syv in (-1, 1):
        for szv in (-1, 1):
            r = sp.Matrix((xm, syv * a / 2, szv * a / 2))
            ur = t + Om.cross(r)
            equations.append(sp.Matrix((1, syv, szv)).dot(ur - c * sp.Matrix((1, 0, 0))))
            equations.append(sp.Matrix((-1, syv, szv)).dot(ur - c * sp.Matrix((1, 0, 0))))
    A, _ = sp.linear_eq_to_matrix(equations, (tx, ty, tz, ox, oy, oz, c))
    rank = int(A.rank())
    null_basis = A.nullspace()
    basis_strings = [[str(x) for x in v] for v in null_basis]
    expected_basis = rank == 5 and len(null_basis) == 2
    add('symbolic_complete_source_motif_leaves_two_rigid_modes', expected_basis,
        {'constraint_rank': rank, 'nullity': len(null_basis), 'basis': basis_strings},
        'rank=5, nullity=2 (axial translation and x-axis rotation)')

    # Decoupled mutation: one line translation plus six rigid bulk motions.
    B0, W0, *_ = combined_extension_matrix(2, coupling=0.0)
    null0, s0 = svd_nullity(np.sqrt(W0)[:, None] * B0)
    add('mutation_zero_source_coupling_restores_disconnected_nullity', null0 == 7,
        null0, 7)

    # Exact finite-dimensional Hamiltonian flow checks on a representative fixture.
    Bd, Wd, bulk_sites_d, nline_d, labels_d = combined_extension_matrix(2)
    Kd = Bd.T @ (Wd[:, None] * Bd)
    masses = np.concatenate((np.ones(nline_d), np.ones(3 * len(bulk_sites_d))))
    Minv = np.diag(1.0 / masses)
    nd = Kd.shape[0]
    Aflow = np.block([[np.zeros((nd, nd)), Minv], [-Kd, np.zeros((nd, nd))]])
    Rtime = np.block([[np.eye(nd), np.zeros((nd, nd))], [np.zeros((nd, nd)), -np.eye(nd)]])
    tflow = 2.75
    Phi = expm(Aflow * tflow)
    rev_resid = float(np.max(np.abs(Rtime @ Phi @ Rtime - expm(-Aflow * tflow))))
    add('finite_flow_time_reversal_identity', rev_resid < 5e-13, rev_resid, '<5e-13')

    z0 = rng.normal(size=2 * nd)
    zt = Phi @ z0
    def energy(z):
        q = z[:nd]; mom = z[nd:]
        return 0.5 * float(mom @ Minv @ mom + q @ Kd @ q)
    E0 = energy(z0); Et = energy(zt)
    energy_rel = abs(Et - E0) / max(abs(E0), 1.0)
    add('finite_flow_energy_conservation_exact_exponential', energy_rel < 5e-13,
        energy_rel, '<5e-13')

    omega_test = 0.73 + 0.11j
    Gdyn = np.linalg.inv(Kd - (omega_test ** 2) * np.diag(masses))
    reciprocity_resid = float(np.max(np.abs(Gdyn - Gdyn.T)))
    add('finite_resolvent_reciprocity_from_symmetric_M_K', reciprocity_resid < 2e-13,
        reciprocity_resid, '<2e-13')

    # Nondimensional scaling covariance.
    ratios = {'m_spine/m_bulk': 1.37, 'k_axis/k_diag': 0.82, 'k_spine/k_diag': 1.91, 'k_couple/k_diag': 0.16}
    p_dimless = model.ModelParameters(a=1.0, m_bulk=1.0, m_spine=ratios['m_spine/m_bulk'],
                                      k_axis=ratios['k_axis/k_diag'], k_diag=1.0,
                                      k_spine=ratios['k_spine/k_diag'], k_couple=ratios['k_couple/k_diag'])
    a_phys, mb_phys, kd_phys = 2.3, 4.7, 7.9
    p_phys = model.ModelParameters(a=a_phys, m_bulk=mb_phys,
                                   m_spine=ratios['m_spine/m_bulk'] * mb_phys,
                                   k_axis=ratios['k_axis/k_diag'] * kd_phys,
                                   k_diag=kd_phys,
                                   k_spine=ratios['k_spine/k_diag'] * kd_phys,
                                   k_couple=ratios['k_couple/k_diag'] * kd_phys)
    max_scale_err = 0.0
    for _ in range(100):
        k = rng.uniform(-math.pi, math.pi, 3)
        D0 = model.bulk_dynamical_matrix(k, p_dimless)
        Dp = model.bulk_dynamical_matrix(k, p_phys) / (kd_phys / mb_phys)
        max_scale_err = max(max_scale_err, float(np.max(np.abs(D0 - Dp))))
        w0 = model.spine_frequency(float(k[0]), p_dimless)
        wp = model.spine_frequency(float(k[0]), p_phys) / math.sqrt(kd_phys / mb_phys)
        max_scale_err = max(max_scale_err, abs(w0 - wp))
    add('nondimensional_scaling_covariance', max_scale_err < 2e-14, max_scale_err, '<2e-14')

    # Every active check must pass.
    all_pass = all(c['passed'] for c in checks)
    result = {
        'schema': 'paper1-p02-independent-model-checks-v1',
        'all_pass': all_pass,
        'check_count': len(checks),
        'checks': checks,
        'interpretation': {
            'model_status': 'No contradiction found in Section-II Hamiltonian, bond inventory, positivity, equilibrium, or units.',
            'nullspace_scope': 'Representative fixtures reproduce the quoted nullities, but the manuscript does not identify the exact finite clusters and producer hashes; that is a presentation/provenance issue, not a demonstrated physics error.',
            'new_exact_deduction': 'If the finite bulk framework is infinitesimally rigid, the scalar line is connected, and one complete positive four-corner source motif is present, the combined stiffness kernel is exactly the common axial translation plus rigid bulk rotation about the line.',
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'all_pass': all_pass, 'check_count': len(checks), 'output': str(OUT)}, indent=2))
    if not all_pass:
        for c in checks:
            if not c['passed']:
                print('FAILED', c['name'], c['actual'], c['expected'])
        raise SystemExit(1)


if __name__ == '__main__':
    main()
