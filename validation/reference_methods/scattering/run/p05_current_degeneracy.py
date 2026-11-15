#!/usr/bin/env python3
from __future__ import annotations

import json
import math

import numpy as np

from p05_common import (
    AXIAL_BONDS,
    FACE_DIAGONAL_BONDS,
    OMEGA,
    P,
    ROOT,
    _unit,
    bulk_dynamical_matrix,
    cross_section_form_factor,
    independent_bulk_matrix,
    independent_dD,
    rel,
    write_json,
)

CONTRACT = json.loads((ROOT / "P05_TEST_CONTRACT.json").read_text(encoding="utf-8"))
TH = CONTRACT["thresholds"]
OUT = ROOT / "results" / "P05_CURRENT_AND_DEGENERACY.json"
BONDS = [(np.asarray(d, float), P.k_axis) for d in AXIAL_BONDS] + [
    (np.asarray(d, float), P.k_diag) for d in FACE_DIAGONAL_BONDS
]


def plane_wave_energy(omega: float, e: np.ndarray) -> float:
    return 0.5 * P.m_bulk * omega**2 * float(np.vdot(e, e).real)


def hf_velocity(k: np.ndarray, e: np.ndarray, omega: float) -> np.ndarray:
    return np.array([float(np.real(np.vdot(e, independent_dD(k, a) @ e))) for a in range(3)]) / (2 * omega)


def microscopic_current(k: np.ndarray, e: np.ndarray, omega: float) -> np.ndarray:
    v0 = -1j * omega * e
    J = np.zeros(3, dtype=float)
    for d, stiffness in BONDS:
        n = _unit(d)
        phase = np.exp(1j * float(k @ d))
        uj = e * phase
        vj = v0 * phase
        force_on_origin = stiffness * n * np.dot(n, uj - e)
        bond_power = -0.25 * np.real(np.vdot(v0 + vj, force_on_origin))
        J += d * bond_power
    return J


def finite_difference_velocity(k: np.ndarray, branch: int, h: float = 2e-6) -> np.ndarray:
    out = []
    for axis in range(3):
        kp = k.copy()
        km = k.copy()
        kp[axis] += h
        km[axis] -= h
        wp = math.sqrt(max(float(np.linalg.eigvalsh(independent_bulk_matrix(kp))[branch]), 0.0))
        wm = math.sqrt(max(float(np.linalg.eigvalsh(independent_bulk_matrix(km))[branch]), 0.0))
        out.append((wp - wm) / (2 * h))
    return np.asarray(out)


def random_current_checks(seed: int = 20260814, target_modes: int = 900) -> dict:
    rng = np.random.default_rng(seed)
    D_errors = []
    current_errors = []
    finite_difference_errors = []
    count = 0
    trials = 0
    while count < target_modes and trials < 100000:
        trials += 1
        k = rng.uniform(-1.35, 1.35, size=3)
        D_errors.append(rel(independent_bulk_matrix(k), bulk_dynamical_matrix(k, P)))
        vals, vecs = np.linalg.eigh(independent_bulk_matrix(k))
        if vals[0] < 1e-7 or np.min(np.diff(vals)) < 2e-5:
            continue
        for branch, w2 in enumerate(vals):
            omega = math.sqrt(float(w2))
            e = vecs[:, branch]
            vg = hf_velocity(k, e, omega)
            J = microscopic_current(k, e, omega)
            target = plane_wave_energy(omega, e) * vg
            current_errors.append(rel(J, target))
            if count % 9 == 0:
                vg_fd = finite_difference_velocity(k, branch)
                finite_difference_errors.append(rel(vg_fd, vg))
            count += 1
            if count >= target_modes:
                break
    if count < target_modes:
        raise RuntimeError("insufficient regular modes")
    return {
        "seed": seed,
        "trials": trials,
        "modes": count,
        "max_D_relative": max(D_errors),
        "max_current_relative": max(current_errors),
        "median_current_relative": float(np.median(current_errors)),
        "p99_current_relative": float(np.quantile(current_errors, 0.99)),
        "finite_difference_samples": len(finite_difference_errors),
        "max_finite_difference_velocity_relative": max(finite_difference_errors),
    }


def core_mode_vector(k: np.ndarray, e: np.ndarray, nseg: int = 6) -> np.ndarray:
    F = cross_section_form_factor(float(k[1]), float(k[2]))
    local = F.conj().T @ e
    return np.concatenate([np.exp(1j * float(k[0]) * s) * local for s in range(nseg)])


def unitary_2(rng: np.random.Generator) -> np.ndarray:
    X = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    Q, R = np.linalg.qr(X)
    phases = np.diag(R)
    phases = np.where(np.abs(phases) > 0, phases / np.abs(phases), 1.0)
    return Q @ np.diag(phases.conj())


def degeneracy_checks(seed: int = 314159, rotations: int = 300) -> dict:
    rng = np.random.default_rng(seed)
    # Along x the lower y/z polarizations form an exact doublet.
    k = np.array([0.73, 0.0, 0.0])
    vals, vecs = np.linalg.eigh(independent_bulk_matrix(k))
    lower_split = float(vals[1] - vals[0])
    E = vecs[:, :2]
    Phi = np.column_stack([core_mode_vector(k, E[:, j]) for j in range(2)])
    measure = Phi @ Phi.conj().T
    projector = E @ E.conj().T
    max_projector = 0.0
    max_measure = 0.0
    max_quadratic = 0.0
    source = rng.normal(size=12) + 1j * rng.normal(size=12)
    base_q = float(np.real(np.vdot(source, measure @ source)))
    for _ in range(rotations):
        U = unitary_2(rng)
        Er = E @ U
        Phir = Phi @ U
        max_projector = max(max_projector, rel(Er @ Er.conj().T, projector))
        Mr = Phir @ Phir.conj().T
        max_measure = max(max_measure, rel(Mr, measure))
        qr = float(np.real(np.vdot(source, Mr @ source)))
        max_quadratic = max(max_quadratic, abs(qr - base_q) / max(abs(base_q), 1e-30))

    # Nondegenerate negative control: unequal weights are not invariant if the basis is rotated without the operator.
    k2 = np.array([0.82, 0.37, -0.24])
    vals2, vecs2 = np.linalg.eigh(independent_bulk_matrix(k2))
    E2 = vecs2[:, :2]
    Phi2 = np.column_stack([core_mode_vector(k2, E2[:, j]) for j in range(2)])
    weights = np.diag([0.4, 1.7])
    base_weighted = Phi2 @ weights @ Phi2.conj().T
    wrong_changes = []
    covariant_errors = []
    for _ in range(rotations):
        U = unitary_2(rng)
        rotated_phi = Phi2 @ U
        wrong = rotated_phi @ weights @ rotated_phi.conj().T
        covariant_weights = U.conj().T @ weights @ U
        covariant = rotated_phi @ covariant_weights @ rotated_phi.conj().T
        wrong_changes.append(rel(wrong, base_weighted))
        covariant_errors.append(rel(covariant, base_weighted))

    # The mutation gate must use a prospectively fixed, deliberately nontrivial
    # basis mixing.  A minimum over random Haar draws is not a valid negative
    # control because the random set can include a matrix arbitrarily close to
    # the identity.  Retain the random distribution only as a diagnostic.
    theta = math.pi / 4
    U_fixed = np.array([[math.cos(theta), -math.sin(theta)],
                        [math.sin(theta),  math.cos(theta)]], dtype=complex)
    fixed_phi = Phi2 @ U_fixed
    fixed_wrong = fixed_phi @ weights @ fixed_phi.conj().T
    fixed_covariant_weights = U_fixed.conj().T @ weights @ U_fixed
    fixed_covariant = fixed_phi @ fixed_covariant_weights @ fixed_phi.conj().T

    return {
        "seed": seed,
        "rotations": rotations,
        "exact_doublet": {
            "k": k.tolist(),
            "eigenvalues": vals.tolist(),
            "lower_split": lower_split,
            "max_projector_relative": max_projector,
            "max_measure_relative": max_measure,
            "max_quadratic_relative": max_quadratic,
        },
        "nondegenerate_weighted_control": {
            "k": k2.tolist(),
            "eigenvalues": vals2.tolist(),
            "fixed_45deg_wrong_weight_change": rel(fixed_wrong, base_weighted),
            "fixed_45deg_covariant_operator_relative": rel(fixed_covariant, base_weighted),
            "minimum_random_wrong_fixed_weight_change_diagnostic": min(wrong_changes),
            "median_random_wrong_fixed_weight_change_diagnostic": float(np.median(wrong_changes)),
            "max_random_covariant_operator_relative": max(covariant_errors),
        },
    }


def branch_state(k: np.ndarray, branch: int) -> tuple[float, np.ndarray]:
    vals, vecs = np.linalg.eigh(independent_bulk_matrix(k))
    return math.sqrt(float(vals[branch])), vecs[:, branch]


def bond_current_at_origin(k: np.ndarray, e: np.ndarray, omega: float, d: np.ndarray, stiffness: float) -> float:
    n = _unit(d)
    v0 = -1j * omega * e
    phase = np.exp(1j * float(k @ d))
    uj = e * phase
    vj = v0 * phase
    force = stiffness * n * np.dot(n, uj - e)
    return float(-0.25 * np.real(np.vdot(v0 + vj, force)))


def coherent_cross_term_control() -> dict:
    # Two independently regular omega=1.2 lower-shell points from the P04 audit.
    k0 = np.array([1.022121111600022, 0.7119903394245296, -0.020431525569168604])
    k1 = np.array([0.8807183236959493, 0.21703145119419856, -0.8514466270114215])
    w0, e0 = branch_state(k0, 0)
    w1, e1 = branch_state(k1, 1)
    if abs(w0 - OMEGA) > 1e-10 or abs(w1 - OMEGA) > 1e-10:
        raise RuntimeError("stored shell witnesses are not on omega=1.2")
    phase = np.exp(1j * 0.371)
    # Sum the complex fields at each endpoint of every representative bond.
    cross_terms = []
    total_currents = []
    modal_sums = []
    for d, stiffness in BONDS:
        n = _unit(d)
        u0_i = e0
        u0_j = e0 * np.exp(1j * float(k0 @ d))
        u1_i = phase * e1
        u1_j = phase * e1 * np.exp(1j * float(k1 @ d))
        vi = -1j * OMEGA * (u0_i + u1_i)
        vj = -1j * OMEGA * (u0_j + u1_j)
        force = stiffness * n * np.dot(n, (u0_j + u1_j) - (u0_i + u1_i))
        total = float(-0.25 * np.real(np.vdot(vi + vj, force)))
        modal = bond_current_at_origin(k0, e0, OMEGA, d, stiffness) + bond_current_at_origin(
            k1, phase * e1, OMEGA, d, stiffness
        )
        total_currents.append(total)
        modal_sums.append(modal)
        cross_terms.append(total - modal)
    return {
        "k0": k0.tolist(),
        "k1": k1.tolist(),
        "relative_cross_term_norm": float(
            np.linalg.norm(cross_terms) / max(np.linalg.norm(total_currents), np.linalg.norm(modal_sums), 1e-30)
        ),
        "maximum_absolute_cross_term": float(np.max(np.abs(cross_terms))),
        "total_representative_bond_currents": total_currents,
        "incoherent_modal_sums": modal_sums,
        "cross_terms": cross_terms,
    }


def main() -> None:
    random_checks = random_current_checks()
    degeneracy = degeneracy_checks()
    coherent = coherent_cross_term_control()
    checks = {
        "independent_D": random_checks["max_D_relative"] < TH["matrix_identity_relative"],
        "microscopic_current_equals_E_vg": random_checks["max_current_relative"] < TH["current_identity_relative"],
        "finite_difference_velocity_control": random_checks["max_finite_difference_velocity_relative"] < TH["finite_difference_velocity_relative"],
        "exact_degenerate_cluster_invariance": max(
            degeneracy["exact_doublet"]["max_projector_relative"],
            degeneracy["exact_doublet"]["max_measure_relative"],
            degeneracy["exact_doublet"]["max_quadratic_relative"],
        )
        < TH["degenerate_cluster_invariance_relative"],
        "nondegenerate_fixed_weight_negative_control": degeneracy["nondegenerate_weighted_control"][
            "fixed_45deg_wrong_weight_change"
        ]
        > TH["wrong_nondegenerate_rotation_minimum_change"],
        "covariant_weighted_operator": max(
            degeneracy["nondegenerate_weighted_control"]["fixed_45deg_covariant_operator_relative"],
            degeneracy["nondegenerate_weighted_control"]["max_random_covariant_operator_relative"],
        )
        < TH["degenerate_cluster_invariance_relative"],
        "local_coherent_cross_term_resolved": coherent["relative_cross_term_norm"] > TH["mutation_minimum_effect"],
    }
    output = {
        "schema": "paper1.p05.current-degeneracy.v1",
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "random_current_checks": random_checks,
        "degeneracy_checks": degeneracy,
        "coherent_local_current_control": coherent,
    }
    write_json(OUT, output)
    print(json.dumps(output, indent=2, sort_keys=True))
    if output["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
