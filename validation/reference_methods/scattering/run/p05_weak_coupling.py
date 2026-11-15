#!/usr/bin/env python3
from __future__ import annotations

import json
import math

import numpy as np
from scipy.optimize import brentq

from p05_common import (
    OMEGA,
    P,
    PROFILE,
    ROOT,
    cross_section_form_factor,
    dressed_objects,
    independent_bulk_matrix,
    independent_dD,
    line_kappa_gamma,
    line_stiffness,
    reduced_static_matrices,
    rel_scalar,
    scattering_summary,
    write_json,
)

CONTRACT = json.loads((ROOT / "P05_TEST_CONTRACT.json").read_text(encoding="utf-8"))
TH = CONTRACT["thresholds"]
OUT = ROOT / "results" / "P05_WEAK_COUPLING.json"
WITNESS_K = np.array([0.6318952450050036, -0.6028746482147381, -0.964106160002423], dtype=float)
BRANCH = 0


def uncoupled_line_caroli_state() -> tuple[np.ndarray, float, complex, float]:
    nq = len(PROFILE) + 1
    kappa, sigma, gamma = line_kappa_gamma(OMEGA)
    KD = line_stiffness(nq)
    SL = np.zeros((nq, nq), dtype=complex)
    SR = np.zeros_like(SL)
    SL[0, 0] = sigma
    SR[-1, -1] = sigma
    A0 = P.m_spine * OMEGA**2 * np.eye(nq) - KD - SL - SR
    eL = np.eye(nq, dtype=complex)[:, 0]
    q_alg = np.linalg.solve(A0, eL)
    chi0 = math.sqrt(gamma) * q_alg
    return chi0, kappa, sigma, gamma


def first_order_source() -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    chi0, _, _, _ = uncoupled_line_caroli_state()
    U1, V1, Kself1 = reduced_static_matrices(PROFILE, P)
    eta1 = V1 @ chi0
    return eta1, U1, V1, float(np.linalg.norm(Kself1))


def gradient_omega2(k: np.ndarray, e: np.ndarray) -> np.ndarray:
    return np.array([float(np.real(np.vdot(e, independent_dD(k, a) @ e))) for a in range(3)])


def production_phi(k: np.ndarray, branch: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eigh(independent_bulk_matrix(k))
    F = cross_section_form_factor(float(k[1]), float(k[2]))
    cm = vecs[:, branch].conj() @ F
    phi = np.concatenate([cm.conj() * np.exp(1j * float(k[0]) * s) for s in range(len(PROFILE))])
    return vals, vecs, phi


def independent_phi(k: np.ndarray, branch: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eigh(independent_bulk_matrix(k))
    # Rebuild the normalized X/R cross-section form factor from the four core corners.
    F = np.zeros((3, 2), dtype=complex)
    corners = [(-1, -1, -1, -1), (-1, 0, -1, 1), (0, -1, 1, -1), (0, 0, 1, 1)]
    for y, z, sy, sz in corners:
        phase = np.exp(-1j * (float(k[1]) * y + float(k[2]) * z))
        F[0, 0] += 0.5 * phase
        F[1, 1] += sy * phase / (2 * math.sqrt(2))
        F[2, 1] += sz * phase / (2 * math.sqrt(2))
    cm = vecs[:, branch].conj() @ F
    phi = np.concatenate([cm.conj() * np.exp(1j * float(k[0]) * s) for s in range(len(PROFILE))])
    return vals, vecs, phi


def witness_record(k: np.ndarray, branch: int, eta1: np.ndarray) -> dict:
    vals, vecs, phi = production_phi(k, branch)
    _, _, phi_i = independent_phi(k, branch)
    overlap = float(abs(np.vdot(phi, eta1)) ** 2)
    overlap_i = float(abs(np.vdot(phi_i, eta1)) ** 2)
    grad = gradient_omega2(k, vecs[:, branch])
    return {
        "branch": branch,
        "k": k.tolist(),
        "omega2_shell_residual": float(vals[branch] - OMEGA**2),
        "gradient_omega2": grad.tolist(),
        "gradient_omega2_norm": float(np.linalg.norm(grad)),
        "internal_lower_gap_omega2": float(vals[1] - vals[0]),
        "upper_cluster_gap_omega2": float(vals[2] - vals[1]),
        "overlap_amplitude_squared_caroli_first_order_source": overlap,
        "independent_overlap_amplitude_squared": overlap_i,
        "overlap_absolute_difference": abs(overlap_i - overlap),
    }


def root_on_direction(direction: np.ndarray, branch: int) -> float:
    d = np.asarray(direction, dtype=float)
    d /= np.linalg.norm(d)
    nz = np.abs(d) > 1e-14
    rmax = float(np.min(np.pi / np.abs(d[nz])))

    def f(r: float) -> float:
        return float(np.linalg.eigvalsh(independent_bulk_matrix(r * d))[branch] - OMEGA**2)

    f0 = f(1e-12)
    f1 = f(rmax * (1 - 1e-12))
    if f0 == 0:
        return 1e-12
    if f0 * f1 < 0:
        return float(brentq(f, 1e-12, rmax * (1 - 1e-12), xtol=2e-14, rtol=2e-14, maxiter=200))
    # Deterministic fallback for an unexpected nonmonotone direction.
    rs = np.linspace(1e-10, rmax * (1 - 1e-10), 1200)
    prev_r = rs[0]
    prev_f = f(prev_r)
    roots = []
    for r in rs[1:]:
        cur_f = f(float(r))
        if prev_f * cur_f < 0:
            roots.append(brentq(f, float(prev_r), float(r), xtol=2e-14, rtol=2e-14, maxiter=200))
        prev_r, prev_f = r, cur_f
    if not roots:
        raise RuntimeError("no shell root found on neighborhood direction")
    target_radius = float(np.linalg.norm(WITNESS_K))
    return float(min(roots, key=lambda x: abs(x - target_radius)))


def witness_neighborhood(eta1: np.ndarray) -> dict:
    d0 = WITNESS_K / np.linalg.norm(WITNESS_K)
    anchor = np.array([1.0, 0.0, 0.0])
    if abs(float(anchor @ d0)) > 0.9:
        anchor = np.array([0.0, 1.0, 0.0])
    t1 = anchor - float(anchor @ d0) * d0
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(d0, t1)
    rows = []
    for da in np.linspace(-0.03, 0.03, 7):
        for db in np.linspace(-0.03, 0.03, 7):
            d = d0 + float(da) * t1 + float(db) * t2
            d /= np.linalg.norm(d)
            r = root_on_direction(d, BRANCH)
            k = r * d
            row = witness_record(k, BRANCH, eta1)
            row["offset_coordinates"] = [float(da), float(db)]
            rows.append(row)
    return {
        "points": len(rows),
        "minimum_overlap": min(r["overlap_amplitude_squared_caroli_first_order_source"] for r in rows),
        "maximum_overlap": max(r["overlap_amplitude_squared_caroli_first_order_source"] for r in rows),
        "minimum_gradient_omega2_norm": min(r["gradient_omega2_norm"] for r in rows),
        "minimum_internal_lower_gap_omega2": min(r["internal_lower_gap_omega2"] for r in rows),
        "minimum_upper_cluster_gap_omega2": min(r["upper_cluster_gap_omega2"] for r in rows),
        "maximum_shell_residual_absolute": max(abs(r["omega2_shell_residual"]) for r in rows),
        "maximum_independent_overlap_difference": max(r["overlap_absolute_difference"] for r in rows),
        "rows": rows,
    }


def weak_sequence(A_qT: float) -> dict:
    s_values = np.array([1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1], dtype=float)
    rows = []
    remainders = []
    singular_values = []
    for s in s_values:
        obj = dressed_objects(s * PROFILE)
        scat = scattering_summary(obj, "left")
        TqT = float(sum(scat["Tbranches"][:2]))
        ratio = TqT / s**2
        relative = abs(ratio - A_qT) / A_qT
        remainder = abs(TqT - s**2 * A_qT)
        remainders.append(remainder)
        singular_values.append(float(np.linalg.svd(np.eye(12) + obj["g0"] @ obj["U"], compute_uv=False)[-1]))
        rows.append(
            {
                "s": float(s),
                "TqT": TqT,
                "TqT_over_s2": ratio,
                "relative_to_AqT": relative,
                "absolute_remainder": remainder,
                "closure_abs": scat["closure_abs"],
                "minimum_dressing_singular_value": singular_values[-1],
            }
        )
    # T=s^2 A + O(s^3), so |T-s^2A| / s^2 should scale as s^1.
    normalized_remainder = np.asarray(remainders[:6]) / s_values[:6] ** 2
    slope = float(np.polyfit(np.log(s_values[:6]), np.log(normalized_remainder), 1)[0])
    return {
        "A_qT": A_qT,
        "rows": rows,
        "maximum_first_four_relative_to_AqT": max(r["relative_to_AqT"] for r in rows[:4]),
        "normalized_remainder_log_slope_first_six": slope,
        "minimum_dressing_singular_value": min(singular_values),
    }


def nodal_negative_control(eta1: np.ndarray) -> dict:
    # The radial X/R source form factor is transverse-nodal at ky=kz=0.
    k = np.array([2 * math.asin(OMEGA / (2 * math.sqrt(P.k_axis / P.m_bulk))), 0.0, 0.0])
    vals, vecs = np.linalg.eigh(independent_bulk_matrix(k))
    overlaps = []
    for branch in (0, 1):
        _, _, phi = production_phi(k, branch)
        overlaps.append(float(abs(np.vdot(phi, eta1)) ** 2))
    return {
        "k": k.tolist(),
        "eigenvalues": vals.tolist(),
        "lower_overlaps": overlaps,
        "maximum_lower_overlap": max(overlaps),
        "note": "This is a transverse form-factor nodal control, not the reference finite-slab witness.",
    }


def main() -> None:
    eta1, _, _, _ = first_order_source()
    raw = json.loads((ROOT / "results" / "P05_CORE_SCATTERING.json").read_text(encoding="utf-8"))
    im_qT = np.load(ROOT / "dependencies" / "phase12_core" / "raw" / "on_shell_N80.npz")
    qT = im_qT["Im_branch_0"][:12, :12] + im_qT["Im_branch_1"][:12, :12]
    A_qT = float(2 * np.real(np.vdot(eta1, qT @ eta1)))
    witness = witness_record(WITNESS_K, BRANCH, eta1)
    neighborhood = witness_neighborhood(eta1)
    sequence = weak_sequence(A_qT)
    nodal = nodal_negative_control(eta1)
    lead_upper = 2 * math.sqrt(P.k_spine / P.m_spine)
    regularity = {
        "scalar_lead_lower_edge_margin": OMEGA,
        "scalar_lead_upper_edge_margin": lead_upper - OMEGA,
        "reference_dressing_minimum_singular_value": float(
            np.linalg.svd(np.eye(12) + dressed_objects(PROFILE)["g0"] @ dressed_objects(PROFILE)["U"], compute_uv=False)[-1]
        ),
        "reference_optical_closure_abs": raw["production_core"]["optical"]["scalar_closure_abs"],
    }

    checks = {
        "witness_on_shell": abs(witness["omega2_shell_residual"]) < TH["weak_shell_residual"],
        "witness_regular_gradient": witness["gradient_omega2_norm"] > TH["weak_gradient_minimum"],
        "witness_spectral_gaps": min(
            witness["internal_lower_gap_omega2"], witness["upper_cluster_gap_omega2"]
        )
        > TH["weak_gap_minimum"],
        "witness_positive_overlap": witness["overlap_amplitude_squared_caroli_first_order_source"] > TH["weak_overlap_minimum"],
        "independent_witness_reimplementation": witness["overlap_absolute_difference"] < TH["matrix_identity_relative"],
        "positive_neighborhood": (
            neighborhood["points"] == 49
            and neighborhood["minimum_overlap"] > TH["weak_neighborhood_minimum_overlap"]
            and neighborhood["minimum_gradient_omega2_norm"] > TH["weak_gradient_minimum"]
            and min(
                neighborhood["minimum_internal_lower_gap_omega2"],
                neighborhood["minimum_upper_cluster_gap_omega2"],
            )
            > TH["weak_gap_minimum"]
            and neighborhood["maximum_shell_residual_absolute"] < TH["weak_shell_residual"]
        ),
        "weak_ratio_convergence": sequence["maximum_first_four_relative_to_AqT"] < TH["weak_ratio_relative_first_four"],
        "weak_remainder_order": TH["weak_remainder_slope_minimum"] < sequence["normalized_remainder_log_slope_first_six"] < TH["weak_remainder_slope_maximum"],
        "regular_away_from_lead_edges_and_poles": min(
            regularity["scalar_lead_lower_edge_margin"],
            regularity["scalar_lead_upper_edge_margin"],
            regularity["reference_dressing_minimum_singular_value"],
            sequence["minimum_dressing_singular_value"],
        )
        > TH["weak_gap_minimum"],
        "nodal_negative_control": nodal["maximum_lower_overlap"] < 1e-24,
    }
    checks = {key: bool(value) for key, value in checks.items()}
    output = {
        "schema": "paper1.p05.weak-coupling.v1",
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "first_order_coefficient": {"A_qT": A_qT},
        "regular_simple_branch_witness": witness,
        "regular_shell_neighborhood": neighborhood,
        "finite_s_sequence": sequence,
        "regularity": regularity,
        "nodal_negative_control": nodal,
    }
    write_json(OUT, output)
    print(json.dumps(output, indent=2, sort_keys=True))
    if output["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
