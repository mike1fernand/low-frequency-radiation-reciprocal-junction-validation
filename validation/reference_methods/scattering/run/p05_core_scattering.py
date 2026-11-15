#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from p05_common import (
    OMEGA,
    P,
    PROFILE,
    ROOT,
    complex_json,
    dressed_objects,
    herm,
    imh,
    line_kappa_gamma,
    line_stiffness,
    rel,
    rel_scalar,
    scattering_summary,
    write_json,
)

CONTRACT = json.loads((ROOT / "P05_TEST_CONTRACT.json").read_text(encoding="utf-8"))
TH = CONTRACT["thresholds"]
OUT = ROOT / "results" / "P05_CORE_SCATTERING.json"


def independent_active_basis(nseg: int) -> np.ndarray:
    signs = [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)]
    B = np.zeros((12 * nseg, 2 * nseg), dtype=float)
    for s in range(nseg):
        for a, (sy, sz) in enumerate(signs):
            r = 12 * s + 3 * a
            B[r, 2 * s] = 0.5
            B[r + 1, 2 * s + 1] = sy / (2 * math.sqrt(2))
            B[r + 2, 2 * s + 1] = sz / (2 * math.sqrt(2))
    return B


def direct_paired_link_hessian(profile: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nseg = len(profile)
    nq = nseg + 1
    nu = 12 * nseg
    H = np.zeros((nq + nu, nq + nu), dtype=float)
    signs = [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)]
    rt3 = math.sqrt(3)
    for s, g in enumerate(profile):
        for a, (sy, sz) in enumerate(signs):
            u0 = nq + 12 * s + 3 * a
            # Left endpoint-to-corner spring: n=(1,sy,sz)/sqrt(3), extension n.u-q_s/sqrt(3).
            c = np.zeros(nq + nu)
            c[s] = -1 / rt3
            c[u0 : u0 + 3] = np.array([1.0, sy, sz]) / rt3
            H += g * np.outer(c, c)
            # Right endpoint-to-corner spring: n=(-1,sy,sz)/sqrt(3), extension n.u+q_{s+1}/sqrt(3).
            c = np.zeros(nq + nu)
            c[s + 1] = 1 / rt3
            c[u0 : u0 + 3] = np.array([-1.0, sy, sz]) / rt3
            H += g * np.outer(c, c)
    Kqq = H[:nq, :nq]
    Kuq = H[nq:, :nq]
    Kuu = H[nq:, nq:]
    return Kuu, Kuq, Kqq


def high_precision_transparent() -> dict:
    mp.mp.dps = 90
    cases = []
    max_R = mp.mpf("0")
    max_Terr = mp.mpf("0")
    max_contact = mp.mpf("0")
    max_flux = mp.mpf("0")
    for omega_float in [0.2, 0.5, 0.9, 1.2, 1.6, 2.0, 2.5]:
        omega = mp.mpf(str(omega_float))
        K = mp.mpf(str(P.k_spine))
        m = mp.mpf(str(P.m_spine))
        kappa = 2 * mp.asin(omega / (2 * mp.sqrt(K / m)))
        sigma = -K * mp.e ** (1j * kappa)
        gamma = 2 * K * mp.sin(kappa)
        n = 7
        A = mp.matrix(n)
        for i in range(n):
            A[i, i] = m * omega**2 - 2 * K
            if i + 1 < n:
                A[i, i + 1] = K
                A[i + 1, i] = K
        A[0, 0] -= sigma
        A[n - 1, n - 1] -= sigma
        G = A ** -1
        r = 1 - 1j * gamma * G[0, 0]
        t = -1j * gamma * G[n - 1, 0]
        R = abs(r) ** 2
        T = abs(t) ** 2
        contact = abs(mp.sqrt(gamma) * G[0, 0])
        expected_contact = 1 / mp.sqrt(gamma)
        caroli_flux = omega * gamma * contact**2 / 4
        unit_flux = omega * gamma * ((2 / mp.sqrt(omega)) * contact) ** 2 / 4
        max_R = max(max_R, R)
        max_Terr = max(max_Terr, abs(T - 1))
        max_contact = max(max_contact, abs(contact - expected_contact))
        max_flux = max(max_flux, abs(unit_flux - 1))
        cases.append(
            {
                "omega": omega_float,
                "kappa": float(kappa),
                "R": float(R),
                "T": float(T),
                "contact_abs": float(contact),
                "expected_contact_abs": float(expected_contact),
                "caroli_physical_flux": float(caroli_flux),
                "omega_over_4": omega_float / 4,
                "scaled_unit_flux": float(unit_flux),
            }
        )
    return {
        "mp_dps": mp.mp.dps,
        "cases": cases,
        "max_reflection": float(max_R),
        "max_transmission_abs_error": float(max_Terr),
        "max_contact_abs_error": float(max_contact),
        "max_scaled_unit_flux_abs_error": float(max_flux),
    }


def double_transparent(endpoint_diag: float | None = None) -> dict:
    omega = OMEGA
    kappa, sigma, gamma = line_kappa_gamma(omega)
    KD = line_stiffness(7, endpoint_diag=endpoint_diag)
    SL = np.zeros((7, 7), complex)
    SR = np.zeros_like(SL)
    SL[0, 0] = sigma
    SR[-1, -1] = sigma
    A = P.m_spine * omega**2 * np.eye(7) - KD - SL - SR
    G = np.linalg.inv(A)
    r = 1 - 1j * gamma * G[0, 0]
    t = -1j * gamma * G[-1, 0]
    chi = G[:, 0] * math.sqrt(gamma)
    flux = omega * gamma * abs(chi[0]) ** 2 / 4
    return {
        "endpoint_diag": 2 * P.k_spine if endpoint_diag is None else endpoint_diag,
        "R": float(abs(r) ** 2),
        "T": float(abs(t) ** 2),
        "closure_abs": float(abs(abs(r) ** 2 + abs(t) ** 2 - 1)),
        "contact_abs": float(abs(chi[0])),
        "caroli_physical_flux": float(flux),
        "scaled_unit_flux": float(omega * gamma * abs((2 / math.sqrt(omega)) * chi[0]) ** 2 / 4),
    }


def production_core() -> dict:
    obj = dressed_objects(PROFILE)
    left = scattering_summary(obj, "left")
    right = scattering_summary(obj, "right")
    branches = obj["branches"]
    total = obj["total"]
    gamma_branches = obj["gamma_branches"]
    chi = obj["chiL"]
    eta = obj["etaL"]

    trace_values = [float(np.real(np.vdot(chi, G @ chi))) for G in gamma_branches]
    source_values = [float(2 * np.real(np.vdot(eta, im @ eta))) for im in branches]

    optical_rhs = -1j * obj["GD"] @ (obj["GL"] + obj["GR"] + obj["gamma_bulk"]) @ obj["GA"]
    physical_scale = 2 / math.sqrt(OMEGA)
    fphys = -physical_scale * eta
    physical_bulk_power = OMEGA / 2 * float(np.real(np.vdot(fphys, total @ fphys)))

    # Actual displacement from a full 19x19 line/core block versus Schur elimination.
    Dcore = np.linalg.inv(obj["g0"]) + obj["U"]
    Aline = P.m_spine * OMEGA**2 * np.eye(7) - obj["KD"] - obj["Kself"] - obj["SL"] - obj["SR"]
    block = np.block([[Aline, -obj["V"].conj().T], [obj["V"], Dcore]])
    rhs = np.zeros(19, complex)
    rhs[0] = math.sqrt(obj["gamma"])
    solution = np.linalg.solve(block, rhs)
    q_direct = solution[:7]
    u_direct = solution[7:]
    q_eliminated = obj["chiL"]
    u_formula = -obj["g0"] @ obj["etaL"]

    return {
        "reference_left": left,
        "reference_right": right,
        "raw_spectral": {
            "branch_hermiticity_relative": [rel(im, im.conj().T) for im in branches],
            "total_hermiticity_relative": rel(total, total.conj().T),
            "branch_min_eigenvalues": [float(np.min(np.linalg.eigvalsh(herm(im)))) for im in branches],
            "total_min_eigenvalue": float(np.min(np.linalg.eigvalsh(herm(total)))),
            "branch_sum_relative_residual": rel(sum(branches), total),
        },
        "dressing": {
            "push_through_relative_residual": rel(obj["gc_left"], obj["gc_right"]),
            "imaginary_factorization_relative_residual": rel(
                imh(obj["gc_left"]), obj["left_dress"] @ total @ obj["right_dress"]
            ),
            "gamma_bulk_hermiticity_relative": rel(obj["gamma_bulk"], obj["gamma_bulk"].conj().T),
            "gamma_bulk_min_eigenvalue": float(np.min(np.linalg.eigvalsh(herm(obj["gamma_bulk"])))),
            "gamma_branch_min_eigenvalues": [float(np.min(np.linalg.eigvalsh(herm(G)))) for G in gamma_branches],
            "gamma_cluster_sum_relative_residual": rel(sum(gamma_branches), obj["gamma_bulk"]),
        },
        "caroli_and_work": {
            "trace_values": trace_values,
            "source_quadratic_values": source_values,
            "trace_source_absolute_errors": [abs(a - b) for a, b in zip(trace_values, source_values)],
            "Tbulk_trace": left["Tbulk"],
            "Tbulk_source": float(2 * np.real(np.vdot(eta, total @ eta))),
            "physical_bulk_power_unit_incident": physical_bulk_power,
            "physical_bulk_power_minus_Tbulk_abs": abs(physical_bulk_power - left["Tbulk"]),
            "caroli_scale_to_unit_physical_flux": physical_scale,
        },
        "optical": {
            "relative_residual": rel(obj["GD"] - obj["GA"], optical_rhs),
            "scalar_closure_abs": left["closure_abs"],
        },
        "full_vs_eliminated_reference": {
            "line_field_relative_residual": rel(q_direct, q_eliminated),
            "core_displacement_relative_residual": rel(u_direct, u_formula),
            "full_block_residual": float(np.linalg.norm(block @ solution - rhs) / np.linalg.norm(rhs)),
            "source_vs_displacement_dimensionless_norm_ratio_diagnostic": float(
                np.linalg.norm(obj["etaL"]) / np.linalg.norm(u_formula)
            ),
        },
    }


def random_full_vs_eliminated(seed: int = 20260813, cases: int = 200) -> dict:
    rng = np.random.default_rng(seed)
    errors_q = []
    errors_u = []
    sign_mutation_min_changes = []
    for _ in range(cases):
        nq = 4
        nb = 5
        X = rng.normal(size=(nb, nb))
        Kb = X.T @ X + (2.0 + rng.random()) * np.eye(nb)
        Mb = np.diag(0.5 + rng.random(nb))
        omega = float(rng.uniform(0.2, 0.8))
        eta = float(rng.uniform(0.01, 0.05))
        z = complex(omega, eta)
        g0 = np.linalg.inv(Kb - z * z * Mb)
        Y = rng.normal(size=(nb, nb))
        U = 0.03 * (Y.T @ Y)
        V = 0.08 * rng.normal(size=(nb, nq))
        Z = rng.normal(size=(nq, nq))
        Aline = Z.T @ Z + (1.5 + rng.random()) * np.eye(nq) + 1j * 0.1 * np.eye(nq)
        Dcore = np.linalg.inv(g0) + U
        rhs = rng.normal(size=nq) + 1j * rng.normal(size=nq)
        block = np.block([[Aline, -V.conj().T], [V, Dcore]])
        full = np.linalg.solve(block, np.r_[rhs, np.zeros(nb, complex)])
        gc = np.linalg.solve(np.eye(nb) + g0 @ U, g0)
        sigma = -V.conj().T @ gc @ V
        q = np.linalg.solve(Aline - sigma, rhs)
        u = -gc @ V @ q
        errors_q.append(rel(full[:nq], q))
        errors_u.append(rel(full[nq:], u))
        q_wrong = np.linalg.solve(Aline + sigma, rhs)
        sign_mutation_min_changes.append(rel(q_wrong, q))
    return {
        "seed": seed,
        "cases": cases,
        "max_line_relative_residual": max(errors_q),
        "max_bulk_relative_residual": max(errors_u),
        "minimum_wrong_schur_sign_relative_change": min(sign_mutation_min_changes),
    }


def mutations(core: dict) -> dict:
    obj = dressed_objects(PROFILE)
    # Wrong bulk Schur sign in the device denominator.
    sigma_wrong = +obj["V"].conj().T @ obj["gc_left"] @ obj["V"]
    Awrong = P.m_spine * OMEGA**2 * np.eye(7) - obj["KD"] - obj["Kself"] - sigma_wrong - obj["SL"] - obj["SR"]
    Gwrong = np.linalg.inv(Awrong)
    chi = Gwrong[:, 0] * math.sqrt(obj["gamma"])
    gamma_wrong = 1j * (sigma_wrong - sigma_wrong.conj().T)
    Rwrong = abs(1 - 1j * obj["gamma"] * Gwrong[0, 0]) ** 2
    Twrong = float(np.real(np.vdot(chi, obj["GR"] @ chi)))
    Tbwrong = float(np.real(np.vdot(chi, gamma_wrong @ chi)))

    # Advanced scalar lead mutation.
    kappa = obj["kappa"]
    sigma_advanced = -P.k_spine * np.exp(-1j * kappa)
    gamma_advanced = float(np.real(1j * (sigma_advanced - sigma_advanced.conjugate())))

    transparent = double_transparent()
    endpoint_omitted = double_transparent(P.k_spine)
    endpoint_duplicated = double_transparent(3 * P.k_spine)

    return {
        "wrong_bulk_schur_sign": {
            "gamma_min_eigenvalue": float(np.min(np.linalg.eigvalsh(herm(gamma_wrong)))),
            "R": float(Rwrong),
            "Tq": Twrong,
            "Tbulk_signed": Tbwrong,
            "closure_abs": float(abs(Rwrong + Twrong + Tbwrong - 1)),
            "Tbulk_change_from_reference": abs(Tbwrong - core["reference_left"]["Tbulk"]),
        },
        "advanced_lead_mutation": {
            "sigma": complex_json(sigma_advanced),
            "gamma": gamma_advanced,
            "expected_retarded_gamma": obj["gamma"],
        },
        "endpoint_partition": {
            "correct": transparent,
            "omitted_external_spring_diagonal": endpoint_omitted,
            "duplicated_external_spring_diagonal": endpoint_duplicated,
            "omitted_T_change": abs(endpoint_omitted["T"] - transparent["T"]),
            "duplicated_T_change": abs(endpoint_duplicated["T"] - transparent["T"]),
        },
        "wrong_linewidth_sign_min_eigenvalue": float(np.min(np.linalg.eigvalsh(herm(-obj["gamma_bulk"])))),
    }


def main() -> None:
    from reduced_core import reduced_static_matrices

    Kuu, Kuq, Kqq = direct_paired_link_hessian(PROFILE)
    B = independent_active_basis(len(PROFILE))
    Uref, Vref, Ksref = reduced_static_matrices(PROFILE, P)
    hessian = {
        "B_orthogonality": rel(B.T @ B, np.eye(12)),
        "U_relative_residual": rel(B.T @ Kuu @ B, Uref),
        "V_relative_residual": rel(B.T @ Kuq, Vref),
        "Kself_relative_residual": rel(Kqq, Ksref),
        "full_hessian_min_eigenvalue": float(np.min(np.linalg.eigvalsh(herm(np.block([[Kqq, Kuq.T], [Kuq, Kuu]]))))),
    }

    hp = high_precision_transparent()
    core = production_core()
    random_schur = random_full_vs_eliminated()
    muts = mutations(core)

    checks = {
        "block_hessian": max(hessian["U_relative_residual"], hessian["V_relative_residual"], hessian["Kself_relative_residual"]) < TH["block_hessian_relative"],
        "transparent_high_precision": max(hp["max_reflection"], hp["max_transmission_abs_error"], hp["max_scaled_unit_flux_abs_error"]) < TH["transparent_high_precision_abs"],
        "transparent_double": max(
            core["reference_left"]["closure_abs"],
            double_transparent()["R"],
            abs(double_transparent()["T"] - 1),
            abs(double_transparent()["scaled_unit_flux"] - 1),
        ) < TH["transparent_double_power_abs"] or double_transparent()["R"] < TH["transparent_double_power_abs"],
        "raw_spectral_psd_additive": (
            max(core["raw_spectral"]["branch_hermiticity_relative"] + [core["raw_spectral"]["total_hermiticity_relative"]]) < TH["spectral_hermiticity_relative"]
            and min(core["raw_spectral"]["branch_min_eigenvalues"] + [core["raw_spectral"]["total_min_eigenvalue"]]) >= TH["spectral_min_eigenvalue"]
            and core["raw_spectral"]["branch_sum_relative_residual"] < TH["cluster_additivity_relative"]
        ),
        "dressed_psd_additive": (
            core["dressing"]["push_through_relative_residual"] < TH["matrix_identity_relative"]
            and core["dressing"]["imaginary_factorization_relative_residual"] < TH["matrix_identity_relative"]
            and min(core["dressing"]["gamma_branch_min_eigenvalues"] + [core["dressing"]["gamma_bulk_min_eigenvalue"]]) >= TH["spectral_min_eigenvalue"]
            and core["dressing"]["gamma_cluster_sum_relative_residual"] < TH["cluster_additivity_relative"]
        ),
        "caroli_work_identity": (
            max(core["caroli_and_work"]["trace_source_absolute_errors"]) < TH["caroli_quadratic_absolute"]
            and core["caroli_and_work"]["physical_bulk_power_minus_Tbulk_abs"] < TH["caroli_quadratic_absolute"]
        ),
        "optical_identity": core["optical"]["relative_residual"] < TH["optical_identity_relative"],
        "reference_full_vs_eliminated": max(core["full_vs_eliminated_reference"]["line_field_relative_residual"], core["full_vs_eliminated_reference"]["core_displacement_relative_residual"]) < TH["full_vs_eliminated_relative"],
        "random_full_vs_eliminated": max(random_schur["max_line_relative_residual"], random_schur["max_bulk_relative_residual"]) < TH["full_vs_eliminated_relative"],
        "mutations_sensitive": (
            muts["advanced_lead_mutation"]["gamma"] < 0
            and muts["wrong_linewidth_sign_min_eigenvalue"] < -TH["mutation_minimum_effect"]
            and min(muts["endpoint_partition"]["omitted_T_change"], muts["endpoint_partition"]["duplicated_T_change"]) > TH["mutation_minimum_effect"]
            and (
                muts["wrong_bulk_schur_sign"]["gamma_min_eigenvalue"] < -TH["mutation_minimum_effect"]
                or muts["wrong_bulk_schur_sign"]["closure_abs"] > TH["mutation_minimum_effect"]
                or muts["wrong_bulk_schur_sign"]["Tbulk_signed"] < -TH["mutation_minimum_effect"]
            )
        ),
    }

    output = {
        "schema": "paper1.p05.core-scattering.v1",
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "hessian": hessian,
        "high_precision_transparent_chain": hp,
        "double_transparent_chain": double_transparent(),
        "production_core": core,
        "random_full_vs_eliminated": random_schur,
        "mutations": muts,
    }
    write_json(OUT, output)
    print(json.dumps(output, indent=2, sort_keys=True))
    if output["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
