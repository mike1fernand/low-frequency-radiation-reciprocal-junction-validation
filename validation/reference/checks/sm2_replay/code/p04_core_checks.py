from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.linalg import norm
from scipy.optimize import brentq, differential_evolution

from p04_common import (
    AXIAL_BONDS,
    FACE_DIAGONAL_BONDS,
    CORE_SIGNS,
    CORE_TRANS,
    KAPPA,
    OMEGA,
    P,
    PROFILE,
    TARGET,
    active_basis,
    aperture,
    bulk_dynamical_matrix,
    bulk_green_realspace,
    d_closed,
    eig,
    finite_core_matrices,
    jsonable,
    polynomial_ql_projector,
    projected_core_green,
    qt_projector,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(20260813)
    checks: dict[str, bool] = {}
    metrics: dict[str, object] = {}

    # Independent closed form versus the canonical bond-sum producer.
    errors = []
    for _ in range(5000):
        k = rng.uniform(-math.pi, math.pi, 3)
        errors.append(norm(bulk_dynamical_matrix(k, P) - d_closed(k), "fro"))
    metrics["D_closed_vs_bond_max_fro"] = float(max(errors))
    checks["D_identity"] = bool(metrics["D_closed_vs_bond_max_fro"] < 1e-12)

    # Long-wave fourth-order matrix remainder and acoustic speeds.
    direction = np.array([0.371, 0.557, 0.743], dtype=float)
    direction /= norm(direction)
    scales = np.logspace(-4, -1, 24)
    residuals = []
    for scale in scales:
        k = scale * direction
        leading = (P.k_diag / P.m_bulk) * (np.dot(k, k) * np.eye(3) + 2.0 * np.outer(k, k))
        residuals.append(norm(d_closed(k) - leading, "fro"))
    slope = float(np.polyfit(np.log(scales[5:]), np.log(residuals[5:]), 1)[0])
    metrics["long_wave_slope"] = slope
    checks["long_wave_order4"] = abs(slope - 4.0) < 0.02
    speed_samples = []
    for scale in [1e-4, 2e-4, 5e-4, 1e-3]:
        speed_samples.append((np.sqrt(eig(scale * direction)[0]) / scale).tolist())
    metrics["speed_samples"] = speed_samples
    checks["speed_limit"] = bool(
        np.max(np.abs(np.array(speed_samples[-1]) - np.array([1.0, 1.0, math.sqrt(3.0)]))) < 2e-6
    )

    # Exact line carrier.
    line = {
        "kappa_a": KAPPA,
        "lambda_over_a": 2.0 * math.pi / KAPPA,
        "vphase_over_cT": OMEGA / KAPPA,
        "vg_over_cT": math.sqrt(P.k_spine / P.m_spine) * math.cos(KAPPA / 2.0),
    }
    metrics["line"] = line
    checks["line_carrier"] = abs(line["kappa_a"] - 0.9594572753515053) < 1e-12

    # Fixed-kappa full transverse-zone minima and maxima.
    fixed_ranges = []
    references = [
        (0.923076923076923, 2.6735611072243266),
        (0.923076923076923, 2.8284271247461903),
        (1.5988161300636134, 3.3388514483400207),
    ]
    for branch in range(3):
        def frequency(x):
            return math.sqrt(max(float(eig(np.array([KAPPA, x[0], x[1]]))[0][branch]), 0.0))

        rmin = differential_evolution(
            frequency,
            [(-math.pi, math.pi), (-math.pi, math.pi)],
            seed=10 + branch,
            tol=1e-11,
            maxiter=180,
            popsize=15,
            polish=True,
            workers=1,
        )
        rmax = differential_evolution(
            lambda x: -frequency(x),
            [(-math.pi, math.pi), (-math.pi, math.pi)],
            seed=20 + branch,
            tol=1e-11,
            maxiter=180,
            popsize=15,
            polish=True,
            workers=1,
        )
        fixed_ranges.append(
            {
                "branch": branch,
                "min": float(rmin.fun),
                "max": float(-rmax.fun),
                "min_kperp": rmin.x,
                "max_kperp": rmax.x,
                "min_optimizer_success": bool(rmin.success),
                "max_optimizer_success": bool(rmax.success),
            }
        )
    metrics["fixed_kappa_ranges"] = fixed_ranges
    range_error = max(
        max(abs(fixed_ranges[b]["min"] - references[b][0]), abs(fixed_ranges[b]["max"] - references[b][1]))
        for b in range(3)
    )
    checks["fixed_kappa_ranges"] = bool(
        range_error < 2e-9
        and all(x["min_optimizer_success"] and x["max_optimizer_success"] for x in fixed_ranges)
    )

    # Radial shell-root census along deterministic directions.
    n_directions = 1200
    idx = np.arange(n_directions)
    z = 1.0 - 2.0 * (idx + 0.5) / n_directions
    phi = (math.pi * (3.0 - math.sqrt(5.0)) * idx) % (2.0 * math.pi)
    rxy = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    directions = np.column_stack([rxy * np.cos(phi), rxy * np.sin(phi), z])
    root_counts: dict[str, list[int]] = {str(branch): [] for branch in range(3)}
    for direction in directions:
        nonzero = np.abs(direction) > 1e-14
        rmax = float(np.min(math.pi / np.abs(direction[nonzero])))
        radii = np.linspace(0.0, rmax, 181)
        values = np.array([eig(radius * direction)[0] - TARGET for radius in radii])
        for branch in range(3):
            root_counts[str(branch)].append(int(np.sum(values[:-1, branch] * values[1:, branch] < 0.0)))
    root_census = {
        key: {"min": min(values), "max": max(values), "directions": len(values)}
        for key, values in root_counts.items()
    }
    metrics["root_census"] = root_census
    checks["single_shell_root_census"] = all(
        item["min"] == 1 and item["max"] == 1 for item in root_census.values()
    )

    # Basis-free polynomial projector for the isolated upper branch.
    projector_errors = []
    operational_errors = []
    idempotence = []
    hermiticity = []
    for _ in range(20000):
        k = rng.uniform(-math.pi, math.pi, 3)
        vals, vecs = eig(k)
        if vals[2] - vals[1] < 1e-8:
            continue
        p_poly = polynomial_ql_projector(d_closed(k), vals)
        p_eig = np.outer(vecs[:, 2], vecs[:, 2])
        error = norm(p_poly - p_eig, "fro")
        projector_errors.append(error)
        if vals[2] - vals[1] > 0.3:
            operational_errors.append(error)
        idempotence.append(norm(p_poly @ p_poly - p_poly, "fro"))
        hermiticity.append(norm(p_poly - p_poly.T, "fro"))
    metrics["polynomial_projector"] = {
        "max_error_all_sampled_BZ": float(max(projector_errors)),
        "mean_error_all_sampled_BZ": float(np.mean(projector_errors)),
        "max_error_gap_above_0p3": float(max(operational_errors)),
        "mean_error_gap_above_0p3": float(np.mean(operational_errors)),
        "max_idempotence": float(max(idempotence)),
        "max_hermiticity": float(max(hermiticity)),
        "samples": len(projector_errors),
    }
    checks["polynomial_projector"] = bool(
        max(operational_errors) < 1e-11 and max(hermiticity) < 1e-11
    )

    # C4v core representation and exact A1 projector.
    corner_positions = [tuple(x) for x in CORE_SIGNS.astype(int)]
    position_to_index = {position: i for i, position in enumerate(corner_positions)}
    operations = []
    for nrot in range(4):
        theta = nrot * math.pi / 2.0
        operations.append(
            np.rint(
                np.array(
                    [[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]]
                )
            ).astype(int)
        )
    reflection = np.array([[1, 0], [0, -1]], dtype=int)
    for nrot in range(4):
        theta = nrot * math.pi / 2.0
        rotation = np.rint(
            np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
        ).astype(int)
        operations.append(rotation @ reflection)
    representations = []
    for op in operations:
        r3 = np.eye(3)
        r3[1:, 1:] = op
        t = np.zeros((12, 12), dtype=float)
        for old, (sy, sz) in enumerate(corner_positions):
            new_position = tuple((op @ np.array([sy, sz])).astype(int))
            new = position_to_index[new_position]
            t[3 * new : 3 * new + 3, 3 * old : 3 * old + 3] = r3
        representations.append(t)
    p_group = sum(representations) / 8.0
    b1 = active_basis(1)
    p_basis = b1 @ b1.conj().T
    c4v = {
        "rank": int(np.linalg.matrix_rank(p_group, tol=1e-12)),
        "group_vs_XR": float(norm(p_group - p_basis, "fro")),
        "idempotence": float(norm(p_group @ p_group - p_group, "fro")),
        "unitarity": float(max(norm(t.T @ t - np.eye(12), "fro") for t in representations)),
        "closure": float(
            max(
                min(norm(ta @ tb - tc, "fro") for tc in representations)
                for ta in representations
                for tb in representations
            )
        ),
    }
    metrics["C4v"] = c4v
    checks["C4v_projector"] = c4v["rank"] == 2 and c4v["group_vs_XR"] < 1e-12

    off_symmetry = []
    for coordinate in [0, 1, 2, 4, 8]:
        e = np.zeros(12)
        e[coordinate] = 1.0
        off_symmetry.append(float(norm((np.eye(12) - p_group) @ e)))
    metrics["off_symmetry_complement_norms"] = off_symmetry
    checks["off_symmetry_control"] = min(off_symmetry) > 0.5

    # Positive control proving the local finite core is not invariant under K_b.
    bonds = [(tuple(d), P.k_axis) for d in AXIAL_BONDS] + [
        (tuple(d), P.k_diag) for d in FACE_DIAGONAL_BONDS
    ]
    core_sites = [(0, int(y), int(z)) for y, z in CORE_TRANS]

    def stiffness_action(kind: str):
        displacement = {site: np.zeros(3) for site in core_sites}
        for site, (sy, sz) in zip(core_sites, CORE_SIGNS):
            displacement[site] = (
                np.array([0.5, 0.0, 0.0])
                if kind == "X"
                else np.array([0.0, sy / (2.0 * math.sqrt(2.0)), sz / (2.0 * math.sqrt(2.0))])
            )
        pairs = set()
        for r in core_sites:
            for d, stiffness in bonds:
                for sign in (-1, 1):
                    s = tuple(r[j] + sign * d[j] for j in range(3))
                    pairs.add((tuple(sorted((r, s))), stiffness))
        ku: dict[tuple[int, int, int], np.ndarray] = {}
        for (r, s), stiffness in pairs:
            vector = np.array(s) - np.array(r)
            nhat = vector / norm(vector)
            force = stiffness * np.outer(nhat, nhat) @ (
                displacement.get(r, np.zeros(3)) - displacement.get(s, np.zeros(3))
            )
            ku[r] = ku.get(r, np.zeros(3)) + force
            ku[s] = ku.get(s, np.zeros(3)) - force
        core_norm = math.sqrt(sum(float(x @ x) for r, x in ku.items() if r in core_sites))
        outside_norm = math.sqrt(sum(float(x @ x) for r, x in ku.items() if r not in core_sites))
        return {
            "core_norm": core_norm,
            "outside_norm": outside_norm,
            "outside_fraction": outside_norm / math.sqrt(core_norm**2 + outside_norm**2),
            "outside_sites": sum(norm(x) > 1e-14 for r, x in ku.items() if r not in core_sites),
        }

    outside_core = {kind: stiffness_action(kind) for kind in ["X", "R"]}
    metrics["outside_core"] = outside_core
    checks["outside_core_positive_control"] = all(
        item["outside_norm"] > 0.1 for item in outside_core.values()
    )

    # Exact 72-to-12 projected-resolvent reduction on two independent FFT grids.
    reduction: dict[str, dict[str, float | int]] = {}
    for n in [24, 32]:
        coords, u, v = finite_core_matrices(PROFILE)
        b = active_basis(6)
        p_core = b @ b.conj().T
        g = projected_core_green(bulk_green_realspace(n), coords)
        g_red = b.conj().T @ g @ b
        u_red = b.conj().T @ u @ b
        v_red = b.conj().T @ v
        g_c_full = np.linalg.solve(np.eye(72) + g @ u, g)
        g_c_red = np.linalg.solve(np.eye(12) + g_red @ u_red, g_red)
        field_errors = []
        self_energy_errors = []
        for _ in range(50):
            q = rng.normal(size=7) + 1j * rng.normal(size=7)
            x_full = g_c_full @ v @ q
            x_red = b @ g_c_red @ v_red @ q
            field_errors.append(norm(x_full - x_red) / max(norm(x_full), 1e-30))
            sigma_full = v.conj().T @ x_full
            sigma_red = v_red.conj().T @ (g_c_red @ v_red @ q)
            self_energy_errors.append(norm(sigma_full - sigma_red) / max(norm(sigma_full), 1e-30))
        # Krylov closure inside the A1 block; this is distinct from direct source rank.
        a_op = g_red @ u_red
        blocks = []
        current = v_red.copy()
        for _ in range(12):
            blocks.append(current)
            current = a_op @ current
        singular_values = np.linalg.svd(np.concatenate(blocks, axis=1), compute_uv=False)
        krylov_rank = int(np.sum(singular_values > 1e-10 * singular_values[0]))
        reduction[str(n)] = {
            "orth": float(norm(b.conj().T @ b - np.eye(12), "fro")),
            "V_leak": float(norm((np.eye(72) - p_core) @ v, "fro") / norm(v, "fro")),
            "U_leak": float(norm((np.eye(72) - p_core) @ u @ b, "fro") / norm(u @ b, "fro")),
            "G_leak": float(norm((np.eye(72) - p_core) @ g @ b, "fro") / norm(g @ b, "fro")),
            "G_commutator": float(norm(g @ p_core - p_core @ g, "fro") / norm(g, "fro")),
            "full_reduced_field_max": float(max(field_errors)),
            "effective_self_energy_max": float(max(self_energy_errors)),
            "direct_source_rank": int(np.linalg.matrix_rank(v_red, tol=1e-12)),
            "krylov_rank": krylov_rank,
            "krylov_smallest_retained_singular_value": float(singular_values[krylov_rank - 1]),
        }
    metrics["reduction"] = reduction
    checks["exact_reduction"] = bool(
        max(
            max(
                float(item[key])
                for key in [
                    "V_leak",
                    "U_leak",
                    "G_leak",
                    "G_commutator",
                    "full_reduced_field_max",
                    "effective_self_energy_max",
                ]
            )
            for item in reduction.values()
        )
        < 5e-13
        and all(item["direct_source_rank"] == 7 for item in reduction.values())
        and all(item["krylov_rank"] == 12 for item in reduction.values())
    )

    # Correct U(2) covariance and deliberately incorrect fixed unequal branch weights.
    max_projector_change = 0.0
    max_covariant_change = 0.0
    max_wrong_fixed_weight_change = 0.0
    for _ in range(1000):
        k = rng.uniform(-2.0, 2.0, 3)
        vals, vecs = eig(k)
        if vals[2] - vals[1] < 0.1:
            continue
        e = vecs[:, :2]
        zmat = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        qmat, rmat = np.linalg.qr(zmat)
        unitary = qmat @ np.diag(np.exp(-1j * np.angle(np.diag(rmat))))
        e_transformed = e @ unitary
        max_projector_change = max(
            max_projector_change,
            float(norm(e_transformed @ e_transformed.conj().T - e @ e.conj().T, "fro")),
        )
        operator = np.array([[1.0, 0.2j], [-0.2j, 2.0]], dtype=complex)
        state = rng.normal(size=2) + 1j * rng.normal(size=2)
        transformed_state = unitary.conj().T @ state
        transformed_operator = unitary.conj().T @ operator @ unitary
        max_covariant_change = max(
            max_covariant_change,
            float(abs(np.vdot(state, operator @ state) - np.vdot(transformed_state, transformed_operator @ transformed_state))),
        )
        field = rng.normal(size=3)
        weights = np.array([1.0, 2.0])
        wrong_change = abs(
            np.sum(weights * np.abs(e.T @ field) ** 2)
            - np.sum(weights * np.abs(e_transformed.conj().T @ field) ** 2)
        )
        max_wrong_fixed_weight_change = max(max_wrong_fixed_weight_change, float(wrong_change))
    metrics["U2"] = {
        "projector_max_change": max_projector_change,
        "covariant_scalar_max_change": max_covariant_change,
        "fixed_unequal_weight_max_change": max_wrong_fixed_weight_change,
    }
    checks["U2_scope"] = bool(
        max_projector_change < 1e-12
        and max_covariant_change < 1e-12
        and max_wrong_fixed_weight_change > 0.1
    )

    # Six-segment aperture spectrum: exact phase-removed closed form and widths.
    a0 = abs(aperture(0.0))
    half_amplitude = brentq(lambda q: abs(aperture(q)) / a0 - 0.5, 0.0, math.pi)
    half_power = brentq(lambda q: (abs(aperture(q)) / a0) ** 2 - 0.5, 0.0, math.pi)
    qgrid = np.linspace(-math.pi, math.pi, 10001)
    analytic = 0.2 * np.cos(2.5 * qgrid) + 0.4 * np.cos(1.5 * qgrid) + 0.4 * np.cos(0.5 * qgrid)
    direct = np.array([(aperture(q) * np.exp(-3j * q)).real for q in qgrid])
    wavelength = 2.0 * math.pi / KAPPA
    aperture_metrics = {
        "amplitude_FWHM_delta_k_a": 2.0 * half_amplitude,
        "power_FWHM_delta_k_a": 2.0 * half_power,
        "first_positive_null_delta_k_a": 2.0 * math.pi / 5.0,
        "central_null_to_null_width_delta_k_a": 4.0 * math.pi / 5.0,
        "interaction_length_over_lambda": 6.0 / wavelength,
        "midpoint_extent_over_lambda": 5.0 / wavelength,
        "analytic_max_abs_error": float(np.max(np.abs(analytic - direct))),
    }
    metrics["aperture"] = aperture_metrics
    checks["aperture_formula"] = bool(
        aperture_metrics["analytic_max_abs_error"] < 1e-14
        and 1.05 < aperture_metrics["power_FWHM_delta_k_a"] < 1.10
    )

    # Exact source norm identity and nodal-set control.
    nodal_errors = []
    for _ in range(10000):
        ky, kz = rng.uniform(-math.pi, math.pi, 2)
        source = np.array(
            [math.sin(ky / 2.0) * math.cos(kz / 2.0), math.cos(ky / 2.0) * math.sin(kz / 2.0)]
        )
        rhs = 0.5 * (1.0 - math.cos(ky) * math.cos(kz))
        nodal_errors.append(abs(float(source @ source) - rhs))
    metrics["source_nodal_identity_max"] = float(max(nodal_errors))
    checks["source_nodal_identity"] = metrics["source_nodal_identity_max"] < 1e-14

    # Counterexample: unweighted qT projection need not equal resonant branch overlap.
    k_counter = np.array([-0.7069989228994551, 1.0255869436172618, 0.015447407579766433])
    vals, vecs = eig(k_counter)
    field = vecs[:, 1]
    p_qt = qt_projector(vecs)
    counterexample = {
        "PT_norm_sq": float(np.vdot(field, p_qt @ field).real),
        "resonant_branch0_overlap_sq": float(abs(np.vdot(vecs[:, 0], field)) ** 2),
        "lower_split": float(vals[1] - vals[0]),
    }
    metrics["projector_overlap_counterexample"] = counterexample
    checks["projector_counterexample"] = bool(
        counterexample["PT_norm_sq"] > 0.999
        and counterexample["resonant_branch0_overlap_sq"] < 1e-25
    )

    status = "passed" if all(checks.values()) else "failed"
    output = {
        "schema": "paper1.p04.core-checks.v1",
        "status": status,
        "summary": {"passed": sum(bool(value) for value in checks.values()), "total": len(checks)},
        "checks": checks,
        "metrics": metrics,
        "producer_contract": {
            "fixed_random_seed": 20260813,
            "deterministic_inputs": True,
            "contains_timestamps": False,
            "external_reproduction_claim": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True, default=jsonable) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "summary": output["summary"]}, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
