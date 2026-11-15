#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from p05_common import (
    ASYM_PROFILE,
    OMEGA,
    P,
    PROFILE,
    ROOT,
    complex_json,
    dressed_objects,
    rel,
    rel_scalar,
    scattering_summary,
    write_json,
)

CONTRACT = json.loads((ROOT / "P05_TEST_CONTRACT.json").read_text(encoding="utf-8"))
TH = CONTRACT["thresholds"]
OUT = ROOT / "results" / "P05_MATCHING_PLANE.json"


def solve_padded(
    profile: np.ndarray,
    left_pad: int,
    right_pad: int,
    incidence: str,
    *,
    wrong_attachment: bool = False,
    defect_incident_interface: bool = False,
) -> dict:
    base = dressed_objects(profile)
    n_active = len(profile) + 1
    n_total = left_pad + n_active + right_pad
    active_start = left_pad
    active_stop = left_pad + n_active

    K = np.zeros((n_total, n_total), dtype=float)
    np.fill_diagonal(K, 2 * P.k_spine)
    for i in range(n_total - 1):
        kval = P.k_spine
        if defect_incident_interface:
            if incidence == "left" and left_pad > 0 and i == left_pad - 1:
                kval *= 1.01
            elif incidence == "right" and right_pad > 0 and i == active_stop - 1:
                kval *= 1.01
        K[i, i + 1] = K[i + 1, i] = -kval
        if kval != P.k_spine:
            delta = kval - P.k_spine
            K[i, i] += delta
            K[i + 1, i + 1] += delta

    Kself = np.zeros((n_total, n_total), dtype=float)
    Kself[active_start:active_stop, active_start:active_stop] = base["Kself"]
    sigma_bulk = np.zeros((n_total, n_total), dtype=complex)
    sigma_bulk[active_start:active_stop, active_start:active_stop] = base["sigma_bulk"]
    gamma_bulk = np.zeros((n_total, n_total), dtype=complex)
    gamma_bulk[active_start:active_stop, active_start:active_stop] = base["gamma_bulk"]
    gamma_branches = []
    for G in base["gamma_branches"]:
        embedded = np.zeros((n_total, n_total), dtype=complex)
        embedded[active_start:active_stop, active_start:active_stop] = G
        gamma_branches.append(embedded)

    left_contact = 0
    right_contact = n_total - 1
    if wrong_attachment:
        if left_pad > 0:
            left_contact = active_start
        if right_pad > 0:
            right_contact = active_stop - 1

    SL = np.zeros((n_total, n_total), dtype=complex)
    SR = np.zeros_like(SL)
    SL[left_contact, left_contact] = base["sigma"]
    SR[right_contact, right_contact] = base["sigma"]
    GL = 1j * (SL - SL.conj().T)
    GR = 1j * (SR - SR.conj().T)

    A = P.m_spine * OMEGA**2 * np.eye(n_total) - K - Kself - sigma_bulk - SL - SR
    GD = np.linalg.inv(A)
    if incidence == "left":
        contact = left_contact
        other = right_contact
        outgoing_gamma = GR
        pad_incident = left_pad
    elif incidence == "right":
        contact = right_contact
        other = left_contact
        outgoing_gamma = GL
        pad_incident = right_pad
    else:
        raise ValueError

    chi = GD[:, contact] * math.sqrt(base["gamma"])
    r = 1 - 1j * base["gamma"] * GD[contact, contact]
    t = -1j * base["gamma"] * GD[other, contact]
    R = float(abs(r) ** 2)
    Tq = float(np.real(np.vdot(chi, outgoing_gamma @ chi)))
    Tbulk = float(np.real(np.vdot(chi, gamma_bulk @ chi)))
    Tbranches = [float(np.real(np.vdot(chi, G @ chi))) for G in gamma_branches]

    return {
        "left_pad": left_pad,
        "right_pad": right_pad,
        "incidence": incidence,
        "R": R,
        "Tq": Tq,
        "Tbulk": Tbulk,
        "Tbranches": Tbranches,
        "closure_abs": float(abs(R + Tq + Tbulk - 1)),
        "r": r,
        "t": t,
        "active_field": chi[active_start:active_stop],
        "incident_pad": pad_incident,
    }


def compare_profile(label: str, profile: np.ndarray) -> dict:
    kappa = dressed_objects(profile)["kappa"]
    baselines = {
        "left": solve_padded(profile, 0, 0, "left"),
        "right": solve_padded(profile, 0, 0, "right"),
    }
    rows = []
    maxima = {
        "R_relative": 0.0,
        "Tq_relative": 0.0,
        "Tbulk_relative": 0.0,
        "branch_relative": 0.0,
        "active_field_deembedded_relative": 0.0,
        "reflection_deembedded_relative": 0.0,
        "transmission_deembedded_relative": 0.0,
        "closure_abs": 0.0,
    }
    for incidence in ("left", "right"):
        b = baselines[incidence]
        for left_pad in range(5):
            for right_pad in range(5):
                row = solve_padded(profile, left_pad, right_pad, incidence)
                incident_pad = left_pad if incidence == "left" else right_pad
                active_deembedded = np.exp(-1j * kappa * incident_pad) * row["active_field"]
                reflection_deembedded = np.exp(-2j * kappa * incident_pad) * row["r"]
                transmission_deembedded = np.exp(-1j * kappa * (left_pad + right_pad)) * row["t"]
                metrics = {
                    "R_relative": rel_scalar(row["R"], b["R"]),
                    "Tq_relative": rel_scalar(row["Tq"], b["Tq"]),
                    "Tbulk_relative": rel_scalar(row["Tbulk"], b["Tbulk"]),
                    "branch_relative": max(rel_scalar(a, c) for a, c in zip(row["Tbranches"], b["Tbranches"])),
                    "active_field_deembedded_relative": rel(active_deembedded, b["active_field"]),
                    "reflection_deembedded_relative": rel(reflection_deembedded, b["r"]),
                    "transmission_deembedded_relative": rel(transmission_deembedded, b["t"]),
                    "closure_abs": row["closure_abs"],
                }
                for key, value in metrics.items():
                    maxima[key] = max(maxima[key], float(value))
                rows.append(
                    {
                        "incidence": incidence,
                        "left_pad": left_pad,
                        "right_pad": right_pad,
                        "R": row["R"],
                        "Tq": row["Tq"],
                        "Tbulk": row["Tbulk"],
                        "Tbranches": row["Tbranches"],
                        "closure_abs": row["closure_abs"],
                        "metrics": metrics,
                    }
                )

    # Deliberately wrong placement: two pristine cells are left in the device but the self-energy is attached at the active boundary.
    wrong = solve_padded(profile, 2, 2, "left", wrong_attachment=True)
    correct = solve_padded(profile, 2, 2, "left")
    wrong_changes = {
        "R_relative_change": rel_scalar(wrong["R"], correct["R"]),
        "Tq_relative_change": rel_scalar(wrong["Tq"], correct["Tq"]),
        "Tbulk_relative_change": rel_scalar(wrong["Tbulk"], correct["Tbulk"]),
    }

    # Physical 1% bond change at the added incident-side/active interface.
    defect = solve_padded(profile, 1, 0, "left", defect_incident_interface=True)
    correct_defect_reference = solve_padded(profile, 1, 0, "left")
    defect_changes = {
        "R_relative_change": rel_scalar(defect["R"], correct_defect_reference["R"]),
        "Tq_relative_change": rel_scalar(defect["Tq"], correct_defect_reference["Tq"]),
        "Tbulk_relative_change": rel_scalar(defect["Tbulk"], correct_defect_reference["Tbulk"]),
    }

    left = baselines["left"]
    right = baselines["right"]
    reciprocity = {
        "transmission_power_relative": rel_scalar(left["Tq"], right["Tq"]),
        "complex_transmission_relative": rel(left["t"], right["t"]),
        "reflection_power_relative": rel_scalar(left["R"], right["R"]),
        "bulk_power_relative": rel_scalar(left["Tbulk"], right["Tbulk"]),
    }

    pass_invariance = (
        max(maxima[k] for k in ["R_relative", "Tq_relative", "Tbulk_relative", "branch_relative"])
        < TH["matching_power_relative"]
        and max(
            maxima[k]
            for k in [
                "active_field_deembedded_relative",
                "reflection_deembedded_relative",
                "transmission_deembedded_relative",
            ]
        )
        < TH["matching_complex_relative"]
        and maxima["closure_abs"] < TH["matching_closure_absolute"]
    )
    pass_negative = (
        max(wrong_changes.values()) > TH["mutation_minimum_effect"]
        and max(defect_changes.values()) > TH["mutation_minimum_effect"]
    )
    pass_reciprocity = (
        reciprocity["transmission_power_relative"] < TH["matching_power_relative"]
        and reciprocity["complex_transmission_relative"] < TH["matching_complex_relative"]
    )

    return {
        "label": label,
        "profile": profile.tolist(),
        "baselines": {
            side: {
                "R": row["R"],
                "Tq": row["Tq"],
                "Tbulk": row["Tbulk"],
                "Tbranches": row["Tbranches"],
                "closure_abs": row["closure_abs"],
                "r": complex_json(row["r"]),
                "t": complex_json(row["t"]),
            }
            for side, row in baselines.items()
        },
        "maxima_over_25_partitions_and_both_directions": maxima,
        "reciprocity": reciprocity,
        "negative_controls": {
            "wrong_attachment": wrong_changes,
            "one_percent_incident_interface_bond_defect": defect_changes,
        },
        "checks": {
            "partition_invariance": pass_invariance,
            "line_transmission_reciprocity": pass_reciprocity,
            "negative_controls_sensitive": pass_negative,
            "asymmetric_bulk_loss_not_forced_equal": (
                True
                if label == "reference_symmetric"
                else reciprocity["bulk_power_relative"] > TH["mutation_minimum_effect"]
            ),
        },
        "rows": rows,
    }


def main() -> None:
    profiles = [
        compare_profile("reference_symmetric", PROFILE),
        compare_profile("asymmetric_profile_control", ASYM_PROFILE),
    ]
    checks = {
        f"{item['label']}::{key}": value
        for item in profiles
        for key, value in item["checks"].items()
    }
    output = {
        "schema": "paper1.p05.matching-plane.v1",
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "profiles": profiles,
        "phase_contract": {
            "active_field": "multiply by exp(-i*kappa*incident_pad)",
            "reflection": "multiply by exp(-2i*kappa*incident_pad)",
            "transmission": "multiply by exp(-i*kappa*(left_pad+right_pad))",
        },
    }
    write_json(OUT, output)
    print(json.dumps(output, indent=2, sort_keys=True))
    if output["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
