#!/usr/bin/env python3
"""Bounded independent exact-data/scalar audit; not a contour reintegration."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    root = args.core.resolve()
    chain = json.loads((root / "tolerance/VALIDATED_PROOF_CHAIN.json").read_text())
    targets = json.loads((root / "r4_certificate/analysis/full_domain_r4/budget/TARGETS.json").read_text())["targets"]
    hashes = {}
    for label, rec in chain["files"].items():
        actual = hashlib.sha256((root / rec["path"]).read_bytes()).hexdigest()
        assert actual == rec["sha256"], label
        hashes[label] = actual

    # Fraction(float) is the exact stored IEEE-754 number, not a decimal surrogate.
    with np.load(root / "tolerance/inputs/REFERENCE_KERNEL.npz", allow_pickle=False) as data:
        G = data["G0"][:12, :12]
        Q = data["Im_qT"][:12, :12]
        assert G.shape == Q.shape == (12, 12)
        for i in range(12):
            for j in range(12):
                re, im = map(F, targets["g"][i][j])
                assert F(float(G[i, j].real)) == re
                assert F(float(G[i, j].imag)) == im
                re, im = map(F, targets["Q"][i][j])
                assert F(float(Q[i, j])) == re and im == 0
        floating_branch_sum_matches = np.array_equal(
            Q, (data["Im_branch_0"] + data["Im_branch_1"])[:12, :12]
        )
        assert floating_branch_sum_matches

    caps = chain["kernel_caps"]
    eG, eQ = F("1/250"), F("1/500")
    assert F(caps["eG_exact_saved_bound"]) <= eG
    assert F(caps["eQ_exact_saved_bound"]) <= eQ
    assert F(caps["Q_recentring_allowance"]) == 0

    # Lead spectral identity, evaluated without approximating sqrt(133).
    omega, k = F("6/5"), F("169/100")
    omega_shift = omega**2 - 2*k
    gamma_squared = F(36*133, 625)
    assert omega_shift == F("-97/50")
    assert 4*k*k-omega_shift**2 == gamma_squared
    assert omega_shift/2 == F("-97/100")
    assert gamma_squared/4 == F(9*133, 625)

    # Recompute advertised point interval independently from outward norm caps.
    h, j, z = F("947/1000"), F("3/5"), F("569/1000")
    al, au = F("4289/5000"), F("8579/10000")
    Z = z/(1-h*eG)
    d = j*eG*Z
    assert h*eG < 1 and al > d
    point_lo = (al-d)**2-2*eQ*Z*Z
    point_hi = (au+d)**2+2*eQ*Z*Z
    pt = chain["nominal_infinite_lattice_interval"]
    assert point_lo == F(pt["lower_rational"])
    assert point_hi == F(pt["upper_rational"])
    assert F(".732166") < point_lo < point_hi < F(".739652")
    # Nominal 0.736 lies inside; the interval is too broad to certify that rounding.
    assert point_lo < F(".736") < point_hi
    assert point_lo < F(".7355") and point_hi > F(".7365")

    # Independent exact propagation of the continuous-box auxiliary caps.
    r, w, v, ell = F(".164"), F(".068"), F(".1648"), F(".0173")
    we, de = F(".0649"), F(".0504")
    h0, j0, dg = F(".994"), F(".6289"), F("1.744")
    cached_floor = al*al-ell-2*au*w*v/(1-r)
    assert cached_floor == F("3634096757/5225000000")
    L = F(".8339")
    assert cached_floor > L*L
    Zc = z+de+we*v/(1-r)
    H = h0+we*dg/(1-r)
    J = j0+w*dg/(1-r)
    assert Zc < F(".633") and H < F("1.13") and J < F(".772")
    Zb = F(".633")/(1-F("1.13")*eG)
    db = F(".772")*eG*Zb
    assert F("1.13")*eG < 1 and L > db
    box_floor = (L-db)**2-2*eQ*Zb*Zb
    box_rec = chain["whole_six_motif_plus_minus5_percent_infinite_lattice_floor"]
    assert box_floor == F(box_rec["lower_rational"]) and box_floor > F(".69")
    report = {
        "status": "PASS",
        "scope": "Independent exact target/hash/lead/scalar binding check; mathematical proof and source review recorded separately.",
        "authenticated_proof_chain_hashes_checked": len(hashes),
        "hashes": hashes,
        "G_exact_complex_entries_checked": 144,
        "Q_exact_real_entries_checked": 144,
        "Q_matches_saved_floating_branch_sum": bool(floating_branch_sum_matches),
        "Q_interpreted_as_exact_target_not_exact_rational_sum_of_individual_branch_arrays": True,
        "kernel_caps_upward_and_zero_recentring_verified": True,
        "exact_real_frequency_lead_polynomial_identity_verified": True,
        "point_lower_rational": str(point_lo),
        "point_upper_rational": str(point_hi),
        "continuous_six_motif_box_lower_rational": str(box_floor),
        "not_performed": ["fresh global contour integration", "new optimization", "new 64-vertex exact replay", "new physical experiment"],
    }
    args.out.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k: report[k] for k in ["status", "authenticated_proof_chain_hashes_checked", "point_lower_rational", "point_upper_rational", "continuous_six_motif_box_lower_rational"]}, indent=2))


if __name__ == "__main__":
    main()
