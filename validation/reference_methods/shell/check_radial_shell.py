#!/usr/bin/env python3
"""Exact algebra and bounded falsification checks for a model-specific shell lemma.

This is not a lattice Green-function integration or a quadrature-error proof.
The exact identities and elementary inequalities are proved in RADIAL_SHELL_PROOF.md.
"""
from __future__ import annotations
import hashlib
import argparse
import importlib.util
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
import sympy as sp
import flint
from flint import arb, fmpq, ctx

HERE = Path(__file__).resolve().parent
REPO = HERE.parent / "baseline/paper1-validation-repository"
DEFAULT_SOURCE = REPO / "validation/current_numerics/analysis/math_resolution/tetra/combined_tetra_diagnostic.py"
EXPECTED_SOURCE_SHA256 = "a04f1469250ce240c8802d1e788aec8a5c2f2c76c47438a88b02cb6b0d9c4818"
ENERGY = 36 / 25
BONDS = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                  [1, 1, 0], [1, -1, 0], [1, 0, 1],
                  [1, 0, -1], [0, 1, 1], [0, 1, -1]], dtype=float)
UNITS = BONDS / np.linalg.norm(BONDS, axis=1)[:, None]


def bulk(k):
    theta = np.asarray(k) @ BONDS.T
    return np.einsum("b,bi,bj->ij", 2 * (1 - np.cos(theta)), UNITS, UNITS)


def radial_derivative(k):
    radius = np.linalg.norm(k)
    theta = np.asarray(k) @ BONDS.T
    return np.einsum("b,bi,bj->ij", 2 * theta * np.sin(theta) / radius, UNITS, UNITS)


def rank_one(k):
    c, s = np.cos(k), np.sin(k)
    d = 2 + c * c - c * (1 + c.sum() - c)
    return 2 * (np.diag(d) + np.outer(s, s))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk-source", type=Path, default=DEFAULT_SOURCE,
                        help="Accepted combined_tetra_diagnostic.py, with stable_scalar.py and centroid_weights.py beside it")
    parser.add_argument("--output", type=Path, default=HERE / "RADIAL_SHELL_CHECKS.json")
    args = parser.parse_args()
    source = args.bulk_source.resolve()
    if hashlib.sha256(source.read_bytes()).hexdigest() != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("The accepted bulk source identity changed; this proof has a fixed model contract")
    # Exact polynomial identities underlying the exterior bound.
    a, b, c, q = sp.symbols("a b c q", real=True)
    da, db = 2 + a*a - a*(1+b+c), 2 + b*b - b*(1+a+c)
    identities = {
        "negative_component": da - sp.Rational(7, 4) - ((a+sp.Rational(1, 2))**2-a*(2+b+c)),
        "positive_exterior_component": da - (2+q*q-3*q) - ((q-a)*(3-a-q)+a*(2-b-c)),
        "exterior_lower_gap": (2+q*q-3*q)-(1-q)-(1-q)**2,
        "positive_other_component": db-(1-q)-((1-b)**2+q*(1-b)+b*(q-a+1-c)),
        "negative_other_component": db-sp.Rational(7, 4)-((b+sp.Rational(1, 2))**2-b*(2+a+c)),
        "isolated_axis_comparison": (3-a*(2+b))-(2-a-b)-(1-a)*(1+b),
        "coordinate_plane_lower_pair_degeneracy":
            (6-2*a*(2+b)-(4-2*a-2*b))*(6-2*b*(2+a)-(4-2*a-2*b))
            -4*(1-a*a)*(1-b*b),
    }
    u, v, w = sp.symbols("u v w", real=True)
    total_u = u + v + w
    d_u, d_v = total_u + u*(2*u-total_u), total_u + v*(2*v-total_u)
    identities.update({
        "gap_u_decomposition": da.subs({a:1-u,b:1-v,c:1-w})-d_u,
        "gap_diagonal_difference": d_u-d_v-(u-v)*(u+v-w),
        "gap_trace_identity": 2*(3*total_u+2*(u*u+v*v+w*w)-total_u*total_u
                                    +2*total_u-(u*u+v*v+w*w))
                              -(10*total_u+2*(u*u+v*v+w*w)-2*total_u*total_u),
    })
    exact = {name: sp.expand(value) == 0 for name, value in identities.items()}
    if not all(exact.values()):
        raise RuntimeError("An exact identity failed")
    # The rank-one equality uses s_i^2=1-c_i^2; verify its three diagonal entries.
    for x, y, z in [(a,b,c), (b,a,c), (c,a,b)]:
        if sp.expand(2*(2+x*x-x*(1+y+z)+1-x*x)-(6-2*x*(1+y+z))) != 0:
            raise RuntimeError("Rank-one diagonal identity failed")

    ctx.prec = 160
    q_arb, e_arb = arb(fmpq(1, 4)), arb(fmpq(36, 25))
    radius_arb = q_arb.acos()
    beta_arb = 2 * radius_arb / arb(15).sqrt()
    slope_arb = 2 * e_arb / arb(45).sqrt()
    if not (fmpq(3,2)-fmpq(36,25) == fmpq(3,50) and beta_arb > 0 and slope_arb > arb(fmpq(429, 1000))):
        raise RuntimeError("Exact-energy or outward numerical constants failed")

    rng = np.random.default_rng(20260920)
    random_points = rng.uniform(-np.pi, np.pi, (2000, 3))
    corners = np.array(np.meshgrid(*([np.array([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])]*3))).reshape(3,-1).T
    points = np.vstack((random_points, corners))
    identity_error = max(float(np.max(np.abs(bulk(k)-rank_one(k)))) for k in points)
    if identity_error > 2e-14:
        raise RuntimeError("Independent bond sum disagrees with rank-one expression")
    sys.path.insert(0, str(source.parent))
    spec = importlib.util.spec_from_file_location("paper1_accepted_bulk_for_shell_check", source)
    accepted = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = accepted
    spec.loader.exec_module(accepted)
    parameters = accepted.Parameters()
    if (parameters.a, parameters.m_bulk, parameters.k_axis, parameters.k_diag) != (1,1,1,1):
        raise RuntimeError("The controlling bare-bulk parameters changed")
    actual_matrices = accepted.bulk_stiffness(points, parameters) / parameters.m_bulk
    producer_error = max(float(np.max(np.abs(actual_matrices[i]-rank_one(k)))) for i,k in enumerate(points))
    if producer_error > 2e-14:
        raise RuntimeError("The controlling producer disagrees with the proved symbol")

    # This sample is a falsification check, not evidence of a global bound.
    outer = points[np.min(np.cos(points), axis=1) <= 0.25]
    exterior_minimum = min(float(np.linalg.eigvalsh(bulk(k))[0]) for k in outer)
    if exterior_minimum < 1.5 - 2e-13:
        raise RuntimeError("The analytic exterior inequality was numerically falsified")

    count = 192
    j = np.arange(count)
    z = 1 - 2*(j+.5)/count
    phi = np.pi*(3-np.sqrt(5))*j
    directions = np.column_stack((np.sqrt(1-z*z)*np.cos(phi), np.sqrt(1-z*z)*np.sin(phi), z))
    special = [np.array(k, float) for k in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1),(-1,1,1),(1,1,1e-8)]]
    special += [np.array([1, epsilon, 0.0]) for epsilon in [1e-8,1e-6,1e-4,.01]]
    special = np.array([k/np.linalg.norm(k) for k in special])
    directions = np.vstack((directions, special))
    R = np.arccos(.25)
    rows = []
    for n in directions:
        upper = R / np.max(np.abs(n))
        for branch in range(3):
            r = brentq(lambda radius: np.linalg.eigvalsh(bulk(radius*n))[branch]-ENERGY,
                       0, upper, xtol=2e-14, rtol=2e-14)
            k = r*n
            values, vectors = np.linalg.eigh(bulk(k))
            derivative = radial_derivative(k)
            loewner_minimum = float(np.linalg.eigvalsh(r*derivative-float(beta_arb.mid())*bulk(k))[0])
            slope = float(vectors[:,branch] @ derivative @ vectors[:,branch])
            rows.append(dict(direction=n.tolist(),branch=branch,radius=r,
                             root_residual=float(values[branch]-ENERGY),
                             radial_slope=slope, highest_branch_gap=float(values[2]-values[1]),
                             radial_Loewner_minimum=loewner_minimum,
                             maximum_component=float(np.max(np.abs(k)))))
    if max(abs(row['root_residual']) for row in rows) > 1e-11:
        raise RuntimeError("Root solve failed its residual check")
    if min(row['radial_slope'] for row in rows) < float(slope_arb.lower())-1e-11:
        raise RuntimeError("The analytic radial-slope bound was numerically falsified")
    if min(row['radial_Loewner_minimum'] for row in rows) < -1e-11:
        raise RuntimeError("The analytic radial matrix inequality was numerically falsified")
    if min(row['highest_branch_gap'] for row in rows) < 18/125-1e-11:
        raise RuntimeError("The analytic highest-branch gap was numerically falsified")
    result = dict(
        schema="paper1.radial-shell-feasibility.v1", status="passed",
        scope="Analytic shell coverage/regularity lemma plus exact algebra and bounded numerical falsification checks; no Green-function error enclosure",
        parameters=dict(lattice_spacing=1,bulk_mass=1,axial_stiffness=1,face_diagonal_stiffness=1,energy="36/25"),
        exact_polynomial_identities=exact,
        outward_constants=dict(precision_bits=160,component_radius=radius_arb.str(35),radial_Loewner_factor=beta_arb.str(35),radial_slope_lower_bound=slope_arb.str(35),exterior_spectral_gap="3/50",highest_branch_gap="18/125"),
        numerical_checks=dict(identity_points=len(points),maximum_bond_rank_identity_error=identity_error,
                              maximum_accepted_producer_rank_identity_error=producer_error,
                              exterior_points=len(outer),sampled_exterior_minimum=exterior_minimum,
                              directions=len(directions),branch_roots=len(rows),
                              maximum_root_residual=max(abs(r['root_residual']) for r in rows),
                              minimum_sampled_radial_slope=min(r['radial_slope'] for r in rows),
                              minimum_sampled_highest_branch_gap=min(r['highest_branch_gap'] for r in rows),
                              minimum_sampled_radial_Loewner_eigenvalue=min(r['radial_Loewner_minimum'] for r in rows)),
        unresolved=["Outward angular/radial integration errors", "Degenerate cluster residues in a complete implementation", "Exact-to-selected Green and qT norm errors"],
        exact_lattice_certificate=False,
        source_sha256={"check_radial_shell.py":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "accepted_bulk_source":hashlib.sha256(source.read_bytes()).hexdigest()},
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,python_flint=flint.__version__),
        rows=rows,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k != "rows"},indent=2))


if __name__ == "__main__":
    main()
