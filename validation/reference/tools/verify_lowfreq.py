"""Check the low-frequency normalization and static force identity in Paper 1, Section V.C.

Run with Python 3 and NumPy. The elementwise midpoint quadrature of the static
bulk is a finite numerical fixture. The analytic results for the infinite
lattice and uniform profile class are proved in SM7A.
"""
import itertools
import json
from math import pi, sqrt

import numpy as np


def check_close(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-12)


def bond_geometry():
    axes = np.eye(3)
    bonds = [*axes]
    for i, j in itertools.combinations(range(3), 2):
        bonds.extend([axes[i] + axes[j], axes[i] - axes[j]])
    bonds = np.asarray(bonds)
    unit = bonds / np.linalg.norm(bonds, axis=1)[:, None]
    return bonds, np.einsum("bi,bj->bij", unit, unit)


def static_compliance(n_segments, grid_order):
    """Bracket-free static cubature on the actual reference lattice and centered core."""
    k1 = -pi + (np.arange(grid_order) + 0.5) * 2 * pi / grid_order
    wavevectors = np.array(np.meshgrid(k1, k1, k1, indexing="ij")).reshape(3, -1).T
    bonds, dyads = bond_geometry()
    dynamical = np.einsum("kb,bij->kij", 2 * (1 - np.cos(wavevectors @ bonds.T)), dyads)
    form = np.zeros((len(wavevectors), 3, 2*n_segments), dtype=complex)
    # Fourier sum all four actual square corners for unit axial and radial patterns.
    for s in range(n_segments):
        for sy, sz in itertools.product((-1, 1), repeat=2):
            phase = np.exp(-1j * (wavevectors @ [s + 0.5, sy/2, sz/2]))
            form[:, 0, 2*s] += phase/2
            form[:, 1, 2*s+1] += phase * sy/sqrt(8)
            form[:, 2, 2*s+1] += phase * sz/sqrt(8)
    static = np.einsum("kia,kib->ab", form.conj(), np.linalg.solve(dynamical, form)) / len(form)
    check_close(static.imag, 0)
    check_close(static.real, static.real.T)
    assert np.linalg.eigvalsh(static.real)[0] > 0
    return static.real


def dressed_static(weights, bulk):
    """Assemble the Hessian from squared extensions before eliminating the core."""
    n = len(weights)
    hess = np.zeros((3*n + 1, 3*n + 1))
    for s, g in enumerate(weights):
        strains = np.zeros((2, 3*n+1))
        strains[0, [s, s+1, n+1+2*s]] = [-1, -1, 1]
        strains[1, [s, s+1, n+2+2*s]] = [-1, 1, sqrt(2)]
        hess += (2*g/3) * strains.T @ strains
    line = np.diff(np.eye(n+1), axis=0)
    stiffness = hess[:n+1, :n+1]
    mixed = hess[n+1:, :n+1]
    core = hess[n+1:, n+1:]
    reduced = 1.69 * line.T @ line + stiffness
    reduced -= mixed.T @ np.linalg.solve(np.linalg.inv(bulk) + core, mixed)
    force = np.linalg.solve(np.eye(2*n) + core @ bulk, mixed)
    return reduced, force


def check_static():
    n = 6
    vertices = np.eye(n)
    profiles = [*vertices, [0.1, 0.2, 0.2, 0.2, 0.2, 0.1],
                np.arange(1, n+1)/sum(range(1, n+1)), [0.5, 0, 0, 0, 0, 0.5]]
    endpoint = np.eye(n+1)[0]
    axial_resultant = np.tile([2, 0], n)
    case_count = 0
    for order in (10, 14):
        bulk = static_compliance(n, order)
        vertex_line = [dressed_static(w, bulk)[0] for w in vertices]
        vertex_gap = min(np.linalg.eigvalsh(line)[0] for line in vertex_line)
        assert vertex_gap > 0
        for weights in profiles:
            line, load = dressed_static(weights, bulk)
            normalized_force = -load @ np.linalg.solve(line, endpoint)
            check_close(axial_resultant @ normalized_force, 1)
            check_close(axial_resultant @ load, -line.sum(axis=0))
            mixture = np.tensordot(weights, vertex_line, axes=1)
            assert np.linalg.eigvalsh(line - mixture)[0] >= -2e-12
            assert np.linalg.eigvalsh(line)[0] >= vertex_gap - 2e-12
            if np.count_nonzero(weights) == 1:
                x = np.arange(n) + 0.5
                mu_x = 2 * x @ normalized_force[::2]
                nu_r = normalized_force[1::2].sum() / sqrt(2)
                nu_xr = x @ normalized_force[1::2] / sqrt(2)
                check_close(nu_xr - mu_x * nu_r, 0)
                assert -2/21 - nu_r**2/3 < 0
            case_count += 1
        # With every motif absent the uniform line translation survives.
        disconnected, _ = dressed_static(np.zeros(n), bulk)
        check_close(disconnected @ np.ones(n+1), 0)
    return case_count


def check_radiation():
    # Calculate the acoustic symbol and its polarization projectors from bonds.
    bonds, dyads = bond_geometry()
    mu, weights = np.polynomial.legendre.leggauss(6)
    angles = 2*pi*np.arange(8)/8
    angular = np.zeros(2)
    for x, w in zip(mu, weights):
        for phi in angles:
            normal = np.array([x, sqrt(1-x*x)*np.cos(phi), sqrt(1-x*x)*np.sin(phi)])
            symbol = np.einsum("b,bij->ij", (bonds @ normal)**2, dyads)
            eigenvalues, modes = np.linalg.eigh(symbol)
            check_close(eigenvalues, [1, 1, 3])
            check_close(modes[:, 2:] @ modes[:, 2:].T, np.outer(normal, normal))
            angular += (2*pi/len(angles)) * w * np.array([sum(modes[0, :2]**2), modes[0, 2]**2])
    check_close(angular, [8*pi/3, 4*pi/3])
    impedance = sqrt(1.69)
    speeds = np.array([1, sqrt(3)])
    # Positive shell spectral weights per frequency in reference units; mass density = 1.
    spectral = pi/(2*pi)**3 * angular/(2*speeds**3)
    coefficients = 4*impedance*spectral
    analytic = np.array([2*impedance/(3*pi), impedance/(3*pi*sqrt(3)**3)])
    check_close(coefficients, analytic)
    # Independent physical work / incident power calculation at leading order.
    omega, amplitude = 0.013, 2.75
    force = 2*impedance*omega*amplitude
    incident = impedance*omega**2*amplitude**2/2
    power = omega*force**2/2 * (omega*spectral)
    check_close(power/(incident*omega**2), analytic)
    eta = coefficients[0]/sum(coefficients)
    return [float(v) for v in coefficients], float(eta)


def main():
    cases = check_static()
    coefficients, share = check_radiation()
    print(json.dumps({
        "status": "PASS",
        "static_profile_grid_cases": cases,
        "static_grid_orders": [10, 14],
        "checked_static_quantities": ["unit axial resultant", "translation identity",
                                     "vertex lower bound", "single-motif force moments",
                                     "uniform line zero mode when disconnected"],
        "C_T": round(coefficients[0], 10), "C_L": round(coefficients[1], 10),
        "emitted_qT_share": round(share, 10),
        "scope": "Numerical static-identity and acoustic-normalization checks. No finite-frequency error or interval certificate."
    }, indent=2))


if __name__ == "__main__":
    main()
