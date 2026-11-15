#!/usr/bin/env python3
"""Independent linear-tetrahedron resolvent and line-to-bulk scattering audit.

This module deliberately does not import the Phase-12 reduced-core, spectral-
histogram, on-shell marching-cubes, or scattering modules.  It reimplements:

* the 3-D central-spring bulk dynamical matrix from the bond geometry;
* the two source-reachable cross-sectional form factors;
* a zero-broadening linear-tetrahedron Brillouin-zone resolvent;
* the finite-slab static coupling matrices;
* the scalar-chain surface self-energy and Caroli/Fisher-Lee scattering.

Inside each tetrahedron the three sorted bulk eigenvalues are linearly
interpolated and the matrix spectral weight is evaluated at the tetrahedron
centroid.  The singular integral is analytic through a third divided
difference of f(x)=1/2 (x-z)^2 log(x-z), for which f'''(x)=1/(x-z).
The method is independent of Lorentzian broadening and spectral histograms.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import argparse
import cmath
import json
import math
import time
from typing import Iterable

import numpy as np
from scipy import linalg
from stable_scalar import stable_tetrahedron_resolvent
from centroid_weights import invariant_centroid_weights


@dataclass(frozen=True)
class Parameters:
    a: float = 1.0
    m_bulk: float = 1.0
    m_spine: float = 1.0
    k_axis: float = 1.0
    k_diag: float = 1.0
    k_spine: float = 1.69

    def validate(self) -> None:
        for key, value in asdict(self).items():
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{key} must be finite and positive")


# Positive bond representatives.  The six face diagonals include both signs
# required to reconstruct the cubic central-spring lattice without double
# counting in Fourier space.
AXIAL = np.array(((1, 0, 0), (0, 1, 0), (0, 0, 1)), dtype=float)
DIAGONAL = np.array(
    ((1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)),
    dtype=float,
)
BONDS = np.concatenate((AXIAL, DIAGONAL), axis=0)
UNITS = BONDS / np.linalg.norm(BONDS, axis=1)[:, None]

# Six tetrahedra around the 000--111 body diagonal.  Vertex numbering is
# x + 2 y + 4 z.
CUBE_VERTICES = np.array(
    ((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
     (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)),
    dtype=int,
)
TETRA = np.array(
    ((0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7),
     (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7)),
    dtype=int,
)


@dataclass
class TetrahedronResult:
    omega: float
    energy_imaginary_regulator: float
    grid: int
    nseg: int
    G0: np.ndarray
    Im_branches: list[np.ndarray]
    Im_longitudinal_projector: np.ndarray
    Im_branch_longitudinal_projectors: list[np.ndarray]
    symmetry_residual: float
    min_imaginary_density_eigenvalue: float
    seconds: float
    diagnostics: dict


@dataclass
class ScatteringResult:
    omega: float
    Rq_amplitude: float
    Rq_flux: float
    Tq: float
    Tbulk: float
    T_branches: list[float]
    TL_modal: float
    TT_modal: float
    TL_projector: float
    TT_projector: float
    branch_longitudinal_powers: list[float]
    branch_transverse_powers: list[float]
    branch_longitudinal_fractions: list[float]
    closure_amplitude: float
    closure_flux: float
    component_closure_modal: float
    component_closure_projector: float
    min_gamma_bulk_eigenvalue: float
    reciprocity_matrix_residual: float


def bulk_stiffness(k: np.ndarray, p: Parameters) -> np.ndarray:
    """Un-normalized 3x3 stiffness K(k), batched over leading axes."""
    p.validate()
    k = np.asarray(k, dtype=float)
    if k.shape[-1] != 3:
        raise ValueError("k must have last dimension 3")
    out = np.zeros(k.shape[:-1] + (3, 3), dtype=float)
    for b, n in zip(BONDS, UNITS):
        ks = p.k_axis if np.count_nonzero(b) == 1 else p.k_diag
        phase = np.einsum("...a,a->...", k, b)
        out += (2.0 * ks * (1.0 - np.cos(phase)))[..., None, None] * np.outer(n, n)
    return out


def active_form_factor(k: np.ndarray) -> np.ndarray:
    """Independent DFT of the normalized equal-axial/radial source patterns.

    Returns F[..., cartesian, active], with active columns X and R.  The four
    cross-sectional nodes use coordinates (y,z)=(-1,-1),(-1,0),(0,-1),(0,0)
    in the phase convention of the exact active reduction.
    """
    k = np.asarray(k, dtype=float)
    if k.shape[-1] != 3:
        raise ValueError("k must have last dimension 3")
    ky = k[..., 1]
    kz = k[..., 2]
    F = np.zeros(k.shape[:-1] + (3, 2), dtype=np.complex128)
    norm_r = 2.0 * math.sqrt(2.0)
    for y in (-1, 0):
        sy = -1.0 if y == -1 else 1.0
        for z in (-1, 0):
            sz = -1.0 if z == -1 else 1.0
            phase = np.exp(-1j * (ky * y + kz * z))
            F[..., 0, 0] += 0.5 * phase
            F[..., 1, 1] += (sy / norm_r) * phase
            F[..., 2, 1] += (sz / norm_r) * phase
    return F


def _f_derivative(x: np.ndarray, z: complex, order: int) -> np.ndarray:
    y = np.asarray(x, dtype=np.complex128) - z
    if order == 0:
        return 0.5 * y * y * np.log(y)
    if order == 1:
        return y * np.log(y) + 0.5 * y
    if order == 2:
        return np.log(y) + 1.5
    if order == 3:
        return 1.0 / y
    raise ValueError("only derivative orders 0..3 are used")


def tetrahedron_resolvent_divided_difference(
    energies: np.ndarray,
    z_energy: complex,
    *,
    confluent_rtol: float = 2.0e-7,
) -> np.ndarray:
    """Reference-simplex integral of 1/(linear energy-z) by divided differences.

    `energies` has shape (...,4).  The return is the integral over the standard
    tetrahedron of volume 1/6.  Multiplication by the cube Jacobian dk**3 gives
    the physical tetrahedron integral.

    Repeated/near-repeated knots use the appropriate confluent derivative
    f^(r)/r!.  This also handles the three acoustic degeneracies at Gamma.
    """
    e = np.sort(np.asarray(energies, dtype=float), axis=-1)
    if e.shape[-1] != 4:
        raise ValueError("energies must have last dimension four")
    dd = _f_derivative(e, z_energy, 0)
    for order in range(1, 4):
        left = e[..., :-order]
        right = e[..., order:]
        den = right - left
        num = dd[..., 1:] - dd[..., :-1]
        scale = 1.0 + np.maximum(np.abs(left), np.abs(right))
        regular = np.abs(den) > confluent_rtol * scale
        nxt = np.empty_like(num)
        nxt[regular] = num[regular] / den[regular]
        if np.any(~regular):
            midpoint = 0.5 * (left[~regular] + right[~regular])
            nxt[~regular] = _f_derivative(midpoint, z_energy, order) / math.factorial(order)
        dd = nxt
    return dd[..., 0]


# The default CLI and imported integrator both use the stable scalar evaluation.
tetrahedron_resolvent_divided_difference = stable_tetrahedron_resolvent


def _assemble_symmetric_toeplitz(positive_blocks: np.ndarray, nseg: int) -> np.ndarray:
    """Assemble the reciprocal complex-symmetric retarded matrix."""
    out = np.empty((2 * nseg, 2 * nseg), dtype=np.complex128)
    for s in range(nseg):
        for t in range(nseg):
            dx = s - t
            block = positive_blocks[dx]  # Both signs were explicitly integrated.
            out[2 * s:2 * s + 2, 2 * t:2 * t + 2] = block
    return out


def integrate_tetrahedron_resolvent(
    omega: float,
    p: Parameters,
    *,
    nseg: int = 6,
    grid: int = 32,
    energy_imaginary_regulator: float = 1.0e-10,
    cube_chunk: int = 4096,
    confluent_rtol: float = 2.0e-7,
    grid_shift: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> TetrahedronResult:
    """Compute the projected retarded bulk Green matrix with no broadening.

    The tiny imaginary number selects the retarded branch of the analytic log;
    it is not a Lorentzian quadrature broadening.  Results should be insensitive
    to it once it is well below the tetrahedron energy resolution.
    """
    p.validate()
    if omega <= 0.0 or grid < 6 or nseg <= 0:
        raise ValueError("invalid omega/grid/nseg")
    if energy_imaginary_regulator <= 0.0:
        raise ValueError("regulator must be positive")
    start = time.time()
    E0 = float(omega * omega)
    zE = complex(E0, energy_imaginary_regulator)
    dk = 2.0 * math.pi / grid
    shift = np.asarray(grid_shift, dtype=float)
    if shift.shape != (3,) or np.any(~np.isfinite(shift)):
        raise ValueError("grid_shift must be a finite 3-vector")
    axes = [-math.pi + dk * (np.arange(grid + 1, dtype=float) + shift[a]) for a in range(3)]

    # Vertex eigenvalues on the periodic cube including the duplicated +pi face.
    KX, KY, KZ = np.meshgrid(axes[0], axes[1], axes[2], indexing="ij")
    k_vertices = np.stack((KX, KY, KZ), axis=-1)
    eigenvalues = np.linalg.eigvalsh(bulk_stiffness(k_vertices, p) / p.m_bulk)
    del KX, KY, KZ, k_vertices

    # Cubes are processed as flattened base indices.
    ii, jj, kk = np.meshgrid(
        np.arange(grid, dtype=int),
        np.arange(grid, dtype=int),
        np.arange(grid, dtype=int),
        indexing="ij",
    )
    bases = np.stack((ii.ravel(), jj.ravel(), kk.ravel()), axis=-1)
    del ii, jj, kk

    positive = np.zeros((2 * nseg - 1, 2, 2), dtype=np.complex128)
    positive_branches = np.zeros((3, 2 * nseg - 1, 2, 2), dtype=np.complex128)
    positive_long = np.zeros((2 * nseg - 1, 2, 2), dtype=np.complex128)
    positive_branch_long = np.zeros((3, 2 * nseg - 1, 2, 2), dtype=np.complex128)
    max_abs_factor = 0.0
    nonfinite_count = 0
    centroid_cluster_counts = {}

    # Physical tetrahedron integral carries dk^3; BZ normalization is (2pi)^-3.
    normalization = (dk ** 3) / (((2.0 * math.pi) ** 3) * p.m_bulk)

    for c0 in range(0, len(bases), cube_chunk):
        base = bases[c0:c0 + cube_chunk]
        base_k = -math.pi + dk * (base.astype(float) + shift[None, :])
        for tet in TETRA:
            offsets = CUBE_VERTICES[tet]  # (4,3)
            # Vertex sorted band energies, shape (ncube,4,3).
            idx = base[:, None, :] + offsets[None, :, :]
            ev = eigenvalues[idx[..., 0], idx[..., 1], idx[..., 2], :]
            centroid_k = base_k + dk * offsets.mean(axis=0)
            vals_c, vec_c = np.linalg.eigh(bulk_stiffness(centroid_k, p) / p.m_bulk)
            F = active_form_factor(centroid_k)
            coeff = np.einsum("...am,...aA->...mA", vec_c.conj(), F, optimize=True)
            principal_k = (centroid_k + math.pi) % (2.0 * math.pi) - math.pi
            knorm = np.linalg.norm(principal_k, axis=-1)
            khat = np.zeros_like(principal_k)
            mask = knorm > 1.0e-14
            khat[mask] = principal_k[mask] / knorm[mask, None]
            longitudinal_fraction = np.abs(
                np.einsum("...a,...am->...m", khat, vec_c, optimize=True)
            ) ** 2

            weights, geometric_weights, tied_counts = invariant_centroid_weights(vals_c, vec_c, F, khat)
            for key, count in tied_counts.items():
                centroid_cluster_counts[key] = centroid_cluster_counts.get(key, 0) + count

            for branch in range(3):
                factor = normalization * tetrahedron_resolvent_divided_difference(
                    ev[..., branch], zE, confluent_rtol=confluent_rtol
                )
                max_abs_factor = max(max_abs_factor, float(np.max(np.abs(factor))))
                finite = np.isfinite(factor.real) & np.isfinite(factor.imag)
                nonfinite_count += int(np.size(finite) - np.count_nonzero(finite))
                if not np.all(finite):
                    factor = np.where(finite, factor, 0.0)
                c = coeff[..., branch, :]
                W = weights[:, branch]
                WL = geometric_weights[:, branch]
                lf = longitudinal_fraction[..., branch]
                for dx in range(-(nseg - 1), nseg):
                    phase = np.exp(1j * centroid_k[:, 0] * dx)
                    block = np.einsum("q,q,qAB->AB", factor, phase, W, optimize=True)
                    positive[dx] += block
                    positive_branches[branch, dx] += block
                    long_block = np.einsum(
                        "q,q,qAB->AB", factor, phase, WL, optimize=True
                    )
                    positive_long[dx] += long_block
                    positive_branch_long[branch, dx] += long_block

    Graw = _assemble_symmetric_toeplitz(positive, nseg)
    branch_raw = [_assemble_symmetric_toeplitz(positive_branches[m], nseg) for m in range(3)]
    long_raw = _assemble_symmetric_toeplitz(positive_long, nseg)
    branch_long_raw = [_assemble_symmetric_toeplitz(positive_branch_long[m], nseg) for m in range(3)]

    # Reciprocity is complex symmetry, not Hermiticity, for a retarded Green matrix.
    symmetry_residual = float(
        np.linalg.norm(Graw - Graw.T) / max(np.linalg.norm(Graw), 1.0e-30)
    )
    G = 0.5 * (Graw + Graw.T)
    branch_mats = [0.5 * (x + x.T) for x in branch_raw]
    long_mat = 0.5 * (long_raw + long_raw.T)
    branch_long_mats = [0.5 * (x + x.T) for x in branch_long_raw]

    Im_branches = [np.imag(x) for x in branch_mats]
    Im_long = np.imag(long_mat)
    Im_branch_long = [np.imag(x) for x in branch_long_mats]
    Im_total = np.imag(G)
    min_eig = float(np.min(np.linalg.eigvalsh(0.5 * (Im_total + Im_total.T))))
    diagnostics = {
        "normalization_per_reference_tetra": normalization,
        "grid_shift": shift.tolist(),
        "centroid_cluster_counts": centroid_cluster_counts,
        "max_abs_scalar_tetra_factor": max_abs_factor,
        "nonfinite_factor_count": nonfinite_count,
        "branch_sum_relative_error": float(
            np.linalg.norm(sum(branch_mats) - G) / max(np.linalg.norm(G), 1.0e-30)
        ),
        "imag_branch_sum_relative_error": float(
            np.linalg.norm(sum(Im_branches) - Im_total) / max(np.linalg.norm(Im_total), 1.0e-30)
        ),
        "raw_G_reciprocity_residual": symmetry_residual,
        "raw_spectral_min_eigenvalue": float(np.linalg.eigvalsh((Graw-Graw.conj().T)/(2j)).min()),
        "raw_branch_spectral_min_eigenvalues": [float(np.linalg.eigvalsh((x-x.conj().T)/(2j)).min()) for x in branch_raw],
        "real_norm": float(np.linalg.norm(np.real(G))),
        "imag_norm": float(np.linalg.norm(Im_total)),
        "trace_real": float(np.trace(np.real(G))),
        "trace_imag": float(np.trace(Im_total)),
    }
    return TetrahedronResult(
        omega=float(omega),
        energy_imaginary_regulator=float(energy_imaginary_regulator),
        grid=int(grid),
        nseg=int(nseg),
        G0=G,
        Im_branches=Im_branches,
        Im_longitudinal_projector=Im_long,
        Im_branch_longitudinal_projectors=Im_branch_long,
        symmetry_residual=symmetry_residual,
        min_imaginary_density_eigenvalue=min_eig,
        seconds=float(time.time() - start),
        diagnostics=diagnostics,
    )


def static_line_core_matrices(g_segments: Iterable[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Independent expansion of the paired diagonal-spring energy."""
    g = np.asarray(tuple(g_segments), dtype=float)
    if g.ndim != 1 or len(g) == 0 or np.any(g < 0.0):
        raise ValueError("invalid coupling profile")
    n = len(g)
    U = np.zeros((2 * n, 2 * n), dtype=float)
    V = np.zeros((2 * n, n + 1), dtype=float)
    Kself = np.zeros((n + 1, n + 1), dtype=float)
    for s, gs in enumerate(g):
        # Four corner nodes, each connected to both segment endpoints.
        U[2 * s, 2 * s] += 2.0 * gs / 3.0
        U[2 * s + 1, 2 * s + 1] += 4.0 * gs / 3.0
        V[2 * s, s] += -2.0 * gs / 3.0
        V[2 * s, s + 1] += -2.0 * gs / 3.0
        radial = 2.0 * math.sqrt(2.0) * gs / 3.0
        V[2 * s + 1, s] += -radial
        V[2 * s + 1, s + 1] += radial
        Kself[s, s] += 4.0 * gs / 3.0
        Kself[s + 1, s + 1] += 4.0 * gs / 3.0
    return U, V, Kself


def surface_self_energy_scalar_chain(omega: float, p: Parameters, eta: float = 1.0e-12) -> complex:
    """Retarded surface self-energy of the semi-infinite scalar chain."""
    z = complex(omega, eta)
    onsite_dynamic = p.m_spine * z * z - 2.0 * p.k_spine
    discriminant = cmath.sqrt(onsite_dynamic * onsite_dynamic - 4.0 * p.k_spine ** 2)
    roots = (0.5 * (onsite_dynamic - discriminant), 0.5 * (onsite_dynamic + discriminant))
    admissible = [r for r in roots if r.imag <= 1.0e-14]
    if admissible:
        return min(admissible, key=abs)
    return min(roots, key=lambda r: r.imag)


def independent_scatter(
    omega: float,
    g_segments: Iterable[float],
    p: Parameters,
    tetra: TetrahedronResult,
    *,
    eta_lead: float = 1.0e-12,
) -> ScatteringResult:
    """Independent NEGF scattering and branch partition."""
    g = np.asarray(tuple(g_segments), dtype=float)
    n = len(g)
    if tetra.nseg != n:
        raise ValueError("tetrahedron result nseg must match coupling segments")
    U, V, Kself = static_line_core_matrices(g)
    G0 = tetra.G0
    # Dyson dressing of the active bulk coordinates by their local static block.
    Gc = linalg.solve(np.eye(2 * n, dtype=complex) + G0 @ U, G0, assume_a="gen")
    Sigma_bulk = -(V.T @ Gc @ V)
    Gamma_bulk = 1j * (Sigma_bulk - Sigma_bulk.conj().T)

    nq = n + 1
    Kq = np.zeros((nq, nq), dtype=float)
    np.fill_diagonal(Kq, 2.0 * p.k_spine)
    for i in range(nq - 1):
        Kq[i, i + 1] = Kq[i + 1, i] = -p.k_spine
    lead_sigma = surface_self_energy_scalar_chain(omega, p, eta_lead)
    SigmaL = np.zeros((nq, nq), dtype=complex)
    SigmaR = np.zeros((nq, nq), dtype=complex)
    SigmaL[0, 0] = lead_sigma
    SigmaR[-1, -1] = lead_sigma
    GammaL = 1j * (SigmaL - SigmaL.conj().T)
    GammaR = 1j * (SigmaR - SigmaR.conj().T)
    z = complex(omega, eta_lead)
    A = (
        p.m_spine * z * z * np.eye(nq)
        - Kq - Kself - Sigma_bulk - SigmaL - SigmaR
    )
    Gq = linalg.inv(A, check_finite=False)
    Gqa = Gq.conj().T
    Tq = float(np.real(np.trace(GammaL @ Gq @ GammaR @ Gqa)))
    Tbulk = float(np.real(np.trace(GammaL @ Gq @ Gamma_bulk @ Gqa)))
    gamma_in = float(np.real(GammaL[0, 0]))
    R_amp = float(abs(1.0 - 1j * gamma_in * Gq[0, 0]) ** 2)
    R_flux = float(1.0 - Tq - Tbulk)

    # Core force generated by a unit-flux left input, then branch spectral power.
    q_column = Gq[:, 0]
    active_force = linalg.solve(np.eye(2 * n) + U @ G0, V @ q_column)

    def branch_power(ImG: np.ndarray) -> float:
        return float(2.0 * gamma_in * np.real(active_force.conj() @ ImG @ active_force))

    branch_T = [branch_power(x) for x in tetra.Im_branches]
    TL_modal = branch_T[2]
    TT_modal = branch_T[0] + branch_T[1]
    branch_long = [branch_power(x) for x in tetra.Im_branch_longitudinal_projectors]
    branch_transverse = [branch_T[i] - branch_long[i] for i in range(3)]
    branch_lf = [branch_long[i] / branch_T[i] if branch_T[i] > 0.0 else 0.0 for i in range(3)]
    TL_proj = branch_power(tetra.Im_longitudinal_projector)
    TT_proj = Tbulk - TL_proj
    return ScatteringResult(
        omega=float(omega),
        Rq_amplitude=R_amp,
        Rq_flux=R_flux,
        Tq=Tq,
        Tbulk=Tbulk,
        T_branches=branch_T,
        TL_modal=TL_modal,
        TT_modal=TT_modal,
        TL_projector=TL_proj,
        TT_projector=TT_proj,
        branch_longitudinal_powers=branch_long,
        branch_transverse_powers=branch_transverse,
        branch_longitudinal_fractions=branch_lf,
        closure_amplitude=float(abs(R_amp + Tq + Tbulk - 1.0)),
        closure_flux=float(abs(R_flux + Tq + Tbulk - 1.0)),
        component_closure_modal=float(abs(TL_modal + TT_modal - Tbulk)),
        component_closure_projector=float(abs(TL_proj + TT_proj - Tbulk)),
        min_gamma_bulk_eigenvalue=float(
            np.min(np.linalg.eigvalsh(0.5 * (Gamma_bulk + Gamma_bulk.conj().T))).real
        ),
        reciprocity_matrix_residual=float(
            np.linalg.norm(Gq - Gq.T) / max(np.linalg.norm(Gq), 1.0e-30)
        ),
    )


def _complex_matrix_to_npz(path: Path, tetra: TetrahedronResult) -> None:
    payload = {
        "G0": tetra.G0,
        "Im_longitudinal_projector": tetra.Im_longitudinal_projector,
        "omega": np.array(tetra.omega),
        "grid": np.array(tetra.grid),
        "regulator": np.array(tetra.energy_imaginary_regulator),
    }
    for i, x in enumerate(tetra.Im_branches):
        payload[f"Im_branch_{i}"] = x
        payload[f"Im_branch_long_{i}"] = tetra.Im_branch_longitudinal_projectors[i]
    np.savez_compressed(path, **payload)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", type=int, default=24)
    parser.add_argument("--omega", type=float, default=1.2)
    parser.add_argument("--nseg", type=int, default=6)
    parser.add_argument("--regulator", type=float, default=1.0e-10)
    parser.add_argument("--chunk", type=int, default=4096)
    parser.add_argument("--shift", nargs=3, type=float, default=(0.0, 0.0, 0.0))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--matrix-output", type=Path)
    args = parser.parse_args()
    p = Parameters()
    g = (0.1, 0.2, 0.2, 0.2, 0.2, 0.1)
    tetra = integrate_tetrahedron_resolvent(
        args.omega, p, nseg=args.nseg, grid=args.grid,
        energy_imaginary_regulator=args.regulator, cube_chunk=args.chunk, grid_shift=tuple(args.shift),
    )
    scatter = independent_scatter(args.omega, g, p, tetra)
    out = {
        "method": "analysis-only scalar/complete-phase/centroid-cluster diagnostic; not convergence certified",
        "tetrahedron": {
            "omega": tetra.omega,
            "grid": tetra.grid,
            "nseg": tetra.nseg,
            "energy_imaginary_regulator": tetra.energy_imaginary_regulator,
            "seconds": tetra.seconds,
            "symmetry_residual": tetra.symmetry_residual,
            "min_imaginary_density_eigenvalue": tetra.min_imaginary_density_eigenvalue,
            "diagnostics": tetra.diagnostics,
        },
        "scattering": asdict(scatter),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    if args.matrix_output:
        args.matrix_output.parent.mkdir(parents=True, exist_ok=True)
        _complex_matrix_to_npz(args.matrix_output, tetra)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
