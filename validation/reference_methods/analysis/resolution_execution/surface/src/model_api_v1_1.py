#!/usr/bin/env python3
"""Canonical minimal model for a 1-D compression spine coupled to a 3-D vector lattice.

Created from scratch in this phase. This is not a recovered or previously existing
'all-axis' package.

Geometry
--------
* Spine nodes: X_i = (i a, 0, 0), scalar axial displacement q_i.
* Bulk nodes: R_{i,j,k} = ((i+1/2)a, (j+1/2)a, (k+1/2)a),
  vector displacement u_{i,j,k}.
* Bulk central springs: axial nearest neighbours plus face diagonals.
* Each spine segment (i,i+1) couples by eight diagonal springs to the four
  nearest bulk nodes at its midpoint: each of the four nodes is connected to
  both segment endpoints.

For k_axis=k_diag the long-wavelength bulk continuum is isotropic with
Lamé coefficients lambda=mu=k_diag (a=m=1 convention), hence
c_T^2=k_diag/m and c_L^2=3*k_diag/m.

The line-to-bulk coupling is reflection symmetric. A positive axial segment
extension drives the four surrounding bulk nodes inward. In the continuum
limit the source couples to the SV-like shear sector but not, at leading order,
to the azimuthal SH sector.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Iterable
import itertools
import math
import numpy as np


@dataclass(frozen=True)
class ModelParameters:
    a: float = 1.0
    m_bulk: float = 1.0
    m_spine: float = 1.0
    k_axis: float = 1.0
    k_diag: float = 1.0
    k_spine: float = 1.69  # c_spine=1.3 for m_spine=1
    k_couple: float = 0.20

    def validate(self) -> None:
        for name, value in asdict(self).items():
            if name == "a":
                if not value > 0:
                    raise ValueError("lattice spacing must be positive")
            elif name == "k_couple":
                if value < 0:
                    raise ValueError("k_couple must be nonnegative")
            elif not value > 0:
                raise ValueError(f"{name} must be positive")

    @property
    def c_spine_longwave(self) -> float:
        return self.a * math.sqrt(self.k_spine / self.m_spine)

    @property
    def c_bulk_T_longwave(self) -> float:
        # exact for k_axis=k_diag; otherwise this returns the [100] T speed
        return self.a * math.sqrt(self.k_diag / self.m_bulk)

    @property
    def c_bulk_L_longwave(self) -> float:
        return self.a * math.sqrt((self.k_axis + 2.0 * self.k_diag) / self.m_bulk)


AXIAL_BONDS: tuple[tuple[int, int, int], ...] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)
FACE_DIAGONAL_BONDS: tuple[tuple[int, int, int], ...] = (
    (1, 1, 0),
    (1, -1, 0),
    (1, 0, 1),
    (1, 0, -1),
    (0, 1, 1),
    (0, 1, -1),
)


def _unit(d: Iterable[float]) -> np.ndarray:
    v = np.asarray(tuple(d), dtype=float)
    return v / np.linalg.norm(v)


def bulk_dynamical_matrix(k: np.ndarray, p: ModelParameters) -> np.ndarray:
    """Return the 3x3 bulk omega^2 matrix at dimensionless Bloch vector k*a.

    The returned matrix is mass-normalized. k is dimensionless in [-pi,pi]^3.
    """
    p.validate()
    k = np.asarray(k, dtype=float)
    if k.shape != (3,):
        raise ValueError("k must have shape (3,)")
    D = np.zeros((3, 3), dtype=float)
    for d in AXIAL_BONDS:
        n = _unit(d)
        phase = float(np.dot(k, d))
        D += 2.0 * p.k_axis * (1.0 - math.cos(phase)) * np.outer(n, n)
    for d in FACE_DIAGONAL_BONDS:
        n = _unit(d)
        phase = float(np.dot(k, d))
        D += 2.0 * p.k_diag * (1.0 - math.cos(phase)) * np.outer(n, n)
    return D / p.m_bulk


def bulk_modes(k: np.ndarray, p: ModelParameters) -> tuple[np.ndarray, np.ndarray]:
    """Return ascending frequencies and column eigenvectors."""
    vals, vecs = np.linalg.eigh(bulk_dynamical_matrix(k, p))
    vals = np.maximum(vals, 0.0)
    return np.sqrt(vals), vecs


def spine_frequency(kx: float, p: ModelParameters) -> float:
    """Uncoupled 1-D spine-chain frequency."""
    return 2.0 * math.sqrt(p.k_spine / p.m_spine) * abs(math.sin(0.5 * kx))


def coupling_vector(k: np.ndarray, p: ModelParameters, normalized_plane_wave: bool = False,
                    transverse_count: int | None = None) -> np.ndarray:
    """Bare line-to-bulk coupling vector for one bulk plane wave.

    The phase convention uses bulk nodes at half-integer coordinates and spine
    nodes at integer x. Overall phase is removed. The returned real vector is

        (8 g / 3) (-c_x c_y c_z,
                    s_x s_y c_z,
                    s_x c_y s_z).

    If normalized_plane_wave=True, divide by sqrt(N_y*N_z), requiring
    transverse_count=N_y*N_z.
    """
    k = np.asarray(k, dtype=float)
    if k.shape != (3,):
        raise ValueError("k must have shape (3,)")
    cx, cy, cz = np.cos(0.5 * k)
    sx, sy, sz = np.sin(0.5 * k)
    v = (8.0 * p.k_couple / 3.0) * np.array(
        [-cx * cy * cz, sx * sy * cz, sx * cy * sz], dtype=float
    )
    if normalized_plane_wave:
        if transverse_count is None or transverse_count <= 0:
            raise ValueError("positive transverse_count required")
        v /= math.sqrt(transverse_count)
    return v


def source_pattern_realspace(kx: float, p: ModelParameters) -> dict[tuple[int, int], np.ndarray]:
    """K_bq entries on the four transverse nodes around one spine segment.

    Indices (-1,-1),(-1,0),(0,-1),(0,0) correspond to positions
    (sigma_y/2,sigma_z/2), sigma_y,sigma_z in {-1,+1}.
    A global complex phase has been removed.
    """
    c = math.cos(0.5 * kx)
    s = math.sin(0.5 * kx)
    out: dict[tuple[int, int], np.ndarray] = {}
    for jy, sy in [(-1, -1.0), (0, 1.0)]:
        for jz, sz in [(-1, -1.0), (0, 1.0)]:
            # complex Hermitian K_bq convention before phase choice
            out[(jy, jz)] = (2.0 * p.k_couple / 3.0) * np.array(
                [-c, 1j * s * sy, 1j * s * sz], dtype=complex
            )
    return out


def fourier_transform_source_pattern(k: np.ndarray, p: ModelParameters) -> np.ndarray:
    """Direct DFT of source_pattern_realspace, for verification."""
    kx, ky, kz = map(float, k)
    total = np.zeros(3, dtype=complex)
    for (jy, jz), v in source_pattern_realspace(kx, p).items():
        y = jy + 0.5
        z = jz + 0.5
        total += np.exp(-1j * (ky * y + kz * z)) * v
    return total


def polarization_basis(k: np.ndarray) -> np.ndarray:
    """Return columns [L, SV, SH] based on the spine axis x.

    L is parallel to k. SH is azimuthal about x, and SV completes the
    right-handed basis. At k_perp=0 an arbitrary y/z basis is used.
    """
    k = np.asarray(k, dtype=float)
    norm = np.linalg.norm(k)
    if norm == 0.0:
        return np.eye(3)
    eL = k / norm
    ex = np.array([1.0, 0.0, 0.0])
    eSH = np.cross(ex, eL)
    n = np.linalg.norm(eSH)
    if n < 1e-14:
        eSH = np.array([0.0, 0.0, 1.0])
    else:
        eSH /= n
    eSV = np.cross(eSH, eL)
    eSV /= np.linalg.norm(eSV)
    return np.column_stack([eL, eSV, eSH])


def assign_modes_L_SV_SH(k: np.ndarray, vecs: np.ndarray) -> tuple[int, int, int]:
    """Assign numerical eigenvectors to L/SV/SH by maximum total overlap."""
    basis = polarization_basis(k)
    overlap = np.abs(basis.T @ vecs) ** 2
    best_score = -1.0
    best_perm: tuple[int, int, int] | None = None
    for perm in itertools.permutations(range(3)):
        score = sum(float(overlap[a, perm[a]]) for a in range(3))
        if score > best_score:
            best_score = score
            best_perm = tuple(int(x) for x in perm)
    assert best_perm is not None
    return best_perm


def branch_couplings(k: np.ndarray, p: ModelParameters) -> dict[str, float]:
    """Squared bare coupling to the three bulk branches."""
    omega, vecs = bulk_modes(k, p)
    indices = assign_modes_L_SV_SH(k, vecs)
    V = coupling_vector(k, p)
    out: dict[str, float] = {}
    for name, idx in zip(("L", "SV", "SH"), indices):
        out[name] = float(abs(np.vdot(vecs[:, idx], V)) ** 2)
        out[f"omega_{name}"] = float(omega[idx])
    return out


def retarded_lorentzian_frequency(omega: float, omega0: np.ndarray, eta: float) -> np.ndarray:
    """Positive-frequency approximation to a retarded spectral delta."""
    return eta / np.pi / ((omega0 - omega) ** 2 + eta**2)


def branch_spectral_density(kx: float, p: ModelParameters, *, n_perp: int = 161,
                            eta: float = 0.012, omega: float | None = None) -> dict[str, float]:
    """Compute bare branch-resolved line-to-bulk spectral density at fixed kx.

    The integral is a uniform Brillouin-zone quadrature over ky,kz. The
    normalization is arbitrary but common to all branches; fractions are
    invariant under the chosen normalization.
    """
    if n_perp < 16:
        raise ValueError("n_perp too small")
    if omega is None:
        omega = spine_frequency(kx, p)
    vals = {"L": 0.0, "SV": 0.0, "SH": 0.0}
    grid = np.linspace(-math.pi, math.pi, n_perp, endpoint=False)
    for ky in grid:
        for kz in grid:
            k = np.array([kx, ky, kz])
            om, vecs = bulk_modes(k, p)
            indices = assign_modes_L_SV_SH(k, vecs)
            V = coupling_vector(k, p)
            weights = retarded_lorentzian_frequency(float(omega), om, eta)
            for name, idx in zip(("L", "SV", "SH"), indices):
                vals[name] += float(abs(np.vdot(vecs[:, idx], V)) ** 2 * weights[idx])
    total = sum(vals.values())
    out = {
        "kx": float(kx),
        "omega": float(omega),
        "eta": float(eta),
        "n_perp": int(n_perp),
        **{f"rho_{key}": float(value) for key, value in vals.items()},
        "rho_total": float(total),
    }
    if total > 0.0:
        out.update({f"fraction_{key}": float(value / total) for key, value in vals.items()})
    return out


def continuum_phase_matching(kx: float, p: ModelParameters) -> dict[str, float | bool]:
    """Long-wavelength continuum phase matching for the bare spine frequency.

    For omega=c_q |kx|, a bulk branch of speed c has transverse wavenumber
    k_perp^2=(c_q^2/c^2-1) kx^2. Positive values are propagating.
    """
    cq = p.c_spine_longwave
    cT = p.c_bulk_T_longwave
    cL = p.c_bulk_L_longwave
    return {
        "c_spine": cq,
        "c_T": cT,
        "c_L": cL,
        "kperp2_T": (cq * cq / (cT * cT) - 1.0) * kx * kx,
        "kperp2_L": (cq * cq / (cL * cL) - 1.0) * kx * kx,
        "T_propagating": cq > cT,
        "L_propagating": cq > cL,
    }


def small_k_continuum_matrix(k: np.ndarray, p: ModelParameters) -> np.ndarray:
    """Quadratic long-wave matrix derived analytically from the bonds."""
    kx, ky, kz = map(float, k)
    ka, kd = p.k_axis, p.k_diag
    D = np.array([
        [(ka + 2*kd)*kx*kx + kd*(ky*ky+kz*kz), 2*kd*kx*ky, 2*kd*kx*kz],
        [2*kd*kx*ky, (ka + 2*kd)*ky*ky + kd*(kx*kx+kz*kz), 2*kd*ky*kz],
        [2*kd*kx*kz, 2*kd*ky*kz, (ka + 2*kd)*kz*kz + kd*(kx*kx+ky*ky)],
    ], dtype=float)
    return (p.a**2 / p.m_bulk) * D
