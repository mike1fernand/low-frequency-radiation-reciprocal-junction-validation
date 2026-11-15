from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from numpy.linalg import norm
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = Path(__file__).resolve().parent / "dependencies"
sys.path.insert(0, str(DEPENDENCIES))

from model_api_v1_1 import (  # noqa: E402
    AXIAL_BONDS,
    FACE_DIAGONAL_BONDS,
    ModelParameters,
    bulk_dynamical_matrix,
)

P = ModelParameters()
OMEGA = 1.2
TARGET = OMEGA**2
PROFILE = np.array([0.1, 0.2, 0.2, 0.2, 0.2, 0.1], dtype=float)
KAPPA = 2.0 * math.asin(OMEGA / (2.0 * math.sqrt(P.k_spine / P.m_spine)))
CORE_TRANS = np.array([[-1, -1], [-1, 0], [0, -1], [0, 0]], dtype=int)
CORE_SIGNS = np.array([[-1.0, -1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, 1.0]], dtype=float)


def d_closed(k: np.ndarray) -> np.ndarray:
    """Independent closed-form mass-normalized bulk dynamical matrix."""
    kx, ky, kz = map(float, k)
    ka, kd, m = P.k_axis, P.k_diag, P.m_bulk
    d = np.zeros((3, 3), dtype=float)
    d[0, 0] = 2 * ka * (1 - np.cos(kx)) + 2 * kd * (
        (1 - np.cos(kx) * np.cos(ky)) + (1 - np.cos(kx) * np.cos(kz))
    )
    d[1, 1] = 2 * ka * (1 - np.cos(ky)) + 2 * kd * (
        (1 - np.cos(kx) * np.cos(ky)) + (1 - np.cos(ky) * np.cos(kz))
    )
    d[2, 2] = 2 * ka * (1 - np.cos(kz)) + 2 * kd * (
        (1 - np.cos(kx) * np.cos(kz)) + (1 - np.cos(ky) * np.cos(kz))
    )
    d[0, 1] = d[1, 0] = 2 * kd * np.sin(kx) * np.sin(ky)
    d[0, 2] = d[2, 0] = 2 * kd * np.sin(kx) * np.sin(kz)
    d[1, 2] = d[2, 1] = 2 * kd * np.sin(ky) * np.sin(kz)
    return d / m


def eig(k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eigh(d_closed(np.asarray(k, dtype=float)))
    return np.maximum(vals, 0.0), vecs


def d_dk(k: np.ndarray, axis: int) -> np.ndarray:
    """Analytic derivative of d_closed with respect to dimensionless k_axis."""
    kx, ky, kz = map(float, k)
    ka, kd, m = P.k_axis, P.k_diag, P.m_bulk
    s = [np.sin(kx), np.sin(ky), np.sin(kz)]
    c = [np.cos(kx), np.cos(ky), np.cos(kz)]
    out = np.zeros((3, 3), dtype=float)
    if axis == 0:
        out[0, 0] = 2 * ka * s[0] + 2 * kd * (s[0] * c[1] + s[0] * c[2])
        out[1, 1] = 2 * kd * s[0] * c[1]
        out[2, 2] = 2 * kd * s[0] * c[2]
        out[0, 1] = out[1, 0] = 2 * kd * c[0] * s[1]
        out[0, 2] = out[2, 0] = 2 * kd * c[0] * s[2]
    elif axis == 1:
        out[1, 1] = 2 * ka * s[1] + 2 * kd * (c[0] * s[1] + s[1] * c[2])
        out[0, 0] = 2 * kd * c[0] * s[1]
        out[2, 2] = 2 * kd * s[1] * c[2]
        out[0, 1] = out[1, 0] = 2 * kd * s[0] * c[1]
        out[1, 2] = out[2, 1] = 2 * kd * c[1] * s[2]
    else:
        out[2, 2] = 2 * ka * s[2] + 2 * kd * (c[0] * s[2] + c[1] * s[2])
        out[0, 0] = 2 * kd * c[0] * s[2]
        out[1, 1] = 2 * kd * c[1] * s[2]
        out[0, 2] = out[2, 0] = 2 * kd * s[0] * c[2]
        out[1, 2] = out[2, 1] = 2 * kd * s[1] * c[2]
    return out / m


def group_velocity(k: np.ndarray, branch: int) -> np.ndarray:
    vals, vecs = eig(k)
    omega = math.sqrt(float(vals[branch]))
    e = vecs[:, branch]
    return np.array([float(e.T @ d_dk(k, a) @ e) / (2.0 * omega) for a in range(3)])


def direction_from_angles(x: np.ndarray) -> np.ndarray:
    theta, phi = map(float, x)
    return np.array(
        [math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)],
        dtype=float,
    )


def root_on_direction(branch: int, direction: np.ndarray) -> float:
    direction = np.asarray(direction, dtype=float)
    direction /= norm(direction)
    nonzero = np.abs(direction) > 1e-14
    rmax = float(np.min(np.pi / np.abs(direction[nonzero])))
    f = lambda r: float(eig(r * direction)[0][branch] - TARGET)
    # The independently verified radial census finds exactly one crossing on each
    # reference shell. Use the direct bracket whenever the Brillouin-zone edge
    # lies above the target; otherwise fall back to a deterministic scan.
    if f(rmax) >= 0.0:
        return float(brentq(f, 0.0, rmax, xtol=1e-13, rtol=1e-13, maxiter=300))
    rs = np.linspace(0.0, rmax, 301)
    fs = np.array([f(r) for r in rs])
    candidates = np.where(fs[:-1] * fs[1:] <= 0.0)[0]
    if candidates.size == 0:
        raise RuntimeError(f"no shell root found for branch {branch}")
    j = int(candidates[0])
    return float(brentq(f, rs[j], rs[j + 1], xtol=1e-13, rtol=1e-13, maxiter=300))


def shell_point(branch: int, angles: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    direction = direction_from_angles(angles)
    radius = root_on_direction(branch, direction)
    k = radius * direction
    vals, vecs = eig(k)
    return k, vals, vecs


def qt_projector(vecs: np.ndarray) -> np.ndarray:
    return np.eye(3) - np.outer(vecs[:, 2], vecs[:, 2])


def polynomial_ql_projector(d: np.ndarray, vals: np.ndarray) -> np.ndarray:
    l0, l1, l2 = map(float, vals)
    eye = np.eye(3)
    return ((d - l0 * eye) @ (d - l1 * eye)) / ((l2 - l0) * (l2 - l1))


def aperture(delta: float) -> complex:
    i = np.arange(PROFILE.size, dtype=float)
    return complex(np.sum(PROFILE * np.exp(1j * delta * (i + 0.5))))


def source_full(k: np.ndarray) -> np.ndarray:
    kx, ky, kz = map(float, k)
    ap = aperture(KAPPA - kx)
    return (8.0 / 3.0) * ap * np.array(
        [
            -math.cos(KAPPA / 2.0) * math.cos(ky / 2.0) * math.cos(kz / 2.0),
            math.sin(KAPPA / 2.0) * math.sin(ky / 2.0) * math.cos(kz / 2.0),
            math.sin(KAPPA / 2.0) * math.cos(ky / 2.0) * math.sin(kz / 2.0),
        ],
        dtype=complex,
    )


def source_radial(k: np.ndarray) -> np.ndarray:
    kx, ky, kz = map(float, k)
    ap = aperture(KAPPA - kx)
    return (8.0 / 3.0) * ap * np.array(
        [
            0.0,
            math.sin(KAPPA / 2.0) * math.sin(ky / 2.0) * math.cos(kz / 2.0),
            math.sin(KAPPA / 2.0) * math.cos(ky / 2.0) * math.sin(kz / 2.0),
        ],
        dtype=complex,
    )


def active_basis(nseg: int) -> np.ndarray:
    b = np.zeros((12 * nseg, 2 * nseg), dtype=complex)
    for s in range(nseg):
        for a, (sy, sz) in enumerate(CORE_SIGNS):
            j = 3 * (4 * s + a)
            b[j, 2 * s] = 0.5
            b[j + 1, 2 * s + 1] = sy / (2.0 * math.sqrt(2.0))
            b[j + 2, 2 * s + 1] = sz / (2.0 * math.sqrt(2.0))
    return b


def finite_core_matrices(g_segments: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nseg = len(g_segments)
    coords = np.empty((4 * nseg, 3), dtype=int)
    u = np.zeros((12 * nseg, 12 * nseg), dtype=complex)
    v = np.zeros((12 * nseg, nseg + 1), dtype=complex)
    for s, g in enumerate(np.asarray(g_segments, dtype=float)):
        for a, ((yy, zz), (sy, sz)) in enumerate(zip(CORE_TRANS, CORE_SIGNS)):
            site = 4 * s + a
            coords[site] = [s, yy, zz]
            sl = slice(3 * site, 3 * site + 3)
            nminus = np.array([1.0, sy, sz]) / math.sqrt(3.0)
            nplus = np.array([-1.0, sy, sz]) / math.sqrt(3.0)
            u[sl, sl] += g * (np.outer(nminus, nminus) + np.outer(nplus, nplus))
            v[sl, s] += g * nminus * (-1.0 / math.sqrt(3.0))
            v[sl, s + 1] += g * nplus * (1.0 / math.sqrt(3.0))
    return coords, u, v


def bulk_green_realspace(n: int, eta: float = 0.03) -> np.ndarray:
    k = 2.0 * math.pi * np.fft.fftfreq(n)
    kx = k[:, None, None]
    ky = k[None, :, None]
    kz = k[None, None, :]
    kgrid = np.zeros((n, n, n, 3, 3), dtype=float)
    bonds = [(d, P.k_axis) for d in AXIAL_BONDS] + [(d, P.k_diag) for d in FACE_DIAGONAL_BONDS]
    for d, stiffness in bonds:
        nhat = np.asarray(d, dtype=float)
        nhat /= norm(nhat)
        phase = kx * d[0] + ky * d[1] + kz * d[2]
        kgrid += 2.0 * stiffness * (1.0 - np.cos(phase))[..., None, None] * np.outer(nhat, nhat)
    z = complex(OMEGA, eta)
    gk = np.linalg.inv(kgrid.astype(complex) - P.m_bulk * z * z * np.eye(3))
    gr = np.empty_like(gk)
    for a in range(3):
        for b in range(3):
            gr[..., a, b] = np.fft.ifftn(gk[..., a, b])
    return gr


def projected_core_green(gr: np.ndarray, coords: np.ndarray) -> np.ndarray:
    shape = np.array(gr.shape[:3])
    m = len(coords)
    g = np.empty((3 * m, 3 * m), dtype=complex)
    for i, ri in enumerate(coords):
        for j, rj in enumerate(coords):
            d = (ri - rj) % shape
            g[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = gr[tuple(d)]
    return g


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, complex):
        return [value.real, value.imag]
    raise TypeError(type(value).__name__)
