#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEP = ROOT / "dependencies" / "phase12_core"
RAW = DEP / "raw"
sys.path.insert(0, str(DEP))

from model_api_v1_1 import (  # noqa: E402
    AXIAL_BONDS,
    FACE_DIAGONAL_BONDS,
    ModelParameters,
    _unit,
    bulk_dynamical_matrix,
)
from reduced_core import (  # noqa: E402
    analytic_q_lead_self_energy,
    cross_section_form_factor,
    reduced_static_matrices,
)

P = ModelParameters()
OMEGA = 1.2
PROFILE = np.array([0.1, 0.2, 0.2, 0.2, 0.2, 0.1], dtype=float)
ASYM_PROFILE = np.array([0.07, 0.18, 0.11, 0.23, 0.16, 0.09], dtype=float)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def herm(a: np.ndarray) -> np.ndarray:
    return (a + a.conj().T) / 2


def imh(a: np.ndarray) -> np.ndarray:
    return (a - a.conj().T) / (2j)


def rel(a: np.ndarray | complex | float, b: np.ndarray | complex | float, floor: float = 1e-300) -> float:
    aa = np.asarray(a)
    bb = np.asarray(b)
    return float(np.linalg.norm(aa - bb) / max(np.linalg.norm(aa), np.linalg.norm(bb), floor))


def rel_scalar(a: float, b: float, floor: float = 1e-30) -> float:
    return float(abs(a - b) / max(abs(a), abs(b), floor))


def complex_json(z: complex) -> list[float]:
    return [float(np.real(z)), float(np.imag(z))]


def load_raw(nseg: int = 6) -> dict[str, Any]:
    pv = np.load(RAW / "pv_mean_M192.npy")[: 2 * nseg, : 2 * nseg]
    shell = np.load(RAW / "on_shell_N80.npz")
    branches = [shell[f"Im_branch_{i}"][: 2 * nseg, : 2 * nseg] for i in range(3)]
    total = shell["Im_total"][: 2 * nseg, : 2 * nseg]
    return {"pv": pv, "branches": branches, "total": total, "g0": pv + 1j * total}


def line_kappa_gamma(omega: float = OMEGA, p: ModelParameters = P) -> tuple[float, complex, float]:
    if not 0 < omega < 2 * math.sqrt(p.k_spine / p.m_spine):
        raise ValueError("omega outside propagating scalar band")
    kappa = 2 * math.asin(omega / (2 * math.sqrt(p.k_spine / p.m_spine)))
    sigma = -p.k_spine * np.exp(1j * kappa)
    gamma = 2 * p.k_spine * math.sin(kappa)
    return float(kappa), complex(sigma), float(gamma)


def line_stiffness(n: int, p: ModelParameters = P, endpoint_diag: float | None = None) -> np.ndarray:
    if n < 1:
        raise ValueError
    diag = 2 * p.k_spine if endpoint_diag is None else endpoint_diag
    K = np.zeros((n, n), dtype=float)
    np.fill_diagonal(K, 2 * p.k_spine)
    if endpoint_diag is not None:
        K[0, 0] = K[-1, -1] = diag
    for i in range(n - 1):
        K[i, i + 1] = K[i + 1, i] = -p.k_spine
    return K


def dressed_objects(profile: np.ndarray = PROFILE, omega: float = OMEGA) -> dict[str, Any]:
    raw = load_raw(len(profile))
    g0 = raw["g0"]
    branches = raw["branches"]
    total = raw["total"]
    U, V, Kself = reduced_static_matrices(profile, P)
    ident = np.eye(2 * len(profile))
    gc_left = np.linalg.solve(ident + g0 @ U, g0)
    gc_right = g0 @ np.linalg.inv(ident + U @ g0)
    sigma_bulk = -V.conj().T @ gc_left @ V
    gamma_bulk = 1j * (sigma_bulk - sigma_bulk.conj().T)
    left_dress = np.linalg.inv(ident + g0.conj().T @ U)
    right_dress = np.linalg.inv(ident + U @ g0)
    gamma_branches = [2 * V.conj().T @ left_dress @ im @ right_dress @ V for im in branches]
    kappa, sigma, gamma = line_kappa_gamma(omega)
    nq = len(profile) + 1
    KD = line_stiffness(nq)
    SL = np.zeros((nq, nq), dtype=complex)
    SR = np.zeros_like(SL)
    SL[0, 0] = sigma
    SR[-1, -1] = sigma
    GL = 1j * (SL - SL.conj().T)
    GR = 1j * (SR - SR.conj().T)
    A = P.m_spine * omega**2 * np.eye(nq) - KD - Kself - sigma_bulk - SL - SR
    GD = np.linalg.inv(A)
    GA = GD.conj().T
    eL = np.eye(nq, dtype=complex)[:, 0]
    eR = np.eye(nq, dtype=complex)[:, -1]
    chiL = GD @ (math.sqrt(gamma) * eL)
    chiR = GD @ (math.sqrt(gamma) * eR)
    etaL = right_dress @ V @ chiL
    etaR = right_dress @ V @ chiR
    return {
        "raw": raw,
        "g0": g0,
        "branches": branches,
        "total": total,
        "U": U,
        "V": V,
        "Kself": Kself,
        "gc_left": gc_left,
        "gc_right": gc_right,
        "sigma_bulk": sigma_bulk,
        "gamma_bulk": gamma_bulk,
        "left_dress": left_dress,
        "right_dress": right_dress,
        "gamma_branches": gamma_branches,
        "kappa": kappa,
        "sigma": sigma,
        "gamma": gamma,
        "KD": KD,
        "SL": SL,
        "SR": SR,
        "GL": GL,
        "GR": GR,
        "GD": GD,
        "GA": GA,
        "chiL": chiL,
        "chiR": chiR,
        "etaL": etaL,
        "etaR": etaR,
    }


def scattering_summary(objects: dict[str, Any], incidence: str = "left") -> dict[str, Any]:
    GD = objects["GD"]
    gamma = objects["gamma"]
    gamma_bulk = objects["gamma_bulk"]
    gamma_branches = objects["gamma_branches"]
    if incidence == "left":
        chi = objects["chiL"]
        other_gamma = objects["GR"]
        contact = 0
        other = -1
    elif incidence == "right":
        chi = objects["chiR"]
        other_gamma = objects["GL"]
        contact = -1
        other = 0
    else:
        raise ValueError
    r = 1 - 1j * gamma * GD[contact, contact]
    t = -1j * gamma * GD[other, contact]
    R = float(abs(r) ** 2)
    T = float(np.real(np.vdot(chi, other_gamma @ chi)))
    Tb = float(np.real(np.vdot(chi, gamma_bulk @ chi)))
    Tbr = [float(np.real(np.vdot(chi, G @ chi))) for G in gamma_branches]
    return {
        "R": R,
        "Tq": T,
        "Tbulk": Tb,
        "Tbranches": Tbr,
        "closure": float(R + T + Tb),
        "closure_abs": float(abs(R + T + Tb - 1)),
        "r": complex_json(r),
        "t": complex_json(t),
    }


def independent_bulk_matrix(k: np.ndarray, p: ModelParameters = P) -> np.ndarray:
    kx, ky, kz = map(float, k)
    ka, kd, mb = p.k_axis, p.k_diag, p.m_bulk
    D = np.zeros((3, 3), dtype=float)
    D[0, 0] = 2 * ka * (1 - math.cos(kx)) + 2 * kd * (
        (1 - math.cos(kx) * math.cos(ky)) + (1 - math.cos(kx) * math.cos(kz))
    )
    D[1, 1] = 2 * ka * (1 - math.cos(ky)) + 2 * kd * (
        (1 - math.cos(kx) * math.cos(ky)) + (1 - math.cos(ky) * math.cos(kz))
    )
    D[2, 2] = 2 * ka * (1 - math.cos(kz)) + 2 * kd * (
        (1 - math.cos(kx) * math.cos(kz)) + (1 - math.cos(ky) * math.cos(kz))
    )
    D[0, 1] = D[1, 0] = 2 * kd * math.sin(kx) * math.sin(ky)
    D[0, 2] = D[2, 0] = 2 * kd * math.sin(kx) * math.sin(kz)
    D[1, 2] = D[2, 1] = 2 * kd * math.sin(ky) * math.sin(kz)
    return D / mb


def independent_dD(k: np.ndarray, axis: int, p: ModelParameters = P) -> np.ndarray:
    kx, ky, kz = map(float, k)
    s = [math.sin(kx), math.sin(ky), math.sin(kz)]
    c = [math.cos(kx), math.cos(ky), math.cos(kz)]
    ka, kd, mb = p.k_axis, p.k_diag, p.m_bulk
    out = np.zeros((3, 3), dtype=float)
    if axis == 0:
        out[0, 0] = 2 * ka * s[0] + 2 * kd * s[0] * (c[1] + c[2])
        out[1, 1] = 2 * kd * s[0] * c[1]
        out[2, 2] = 2 * kd * s[0] * c[2]
        out[0, 1] = out[1, 0] = 2 * kd * c[0] * s[1]
        out[0, 2] = out[2, 0] = 2 * kd * c[0] * s[2]
    elif axis == 1:
        out[1, 1] = 2 * ka * s[1] + 2 * kd * s[1] * (c[0] + c[2])
        out[0, 0] = 2 * kd * c[0] * s[1]
        out[2, 2] = 2 * kd * s[1] * c[2]
        out[0, 1] = out[1, 0] = 2 * kd * s[0] * c[1]
        out[1, 2] = out[2, 1] = 2 * kd * c[1] * s[2]
    elif axis == 2:
        out[2, 2] = 2 * ka * s[2] + 2 * kd * s[2] * (c[0] + c[1])
        out[0, 0] = 2 * kd * c[0] * s[2]
        out[1, 1] = 2 * kd * c[1] * s[2]
        out[0, 2] = out[2, 0] = 2 * kd * s[0] * c[2]
        out[1, 2] = out[2, 1] = 2 * kd * s[1] * c[2]
    else:
        raise ValueError
    return out / mb


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")
