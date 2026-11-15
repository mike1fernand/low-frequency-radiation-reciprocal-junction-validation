#!/usr/bin/env python3
"""Target-blind global exact-cluster time-domain producer for Paper 1 P1-S.

This script evolves the original full real-space Hamiltonian with exactly zero
initial bulk displacement and velocity. It records only direct-dynamics
quantities until the target-blind reporting window is selected. The
frequency-domain packet prediction is intentionally absent from this producer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np


def symmetric_eigenvalues_3x3(a: np.ndarray) -> np.ndarray:
    """Vectorized ordered eigenvalues for real symmetric 3x3 matrices."""
    a00, a01, a02 = a[..., 0, 0], a[..., 0, 1], a[..., 0, 2]
    a11, a12, a22 = a[..., 1, 1], a[..., 1, 2], a[..., 2, 2]
    q = (a00 + a11 + a22) / 3.0
    p2 = ((a00-q)**2 + (a11-q)**2 + (a22-q)**2 + 2.0*(a01*a01 + a02*a02 + a12*a12)) / 6.0
    p = np.sqrt(np.maximum(p2, 0.0))
    nz = p > 1.0e-30
    b00 = np.divide(a00-q, p, out=np.zeros_like(q), where=nz)
    b01 = np.divide(a01, p, out=np.zeros_like(q), where=nz)
    b02 = np.divide(a02, p, out=np.zeros_like(q), where=nz)
    b11 = np.divide(a11-q, p, out=np.zeros_like(q), where=nz)
    b12 = np.divide(a12, p, out=np.zeros_like(q), where=nz)
    b22 = np.divide(a22-q, p, out=np.zeros_like(q), where=nz)
    det_b = b00*(b11*b22-b12*b12) - b01*(b01*b22-b12*b02) + b02*(b01*b12-b11*b02)
    phi = np.arccos(np.clip(det_b/2.0, -1.0, 1.0)) / 3.0
    lam_max = q + 2.0*p*np.cos(phi)
    lam_min = q + 2.0*p*np.cos(phi + 2.0*np.pi/3.0)
    lam_mid = 3.0*q - lam_min - lam_max
    lam_min = np.where(nz, lam_min, q)
    lam_mid = np.where(nz, lam_mid, q)
    lam_max = np.where(nz, lam_max, q)
    return np.stack((lam_min, lam_mid, lam_max), axis=-1)


def upper_eigenvector_3x3(a: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """Vectorized normalized eigenvector of the isolated upper branch."""
    m = a - lam[..., None, None] * np.eye(3)
    r0, r1, r2 = m[..., 0, :], m[..., 1, :], m[..., 2, :]
    c01, c02, c12 = np.cross(r0, r1), np.cross(r0, r2), np.cross(r1, r2)
    candidates = np.stack((c01, c02, c12), axis=-2)
    norms2 = np.stack((np.sum(c01*c01, axis=-1), np.sum(c02*c02, axis=-1), np.sum(c12*c12, axis=-1)), axis=-1)
    pick = np.argmax(norms2, axis=-1)
    v = np.take_along_axis(candidates, pick[..., None, None], axis=-2)[..., 0, :]
    norm = np.linalg.norm(v, axis=-1)
    fallback = np.zeros_like(v)
    fallback[..., 2] = 1.0  # matches the conventional k=0 ordering used by eigh
    return np.divide(v, norm[..., None], out=fallback, where=norm[..., None] > 1.0e-13)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase12-root", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--nx", type=int, default=120)
    ap.add_argument("--ny", type=int, default=56)
    ap.add_argument("--nz", type=int, default=56)
    ap.add_argument("--dt", type=float, default=0.04)
    ap.add_argument("--tmax", type=float, default=57.0)
    ap.add_argument("--sample-start", type=float, default=48.0)
    ap.add_argument("--sample-step", type=float, default=1.0)
    ap.add_argument("--trace-step", type=float, default=1.0)
    ap.add_argument("--omega", type=float, default=1.2)
    ap.add_argument("--sigma", type=float, default=7.0)
    ap.add_argument("--x0", type=float, default=18.0)
    ap.add_argument("--slab-start", type=int, default=60)
    ap.add_argument("--amplitude", type=float, default=0.002)
    ap.add_argument("--resume-from")
    ap.add_argument("--checkpoint-out")
    ap.add_argument("--segment-end", type=float)
    args = ap.parse_args()

    root = Path(args.phase12_root).resolve()
    sys.path.insert(0, str(root))
    from model_api_v1_1 import ModelParameters
    from open_boundary_time_domain import (
        carrier_k,
        exact_right_velocity,
        forces_open,
        coupling_current,
        subsystem_energies,
    )
    from reduced_core import bulk_stiffness_batch

    p = ModelParameters()
    nx, ny, nz = args.nx, args.ny, args.nz
    g_segments = np.array([0.1, 0.2, 0.2, 0.2, 0.2, 0.1], dtype=float)
    gprof = np.zeros(nx - 1, dtype=float)
    gprof[args.slab_start : args.slab_start + len(g_segments)] = g_segments
    scientific_parameters = {
        "nx": nx, "ny": ny, "nz": nz, "dt": args.dt, "tmax": args.tmax,
        "sample_start": args.sample_start, "sample_step": args.sample_step, "trace_step": args.trace_step,
        "omega": args.omega, "sigma": args.sigma, "x0": args.x0, "slab_start": args.slab_start,
        "amplitude": args.amplitude, "carrier_k": carrier_k(args.omega, p), "g_segments": g_segments.tolist(),
    }
    fingerprint = hashlib.sha256(json.dumps(scientific_parameters, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    if args.resume_from:
        cp = np.load(args.resume_from, allow_pickle=False)
        meta = json.loads(str(cp["metadata_json"].item()))
        if meta["parameter_fingerprint"] != fingerprint:
            raise ValueError("Checkpoint parameters do not match this run")
        q, vq, u, vu = cp["q"], cp["vq"], cp["u"], cp["vu"]
        fq, fu = cp["fq"], cp["fu"]
        start_step = int(meta["step"])
        E0, Eb0, Ec0, total0 = (float(meta[x]) for x in ("E0", "Eb0", "Ec0", "total0"))
        jprev = float(meta["jprev"])
        cumulative_current = float(meta["cumulative_current"])
        initial = meta["initial"]
        trace = meta["trace"]
        snapshots = meta["snapshots"]
        elapsed_accum = float(meta.get("elapsed_accum", 0.0))
    else:
        x = np.arange(nx, dtype=float)
        k0 = scientific_parameters["carrier_k"]
        q = args.amplitude * np.exp(-0.5 * ((x - args.x0) / args.sigma) ** 2) * np.cos(k0 * (x - args.x0))
        vq = exact_right_velocity(q, p)
        u = np.zeros((nx, ny, nz, 3), dtype=float)
        vu = np.zeros_like(u)
        initial_bulk_nonzero_count = int(np.count_nonzero(u) + np.count_nonzero(vu))
        initial_bulk_displacement_norm = float(np.linalg.norm(u))
        initial_bulk_velocity_norm = float(np.linalg.norm(vu))
        fq = np.zeros_like(q)
        fu = np.zeros_like(u)
        forces_open(q, u, gprof, p.k_spine, p.k_axis, p.k_diag, fq, fu)
        E0, Eb0, Ec0 = subsystem_energies(
            q, vq, u, vu, gprof, p.m_spine, p.m_bulk, p.k_spine, p.k_axis, p.k_diag
        )
        total0 = E0 + Eb0 + Ec0
        jprev = coupling_current(q, vq, u, vu, gprof)
        cumulative_current = 0.0
        start_step = 0
        trace, snapshots = [], []
        elapsed_accum = 0.0
        initial = {
            "line_energy_E0": E0, "bulk_energy": Eb0, "coupling_energy": Ec0,
            "bulk_exactly_zero": bool(initial_bulk_nonzero_count == 0),
            "bulk_nonzero_coordinate_count": initial_bulk_nonzero_count,
            "bulk_displacement_norm": initial_bulk_displacement_norm,
            "bulk_velocity_norm": initial_bulk_velocity_norm,
            "Ec0_over_E0": Ec0 / E0,
        }

    sample_times = np.arange(args.sample_start, args.tmax + 0.5 * args.sample_step, args.sample_step)
    sample_steps = {int(round(t / args.dt)): float(t) for t in sample_times}
    trace_stride = max(1, int(round(args.trace_step / args.dt)))
    # Precompute periodic k-grid and dynamical matrices for exact branch energies.
    kx = 2 * np.pi * np.fft.fftfreq(nx)[:, None, None]
    ky = 2 * np.pi * np.fft.fftfreq(ny)[None, :, None]
    kz = 2 * np.pi * np.fft.fftfreq(nz)[None, None, :]
    dyn = bulk_stiffness_batch(kx, ky, kz, p) / p.m_bulk
    eigvals = symmetric_eigenvalues_3x3(dyn)
    lam_upper = eigvals[..., 2]
    upper_vec = upper_eigenvector_3x3(dyn, lam_upper)
    N = nx * ny * nz

    def cluster_energies() -> tuple[float, float, float]:
        """Return invariant lower-two, upper-one and total periodic bulk energies."""
        uh = np.fft.fftn(u, axes=(0, 1, 2))
        vh = np.fft.fftn(vu, axes=(0, 1, 2))
        au = np.einsum("...a,...a->...", upper_vec, uh, optimize=True)
        av = np.einsum("...a,...a->...", upper_vec, vh, optimize=True)
        e_upper_total = float((0.5 * p.m_bulk * (np.abs(av)**2 + lam_upper*np.abs(au)**2) / N).sum())
        kinetic = float(np.sum(np.abs(vh)**2))
        potential = float(np.real(np.einsum("...a,...ab,...b->", uh.conj(), dyn, uh, optimize=True)))
        e_total = 0.5 * p.m_bulk * (kinetic + potential) / N
        e_lower = e_total - e_upper_total
        return e_lower, e_upper_total, e_total

    def boundary_shell_kinetic_fraction(width: int = 2) -> float:
        kinetic = np.sum(vu * vu, axis=-1)
        mask = np.zeros((nx, ny, nz), dtype=bool)
        # Transverse boundaries are the relevant earliest-return surfaces.
        mask[:, :width, :] = True
        mask[:, -width:, :] = True
        mask[:, :, :width] = True
        mask[:, :, -width:] = True
        den = float(kinetic.sum())
        return float(kinetic[mask].sum() / den) if den else 0.0

    steps = int(round(args.tmax / args.dt))
    segment_end = args.tmax if args.segment_end is None else args.segment_end
    end_step = int(round(segment_end / args.dt))
    if not (start_step < end_step <= steps):
        raise ValueError(f"Require start_step < end_step <= total steps; got {start_step}, {end_step}, {steps}")
    started = time.time()
    for step in range(start_step + 1, end_step + 1):
        vq += 0.5 * args.dt * fq / p.m_spine
        vu += 0.5 * args.dt * fu / p.m_bulk
        q += args.dt * vq
        u += args.dt * vu
        forces_open(q, u, gprof, p.k_spine, p.k_axis, p.k_diag, fq, fu)
        vq += 0.5 * args.dt * fq / p.m_spine
        vu += 0.5 * args.dt * fu / p.m_bulk

        jnow = coupling_current(q, vq, u, vu, gprof)
        cumulative_current += 0.5 * args.dt * (jprev + jnow)
        jprev = jnow

        do_trace = (step % trace_stride == 0 or step == steps)
        do_snapshot = step in sample_steps
        if do_trace or do_snapshot:
            Eq, Eb, Ec = subsystem_energies(
                q, vq, u, vu, gprof, p.m_spine, p.m_bulk, p.k_spine, p.k_axis, p.k_diag
            )
            common = {
                "T_current": cumulative_current / E0,
                "T_coupling": Ec / E0,
                "T_bulk_open": Eb / E0,
                "total_energy_relative_drift": (Eq + Eb + Ec - total0) / total0,
                "instantaneous_coupling_current_over_E0": jnow / E0,
            }
            if do_trace:
                trace.append({"time": step * args.dt, **common})
            if do_snapshot:
                e_lower, e_upper_energy, tmodal = cluster_energies()
                snapshots.append(
                    {
                        "time": sample_steps[step],
                        **common,
                        "cluster_energy_over_E0": {"lower_two": e_lower / E0, "upper_one": e_upper_energy / E0},
                        "TT_modal": e_lower / E0,
                        "TL_modal": e_upper_energy / E0,
                        "Tbulk_modal": tmodal / E0,
                        "periodic_open_bulk_relative_difference": abs(tmodal - Eb) / Eb if Eb else 0.0,
                        "boundary_shell_kinetic_fraction_width2": boundary_shell_kinetic_fraction(2),
                    }
                )
                print(f"snapshot t={sample_steps[step]:.1f} TT={snapshots[-1]['TT_modal']:.9f}", flush=True)

    segment_elapsed = time.time() - started
    elapsed_total = elapsed_accum + segment_elapsed
    if end_step < steps:
        if not args.checkpoint_out:
            raise ValueError("A nonfinal segment requires --checkpoint-out")
        checkpoint_meta = {
            "parameter_fingerprint": fingerprint, "scientific_parameters": scientific_parameters,
            "step": end_step, "E0": E0, "Eb0": Eb0, "Ec0": Ec0, "total0": total0,
            "jprev": jprev, "cumulative_current": cumulative_current, "initial": initial,
            "trace": trace, "snapshots": snapshots, "elapsed_accum": elapsed_total,
        }
        np.savez(
            args.checkpoint_out, q=q, vq=vq, u=u, vu=vu, fq=fq, fu=fu,
            metadata_json=np.array(json.dumps(checkpoint_meta, sort_keys=True)),
        )
        print(json.dumps({"checkpoint": args.checkpoint_out, "completed_time": end_step*args.dt, "segment_elapsed_seconds": segment_elapsed}, indent=2))
        return

    out = {
        "method": "full Phase-12 Hamiltonian; open real-space velocity-Verlet evolution; exact invariant lower-two/upper-one periodic Bloch-cluster snapshot energies via analytic symmetric-3x3 spectral projector",
        "target_blind": True,
        "frequency_domain_prediction_used_during_production": False,
        "parameters": scientific_parameters | {"phase12_dependency": Path(args.phase12_root).name},
        "initial": initial,
        "trace": trace,
        "snapshots": snapshots,
        "elapsed_seconds": elapsed_total,
        "segmented_execution": bool(args.resume_from),
    }
    Path(args.output).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"output": args.output, "elapsed_seconds": out["elapsed_seconds"], "snapshots": len(snapshots)}, indent=2))


if __name__ == "__main__":
    main()
