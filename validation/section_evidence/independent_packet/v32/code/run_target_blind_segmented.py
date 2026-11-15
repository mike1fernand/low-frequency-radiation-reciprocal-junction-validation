#!/usr/bin/env python3
"""Portable segmented driver for the P1-S target-blind direct-dynamics producer.

Each segment runs in a fresh Python process. This avoids long single-process
execution limits without changing a single velocity-Verlet step. Checkpointed
and unsegmented trajectories have been verified bit-for-bit on a control case.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase12-root', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--temp-dir')
    ap.add_argument('--segment-duration', type=float, default=5.0)
    ap.add_argument('--keep-checkpoints', action='store_true')
    ap.add_argument('--nx', type=int, default=120)
    ap.add_argument('--ny', type=int, default=56)
    ap.add_argument('--nz', type=int, default=56)
    ap.add_argument('--dt', type=float, default=0.04)
    ap.add_argument('--tmax', type=float, default=57.0)
    ap.add_argument('--sample-start', type=float, default=48.0)
    ap.add_argument('--sample-step', type=float, default=1.0)
    ap.add_argument('--trace-step', type=float, default=1.0)
    ap.add_argument('--omega', type=float, default=1.2)
    ap.add_argument('--sigma', type=float, default=7.0)
    ap.add_argument('--x0', type=float, default=18.0)
    ap.add_argument('--slab-start', type=int, default=60)
    ap.add_argument('--amplitude', type=float, default=0.002)
    args = ap.parse_args()

    worker = Path(__file__).with_name('run_target_blind_global_td.py')
    phase12 = Path(args.phase12_root)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    base = [
        sys.executable, str(worker), '--phase12-root', str(phase12), '--output', str(output),
        '--nx', str(args.nx), '--ny', str(args.ny), '--nz', str(args.nz),
        '--dt', str(args.dt), '--tmax', str(args.tmax), '--sample-start', str(args.sample_start),
        '--sample-step', str(args.sample_step), '--trace-step', str(args.trace_step),
        '--omega', str(args.omega), '--sigma', str(args.sigma), '--x0', str(args.x0),
        '--slab-start', str(args.slab_start), '--amplitude', str(args.amplitude),
    ]
    ends = []
    t = args.segment_duration
    while t < args.tmax - 0.5 * args.dt:
        ends.append(round(t, 12))
        t += args.segment_duration
    ends.append(args.tmax)

    owned_temp = args.temp_dir is None
    tmp_context = tempfile.TemporaryDirectory(prefix='p1s_td_') if owned_temp else None
    temp_dir = Path(tmp_context.name if tmp_context else args.temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    resume: Path | None = None
    checkpoints: list[Path] = []
    try:
        for idx, end in enumerate(ends):
            cmd = list(base)
            if resume is not None:
                cmd += ['--resume-from', str(resume)]
            final = math.isclose(end, args.tmax, rel_tol=0.0, abs_tol=0.5 * args.dt)
            if not final:
                checkpoint = temp_dir / f'checkpoint_{idx+1:02d}_t{end:g}.npz'
                cmd += ['--segment-end', str(end), '--checkpoint-out', str(checkpoint)]
                checkpoints.append(checkpoint)
            subprocess.run(cmd, check=True)
            if not final:
                resume = checkpoint
        print(f'completed segmented target-blind trajectory: {output}')
    finally:
        if not args.keep_checkpoints:
            for p in checkpoints:
                p.unlink(missing_ok=True)
        if tmp_context is not None:
            tmp_context.cleanup()


if __name__ == '__main__':
    main()
