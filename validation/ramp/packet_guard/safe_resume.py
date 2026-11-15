#!/usr/bin/env python3
"""Validate complete packet histories before checkpoint continuation.

This additive wrapper leaves the archived engine unchanged.  It checks the
checkpoint/model binding, complete trace sizes, monitor/time grids, and saved
endpoints.  A checkpoint alone is deliberately insufficient: the standard
analyzer requires the preceding history, beginning at step zero.

Dry-run is the default. With --execute, the pinned engine source is freshly
compiled and only then invoked. The prior configuration is backed up exactly;
successful continuation updates its final target fields for the analyzer. Final
targets at/after tau48 must be on the 0.2 snapshot grid. The wrapper refuses a
stale checkpoint, missing history, an existing receipt/lock, or a different source/model. It does not
promise byte-for-byte agreement across platforms, transaction rollback after an
interrupted engine, or detect arbitrary corruption of interior trace values.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

SOURCE_SHA256 = "661925e394cb680600e0b0f98f6019216bff063c71900e917eccb27a22dfc542"
BONDS = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0),
         (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
PROFILE = tuple((i + 1) / 7 for i in range(6))
CONSTANTS = (1.0, 1.0, 1.0, 1.69)
HEADER_SIZE = 36 + 48 + 108 + 32 + 24
TRACE_NAMES = ("line.bin", "corner.bin", "snap.bin", "snap_times.bin", "monitor.csv")


class InvalidResume(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidResume(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_float(a, b):
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= 2e-12 * max(1.0, abs(a), abs(b))


def checkpoint_info(path):
    data = path.read_bytes()
    require(len(data) >= HEADER_SIZE, "truncated checkpoint header")
    magic, step, n, h, length, dt = struct.unpack_from("<8sqiiid", data)
    require(magic == b"P1PKT01\0", "wrong checkpoint magic")
    require(step >= 0 and n == 120 and h >= 1, "invalid checkpoint dimensions/step")
    require(length == n + 3 * n * h * h, "checkpoint degrees of freedom disagree")
    require(len(data) == HEADER_SIZE + 16 * length, "checkpoint byte count disagrees")
    require(math.isfinite(dt) and dt > 0, "invalid checkpoint timestep")
    g = struct.unpack_from("<6d", data, 36)
    b = tuple(zip(*[iter(struct.unpack_from("<27i", data, 84))] * 3))
    constants = struct.unpack_from("<4d", data, 192)
    cumulative, last_current, initial_energy = struct.unpack_from("<3d", data, 224)
    require(g == PROFILE and b == tuple(BONDS) and constants == CONSTANTS,
            "checkpoint is not the pinned six-motif/48-spring model")
    require(all(math.isfinite(x) for x in (cumulative, last_current, initial_energy)) and initial_energy > 0,
            "invalid checkpoint current/initial-energy metadata")
    # No state reconstruction or diagonal averaging: retain every x/v byte.
    states = data[HEADER_SIZE:]
    require(all(math.isfinite(x[0]) for x in struct.iter_unpack("<d", states)),
            "nonfinite checkpoint state")
    return dict(step=step, N=n, H=h, L=length, dt=dt, profile=g, bonds=b,
                constants=constants, cumulative=cumulative, last_current=last_current,
                initial_energy=initial_energy, data=data, path=str(path))


def tail_bytes(path, count):
    with path.open("rb") as stream:
        stream.seek(-count, os.SEEK_END)
        return stream.read(count)


def validate(packet_dir, run_dir, checkpoint_path, target_step, every):
    require(sys.byteorder == "little" and struct.calcsize("d") == 8 and struct.calcsize("i") == 4,
            "archived native binary layout requires little-endian float64/int32 runtime")
    require(type(target_step) is int and target_step >= 0 and type(every) is int and every > 0,
            "target step and checkpoint cadence must be positive integral parameters")
    source = packet_dir / "evolve_candidate.cpp"
    require(source.is_file() and sha256(source) == SOURCE_SHA256, "engine source hash is not the pinned source")
    require(run_dir.is_dir(), "run directory is absent")
    cp = checkpoint_info(checkpoint_path)
    n, h, length, dt, step = (cp[k] for k in ("N", "H", "L", "dt", "step"))
    require(target_step > step, "target must be later than checkpoint")
    config_path = run_dir / "configuration.json"
    require(config_path.is_file(), "run configuration.json is required")
    config = json.loads(config_path.read_text())
    require(config.get("N_line") == n and config.get("bulk_shape") == [n, 2*h, 2*h, 3]
            and config.get("quarter_shape") == [n, h, h, 3], "configuration geometry disagrees with checkpoint")
    require(config.get("dt") == dt and tuple(config.get("profile", [])) == PROFILE,
            "configuration timestep/profile disagrees with checkpoint")
    require(config.get("B_plus") == [list(b) for b in BONDS]
            and config.get("contact_segments") == list(range(60, 66)), "configuration bonds/contact geometry disagree")
    require(config.get("masses") == 1.0 and config.get("scalar_line_stiffness") == 1.69
            and config.get("bulk_axis_and_face_stiffness") == 1.0, "configuration material constants disagree")
    require(config.get("save_cadence") == .2 and config.get("snapshot_start") == 48.0,
            "configuration snapshot timing disagrees with engine")
    require(config.get("source_sha256") == SOURCE_SHA256, "configuration source binding disagrees")
    cadence = round(.2 / dt)
    require(cadence >= 1 and abs(cadence * dt - .2) <= 1e-12, "timestep must divide saved cadence 0.2")
    require(target_step * dt < 48.0 - 1e-12 or target_step % cadence == 0,
            "analyzer-ready target at/after tau48 must align with saved cadence 0.2")
    for name in TRACE_NAMES:
        require((run_dir / name).is_file(), f"required preceding history missing: {name}")
    # line/corner are recorded at *every* integral step, including zero.
    sizes = {name: (run_dir / name).stat().st_size for name in TRACE_NAMES}
    require(sizes["line.bin"] == (step + 1) * 2 * n * 8,
            "line history byte count does not end exactly at checkpoint")
    require(sizes["corner.bin"] == (step + 1) * 6 * 2 * 3 * 8,
            "corner history byte count does not end exactly at checkpoint")
    expected_steps = [i for i in range(0, step + 1, cadence) if i * dt >= 48.0 - 1e-12]
    require(sizes["snap.bin"] == len(expected_steps) * 2 * length * 8,
            "bulk snapshot byte count disagrees with full preceding history")
    require(sizes["snap_times.bin"] == len(expected_steps) * 8,
            "snapshot time byte count disagrees with checkpoint")
    times = [x[0] for x in struct.iter_unpack("<d", (run_dir / "snap_times.bin").read_bytes())]
    require(all(same_float(t, i * dt) for t, i in zip(times, expected_steps)),
            "snapshot timestamps disagree with required complete time grid")
    data = cp["data"]; xoff = HEADER_SIZE; voff = HEADER_SIZE + 8 * length
    expected_line = data[xoff:xoff + n*8] + data[voff:voff + n*8]
    require(tail_bytes(run_dir / "line.bin", 2*n*8) == expected_line,
            "line endpoint does not equal complete checkpoint state")
    corner = bytearray()
    for i in range(60, 66):
        offset = (n + 3 * i * h * h) * 8
        corner += data[xoff + offset:xoff + offset + 24]
        corner += data[voff + offset:voff + offset + 24]
    require(tail_bytes(run_dir / "corner.bin", len(corner)) == bytes(corner),
            "corner endpoint does not equal complete checkpoint state")
    if expected_steps and expected_steps[-1] == step:
        require(tail_bytes(run_dir / "snap.bin", 16 * length) == data[HEADER_SIZE:],
                "bulk snapshot endpoint does not equal complete checkpoint state")
    with (run_dir / "monitor.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames == ["step", "time", "energy", "initial_energy", "relative_drift", "current", "cumulative_current"],
                "monitor header disagrees with engine")
        rows = list(reader)
    require(rows, "empty monitor history")
    monitor_steps = []
    for row in rows:
        try:
            i = int(row["step"]); values = {k: float(v) for k, v in row.items() if k != "step"}
        except (TypeError, ValueError):
            raise InvalidResume("invalid monitor row")
        require(all(math.isfinite(v) for v in values.values()), "nonfinite monitor row")
        require(same_float(values["time"], i * dt) and values["initial_energy"] == cp["initial_energy"],
                "monitor time/initial energy disagrees with checkpoint")
        monitor_steps.append(i)
    require(monitor_steps[0] == 0 and monitor_steps[-1] == step
            and all(a < b for a, b in zip(monitor_steps, monitor_steps[1:])),
            "monitor does not span step zero through checkpoint exactly once")
    # Earlier stage endpoints may produce legitimate extra monitor rows.
    require(set(range(0, step + 1, every)).issubset(monitor_steps),
            "monitor misses required checkpoints for requested cadence; retain original cadence")
    last = rows[-1]
    require(float(last["current"]) == cp["last_current"]
            and float(last["cumulative_current"]) == cp["cumulative"],
            "monitor current endpoint disagrees with checkpoint")
    return dict(status="PASS", source_sha256=SOURCE_SHA256, configuration_sha256=sha256(config_path),
                checkpoint_sha256=sha256(checkpoint_path), checkpoint_step=step, target_step=target_step,
                N=n, H=h, L=length, dt=dt, monitor_every=every, trace_bytes=sizes,
                snapshot_count=len(expected_steps), monitor_count=len(rows),
                scope="Complete state, geometry/runtime binding, complete history lengths/time grids and endpoint consistency; interior trace corruption is not excluded.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-dir", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--target-step", type=int, required=True)
    parser.add_argument("--every", type=int, required=True, help="retain the original monitor/checkpoint cadence in steps")
    parser.add_argument("--execute", action="store_true", help="compile pinned source and append after validation")
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    packet_dir, run_dir, checkpoint_path = [p.resolve() for p in (args.packet_dir, args.run_dir, args.checkpoint)]
    if not args.execute:
        print(json.dumps(validate(packet_dir, run_dir, checkpoint_path, args.target_step, args.every), indent=2))
        return
    lock = run_dir / ".safe_resume.lock"
    receipt = run_dir / f"SAFE_RESUME_TO_{args.target_step}.json"
    backup = run_dir / f"CONFIGURATION_BEFORE_SAFE_RESUME_TO_{args.target_step}.json"
    temporary_config = run_dir / ".configuration.safe_resume.tmp"
    require(not receipt.exists(), "continuation receipt already exists; do not overwrite")
    require(not backup.exists(), "configuration backup already exists; do not overwrite")
    require(not temporary_config.exists(), "temporary configuration exists; diagnose prior interrupted write")
    # Atomic lock protects this wrapper against another concurrent invocation.
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.close(fd)
    result = None
    try:
        result = validate(packet_dir, run_dir, checkpoint_path, args.target_step, args.every)
        configuration_before = (run_dir / "configuration.json").read_bytes()
        with backup.open("xb") as stream:
            stream.write(configuration_before)
        result["configuration_before_backup"] = backup.name
        with tempfile.TemporaryDirectory(prefix="p1_safe_resume_") as temporary:
            engine = Path(temporary) / "evolve_candidate"
            compile_cmd = [args.compiler, "-std=c++17", "-O3", "-fopenmp", str(packet_dir / "evolve_candidate.cpp"), "-o", str(engine)]
            result["compile_command"] = compile_cmd
            compile_run = subprocess.run(compile_cmd, capture_output=True, text=True)
            result.update(compile_command=compile_cmd, compile_stdout=compile_run.stdout,
                          compile_stderr=compile_run.stderr, compile_returncode=compile_run.returncode)
            require(compile_run.returncode == 0, "compilation failed; see saved failure receipt")
            # Recheck immediately before an engine capable of appending starts.
            validate(packet_dir, run_dir, checkpoint_path, args.target_step, args.every)
            command = [str(engine), str(run_dir), str(args.target_step), format(result["dt"], ".17g"),
                       str(result["H"]), *[format(g, ".17g") for g in PROFILE], str(args.every), str(checkpoint_path)]
            env = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1")
            run = subprocess.run(command, capture_output=True, text=True, env=env)
            result.update(engine_sha256=sha256(engine), compile_command=compile_cmd, engine_stdout=run.stdout,
                          engine_stderr=run.stderr, engine_returncode=run.returncode)
            require(run.returncode == 0, "engine continuation failed; preserve partial output and diagnose before retrying")
        final_checkpoint = run_dir / f"checkpoint_{args.target_step}.bin"
        result["postvalidation"] = validate(packet_dir, run_dir, final_checkpoint,
                                            args.target_step + round(.2 / result["dt"]), args.every)
        result["postvalidation"].pop("target_step")  # Only validate the achieved state/history here.
        configuration_after = json.loads(configuration_before)
        configuration_after.update(target_final_step=args.target_step, target_time=args.target_step * result["dt"],
                                   binary_sha256=result["engine_sha256"])
        with temporary_config.open("x") as stream:
            stream.write(json.dumps(configuration_after, indent=2) + "\n")
        temporary_config.replace(run_dir / "configuration.json")
        result["configuration_after_sha256"] = sha256(run_dir / "configuration.json")
        receipt.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
    except BaseException as exc:
        failure = dict(result or {})
        failure.update(status="FAILED", reason=str(exc),
                       recovery="Preserve all partial histories and this receipt; diagnose before any retry. No automatic rollback is attempted.")
        receipt.write_text(json.dumps(failure, indent=2) + "\n")
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    try:
        main()
    except (InvalidResume, OSError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(json.dumps({"status": "REJECTED", "reason": str(exc)}), file=sys.stderr)
        sys.exit(2)
