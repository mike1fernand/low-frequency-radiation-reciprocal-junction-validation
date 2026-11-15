#!/usr/bin/env python3
"""Relocate frozen fresh-kernel checks without changing their scientific code.

Install beside the original fresh_kernel/SHA256SUMS.json. All old files remain
unchanged. Each command stages and runs one declared check in a new directory.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

MANIFEST_HASH = "e551893c6b73e200c0f4aa4262a7ac6f21f3314d886f69f3bd8bb71c5dc6d570"
DRIVER_HASHES = {
    "run_fresh_kernel.py": "4f9132bc190ccf64f07a8060eb08ee5f2bcbae1afa629bc0b50bb0903a9ec9c3",
    "run_fresh_tetra_N48.py": "7b8fef40d34eb9903bfa2906ce442295de0ac41009cc04beddc0a5e3f292d12a",
    "check_independent_tetra.py": "232d8edd91397e99f7fe8dc90d74b7607056e35d3b0ec44f32432f1d728a3b86",
}
OLD_PACKAGE = "PACKAGE=Path('/workspace/scratch/044aa7b17d80/p1_implemented_2026_11_08/revision/paper1-significance')"
OLD_R8 = "R8=Path('/workspace/scratch/044aa7b17d80/p1_736_validation_2026_11_08/r8_selected/paper1-validation-integrated-r8/legacy_rc4/validation/current_numerics')"
SURFACE_INPUTS = {
    "recheck/purity_inputs/model_api_v1_1.py": "inputs/model_api_v1_1.py",
    "recheck/purity_inputs/reduced_core.py": "inputs/reduced_core.py",
    "tradeoff_validation/inputs/refined_surface_pv/surface/src/on_shell_basis.py": "inputs/on_shell_basis.py",
    "tradeoff_validation/inputs/refined_surface_pv/surface/src/spectral_histogram.py": "inputs/spectral_histogram.py",
    "tradeoff_validation/inputs/refined_surface_pv/selected/REFERENCE_KERNEL.npz": "inputs/REFERENCE_KERNEL.npz",
    "tradeoff_validation/inputs/refined_surface_pv/selected/REFERENCE_KERNEL.json": "inputs/REFERENCE_KERNEL.json",
}
TETRA_INPUTS = {
    "analysis/resolution_execution/tetra/multifrequency_tetra.py": "fresh_tetra_N48/sources/analysis/resolution_execution/tetra/multifrequency_tetra.py",
    "analysis/math_resolution/tetra/combined_tetra_diagnostic.py": "fresh_tetra_N48/sources/analysis/math_resolution/tetra/combined_tetra_diagnostic.py",
    "analysis/math_resolution/tetra/stable_scalar.py": "fresh_tetra_N48/sources/analysis/math_resolution/tetra/stable_scalar.py",
    "analysis/math_resolution/tetra/centroid_weights.py": "fresh_tetra_N48/sources/analysis/math_resolution/tetra/centroid_weights.py",
}
FINE_INPUTS = {
    "analysis/math_resolution/tetra/" + name: "independent_tetra/inputs/" + name
    for name in ("combined_tetra_diagnostic.py", "stable_scalar.py", "centroid_weights.py")
}
for case, stem in (
    ("reference_N96", "omega_001"), ("reference_N128", "omega_001"),
    ("reference_N128_halfshift", "omega_001"), ("reference_N160", "omega_001"),
    ("reference_N192", "omega_000"),
):
    for name in ("CONFIG.json", stem + ".json", stem + ".npz"):
        FINE_INPUTS["analysis/resolution_execution/tetra/" + case + "/" + name] = (
            "independent_tetra/inputs/" + case + "/" + name
        )
MODES = {
    "pilot": ["run_fresh_kernel.py", "--phase", "pilot"],
    "refine": ["run_fresh_kernel.py", "--phase", "refine"],
    "independent-tetra": ["check_independent_tetra.py"],
    "fresh-tetra-n48": ["run_fresh_tetra_N48.py"],
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def verify_source(source: Path) -> dict[str, str]:
    manifest_path = source / "SHA256SUMS.json"
    if digest(manifest_path) != MANIFEST_HASH:
        raise ValueError("Frozen SHA256SUMS.json does not match the audited release.")
    manifest = json.loads(manifest_path.read_text())
    for name, expected in manifest.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Unsafe frozen manifest path: " + name)
        path = source / relative
        if not path.is_file() or digest(path) != expected:
            raise ValueError("Missing or changed frozen file: " + name)
    for name, expected in DRIVER_HASHES.items():
        if manifest.get(name) != expected:
            raise ValueError("Frozen driver hash mismatch: " + name)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=MODES, required=True)
    parser.add_argument("--output", type=Path, required=True,
                        help="New output directory; must not already exist.")
    parser.add_argument("--stage-only", action="store_true",
                        help="Stage verified copies, but do not run numerical code.")
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    output = args.output.resolve()
    if output.exists():
        parser.error("Output already exists. Choose a new directory; old results are preserved.")
    if output.is_relative_to(source) or source.is_relative_to(output):
        parser.error("Choose an output directory outside the frozen source tree.")
    manifest = verify_source(source)
    dependencies = {}
    if not args.stage_only:
        required = (("numpy", "scipy", "skimage", "lazy_loader")
                    if args.mode in ("pilot", "refine")
                    else ("numpy", "scipy", "mpmath"))
        for name in required:
            try:
                module = importlib.import_module(name)
            except ImportError as error:
                parser.error("Missing scientific dependency " + name + ": " + str(error)
                             + ". Install requirements-fresh-pinned.txt in the active environment.")
            dependencies[name] = {"version": getattr(module, "__version__", "unknown"),
                                  "module_file": getattr(module, "__file__", None)}
    # Prepare every bundled upstream input layout. This is copying only, and all
    # source and destination hashes are checked before any driver is executed.
    staged = output / "fresh_kernel"
    upstream_package = output / "upstream" / "paper1-significance"
    upstream_r8 = output / "upstream" / "current_numerics"
    staged.mkdir(parents=True)
    copied = []
    for target_root, mapping in ((upstream_package, SURFACE_INPUTS),
                                 (upstream_r8, TETRA_INPUTS),
                                 (upstream_r8, FINE_INPUTS)):
        for target_name, source_name in mapping.items():
            target = target_root / target_name
            expected = manifest[source_name]
            target.parent.mkdir(parents=True, exist_ok=True)
            # Two input layouts share three scientific modules; reject disagreement.
            if target.exists() and digest(target) != expected:
                raise ValueError("Conflicting bundled upstream input: " + str(target))
            shutil.copyfile(source / source_name, target)
            if digest(target) != expected:
                raise ValueError("Staged input copy mismatch: " + str(target))
            copied.append({"source": source_name, "target": str(target.relative_to(output)),
                           "sha256": expected})
    drivers = []
    for name in DRIVER_HASHES:
        original = (source / name).read_text()
        old, new = ((OLD_PACKAGE, "PACKAGE=Path(" + repr(str(upstream_package)) + ")")
                    if name == "run_fresh_kernel.py"
                    else (OLD_R8, "R8=Path(" + repr(str(upstream_r8)) + ")"))
        if original.count(old) != 1:
            raise ValueError("Expected exactly one frozen input-root assignment: " + name)
        adapted = original.replace(old, new, 1)
        if adapted.replace(new, old, 1) != original:
            raise ValueError("Unexpected non-I/O driver modification: " + name)
        target = staged / name
        target.write_text(adapted)
        drivers.append({"name": name, "original_sha256": DRIVER_HASHES[name],
                        "staged_sha256": digest(target), "old_assignment": old,
                        "new_assignment": new, "changed_lines": 1,
                        "scientific_code_unchanged": True})
    command = [sys.executable, str(staged / MODES[args.mode][0]), *MODES[args.mode][1:]]
    record = {
        "status": "STAGED", "mode": args.mode, "original_source": str(source),
        "frozen_manifest_sha256": MANIFEST_HASH, "verified_frozen_files": len(manifest),
        "drivers": drivers, "copied_inputs": copied, "command": command,
        "runtime": {"python": sys.version, "platform": platform.platform()},
        "scientific_dependencies": dependencies,
        "thread_environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
        "python_guard_environment": {"PYTHONOPTIMIZE": "0", "PYTHONDONTWRITEBYTECODE": "1"},
        "scope": "I/O-only relocation; fixed archived checks; no parameter search or stronger physical claim.",
    }
    write_json(output / "RELOCATION_MANIFEST.json", record)
    print("Verified and staged:", output, flush=True)
    if args.stage_only:
        return 0
    env = os.environ.copy()
    env.update(record["thread_environment"])
    env.update(record["python_guard_environment"])
    start = time.perf_counter()
    with (output / "EXECUTION.log").open("w") as log:
        completed = subprocess.run(command, cwd=staged, env=env,
                                   stdout=log, stderr=subprocess.STDOUT)
    # Verifying after the child exits also covers unsuccessful runs.
    verify_source(source)
    record.update({"status": "PASS" if completed.returncode == 0 else "FAIL",
                   "returncode": completed.returncode,
                   "seconds": time.perf_counter() - start,
                   "frozen_source_unchanged_after_execution": True,
                   "log": "EXECUTION.log"})
    write_json(output / "RELOCATION_MANIFEST.json", record)
    print(record["status"], "-- execution log:", output / "EXECUTION.log", flush=True)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
