#!/usr/bin/env python3
"""Cheap controls for the additive resume guard, without a new packet run."""
from pathlib import Path
import argparse
import importlib.util
import json
import shutil
import struct
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("safe_resume", ROOT / "safe_resume.py")
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet_dir.resolve()
    output = args.output.resolve()
    assert not output.exists(), "keep test receipts immutable"
    config = json.loads((packet / "PREDECLARED_PACKET_CONFIGURATION.json").read_text())
    config.update(bulk_shape=[120, 4, 4, 3], quarter_shape=[120, 2, 2, 3], target_final_step=20, target_time=.8)
    results = []
    with tempfile.TemporaryDirectory(prefix="p1_resume_guard_test_") as temporary:
        root = Path(temporary)
        valid = root / "valid"
        fixture = ROOT / "fixtures" / "resume_control_continuous"
        shutil.copytree(fixture, valid)
        (valid / "configuration.json").write_text(json.dumps(config))
        receipt = guard.validate(packet, valid, valid / "checkpoint_20.bin", 21, 10)
        results.append(dict(test="existing_complete_history_acceptance", status="PASS", validation=receipt))

        def negative(name, change, checkpoint="checkpoint_20.bin"):
            case = root / name
            shutil.copytree(valid, case)
            change(case)
            try:
                guard.validate(packet, case, case / checkpoint, 21, 10)
            except guard.InvalidResume as exc:
                results.append(dict(test=name, status="PASS", expected_rejection=str(exc)))
            else:
                raise AssertionError(f"negative fixture accepted: {name}")

        negative("stale_checkpoint_into_later_history", lambda p: None, "checkpoint_10.bin")
        negative("truncated_line", lambda p: (p / "line.bin").write_bytes((p / "line.bin").read_bytes()[:-8]))
        negative("missing_prior_corner_history", lambda p: (p / "corner.bin").unlink())

        def checkpoint_only(p):
            for name in guard.TRACE_NAMES:
                (p / name).unlink()
        negative("complete_checkpoint_without_prior_histories", checkpoint_only)

        def wrong_endpoint(p):
            data = bytearray((p / "corner.bin").read_bytes())
            data[-8:] = struct.pack("<d", 999.)
            (p / "corner.bin").write_bytes(data)
        negative("same_size_wrong_corner_endpoint", wrong_endpoint)

        def wrong_metadata(p):
            data = bytearray((p / "checkpoint_20.bin").read_bytes())
            struct.pack_into("<d", data, 28, .02)
            (p / "checkpoint_20.bin").write_bytes(data)
        negative("checkpoint_configuration_timestep_mismatch", wrong_metadata)

        def extra_snapshot(p):
            (p / "snap_times.bin").write_bytes(struct.pack("<d", 48.))
        negative("inconsistent_snapshot_time_history", extra_snapshot)

        def late_monitor(p):
            text = (p / "monitor.csv").read_text()
            (p / "monitor.csv").write_text(text + text.splitlines()[-1] + "\n")
        negative("duplicate_monitor_endpoint", late_monitor)

        # Synthetic structural fixture reaches the first saved bulk frame at
        # tau=48. It tests file-layout guards only, never a physical trajectory.
        snapshot = root / "synthetic_snapshot_layout"
        shutil.copytree(valid, snapshot)
        checkpoint_data = bytearray((snapshot / "checkpoint_20.bin").read_bytes())
        struct.pack_into("<q", checkpoint_data, 8, 1200)
        (snapshot / "checkpoint_1200.bin").write_bytes(checkpoint_data)
        last_line = (snapshot / "line.bin").read_bytes()[-1920:]
        last_corner = (snapshot / "corner.bin").read_bytes()[-288:]
        (snapshot / "line.bin").write_bytes(last_line * 1201)
        (snapshot / "corner.bin").write_bytes(last_corner * 1201)
        (snapshot / "snap.bin").write_bytes(checkpoint_data[guard.HEADER_SIZE:])
        (snapshot / "snap_times.bin").write_bytes(struct.pack("<d", 48.))
        lines = (snapshot / "monitor.csv").read_text().splitlines()
        last_fields = lines[-1].split(",")
        monitor = [lines[0]]
        for step in range(0, 1201, 300):
            monitor.append(",".join([str(step), repr(step * .04), *last_fields[2:]]))
        (snapshot / "monitor.csv").write_text("\n".join(monitor) + "\n")
        guard.validate(packet, snapshot, snapshot / "checkpoint_1200.bin", 1205, 300)
        results.append(dict(test="synthetic_aligned_bulk_snapshot_layout", status="PASS",
                            scope="Structural fixture, not an evolved state history."))
        for name, filename, change in [
            ("same_size_wrong_bulk_snapshot_endpoint", "snap.bin", lambda b: b[:-8] + struct.pack("<d", 999.)),
            ("same_size_wrong_snapshot_timestamp", "snap_times.bin", lambda b: struct.pack("<d", 48.2)),
        ]:
            case = root / name
            shutil.copytree(snapshot, case)
            path = case / filename
            path.write_bytes(change(path.read_bytes()))
            try:
                guard.validate(packet, case, case / "checkpoint_1200.bin", 1205, 300)
            except guard.InvalidResume as exc:
                results.append(dict(test=name, status="PASS", expected_rejection=str(exc)))
            else:
                raise AssertionError(f"negative fixture accepted: {name}")

        try:
            guard.validate(packet, snapshot, snapshot / "checkpoint_1200.bin", 1203, 300)
        except guard.InvalidResume as exc:
            results.append(dict(test="off_snapshot_grid_production_target", status="PASS", expected_rejection=str(exc)))
        else:
            raise AssertionError("off-grid production target accepted")

        failed_compile = root / "failed_compile"
        shutil.copytree(valid, failed_compile)
        before_failure = {name: guard.sha256(failed_compile / name) for name in guard.TRACE_NAMES}
        failure_command = [sys.executable, str(ROOT / "safe_resume.py"), "--packet-dir", str(packet),
                           "--run-dir", str(failed_compile), "--checkpoint", str(failed_compile / "checkpoint_20.bin"),
                           "--target-step", "21", "--every", "10", "--execute", "--compiler", "__p1_missing_cpp_compiler_for_test__"]
        failure_run = subprocess.run(failure_command, capture_output=True, text=True)
        assert failure_run.returncode == 2
        failure_receipt = json.loads((failed_compile / "SAFE_RESUME_TO_21.json").read_text())
        assert failure_receipt["status"] == "FAILED" and failure_receipt["compile_command"][0] == "__p1_missing_cpp_compiler_for_test__"
        assert before_failure == {name: guard.sha256(failed_compile / name) for name in guard.TRACE_NAMES}
        results.append(dict(test="compilation_failure_preserves_receipt_and_histories", status="PASS",
                            expected_rejection=failure_receipt["reason"]))

        # One extra small-box step exercises the valid append path; this is not
        # a fresh physical packet simulation or a new accuracy requirement.
        command = [sys.executable, str(ROOT / "safe_resume.py"), "--packet-dir", str(packet),
                   "--run-dir", str(valid), "--checkpoint", str(valid / "checkpoint_20.bin"),
                   "--target-step", "21", "--every", "10", "--execute"]
        run = subprocess.run(command, capture_output=True, text=True)
        assert run.returncode == 0, run.stderr + run.stdout
        execute = json.loads(run.stdout)
        config_after = json.loads((valid / "configuration.json").read_text())
        assert config_after["target_final_step"] == 21 and config_after["target_time"] == 21 * .04
        assert config_after["binary_sha256"] == execute["engine_sha256"]
        assert (valid / execute["configuration_before_backup"]).read_bytes() == json.dumps(config).encode()
        reference = fixture / "checkpoint_20.bin"
        assert guard.sha256(reference) == receipt["checkpoint_sha256"], "archived checkpoint was altered"
        results.append(dict(test="valid_one_step_staged_append", status="PASS", result=execute))
    result = dict(status="PASS", tests=results, test_count=len(results),
                  guard_source_sha256=guard.sha256(ROOT / "safe_resume.py"),
                  test_source_sha256=guard.sha256(ROOT / "test_safe_resume.py"),
                  fixture_provenance_sha256=guard.sha256(ROOT / "fixtures" / "FIXTURE_PROVENANCE.json"),
                  scope="Twelve expected rejections/failure controls, two structural positive checks and one H=2 step from the preserved state, including exact prior-configuration backup and final-target update. Synthetic snapshot fixtures check layout only. No long dynamics; original science files unchanged.")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "test_count": len(results), "output": str(output)}))


if __name__ == "__main__":
    main()
