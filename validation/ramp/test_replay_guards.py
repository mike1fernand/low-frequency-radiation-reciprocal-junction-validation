#!/usr/bin/env python3
"""Test replay-runner control flow with subprocess mocks, without physics runs.

Verification subprocesses are deliberately mocked PASS. These tests exercise
error receipts, selected-window checks and isolation; they do not authenticate
scientific payloads or replace verify_repository.py or a numerical replay.
"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("frozen_replay_runner", ROOT / "replay.py")
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)

PACKET_CASES = {
    "candidate_dt004_box80": "CANDIDATE_PACKET_BOX80_RESULT.json",
    "control_dt002_box80": "CONTROL_DT002_BOX80_RESULT.json",
    "control_dt004_box96": "CONTROL_DT004_BOX96_RESULT.json",
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def invoke(argv, fake_subprocess):
    output = StringIO()
    with patch.object(sys, "argv", [str(ROOT / "replay.py"), *argv]), \
            patch.object(replay.subprocess, "run", side_effect=fake_subprocess), redirect_stdout(output):
        code = replay.main()
    return code, output.getvalue()


def controls():
    records = []
    with tempfile.TemporaryDirectory(prefix="p1_replay_guard_controls_") as temporary:
        base = Path(temporary)

        def child_case(name, behavior, expect_status, postcheck_failure=False):
            out = base / name
            verification_count = 0

            def fake(argv, **kwargs):
                nonlocal verification_count
                if any(str(a).endswith("verify_repository.py") for a in argv):
                    verification_count += 1
                    failed = postcheck_failure and verification_count == 2
                    return subprocess.CompletedProcess(argv, 1 if failed else 0,
                                                       "FAKE_VERIFICATION_FAILURE" if failed else "FAKE_VERIFICATION_PASS", "")
                if behavior == "interrupt":
                    raise KeyboardInterrupt("fixture interruption")
                return subprocess.CompletedProcess(argv, 7 if behavior == "fail" else 0,
                                                   "MOCK_CHILD_STDOUT", "MOCK_CHILD_STDERR")

            code, stdout = invoke(["response", "--output-dir", str(out)], fake)
            receipt = json.loads((out / "REPLAY_RECEIPT.json").read_text())
            require(receipt["status"] == expect_status, f"wrong status: {name}")
            require(code == (0 if expect_status == "PASS" else 1), f"wrong return code: {name}")
            require(json.loads(stdout)["status"] == expect_status, f"stdout/receipt mismatch: {name}")
            require(verification_count == 2, f"pre/post checks not both invoked: {name}")
            if behavior == "interrupt":
                require("KeyboardInterrupt" in receipt["error"], "interruption was not retained")
                require(not receipt["commands"], "interrupted child falsely recorded as completed")
            else:
                require((out / "response.stdout").read_text() == "MOCK_CHILD_STDOUT", "child stdout lost")
                require((out / "response.stderr").read_text() == "MOCK_CHILD_STDERR", "child stderr lost")
                require(len(receipt["commands"]) == 1, "child command not recorded")
            records.append(dict(test=name, status="PASS", expected_receipt_status=expect_status))

        child_case("successful_child", "success", "PASS")
        child_case("failed_subprocess_preserves_logs", "fail", "FAIL")
        child_case("keyboard_interrupt_never_reports_pass", "interrupt", "FAIL")
        child_case("failed_postcheck_overrides_success", "success", "FAIL", postcheck_failure=True)

        for name, out in [("existing_output_rejected", base / "already_exists"),
                          ("inside_frozen_source_rejected", ROOT / "__never_create_replay_guard_test_output__")]:
            if name == "existing_output_rejected":
                out.mkdir()
            else:
                require(not out.exists(), "inside-source fixture name already exists")

            def should_not_run(*args, **kwargs):
                raise AssertionError("subprocess invoked before output-isolation rejection")

            try:
                invoke(["response", "--output-dir", str(out)], should_not_run)
            except RuntimeError as exc:
                require("NEW output" in str(exc) or "outside the frozen" in str(exc), "unexpected refusal")
            else:
                raise AssertionError(f"unsafe output accepted: {name}")
            if name == "inside_frozen_source_rejected":
                require(not out.exists(), "refused inside-source path was created")
            records.append(dict(test=name, status="PASS"))

        for bad_case in [None, *PACKET_CASES]:
            name = "all_packet_controls_pass" if bad_case is None else f"false_predicate_rejected_{bad_case}"
            out = base / name
            analyses_seen = []
            evolution_arguments = []

            def packet_fake(argv, **kwargs):
                strings = [str(a) for a in argv]
                if any(a.endswith("analyze_candidate.py") for a in strings):
                    case = Path(strings[-1]).name
                    analyses_seen.append(case)
                    data = json.loads((ROOT / "packet" / PACKET_CASES[case]).read_text())
                    data = copy.deepcopy(data)
                    # Keep PASS, window and numeric means identical. The boolean
                    # must independently prevent a false replay PASS.
                    data["frozen_window_all_predicates_pass"] = case != bad_case
                    return subprocess.CompletedProcess(argv, 0, json.dumps(data), "")
                if any(a.endswith("run_control.py") for a in strings):
                    evolution_arguments.append(strings)
                return subprocess.CompletedProcess(argv, 0, "FAKE_CHILD_OR_VERIFICATION_PASS", "")

            code, stdout = invoke(["packet", "--case", "all", "--output-dir", str(out)], packet_fake)
            receipt = json.loads((out / "REPLAY_RECEIPT.json").read_text())
            expected = "PASS" if bad_case is None else "FAIL"
            require(receipt["status"] == expected and code == (0 if expected == "PASS" else 1), name)
            require(json.loads(stdout)["status"] == expected, "packet stdout/receipt mismatch")
            if bad_case is None:
                require(analyses_seen == list(PACKET_CASES), "--case all omitted a packet case")
                require(len(receipt["task_result"]["comparisons"]) == 3, "not all means compared")
                require("--frozen-window" not in evolution_arguments[0], "primary selector was frozen after the fact")
                for argv in evolution_arguments[1:]:
                    i = argv.index("--frozen-window")
                    require(argv[i+1:i+3] == ["57.4", "64.2"], "control window changed")
            else:
                require(analyses_seen[-1] == bad_case, "wrong packet case caused rejection")
                require("selected-window packet predicate failed" in receipt["error"], "boolean predicate not load-bearing")
            records.append(dict(test=name, status="PASS", expected_receipt_status=expected,
                                analyzed_cases=analyses_seen))

        out = base / "report_inputs_are_staged"

        def report_fake(argv, **kwargs):
            return subprocess.CompletedProcess(argv, 0, "FAKE_CHILD_OR_VERIFICATION_PASS", "")

        code, _ = invoke(["report", "--output-dir", str(out)], report_fake)
        receipt = json.loads((out / "REPLAY_RECEIPT.json").read_text())
        require(code == 0 and receipt["status"] == "PASS", "report staging rejected")
        for relative in ["tolerance/VALIDATED_PROOF_CHAIN.json", "fresh_kernel/SUMMARY.json",
                         "packet_spectrum/PACKET_SPECTRUM.json",
                         *["packet/" + filename for filename in PACKET_CASES.values()],
                         "assets/fonts/DejaVuSans.ttf", "assets/fonts/DejaVuSans-Bold.ttf"]:
            require((out / "report" / relative).is_file(), "report input not staged: " + relative)
        staged = (out / "report" / "make_validation_report.py").read_text()
        require("'/usr/share/fonts/truetype/dejavu/" not in staged, "system font dependency retained")
        require(len(receipt["task_result"]["io_adapter"]["exact_substitutions"]) == 2, "font substitutions not exactly recorded")
        records.append(dict(test="report_inputs_and_exact_font_adapter_are_staged", status="PASS",
                            scope="Checks staging only; does not render the report."))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New external test receipt path")
    args = parser.parse_args()
    output = args.output.resolve()
    require(not output.exists(), "do not overwrite a test receipt")
    require(ROOT not in output.parents, "test receipt must be outside the frozen repository")
    replay_hash = digest(ROOT / "replay.py")
    tests = controls()
    require(replay_hash == digest(ROOT / "replay.py"), "runner source changed during controls")
    result = dict(status="PASS", tests=tests, test_count=len(tests),
                  replay_source_sha256=replay_hash, test_source_sha256=digest(Path(__file__)),
                  scope="Standard-library subprocess mocks verify runner control flow only. Verification is mocked PASS; no scientific calculation, certificate replay or payload authentication is performed.")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "test_count": len(tests), "output": str(output)}))


if __name__ == "__main__":
    main()
