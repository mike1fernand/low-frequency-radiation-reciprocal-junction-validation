"""Check the archived SM2 finite-grid result and rerun its original P04 producer."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
FIELDS = (
    "V_leak",
    "U_leak",
    "G_leak",
    "G_commutator",
    "full_reduced_field_max",
    "effective_self_energy_max",
)
CAP = 5e-16


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_metrics(path):
    data = json.loads(path.read_text())
    require(data["status"] == "passed" and data["summary"] == {"passed": 15, "total": 15},
            "The original P04 suite did not pass all 15 checks")
    raw = data["metrics"]["reduction"]["32"]
    result = {key: raw[key] for key in (*FIELDS, "orth")}
    require(all(math.isfinite(v) and v >= 0 for v in result.values()),
            "A reduction metric is negative or nonfinite")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-only", action="store_true",
                        help="check the original records without rerunning the producer")
    args = parser.parse_args()
    provenance = json.loads((HERE / "SOURCE.json").read_text())
    for name, expected in provenance["files"].items():
        require(sha256(HERE / name) == expected, f"Archive member digest mismatch: {name}")

    core = HERE / "records/P04_CORE_CHECKS.json"
    archived = read_metrics(core)
    require(all(value < CAP for value in archived.values()), "Archived SM2 threshold failed")
    independent = json.loads((HERE / "records/ARCHIVED_SECTIONIV_INDEPENDENT_CHECKS_v1.json").read_text())
    green_match = independent["reduction"]["G_match"]
    require(green_match == 1.6311799907932728e-15,
            "Unexpected independent Green-matrix comparison")
    require(all(archived[key] == independent["reduction"][key]
                for key in ("orth", "V_leak", "U_leak", "G_leak", "G_commutator")),
            "Shared diagnostics in the independent record disagree")
    report = {"archived_N32_six_relative_residuals": {key: archived[key] for key in FIELDS},
              "archived_basis_orthonormality_absolute_Frobenius_error": archived["orth"],
              "archived_independent_G_match": green_match,
              "fresh_P04": None}

    if not args.records_only:
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "P04_CORE_CHECKS.json"
            env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
            subprocess.run([sys.executable, str(HERE / "code/p04_core_checks.py"), "--output", str(out)],
                           env=env, stdout=subprocess.DEVNULL, check=True)
            fresh = read_metrics(out)
            # The manuscript reports the archived run. Last bits may differ across
            # floating-point libraries. The original producer checks its own margins.
            report["fresh_P04"] = {
                "passed_checks": 15,
                "total_checks": 15,
                "output_byte_identical_to_original": sha256(out) == sha256(core),
                "six_relative_residuals": {key: fresh[key] for key in FIELDS},
                "basis_orthonormality_absolute_Frobenius_error": fresh["orth"],
                "all_seven_strictly_below_archived_display_bound": all(x < CAP for x in fresh.values()),
            }
    print(json.dumps({"status": "PASS", **report}, indent=2))


if __name__ == "__main__":
    main()
