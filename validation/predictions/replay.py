"""Stage original calculation scripts outside the release and check selected results."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compare(value, reference, label):
    if isinstance(reference, dict):
        require(isinstance(value, dict) and value.keys() == reference.keys(), label)
        for key in reference:
            compare(value[key], reference[key], f"{label}.{key}")
    elif isinstance(reference, list):
        require(isinstance(value, list) and len(value) == len(reference), label)
        for i, (a, b) in enumerate(zip(value, reference)):
            compare(a, b, f"{label}[{i}]")
    elif isinstance(reference, float):
        power_field = ".powers[" in label or ".leading_prediction[" in label
        absolute_tolerance = 1e-27 if power_field else 1e-11
        require(isinstance(value, (int, float)) and math.isfinite(value)
                and math.isclose(value, reference, rel_tol=1e-8, abs_tol=absolute_tolerance), label)
    else:
        require(value == reference, label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=("predictions", "upstream", "fresh-cartesian", "fresh-gb"))
    parser.add_argument("output", type=Path, help="New external directory for working files, logs and receipt")
    args = parser.parse_args()
    output = args.output.resolve()
    require(output != ROOT and ROOT not in output.parents and not output.exists(), "Use a new external directory")
    records = read(ROOT / "SOURCE_FILES.json")["files"]
    for name, record in records.items():
        require(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == record["sha256"], "Source changed: " + name)
    output.mkdir(parents=True)
    work = output / "work"
    for name in records:
        target = work / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONOPTIMIZE="0",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    executed = []
    verified = []

    def run(script, *arguments):
        command = [str(work / script), *map(str, arguments)]
        log = output / f"{len(executed) + 1:02d}_{Path(script).stem}.log"
        with log.open("w", encoding="utf-8") as stream:
            subprocess.run([sys.executable, "-B", *command], cwd=work, env=env,
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
        executed.append([script, *arguments])

    def compare_output(script, name, fields, *arguments):
        generated = work / name
        generated.unlink()
        run(script, *arguments)
        result = read(generated)
        expected = read(ROOT / name)
        for field in fields:
            compare(result[field], expected[field], f"{name}.{field}")
            verified.append(f"{name}.{field}")
        return result

    if args.profile == "predictions":
        compare_output("lowfreq_audit/audit_lowfreq.py", "lowfreq_audit/audit_results.json",
                       ["independent_point_force_coefficients", "cache_reproduction", "nonnegative_sparse_profile_checks"])
        compare_output("controls_audit/interference_reaudit.py", "controls_audit/interference_reaudit.json",
                       ["replayed_31_rows", "single_motif_checks"])
        compare_output("controls_audit/bandwidth/independent_reaudit.py", "controls_audit/bandwidth/INDEPENDENT_CHECKS.json",
                       ["surface_count", "tetra_count", "minimum_surface", "minimum_tetra", "all_declared_input_hashes_match"])
        compare_output("angular_audit/check_angular_independent.py", "angular_audit/INDEPENDENT_ANGULAR_CHECK.json",
                       ["cones", "local_coefficients"])
        compare_output("angular_audit/make_tradeoffs.py", "angular_audit/COLLECTION_TRADEOFFS.json", ["rows"])
        compare_output("collection_extension/cone_limits.py", "collection_extension/cone_w0p00312500_n24.json",
                       ["omega", "order", "source_input", "rows"], "--omega", ".003125", "--order", "24")
    elif args.profile == "fresh-gb":
        name = "gb/fresh_reference_results.json"
        (work / name).unlink()
        for pattern in ("fresh_G_*.npz", "fresh_residue_*.npz"):
            for kernel in (work / "gb").glob(pattern):
                kernel.unlink()
        run("gb/fresh_reference_check.py")
        expected, generated = read(ROOT / name), read(work / name)
        require(len(generated) == len(expected) == 3, "Expected three regulator-grid results")
        for row, old in zip(generated, expected):
            compare(row["Nperp"], old["Nperp"], "Transverse grid")
            compare(row["results"], old["results"], "Extrapolated results at " + str(old["Nperp"]))
            verified.append(f"reference regulator grid {row['Nperp']}")
    else:
        if args.profile == "fresh-cartesian":
            shutil.rmtree(work / "upstream/cartesian_quadrant/N384")
            for start in range(0, 384, 24):
                run("upstream/reconstruct_cartesian.py", 384, start, start + 24)
        name = "upstream/cartesian_quadrant/N384/results.json"
        (work / name).unlink(missing_ok=True)
        run("upstream/reconstruct_cartesian.py", 384, "combine")
        compare(read(work / name)["zero_width_fits_in_epsilon"], read(ROOT / name)["zero_width_fits_in_epsilon"], name)
        verified.append("three zero-width fits using all 16 Cartesian chunks")
        if args.profile == "upstream":
            compare_output("upstream/reconstruct_reflection_check.py", "upstream/reflection_check_results.json", ["rows"])
            r = compare_output("upstream/reconstruct_root_checks.py", "upstream/root_checks_results.json",
                               ["all_root_counts_match", "difficult_count", "random_count", "precision_decimal_digits"])
            require(len(r["records"]) == 80 and r["maximum_shell_eigenvalue_residual"] < 1.1e-7, "Sampled root residual")
            compare_output("upstream/reconstruct_realaxis.py", "upstream/realaxis/w2.90_N768_result.json",
                           ["T_minus", "T_plus", "simultaneous_T_minus", "simultaneous_T_plus"], 2.9, 768)
    receipt = {"status": "PASS", "profile": args.profile, "executed": executed, "verified": verified}
    (output / "REPLAY.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"PASS {args.profile}, {len(executed)} producer invocations. See {output / 'REPLAY.json'}")


if __name__ == "__main__":
    main()
