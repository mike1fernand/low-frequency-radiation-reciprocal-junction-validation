"""Rerun the November significance calculations in an external working copy."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT.parents[1]) or output in ROOT.parents:
        parser.error("Choose a new directory outside the repository.")
    output.mkdir(parents=True, exist_ok=False)
    work = output / "working-copy"
    shutil.copytree(ROOT, work)
    checks = json.loads((work / "tasks.json").read_text())
    env = dict(os.environ, PYTHONOPTIMIZE="0", PYTHONDONTWRITEBYTECODE="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    completed = []
    for task in checks:
        script = work / task["script"]
        result = script.with_name(task["output"])
        expected = json.loads(result.read_text())
        result.unlink()
        cmd = [sys.executable, "-B", str(script)]
        if task["output_argument"]:
            cmd.extend(["--output", str(result)])
        with (output / (script.stem + ".log")).open("w") as log:
            subprocess.run(cmd, cwd=script.parent, env=env, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        actual = json.loads(result.read_text())
        if not actual or actual.get("status") != expected.get("status"):
            raise ValueError(f"Unexpected result status from {script.name}")
        if script.stem == "check_two_profiles":
            if [s["name"] for s in actual["states"]] != [s["name"] for s in expected["states"]]:
                raise ValueError("Profile identities differ from the archived result")
            for a, b in zip(actual["states"], expected["states"]):
                for field in ["Rq", "Tline", "TqT", "TqL", "bulk_qT_share"]:
                    if not math.isclose(a[field], b[field], abs_tol=1e-11, rel_tol=0):
                        raise ValueError(f"{a['name']} {field} differs from the archived result")
            if actual["direction_survives"] is not True:
                raise ValueError("Directional tradeoff failed")
        completed.append(task["script"])
        print("PASS", task["script"], flush=True)
    (output / "replay.json").write_text(json.dumps({"completed": completed,
        "tests": "Original producer checks. Two-profile comparison additionally checked against archived power fractions."}, indent=2) + "\n")


if __name__ == "__main__":
    main()
