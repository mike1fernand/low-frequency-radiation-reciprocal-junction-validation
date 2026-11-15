"""Compare staged real-axis and Cartesian estimates of the 24-motif upstream control."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[3]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def stage_file(source, destination, expected):
    if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Unexpected digest for {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def compare(actual, expected, label):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(f"Unexpected keys in {label}")
        for key, value in expected.items():
            compare(actual[key], value, f"{label}/{key}")
    elif not math.isclose(actual, expected, rel_tol=2e-6, abs_tol=1e-8):
        raise ValueError(f"{label}: computed {actual}, expected {expected}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("records", "fresh"))
    parser.add_argument("output", type=Path, help="New external working directory")
    parser.add_argument("--workers", type=int, choices=range(1, 5), default=1)
    args = parser.parse_args()
    out = args.output.resolve()
    if out == REPOSITORY or REPOSITORY in out.parents or out.exists():
        parser.error("Choose a new output directory outside the repository.")
    provenance = read(HERE / "SOURCES.json")
    for name, digest in provenance["original_sources"].items():
        if name != "reconstruct_cartesian.py":  # That historical producer uses frequency 2.9.
            stage_file(HERE.parent / name, out / name, digest)
    stage_file(HERE / "reconstruct_cartesian.py", out / "reconstruct_cartesian.py",
               provenance["changed_source"]["sha256"])
    shutil.copyfile(HERE / "check_results.py", out / "check_results.py")
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1",
               PYTHONOPTIMIZE="0", PYTHONDONTWRITEBYTECODE="1")

    def run(script, *arguments):
        name = "-".join([Path(script).stem, *(str(arg) for arg in arguments)])
        with (out / f"{name}.log").open("w", encoding="utf-8") as log:
            subprocess.run([sys.executable, "-B", str(out / script), *map(str, arguments)],
                           cwd=out, env=env, check=True, stdout=log, stderr=subprocess.STDOUT)

    if args.mode == "records":
        for directory in ("cartesian_quadrant", "realaxis"):
            shutil.copytree(HERE / "records" / directory, out / directory)
    else:
        for n in (256, 384, 512):
            step = n // 16
            def chunk(start):
                run("reconstruct_cartesian.py", n, start, start + step)
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                list(pool.map(chunk, range(0, n, step)))
            run("reconstruct_cartesian.py", n, "combine")
        for n in (512, 1024, 1536, 2048):
            run("reconstruct_realaxis.py", 2.6, n)
        run("reconstruct_realaxis.py", 2.9, 768)
    run("check_results.py")
    actual = read(out / "method-comparison.json")
    compare(actual, read(HERE / "records/method-comparison.json"), "method-comparison")
    print(f"PASS. {args.mode} comparison: {out / 'method-comparison.json'}")
    print("2.6 archived shell and stored real part:", actual["historical_control"])
    print("2.6 Cartesian 512 cubed, cubic fit to finite-width powers:",
          actual["combined_finite_powers"]["512_3"])
    print("2.9 real-axis Nperp=768:", actual["realaxis_all"]["2.90_N768"])


if __name__ == "__main__":
    main()
