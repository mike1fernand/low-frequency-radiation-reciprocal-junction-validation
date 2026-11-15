"""Verify the supplied release, optionally checking all its staged Git files."""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
MANIFEST = "REPOSITORY_SHA256SUMS.json"


def check_index(paths):
    """Check actual staged blob identities, including in subdirectory installs."""
    try:
        listing = subprocess.check_output(
            ["git", "ls-files", "--stage", "-z", "--", "."], cwd=ROOT
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError("--git-index requires Git and a local checkout.") from exc
    entries = {}
    for record in listing.split(b"\0"):
        if record:
            header, name = record.split(b"\t", 1)
            mode, oid, stage = header.split()
            if stage == b"0":
                entries[os.fsdecode(name)] = (mode, oid.decode("ascii"))
    missing = sorted(paths - entries.keys())
    if missing:
        raise RuntimeError(
            "Missing from the Git index (stage the release and check ignore rules):\n"
            + "\n".join(missing)
        )
    for name in sorted(paths):
        mode, oid = entries[name]
        if mode not in (b"100644", b"100755"):
            raise RuntimeError("Expected a regular file in the Git index: " + name)
        path = ROOT / name
        algorithm = {40: "sha1", 64: "sha256"}.get(len(oid))
        if algorithm is None:
            raise RuntimeError("Unsupported Git object ID format for " + name)
        digest = hashlib.new(algorithm)
        digest.update(f"blob {path.stat().st_size}\0".encode("ascii"))
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
        if digest.hexdigest() != oid:
            raise RuntimeError(
                "Staged bytes differ from the verified working file: " + name
            )
    print(f"PASS. All {len(paths)} release files have identical bytes in the Git index.")


def check_ramp_table():
    """Compare the printed SM5D table's frequency/incident-power pairs with its data."""
    spectrum_path = ROOT / "validation/ramp/packet_spectrum/PACKET_SPECTRUM.json"
    spectrum = json.loads(spectrum_path.read_text(encoding="utf-8"))
    family, = (family for family in spectrum["families"] if family["label"] == "surface_N160")
    rows = family["rows"]
    expected = [(format(row["omega"], ".4f"), format(100 * row["states"]["ramp"]["qT"], ".1f")) for row in rows]

    source = (ROOT / "manuscript/source/paper1_SM5D.tex").read_text(encoding="utf-8")
    tables = re.findall(r"\\begin\{table\}.*?\\end\{table\}", source, re.DOTALL)
    table, = (table for table in tables if r"\label{tab:sm-ramp-frequency-samples}" in table)
    body = table.split(r"\midrule", 1)[1].split(r"\bottomrule", 1)[0]
    left, right = [], []
    for line in body.split(r"\\"):
        if line.strip():
            frequency_a, percent_a, frequency_b, percent_b = (cell.strip() for cell in line.split("&"))
            left.append((frequency_a, percent_a))
            right.append((frequency_b, percent_b))
    if len(expected) != 14 or left + right != expected:
        raise RuntimeError("The printed SM5D ramp table differs from the archived surface_N160 pairs.")
    print("PASS. All fourteen SM5D frequency/conversion pairs match the archived data at printed precision.")


COMPACT_REFERENCE = (
    "verify_incident_power.py", "verify_lowfreq.py", "verify_scattering.py",
    "verify_vi_results.py", "verify_viii_clarifications.py", "verify_angles_controls.py",
    "verify_clarity_observables.py", "verify_robustness.py", "verify_line_mass.py",
    "verify_scaling_and_projections.py", "verify_packet_comparison.py",
    "verify_packet_reference.py", "verify_ramp_packet.py", "verify_packet_spectrum.py",
)


def run_check(script, env, *arguments):
    print("Running " + script, flush=True)
    subprocess.run(
        [sys.executable, "-B", str(ROOT / script), *arguments],
        cwd=ROOT, env=env, check=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--git-index", action="store_true",
        help="also verify the staged copy of every release file after git add",
    )
    parser.add_argument(
        "--compact", action="store_true",
        help="also run reference response, saved packet, and independent analytic checks",
    )
    args = parser.parse_args()
    files = json.loads((ROOT / MANIFEST).read_text(encoding="utf-8"))["files"]
    for name, expected in files.items():
        path = ROOT / name
        if not path.is_file():
            raise RuntimeError("Missing packaged file: " + name)
        with path.open("rb") as source:
            if hashlib.file_digest(source, "sha256").hexdigest() != expected:
                raise RuntimeError("Changed packaged file: " + name)
    binding = json.loads((ROOT / "MANUSCRIPT_BINDING.json").read_text(encoding="utf-8"))
    for name, expected in binding["sources"].items():
        path = ROOT / "manuscript" / "source" / name
        with path.open("rb") as source:
            if hashlib.file_digest(source, "sha256").hexdigest() != expected:
                raise RuntimeError("Manuscript source changed: " + name)
    check_ramp_table()
    if args.git_index:
        check_index(set(files) | {MANIFEST})

    env = {
        **os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONOPTIMIZE": "0",
        "PYTHONUTF8": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    }
    table_check = "validation/referee_checks/verify_contact_table.py"
    run_check(table_check, env, *(["--replay"] if args.compact else []))
    for script in (
        "validation/reference/tools/verify_package.py",
        "validation/ramp/verify_repository.py",
        "validation/referee_checks/verify_reference_shell.py",
    ):
        run_check(script, env)
    print("PASS. Packaged bytes, current manuscript source binding, validation inventories and reference-shell localization.", flush=True)

    if args.compact:
        for script in COMPACT_REFERENCE:
            run_check("validation/reference/tools/" + script, env)
        run_check(
            "validation/reference/tools/verify_robustness_bounds.py", env,
            "--campaign", str(ROOT / "validation/ramp/r4_certificate/analysis/full_domain_r4"),
        )
        for script in ("reference_shell_algebra.py", "independent_shell.py", "check_force_moment.py", "verify_profile_reversal.py"):
            run_check("validation/referee_checks/" + script, env)
        with tempfile.TemporaryDirectory(prefix="paper1-directional-") as directory:
            run_check("validation/predictions/upstream/control_26/replay.py", env,
                      "records", str(Path(directory) / "replay"))
        print("PASS. Compact stationary, profile reversal, directional-control, saved-packet, bound and analytic checks.")
    print("RUNNING_VALIDATIONS.md lists the separate commands for fresh dynamics and contour-certificate replay.")


if __name__ == "__main__":
    main()
