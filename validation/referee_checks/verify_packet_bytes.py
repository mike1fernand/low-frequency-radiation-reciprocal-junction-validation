"""Compare regenerated primary packet files to their recorded SHA-256 digests."""
import argparse
import hashlib
import json
from pathlib import Path

VALIDATION = Path(__file__).resolve().parents[1]
RAMP_NAMES = (
    "line.bin", "corner.bin", "snap_times.bin", "snap.bin",
    "checkpoint_0.bin", "checkpoint_1650.bin",
)
REFERENCE_NAMES = {
    "trajectory.csv": "rows", "snapshots.bin": "snapshots",
    "final.bin": "final_checkpoint",
}


def expected_digests(case):
    if case == "reference":
        receipt = VALIDATION / "reference/checks/REFERENCE_FULL_REPLAY.json"
        digests = json.loads(receipt.read_text(encoding="utf-8"))["raw_file_sha256"]
        return {name: digests[key] for name, key in REFERENCE_NAMES.items()}
    manifest = VALIDATION / "ramp/RAW_DELIVERY_SCOPE.json"
    rows = json.loads(manifest.read_text(encoding="utf-8"))["files"]
    digests = {entry["path"]: entry["sha256"] for entry in rows}
    return {name: digests["candidate_dt004_box80/" + name] for name in RAMP_NAMES}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("reference", "ramp"))
    parser.add_argument("directory", type=Path, help="directory containing the raw files")
    args = parser.parse_args()
    expected = expected_digests(args.case)
    for name, digest in expected.items():
        path = args.directory / name
        if not path.is_file():
            parser.exit(1, f"Missing raw packet file: {path}\n")
        with path.open("rb") as source:
            actual = hashlib.file_digest(source, "sha256").hexdigest()
        if actual != digest:
            parser.exit(1, f"SHA-256 mismatch in {path}\nExpected {digest}\nActual   {actual}\n")
        print(f"{actual}  {name}")
    print(f"PASS. All {len(expected)} {args.case} primary packet files match the archived bytes.")


if __name__ == "__main__":
    main()
