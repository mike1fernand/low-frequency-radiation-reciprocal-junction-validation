#!/usr/bin/env python3
"""Check exact saved inventory and frozen scientific input bindings.

This is a fast saved-evidence check, not numerical proof re-execution. Use
replay.py for separate proof and numerical modes. --files-only is intended for
integrity fixtures and makes no scientific claim. Checks remain active with -O.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys

MANIFEST_NAME = "PACKAGE_SHA256SUMS.json"
R4_MANIFEST_PATH = "r4_certificate/MANIFEST.json"
R4_MANIFEST_SHA256 = "e9acff005ca2620b8f435bbf7243226c3075fed6077699b623ee2b6b482a582f"
R4_PAYLOAD_COUNT = 682  # These plus the original pinned manifest = 683 files.
CHAIN_PATH = "tolerance/VALIDATED_PROOF_CHAIN.json"
CHAIN_SHA256 = "2e6479af3c223dd7568be5adc06f48d3a16f98675b715dee4306b72e0296a6a1"
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")

class VerificationError(Exception):
    """A failed acceptance condition, never disabled by Python -O."""

def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)

def no_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key: " + repr(key))
        result[key] = value
    return result

def load_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=no_duplicate_keys)
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise VerificationError("Invalid JSON in " + str(path) + ": " + str(error)) from error

def safe_name(name: object) -> str:
    require(isinstance(name, str), "Manifest paths must be strings.")
    require(bool(name) and "\\" not in name and not any(ord(c) < 32 for c in name),
            "Invalid manifest path: " + repr(name))
    path = PurePosixPath(name)
    require(not path.is_absolute() and path.as_posix() == name and name != "."
            and all(part not in ("", ".", "..") for part in name.split("/"))
            and re.match(r"^[A-Za-z]:", name) is None,
            "Noncanonical, absolute, or traversing manifest path: " + repr(name))
    return name

def checked_manifest(value: object) -> dict[str, str]:
    require(isinstance(value, dict) and bool(value), "Manifest must be a nonempty object.")
    result = {}
    for name, digest in value.items():
        name = safe_name(name)
        require(isinstance(digest, str) and SHA256_PATTERN.fullmatch(digest) is not None,
                "Invalid SHA-256 digest for " + name)
        result[name] = digest
    return result

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def inventory(root: Path) -> set[str]:
    """Reject symlinks, including directory links, without following them."""
    result = set()
    def walk_error(error: OSError) -> None:
        raise VerificationError("Cannot enumerate repository: " + str(error))
    for directory, dirs, files in os.walk(root, followlinks=False, onerror=walk_error):
        for name in dirs + files:
            path = Path(directory) / name
            mode = path.lstat().st_mode
            relative = path.relative_to(root).as_posix()
            require(not stat.S_ISLNK(mode), "Symlink rejected: " + relative)
            require(stat.S_ISDIR(mode) or stat.S_ISREG(mode),
                    "Nonregular filesystem entry rejected: " + relative)
            if stat.S_ISREG(mode):
                result.add(safe_name(relative))
    return result

def verify_files(root: Path) -> tuple[dict[str, str], set[str]]:
    require(root.is_dir(), "Repository root is not a directory: " + str(root))
    actual = inventory(root)
    require(MANIFEST_NAME in actual, "Root inventory manifest is missing.")
    manifest = checked_manifest(load_json(root / MANIFEST_NAME))
    require(MANIFEST_NAME not in manifest,
            "Root manifest must not list itself; it is the sole inventory exclusion.")
    expected = set(manifest) | {MANIFEST_NAME}
    missing, extra = sorted(expected - actual), sorted(actual - expected)
    require(not missing and not extra,
            "Inventory mismatch: " + json.dumps({"missing": missing, "extra": extra}))
    for name, expected_digest in manifest.items():
        require(sha256(root / name) == expected_digest, "Payload SHA-256 mismatch: " + name)
    return manifest, actual

def verify_r4(root: Path, manifest: dict[str, str], actual: set[str]) -> dict[str, object]:
    require(manifest.get(R4_MANIFEST_PATH) == R4_MANIFEST_SHA256,
            "Original full R4 manifest is missing or changed.")
    saved = load_json(root / R4_MANIFEST_PATH)
    require(isinstance(saved, dict) and saved.get("schema") == "paper1-r4-export-manifest-v1",
            "Unexpected original R4 manifest schema.")
    files = checked_manifest(saved.get("files"))
    require(len(files) == R4_PAYLOAD_COUNT and "MANIFEST.json" not in files,
            "Original R4 manifest must cover exactly 682 payloads.")
    expected = {"r4_certificate/" + name for name in files} | {R4_MANIFEST_PATH}
    r4_actual = {name for name in actual if name.startswith("r4_certificate/")}
    require(r4_actual == expected, "Full original R4 file inventory is not preserved.")
    for name, digest in files.items():
        require(manifest.get("r4_certificate/" + name) == digest,
                "Original R4 payload changed: " + name)
    return {"regular_files": len(expected), "pinned_manifest_sha256": R4_MANIFEST_SHA256,
            "status": "PASS_ORIGINAL_R4_PRESERVED"}

def verify_bindings(root: Path, manifest: dict[str, str]) -> dict[str, object]:
    require(manifest.get(CHAIN_PATH) == CHAIN_SHA256, "Frozen proof-chain record changed.")
    chain = load_json(root / CHAIN_PATH)
    require(isinstance(chain, dict) and
            chain.get("status") == "PASS_RIGOROUS_INFINITE_LATTICE_CANDIDATE_RESPONSE_CERTIFICATES",
            "Saved proof-chain status is not the frozen passing status.")
    bindings = chain.get("files")
    require(isinstance(bindings, dict) and len(bindings) == 19,
            "Expected all 19 frozen proof-chain bindings.")
    bound_paths = set()
    for name, binding in bindings.items():
        require(isinstance(binding, dict), "Malformed proof binding: " + str(name))
        path = safe_name(binding.get("path"))
        digest = binding.get("sha256")
        require(isinstance(digest, str) and SHA256_PATTERN.fullmatch(digest) is not None,
                "Invalid proof-binding digest: " + str(name))
        require(path not in bound_paths, "Duplicate proof-binding path: " + path)
        bound_paths.add(path)
        require(manifest.get(path) == digest, "Proof binding mismatch: " + str(name))
    point = chain["nominal_infinite_lattice_interval"]
    box = chain["whole_six_motif_plus_minus5_percent_infinite_lattice_floor"]
    lower, upper = Fraction(point["lower_rational"]), Fraction(point["upper_rational"])
    floor = Fraction(box["lower_rational"])
    require(Fraction(732, 1000) < lower <= upper < Fraction(740, 1000),
            "Saved nominal interval does not imply 0.732 < T_qT < 0.740.")
    require(Fraction(69, 100) < floor <= 1,
            "Saved whole-motif box floor does not imply T_qT > 0.69.")
    upstream = chain["upstream_status"]
    for key in ("saved_ledger_reproduced_exactly", "physical_controls_replayed",
                "original_reference_bridge_replayed",
                "combined_replay_with_physical_and_bridge_status_certified"):
        require(upstream.get(key) is True, "Missing passing upstream saved status: " + key)
    require(upstream.get("saved_R4_contour_jobs_replayed") == 222,
            "Unexpected saved upstream replay count.")
    require(upstream.get("fresh_global_contour_integration") is False,
            "Saved-job check must not be labelled fresh global integration.")
    caps = chain["kernel_caps"]
    require(caps.get("actual_cached_G_target_match") is True
            and caps.get("actual_cached_Q_target_match") is True
            and Fraction(caps["Q_recentring_allowance"]) == 0,
            "Saved kernel-center correspondence is not preserved.")
    require(Fraction(caps["rounded_up_eG"]) == Fraction(1, 250)
            and Fraction(caps["rounded_up_eQ"]) == Fraction(1, 500),
            "Unexpected saved kernel-error caps.")
    point_receipt = load_json(root / "tolerance/POINT_RATIONAL_CERTIFICATE.json")
    box_receipt = load_json(root / "tolerance/BOX_RATIONAL_CERTIFICATE.json")
    require(point_receipt.get("status") == box_receipt.get("status") == "PASS",
            "Saved response certificate status failed.")
    for key in ("eG", "eQ", "lower_rational", "upper_rational"):
        require(point_receipt["tighter_conditional_interval"][key] == point[key],
                "Point certificate/chain disagreement: " + key)
    for key in ("eG", "eQ", "lower_rational"):
        require(box_receipt["tight_kernel_box_lower"][key] == box[key],
                "Box certificate/chain disagreement: " + key)
    require(box_receipt.get("vertex_count") == 64
            and box_receipt.get("uses_transmission_vertex_extrema_assumption") is False,
            "Saved whole-motif certificate method is not preserved.")
    return {"proof_chain_bindings": len(bindings),
            "nominal_interval_rational": [str(lower), str(upper)],
            "whole_motif_plus_minus5_percent_floor_rational": str(floor),
            "saved_upstream_jobs": 222}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent,
                        help="Repository or integrity-fixture directory.")
    parser.add_argument("--files-only", action="store_true",
                        help="Check only exact file inventory and hashes; no scientific pass.")
    args = parser.parse_args()
    try:
        require(not args.root.is_symlink(), "Repository root must not be a symlink.")
        root = args.root.resolve()
        manifest, actual = verify_files(root)
        result = {"status": "PASS_FILE_INTEGRITY_ONLY" if args.files_only
                  else "PASS_SAVED_BYTES_AND_SCIENTIFIC_BINDINGS",
                  "payloads": len(manifest), "inventory_exact": True,
                  "numerical_proofs_reexecuted": False,
                  "scientific_bindings_checked": not args.files_only,
                  "scope": "Saved bytes and bindings only; use replay.py for proof re-execution."}
        if args.files_only:
            result["scope"] = "File-integrity check only; no scientific certificate or result verified."
        else:
            result["full_r4_preservation"] = verify_r4(root, manifest, actual)
            result.update(verify_bindings(root, manifest))
        print(json.dumps(result, indent=2))
        return 0
    except (VerificationError, OSError, ValueError, TypeError, KeyError) as error:
        print(json.dumps({"status": "FAIL", "error": str(error)}, indent=2), file=sys.stderr)
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
