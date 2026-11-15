#!/usr/bin/env python3
"""Small adverse fixtures for verify_repository.py; never edits the repository."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

VERIFIER = Path(__file__).resolve().with_name("verify_repository.py")
MANIFEST = "PACKAGE_SHA256SUMS.json"

class IntegrityVerifierTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="p1-integrity-fixture-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / "one.txt").write_bytes(b"one\n")
        (self.root / "nested").mkdir()
        (self.root / "nested/two.bin").write_bytes(bytes(range(64)))
        self.manifest = {p.relative_to(self.root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in self.root.rglob("*") if p.is_file()}
        self.save_manifest()

    def save_manifest(self):
        (self.root / MANIFEST).write_text(json.dumps(self.manifest), encoding="utf-8")

    def run_verifier(self, optimized=False):
        command = [sys.executable]
        if optimized:
            command.append("-O")
        command += ["-B", str(VERIFIER), "--root", str(self.root), "--files-only"]
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        return subprocess.run(command, text=True, capture_output=True, env=env)

    def rejected(self, optimized=False):
        result = self.run_verifier(optimized)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertEqual(json.loads(result.stderr)["status"], "FAIL")

    def test_valid_inventory_is_only_a_file_integrity_pass(self):
        result = self.run_verifier()
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt = json.loads(result.stdout)
        self.assertEqual(receipt["status"], "PASS_FILE_INTEGRITY_ONLY")
        self.assertFalse(receipt["scientific_bindings_checked"])
        self.assertFalse(receipt["numerical_proofs_reexecuted"])

    def test_valid_inventory_with_optimized_python(self):
        self.assertEqual(self.run_verifier(optimized=True).returncode, 0)

    def test_manifest_omission(self):
        del self.manifest["one.txt"]
        self.save_manifest()
        self.rejected()

    def test_extra_file(self):
        (self.root / "extra.txt").write_text("untracked")
        self.rejected()

    def test_missing_file(self):
        (self.root / "one.txt").unlink()
        self.rejected()

    def test_payload_mutation(self):
        (self.root / "one.txt").write_text("changed")
        self.rejected()

    def test_payload_mutation_with_optimized_python(self):
        (self.root / "one.txt").write_text("changed")
        self.rejected(optimized=True)

    def test_file_symlink(self):
        (self.root / "one.txt").unlink()
        (self.root / "one.txt").symlink_to(self.root / "nested/two.bin")
        self.rejected()

    def test_directory_symlink(self):
        (self.root / "linked").symlink_to(self.root / "nested", target_is_directory=True)
        self.rejected()

    def test_traversal_path(self):
        self.manifest["../outside.txt"] = "0" * 64
        self.save_manifest()
        self.rejected()

    def test_absolute_path(self):
        self.manifest["/absolute.txt"] = "0" * 64
        self.save_manifest()
        self.rejected()

    def test_noncanonical_paths(self):
        for name in ("./one.txt", "nested//two.bin", "nested/../one.txt", "one.txt/",
                     "nested\\two.bin", "C:/one.txt"):
            with self.subTest(path=name):
                self.manifest[name] = "0" * 64
                self.save_manifest()
                self.rejected()
                del self.manifest[name]

    def test_duplicate_json_keys(self):
        digest = self.manifest["one.txt"]
        (self.root / MANIFEST).write_text('{"one.txt":"' + digest + '","one.txt":"' + digest + '"}')
        self.rejected()

    def test_empty_manifest(self):
        (self.root / MANIFEST).write_text("{}")
        self.rejected()

    def test_bad_digest(self):
        self.manifest["one.txt"] = "not-a-sha256"
        self.save_manifest()
        self.rejected()

    def test_manifest_self_listing(self):
        self.manifest[MANIFEST] = "0" * 64
        self.save_manifest()
        self.rejected()

if __name__ == "__main__":
    unittest.main(verbosity=2)
