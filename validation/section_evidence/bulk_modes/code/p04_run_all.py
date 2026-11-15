from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def run(command: list[str], cwd: Path) -> None:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    result = subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    code = root / "code"
    records = root / "records"
    records.mkdir(parents=True, exist_ok=True)
    run([sys.executable, str(code / "p04_core_checks.py"), "--output", str(records / "P04_CORE_CHECKS.json")], root)
    for branch in range(3):
        run(
            [
                sys.executable,
                str(code / "p04_shell_check.py"),
                "--branch",
                str(branch),
                "--output",
                str(records / f"P04_SHELL_BRANCH_{branch}.json"),
            ],
            root,
        )
    run(
        [
            sys.executable,
            str(code / "p04_merge_results.py"),
            "--records",
            str(records),
            "--output",
            str(records / "P04_FRESH_AUDIT_RESULTS.json"),
        ],
        root,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
