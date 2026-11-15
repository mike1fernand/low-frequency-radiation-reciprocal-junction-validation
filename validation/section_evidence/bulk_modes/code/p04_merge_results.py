from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = [args.records / "P04_CORE_CHECKS.json"] + [
        args.records / f"P04_SHELL_BRANCH_{branch}.json" for branch in range(3)
    ]
    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    status = "passed" if all(payload.get("status") == "passed" for payload in payloads) else "failed"
    output = {
        "schema": "paper1.p04.bulk-modes-symmetry-audit.v1",
        "status": status,
        "parts": [
            {"file": path.name, "sha256": sha256(path), "status": payload.get("status")}
            for path, payload in zip(paths, payloads)
        ],
        "core": payloads[0],
        "shells": payloads[1:],
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "parts": output["parts"]}, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
