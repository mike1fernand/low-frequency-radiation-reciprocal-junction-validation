"""Check each delivered payload against the package manifest."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
manifest = root / 'PACKAGE_SHA256SUMS.json'
expected = json.loads(manifest.read_bytes())['files']
actual = {}
for path in root.rglob('*'):
    relative = path.relative_to(root)
    if path.is_file() and path != manifest and not {'.git', '__pycache__'}.intersection(relative.parts):
        actual[relative.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
if actual != expected:
    added_or_missing = sorted(actual.keys() ^ expected.keys())
    changed = sorted(key for key in actual.keys() & expected.keys() if actual[key] != expected[key])
    raise SystemExit(f'FAIL. Added or missing payloads: {added_or_missing}. Changed: {changed}.')
print(json.dumps({'status': 'PASS', 'payloads': len(actual), 'scope': 'Delivered payload inventory and file hashes'}, indent=2))
