"""Isolated interpreter adapter to the unchanged R3 evidence auditor."""
from pathlib import Path
import argparse
import hashlib
import json
import sys


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();sys.path.insert(0,str(root))
    from complete_ledger import read_evidence
    from matrix.matrix_budget import json_exact
    path=root/'EVIDENCE_MANIFEST.json'
    manifest=json.loads(path.read_text())
    records=[read_evidence(root,x) for x in manifest['evidence']]
    print(json.dumps(json_exact({'manifest_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                                'records':records})))


if __name__=='__main__':main()
