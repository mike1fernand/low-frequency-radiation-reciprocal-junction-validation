#!/usr/bin/env python3
"""Verify and materialize the unchanged nested R3 evidence archive."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile
import zipfile

ARCHIVE_NAME = 'PAPER1_COMPLETE_DOMAIN_VALIDATION_R3.zip'
ARCHIVE_SHA256 = 'd1b3dcac19297522b6e271711a956fce24f681662c7c978c25fb679289039c88'
BASELINE_NAME = 'paper1-complete-domain-validation-r3'
PACKAGE_ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    hasher=hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda:source.read(1024*1024),b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def safe_relative(name):
    path=PurePosixPath(name)
    if not name or '\\' in name or path.is_absolute() or '..' in path.parts:
        raise ValueError('unsafe relative path: '+repr(name))
    return path


def verify_manifest(root, manifest_path):
    root=Path(root).resolve()
    manifest=json.loads(Path(manifest_path).read_text())
    files=manifest.get('files')
    if not isinstance(files,dict) or not files:
        raise ValueError('manifest requires a nonempty files hash mapping')
    for name,expected in files.items():
        relative=safe_relative(name)
        target=(root/relative).resolve()
        if root not in target.parents or not target.is_file() or sha(target)!=expected:
            raise ValueError('manifest mismatch: '+name)
    return {'manifest_sha256':sha(manifest_path),'verified_payloads':len(files)}


def bootstrap(package_root=PACKAGE_ROOT, archive=None):
    package_root=Path(package_root).resolve()
    archive=Path(archive) if archive is not None else package_root/'baseline'/ARCHIVE_NAME
    if sha(archive)!=ARCHIVE_SHA256:
        raise ValueError('inherited archive SHA-256 mismatch')
    destination=package_root/'output_r3'
    baseline=destination/BASELINE_NAME
    with zipfile.ZipFile(archive) as z:
        members=z.infolist()
        names=[x.filename for x in members]
        if len(names)!=len(set(names)):
            raise ValueError('duplicate archive member')
        for member in members:
            relative=safe_relative(member.filename)
            if not relative.parts or relative.parts[0]!=BASELINE_NAME:
                raise ValueError('unexpected inherited archive prefix')
            if stat.S_ISLNK(member.external_attr>>16):
                raise ValueError('symlink in inherited archive')
        if not baseline.exists():
            destination.mkdir(parents=True,exist_ok=True)
            # Extraction is checkpointed as one directory rename. A failed
            # attempt cannot leave an apparently complete baseline directory.
            temporary=Path(tempfile.mkdtemp(prefix='.r3-extract-',dir=destination))
            try:
                z.extractall(temporary)
                verify_manifest(temporary/BASELINE_NAME,temporary/BASELINE_NAME/'MANIFEST.json')
                os.replace(temporary/BASELINE_NAME,baseline)
            finally:
                shutil.rmtree(temporary,ignore_errors=True)
    if baseline.is_symlink() or baseline.resolve().parent!=destination.resolve():
        raise ValueError('baseline directory escapes expected location')
    verification=verify_manifest(baseline,baseline/'MANIFEST.json')
    return {'status':'passed','archive_sha256':ARCHIVE_SHA256,
            'archive_bytes':archive.stat().st_size,'baseline':str(baseline),
            **verification}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package-root',type=Path,default=PACKAGE_ROOT)
    p.add_argument('--archive',type=Path)
    a=p.parse_args()
    print(json.dumps(bootstrap(a.package_root,a.archive),indent=2))
