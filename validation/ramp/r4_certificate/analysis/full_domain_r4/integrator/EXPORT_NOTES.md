# Portable export strategy and executable replay instructions

These are packaging instructions for the final R4 assembler. The frozen
integrator producer is unchanged. The export itself and its final cold replay
must be performed after the live campaign is paused or finished so that saved
states, aggregate results, and the final manifest describe one snapshot.

## Layout that preserves current relative paths

Use this outer package layout:

```text
README.md
MANIFEST.json
requirements.txt
baseline/PAPER1_COMPLETE_DOMAIN_VALIDATION_R3.zip
analysis/full_domain_r4/...
```

Include the R4 sources, settings, complete and explicitly partial checkpoints,
CAMPAIGN.json, actual run records, exact budget inputs/results, review records,
and reports. Exclude virtual environments, __pycache__, *.pyc, temporary files,
and any unrelated upload credentials/receipts. Include one original R3 ZIP;
do not also package the expanded R3 evidence tree. Generated replay outputs
belong outside the sealed outer package and must not silently replace saved
evidence. An outer manifest should hash every shipped payload except itself,
including the untouched nested archive. A verification report created after
sealing should sit beside the archive rather than require a recursive manifest.

The checked inherited archive has:

| Property | Verified value |
|---|---|
| Name | PAPER1_COMPLETE_DOMAIN_VALIDATION_R3.zip |
| Bytes | 53,490,635 |
| SHA-256 | d1b3dcac19297522b6e271711a956fce24f681662c7c978c25fb679289039c88 |
| File members | 235 |
| Uncompressed file bytes | 342,354,460 |
| Sole member prefix | paper1-complete-domain-validation-r3 |

This ZIP already contains the inherited R2 evidence and selected repository
provenance required by R3. No separate 295 MB validation repository archive
is needed for the R4 coverage, source, moment, or matrix-budget checks. This
does not assert that R4 replaces the paper's entire original validation repo.

## Identified path and mutation assumptions

| Component | Path behavior | Export action |
|---|---|---|
| build_campaign.py | R3 is `ROOT.parents[1]/output_r3/paper1-complete-domain-validation-r3`; writes CAMPAIGN.json and settings below its own directory | Preserve `analysis/full_domain_r4` depth and extract R3 beneath the outer root's `output_r3`; run regeneration in a replay copy |
| CAMPAIGN.json | Baseline field is `../../output_r3/paper1-complete-domain-validation-r3`; state/settings names are relative to R4 root | Preserve directory layout and relative names |
| run_campaign.py | Invokes its own interpreter via `sys.executable`; writes states, RUN_LOG.jsonl, LAST_RUN.json beneath its own root | Invoke the intended environment's interpreter directly; run in an isolated replay copy |
| integrator/run_stage.py | Settings and output supplied explicitly | Safe for a small physical replay to an external output path |
| compact_controller.py | Source filenames are beside the module; source identity hashes content, not absolute paths | Copy these files byte-for-byte; moving the tree does not alter their identity |
| compact controller resume | Requires identical source, settings, Python `sys.version`, python-flint version and endpoint certificate | If runtime differs, do a fresh run in an empty output location; do not edit checkpoint identities |
| budget/matrix_budget.py and budget_ledger.py | Target defaults beside script; explicit input/output CLI | Ship TARGETS.json and use explicit external output |
| review/check_new_arithmetic.py | Finds R4 root from its own path; writes review/NEW_ARITHMETIC_CHECKS.json | Execute in the replay copy to preserve original evidence |
| inherited R3 verify.py | Reads baseline root and accepts explicit external replay directory | It can replay the full inherited suite after extraction; no original production simulation is repeated |

`run_campaign.py` skips states marked complete during scheduling. Its skip
decision is not the source/coverage/aggregate audit: the final aggregation
adapter must audit every selected state before accepting it. Scheduler
progress percentages are operational diagnostics, not proof of full-domain
completion.

The final budget adapter CLI is `budget/campaign_ledger.py --root <R4root>
--output <fresh-ledger.json>`, with optional `--manifest`. It takes the
baseline from CAMPAIGN.json and calls inherited_reader.py in a subprocess
to avoid R3/R4 module-name collisions. It does not write to the R3 baseline.
The portable default verifier is `export/verify_saved.py`; it audits saved
evidence and recomputes the ledger without repeating global quadrature.

## Safe inherited-archive bootstrap

Run from the outer package root. This exact command verifies the archive hash,
rejects duplicate/absolute/traversal/symlink entries, extracts only to a new
baseline directory, and checks the extracted R3 manifest. It permits an
existing extracted baseline only after rechecking the entire R3 manifest.

```sh
python - <<'PY'
import hashlib, json, stat, zipfile
from pathlib import Path, PurePosixPath

root = Path.cwd()
archive = root/'baseline/PAPER1_COMPLETE_DOMAIN_VALIDATION_R3.zip'
expected = 'd1b3dcac19297522b6e271711a956fce24f681662c7c978c25fb679289039c88'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == expected
destination = root/'output_r3'
baseline = destination/'paper1-complete-domain-validation-r3'
with zipfile.ZipFile(archive) as z:
    members = z.infolist()
    names = [m.filename for m in members]
    assert len(names) == len(set(names)), 'duplicate ZIP member'
    for member in members:
        p = PurePosixPath(member.filename)
        assert not p.is_absolute() and '..' not in p.parts
        assert p.parts and p.parts[0] == 'paper1-complete-domain-validation-r3'
        assert '\\' not in member.filename
        assert not stat.S_ISLNK(member.external_attr >> 16), 'symlink ZIP member'
    if not baseline.exists():
        destination.mkdir(exist_ok=True)
        z.extractall(destination)
manifest = json.loads((baseline/'MANIFEST.json').read_text())
for name, sha in manifest['files'].items():
    p = PurePosixPath(name)
    assert not p.is_absolute() and '..' not in p.parts
    assert hashlib.sha256((baseline/name).read_bytes()).hexdigest() == sha, name
print('Inherited archive and all manifested R3 payloads verified.')
PY
```

If extraction is interrupted, remove only that generated baseline directory
and rerun bootstrap; do not accept a partial tree or overwrite any source
archive. The final exporter should perform bootstrap inside its temporary
cold-replay directory, not inside the sealed source package.

## Environment and bounded replay commands

The producing runtime recorded Python 3.12.14, numpy 2.3.5, python-flint 0.8.0,
FLINT 3.3.1, sympy 1.14.0 and mpmath 1.3.0. A pinned requirements file should
contain the existing R3 requirements. Read-only audit records both producer and audit Python builds, requires the
reviewed python-flint/FLINT versions and unchanged sources, and rechecks the
certificates, exact aggregates and coverage. Integration resume separately
retains strict producer-runtime identity. Creating an environment is separate
from validating compatibility with either operation. Use its Python
entry path directly; resolving a virtual environment's Python symlink can
accidentally invoke the system interpreter without that environment's packages.

From an isolated copy of the outer package after successful bootstrap:

```sh
python -m venv .venv
.venv/bin/python -m pip install -r output_r3/paper1-complete-domain-validation-r3/requirements.txt
mkdir replay
.venv/bin/python analysis/full_domain_r4/budget/check_budget_ledger.py --output replay/budget_checks.json
.venv/bin/python analysis/full_domain_r4/integrator/check_compact.py --reference output_r3/paper1-complete-domain-validation-r3/quadrature/qt_chart1_body_preconditioned_gl5.json --output replay/compact_checks.json
.venv/bin/python analysis/full_domain_r4/review/check_new_arithmetic.py
```

The compact checker performs three fresh body-cell integrations (interrupted
and direct scalar runs, then heterogeneous vector targets), each split into
bounded calls if needed. Individual stage limits are not a claim that the
whole checking command will finish within one stage under heavy system load.
Check its explicit success JSON and exit status; a process launch is not a
successful validation. The checkpoint identity of these *fresh* runs uses
the current runtime and does not rewrite inherited producer identities.

An independent small fresh physical replay uses a new output file:

```sh
.venv/bin/python analysis/full_domain_r4/integrator/run_stage.py --settings analysis/full_domain_r4/integrator/BODY_SAMPLE_SETTINGS.json --output replay/fresh_body.json --max-seconds 20
```

Repeat that command if its status is partial. Do not label a partial result a
replayed whole cell. Compare its entire exact outward moment intervals with
the saved result; timing and stage hashes need not match. The existing compact
checker already tests whole-cell overlap and exact equality across interrupted
versus uninterrupted runs under one runtime.

For the optional full inherited suite:

```sh
.venv/bin/python output_r3/paper1-complete-domain-validation-r3/verify.py --output replay/inherited_r3 --physical-replay
```

The R3 verifier uses separately bounded subprocesses. Its historical fixed
per-stage wall timeout may need an explicitly documented rerun on a less
loaded system; a timeout is not a physics failure or proof of completion.
Do not edit its sealed source to conceal an unsuccessful replay.

## Campaign continuation and final acceptance

Copy the sealed package to an external work directory before continuing.
Preserve the stated layout and all inherited checkpoints. Run:

```sh
.venv/bin/python analysis/full_domain_r4/run_campaign.py --workers 4 --stage-seconds 20 --max-seconds 300
```

If the exact checkpoint producer runtime is unavailable, create a fresh work
copy with an empty `campaign/states` directory and no old run logs. Preserve
the original states in the sealed package; this new run is a replication,
not continuation of the old source/runtime chain. Do not run build_campaign
blindly after tolerance refinements: it regenerates its original exploratory
settings and may supersede refined settings. Use the final saved campaign
and settings as the authoritative restart inputs.

Before claiming the export passes, the final cold replay must verify the
outer manifest, inherited archive and its manifest, all selected compact
source/settings/certificate/hash/coverage/endpoint sums, exact physical
coverage including symmetry multiplicities and Green exterior, final
moment-to-matrix assembly, cache-center displacement, both scaled budgets,
and the separate PSD/bridge premises required for any exact-lattice floor.
If a budget remains unresolved, a correctly reproduced unresolved result is
a valid replay outcome but is not a certified floor.

No cold export replay or final full-domain campaign is claimed by this note.
This file documents the inspected interfaces, verified nested-archive
metadata, executable commands, and the checks the final assembler must run.

Read-only audit and resume are separate operations. Export verification now
permits a different Python build for the read-only audit, records producer and
auditor builds in VERIFICATION.json, and retains reviewed python-flint/FLINT
requirements and all exact checks. The mathematical ledger is runtime-neutral;
its complete JSON, including all state/source/target identities, is compared
without excluding fields. Strict numerical-checkpoint resume is unchanged.
This interface change has not itself established a second-platform replay.
