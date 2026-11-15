# Saved evidence audit and continuation

The package contains `analysis/full_domain_r4/` and one unchanged inherited
archive at `baseline/PAPER1_COMPLETE_DOMAIN_VALIDATION_R3.zip`. The bootstrap
checks that archive's SHA-256, extracts it atomically to `output_r3/`, and
verifies its internal manifest. The outer manifest must contain every shipped
payload, including the nested archive. Expanded R3 files need not be duplicated
in the export.

Use Python 3.12.14 with numpy 2.3.5, sympy 1.14.0, mpmath 1.3.0, and
python-flint 0.8.0. The producing FLINT library was 3.3.1. Read-only evidence
auditing records both the producing and auditing Python builds and permits a
different Python build if all source, reviewed python-flint/FLINT, certificate,
coverage and exact aggregate checks pass. Continuing an integration checkpoint
retains the stricter full Python `sys.version` identity requirement. Its exact
original string is present in every checkpoint. A fresh numerical run can use
a separately recorded runtime, but must write fresh state rather than editing
old producer identities. Audit-runtime metadata is kept in VERIFICATION.json;
the complete mathematical ledger, including all input identities, is compared
without excluding any fields.

From the outer package root, create an environment if needed:

```sh
python -m venv .venv
.venv/bin/python -m pip install -r analysis/full_domain_r4/export/requirements.txt
```

The default verification audits saved evidence and recomputes the global
moment and matrix-error ledger. It does not repeat the full quadrature:

```sh
.venv/bin/python analysis/full_domain_r4/export/verify_saved.py --output ../paper1_r4_replay
```

Output must be a new or empty directory outside the sealed package. Add
`--expected-ledger analysis/full_domain_r4/budget/FINAL_CAMPAIGN_LEDGER.json`
only if that exact selected final filename is present in the export; replace
it with the actual selected ledger path documented in the top-level README.
This option requires equality of the complete audited JSON result.

Add `--physical-check` for the compact controller's fresh whole-cell scalar
and vector integrations, restart checks, and corruption controls. Each
subprocess is timeout-bounded; inspect `VERIFICATION.json` for completion.
A passed verification with unresolved arithmetic budgets reproduces an
unresolved result and does not certify an exact-lattice floor. The floor's
PSD and bridge gates are separate.

When the shipped complete-domain ledger passes both exact budgets with zero
fallback area, audit the evidence and replay the original bridge in one command:

```sh
.venv/bin/python analysis/full_domain_r4/export/verify_saved.py --output ../paper1_r4_bridge_replay --bridge
```

The optional bridge runs in the external output directory only after the
audited ledger and any requested saved-ledger comparison pass. The verifier
binds the returned original proof to that same global ledger, source, selected
targets and exact-Q PSD proof, and checks the strict rational floor contract.
Only successful completion sets `exact_lattice_floor_certified` to true.
Failed prerequisites, subprocesses or malformed results leave a terminal
failed verification record. Without `--bridge`, passing arithmetic budgets
remain separate from a certified floor. This option does not repeat the global
quadrature campaign.

For bounded continuation, copy the package to a working directory, bootstrap
the baseline, and invoke the same environment that produced the checkpoints:

```sh
.venv/bin/python analysis/full_domain_r4/export/bootstrap.py
.venv/bin/python analysis/full_domain_r4/run_campaign.py --workers 4 --stage-seconds 20 --max-seconds 300
```

Repeat the scheduler command to continue. It stops launching new work after
the wall budget and retains atomic per-tile checkpoints. Existing complete
states are skipped operationally and must still pass the final ledger audit.
Do not regenerate the initial campaign after settings were refined: use the
saved settings and manifest. If the original runtime is unavailable, retain
the sealed states and start a fresh replication in a separate work copy with
an empty `campaign/states` directory and fresh run logs.

The derivation, compact arithmetic and coverage proof are in
`integrator/COMPACTION_PROOF.md`. Packaging details and inspected path
assumptions are in `integrator/EXPORT_NOTES.md`.
