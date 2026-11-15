# Paper 1 complete-domain validation R4

The complete-domain enclosures and original perturbation bridge certify **T_qT > 1037/10000 = 0.1037** at omega = 1.2 in the declared exact infinite lattice, throughout the six independent complete-motif coupling-factor intervals [0.95, 1.05]. The bare media and geometry are fixed; each segment's eight links vary together. The nominal approximately 12.6% monochromatic coefficient and separate packet result retain their numerical status.

All 222 production jobs completed with zero unresolved fallback area. Certified scaled operator-error upper bounds, rounded upward, are g <= 0.003903452 < 0.020 and Q <= 0.001966909 < 0.005. The simpler direct g bound also passes. See the exact rational ledger and linked bridge below for the controlling quantities and proof premises.

This is a versioned validation addition with a revised manuscript candidate. It includes the unchanged sealed R3 archive and does not overwrite the recovered original repository. There has been no GitHub or Zenodo publication. The earlier finite-matrix producer intentionally retains its conditional conclusion; the linked R4 bridge discharges its integration-error premises.

## Read and reproduce

All paths below are relative to this archive's root. Start with:

- `analysis/full_domain_r4/RESULT.md`: result, scope, corrected issues, and unperformed checks.
- `analysis/full_domain_r4/PLAN.md`: reviewed implementation plan and checking gates.
- `analysis/full_domain_r4/budget/FINAL_CAMPAIGN_LEDGER.json`: exact complete-domain matrix budgets.
- `analysis/full_domain_r4/bridge_replay/EXACT_LATTICE_BRIDGE.json`: original proof replay and input bindings.
- `analysis/full_domain_r4/review/FINAL_REVIEW.md`: independent targeted scientific review.
- `analysis/full_domain_r4/manuscript/candidate/`: revised TeX sources and figures.
- `analysis/full_domain_r4/manuscript/build/candidate/`: final main and supplement PDFs.
- `analysis/full_domain_r4/manuscript/PROVENANCE.json` and `APPLIED_PATCH.diff`: parent/candidate identities and exact source changes.
- `analysis/full_domain_r4/manuscript/BUILD_RECORDS.json` and `qa/QA_RESULTS.json`: build and visual review evidence; main 18 pages, supplement 11, unchanged.
- `analysis/full_domain_r4/export/README.md`: runtime, evidence auditing and continuation details.

Create a Python 3.12 environment and install the pinned requirements (reviewed producer: Python 3.12.14, python-flint 0.8.0, FLINT 3.3.1):

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r analysis/full_domain_r4/export/requirements.txt
.venv/bin/python -B analysis/full_domain_r4/export/verify_saved.py --output ../paper1_r4_replay --expected-ledger analysis/full_domain_r4/budget/FINAL_CAMPAIGN_LEDGER.json --physical-check --bridge
```

The output directory must be new or empty and outside the package. The verifier checks the outer manifest, bootstraps and verifies all inherited R3 payloads, audits every saved production state, reconstructs the complete matrix ledger, requires exact equality to the selected final ledger, runs fresh whole-cell/restart/corruption controls, and replays the original bridge. Each computational subprocess is bounded to at most 55 seconds by default. A successful run must report `status: passed`, `saved_ledger_reproduced_exactly: true`, and `exact_lattice_floor_certified: true` in the external `VERIFICATION.json`.

This audit does not repeat all 449,698 quadrature leaves. The full production campaign was executed once. Campaign settings, compact states, producer code, exact coverage, logs and dependencies are included for fresh reproduction or documented continuation. Use a working copy for any new computation; never change old producer identities. Historical interim ledgers and proposal records are preserved as history; the final ledger, bridge and applied-patch provenance named above control this release.

The included REVTeX tree supports manuscript compilation with an installed TeX distribution, `pdflatex`, `kpsewhich`, and Poppler tools. Rebuilding can change PDF metadata and byte hashes; compare sources, page content and build diagnostics. Perform builds in a working copy:

```sh
python analysis/full_domain_r4/build_manuscript.py --include-parent
```

## Scope and delivery checks

The certificate is for the specified reference-frequency harmonic model and motifwise coupling box. It is not a certified packet average, arbitrary bondwise perturbation result, whole-band theorem, or experiment. A stronger finite-window/asymptotic identity and a rigorous band extension remain future work. No second-platform or independent-person execution was performed.

The companion `PAPER1_R4_EXPORT_VERIFICATION.json` identifies the exact ZIP hash and records the actual fresh-extraction replay and PDF identity checks performed for delivery. It is separate to avoid making the sealed archive hash depend on a report containing that same hash. The companion continuation report records the completed delivery outcome. The intermediate `WORK_IN_PROGRESS.json` inside the sealed source tree records the pre-export checkpoint only.
