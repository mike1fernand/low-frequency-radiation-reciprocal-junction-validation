# Paper One: full validation repository for the approximately 73.6% candidate

**Public edition.** The R4 certificate capsule and all 19 proof-chain bindings retain their historical hashes. The current public package inventory includes the scientific parameter-comparison code, packet-producer sources and supporting data. It omits private editorial records and links. For exact differences from the November 8 archive, see `../../PUBLIC_SOURCE_PROVENANCE.json`. Dated preservation manifests and delivery receipts describe the earlier archive. Use `verify_repository.py` with the adjacent `PACKAGE_SHA256SUMS.json` for this edition.


Original scientific release dated 8 November 2026. This self-contained replay release supersedes the earlier compact-core-plus-two-source-parts download instructions. No separately downloaded baseline archive is needed. The original scientific results and manuscript are unchanged.

## What the evidence establishes

The passive harmonic Paper One contact has six complete eight-spring motifs with per-link coefficients `g = (1,2,3,4,5,6)/7`, at structure-time frequency `omega = 6/5`, with scalar-line stiffness `169/100` and the same unit-mass braced central-spring bulk.

| Statement | Result | Scope |
|---|---:|---|
| Nominal numerical carrier estimate | 0.7358742883948081 | Approximate value, about 73.6%. |
| Rigorous nominal containing interval | 73.2166%–73.9652% | Outward display; does not certify one-decimal rounding. |
| Continuous six-motif ±5% floor | > 69% | Eight links within each motif vary together. |
| Direct broadband primary packet mean | 0.7164564184629316 | Finite open evolution; invariant lower-branch allocation midpoint. |
| Corresponding stationary packet prediction | 0.7163064862075107 | Finite-spectrum weighting, open-line energy normalization. |

The rigorous response certificates are tied to the authenticated real-frequency kernel bounds. Saved baseline replay, independent fresh quadrature and packet evolution have separate roles. No claim of a global optimum, bandwidth guarantee, 48 independent spring errors, electromagnetic identification, experimental efficiency or microscopic DST derivation is made here. Read `AUDIT_AND_CORRECTIONS.md` for the audit, corrections and limits.

## First check after extracting

Requires Python 3.12 for the recorded environment; the integrity checker itself uses only the standard library. Run from this extracted directory:

```sh
python verify_repository.py
python test_integrity_verifier.py
```

The first command checks exact file inventory, hashes and saved proof bindings. It does not re-execute the numerical proofs. It uses explicit errors and remains effective under `python -O`. The negative tests create temporary fixtures rather than changing the evidence. A manifest detects accidental corruption; it is not an external digital signature.

Use the root runner for replay. Every output path below must be a NEW directory outside this repository. It checks frozen evidence before and after, runs historical code with optimization disabled, and leaves work, outputs and failure logs in the external directory. Do not run legacy producers directly in the frozen evidence tree.

## Environments and replay commands

Create separate Python environments as needed. Package installation requires your ordinary package source; no scientific input download is required. C++ packet evolution and guard execution additionally require `g++` with C++17 and OpenMP. Tested native binary layout is little-endian float64/int32 Linux. Numerical bitwise identity across platforms is not promised.

| Task | Requirements file | Command, replacing `/path/to/new-*` |
|---|---|---|
| Exact-rational candidate point and box proofs | `environments/proof-requirements.txt` | `python replay.py response --output-dir /path/to/new-response` |
| Full saved R4 certificate replay, physical checks and bridge | `environments/proof-requirements.txt` | `python replay.py baseline --output-dir /path/to/new-baseline` |
| Fourteen-node packet spectrum, five saved kernel families | `environments/spectrum-requirements.txt` | `python replay.py spectrum --output-dir /path/to/new-spectrum` |
| Relocated fresh N48 quadrature pilot | `environments/fresh-requirements.txt` | `python replay.py quadrature --mode pilot --output-dir /path/to/new-pilot` |
| Primary packet evolution and analysis | `environments/packet-requirements.txt` | `python replay.py packet --case base --output-dir /path/to/new-packet` |
| Primary packet plus both frozen-window controls | `environments/packet-requirements.txt` | `python replay.py packet --case all --output-dir /path/to/new-packet-all` |
| Safe-resume guard adverse tests and small-box continuation | Python standard library and g++/OpenMP | `python replay.py resume-tests --output-dir /path/to/new-guard-tests` |
| Dated scientific report rebuild with bundled fonts | `environments/report-requirements.txt` | `python replay.py report --output-dir /path/to/new-report` |

For example, create a proof environment with `python -m venv /path/to/proof-env`, then install its requirements using `/path/to/proof-env/bin/python -m pip install -r environments/proof-requirements.txt`. Invoke the runner with that environment's Python. The fresh environment uses mpmath 1.4.1; the baseline proof uses mpmath 1.3.0. They are intentionally separate.

The baseline replay reproduces 222 saved compact states and the ledger, then performs fresh compact physical/restart/corruption checks and the bridge. It does not reintegrate the 449,698 contour leaves. The inherited upstream archive and all 683 R4 files are bundled. The new point and continuous-box response verifiers use exact rational arithmetic around the same authenticated target arrays.

Fresh quadrature also offers `--mode refine`, `--mode independent-tetra`, and `--mode fresh-tetra-n48`; each requires its own new output directory. `fresh_kernel/REPRODUCE_FRESH.md` gives their exact scope. These are corroborating numerical cases, not rigorous error enclosures or a profile search.

The full three-case packet replay writes approximately 3 GB of generated traces and checkpoints plus temporary analysis data; allow several GB of disk and memory for Fourier postprocessing. Two-thread producer runs were tens of seconds each in the recorded environment; postprocessing adds time. Cases run sequentially. For staged continuation after regenerating histories, use `packet_guard/safe_resume.py` as documented in `packet_guard/README.md`. A checkpoint alone is insufficient for the standard historical analyzer.

## Contents and authority

| Path | Contents |
|---|---|
| `tolerance/` | Exact-rational point/box derivations, verifiers, target kernel, binding chain and candidate certificates. |
| `r4_certificate/` | All 683 upstream saved-certificate files and inherited R3 archive. |
| `baseline_audit/`, `bound_review/` | Earlier full-baseline and independent response reviews. |
| `fresh_kernel/` | Frozen independent quadrature inputs/results plus additive portable adapter. |
| `packet_spectrum/` | All fourteen frequency nodes and five saved kernel families, reference regressions, states and results. |
| `packet/` | Dynamics and analysis code, exact preparation, configurations, recorded results, raw inventory, primary τ66 checkpoint and original producer source. |
| `packet_review/`, `packet_guard/` | Independent packet review, 15 safe-resume controls and small self-contained history fixture. |
| `audits/` | Additional mathematical, packaging and portability audits. |
| `context/` | Reproducible numerical parts of the earlier coupling and selected-band applicability studies; see `context/README.md`. |
| `assets/fonts/`, `environments/` | Licensed report fonts and separate pinned requirements. |
| `history/before_full_audit/` | Original root documentation/checker/manifest/state, preserved byte-for-byte. |

Current authority is this README, `AUDIT_AND_CORRECTIONS.md`, `CURRENT_STATE.json` and `FULL_RELEASE_CHECKS.json`, together with the bound scientific source files. The PDF remains the original dated four-page scientific report. Its earlier core/source-parts delivery paragraph, root `COLD_PACKAGE_REPLAY.json`, `PACKAGING_CHECK.json` and `SOURCE_PARTS.json` are historical provenance, not instructions for the present release. The new full-release replay receipt supersedes their packaging scope.

Within packet results, `CURRENT_CADENCE_CHECKS.json` and `CADENCE_FIELD_CORRECTION.json` supersede legacy secondary cadence means. Primary means are unchanged. The initial failed 56-box attempt is retained. The 80-box adaptation was diagnostic-informed and coefficient-blind; only the subsequent control interval was frozen in advance. The lower-branch stored-energy diagnostic is not a local radiation-flux measurement.

## Generated raw files and preservation

All public replay entry points use the source and inputs supplied here. Of the 58 raw packet files inventoried (3,162,169,048 bytes total), it includes the primary 9,218,168-byte checkpoint and omits the other 57 regenerable files (3,152,950,880 bytes). `packet/RAW_DATA_MANIFEST.json` records their paths, sizes and SHA-256 values. All code, preparation and configurations needed to regenerate the scientific runs are included. One complete primary continuation checkpoint at τ66 is included, but preceding large histories are omitted. The guard correctly rejects that compact directory for history-based continuation until those histories are regenerated.

The historical `SOURCE_PRESERVATION.json` records how the November 8 release preserved its precursor. Changes specific to this public edition are listed in `../../PUBLIC_SOURCE_PROVENANCE.json`. `PACKAGE_SHA256SUMS.json` covers every regular file except itself. External output folders and installed dependencies are not part of the frozen repository. No manuscript or physics-model change was needed by this audit.
