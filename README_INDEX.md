# Guide to the README files

Start with [README.md](README.md), then [RUNNING_VALIDATIONS.md](RUNNING_VALIDATIONS.md). You do not need to read every specialist guide before running the basic check. [WINDOWS_GUIDE.md](WINDOWS_GUIDE.md) covers extraction/setup; [GITHUB_GUIDE.md](GITHUB_GUIDE.md) covers complete uploads.

**Repository root** means the extracted folder containing `verify.py`. It is `C:\p1\p1val` for the recommended Windows extraction. A guide's command directory matters: several inherited guides use paths relative to their own folders, whereas the root run guide uses complete repository-relative script paths.

The table covers all 27 README files inherited from the 13 November release, including the two third-party notices. Historical guides and sealed capsules are retained for provenance; the current root guides control navigation, software setup and external output locations.

## Current entry points and specialist evidence

| README | Use and command directory |
|---|---|
| [Repository](README.md) | First-time entry point; run from repository root. |
| [Manuscript snapshot](manuscript/README.md) | Optional build of the bound 13 November source snapshot. Updated command blocks run from repository root and write to a sibling directory. |
| [Artwork](validation/artwork/README.md) | Optional figure reproduction. Run from repository root in the reference environment; use the extra dependency listed there. |
| [Low-frequency panel](validation/artwork/low_frequency_ratio/README.md) | Optional low-frequency figure replay; command blocks run from repository root. |
| [Referee checks](validation/referee_checks/README.md) | Independent table, algebra, shell, force-moment and profile-order checks. The root entry commands use repository-relative paths. |
| [Reference suite](validation/reference/README.md) | Inventory and smaller numerical checks; begin with the root run guide. The code path in the original supplement is explained here. |
| [SM2 replay](validation/reference/checks/sm2_replay/README.md) | Saved shell/core records and optional new shell searches. The updated commands run from repository root; numerical searches require the reference environment. |
| [Earlier numerical methods](validation/reference_methods/README.md) | Reference-taper three-method, rigidity, scattering, continuum, checkpoint and clock fixtures. Its commands run from `validation/reference_methods/`, using new external output directories. |
| [Section evidence](validation/section_evidence/README.md) | Recovered mechanics, shell and independent packet evidence. Paths start in `validation/section_evidence/`. Run old producers only in an external working copy as explained there. |
| [Prediction and upstream tests](validation/predictions/README.md) | Distinguishes retained-kernel replays from fresh directional integrations. Entry commands run from repository root, with new external outputs. |
| [Directional control at 2.6](validation/predictions/upstream/control_26/README.md) | Archived-record replay versus fresh bulk integrations. Commands run from repository root. The 24-motif directional experiment is separate from the six-motif ramp. |
| [Enhancement calculations](validation/significance/README.md) | Fixed-total comparisons and other November fixtures. The driver command runs from repository root and copies its inputs to a new external directory. |
| [Ramp validation](validation/ramp/README.md) | Scientific scope, proof chain, fresh quadrature and packet distinctions. Short local commands assume `validation/ramp/`; prefer the full paths and separate environments in the root run guide. |

## Retained ramp records and sealed capsules

The following documents contain useful detail but are not first-time setup instructions. Their dated statements describe what was run then. A historical “new” result is not evidence that this documentation update reran that result.

| README | Scope and navigation caution |
|---|---|
| [Earlier applicability checks](validation/ramp/context/README.md) | Relative paths start in this context folder. Public replay wrappers copy inputs to external output directories. |
| [Coupling-scale probe](validation/ramp/context/coupling-study/numerical_probe/README.md) | Historical eight-point scale probe on saved kernels; not a global optimum. Use the public context replay wrapper rather than modifying its stored outputs. |
| [Profile comparison](validation/ramp/context/coupling-study/numerical_probe/profile_comparison/README.md) | Historical fixed-profile/fixed-strength survey; relative source paths belong to that original workspace. Current reproducible carrier comparisons are also in the root compact suite. |
| [Fresh-kernel record](validation/ramp/fresh_kernel/README.md) | Distinguishes originally fresh grids from recovered fine-grid kernels. To run new quadrature, follow the root guide and the linked relocation instructions. |
| [Ramp packet](validation/ramp/packet/README.md) | Historical broadband evolution, selection rule and controls. Use the root packet wrapper to regenerate raw histories externally. |
| [Recovered packet device](validation/ramp/packet/recovered/p1_packet_reproduction/device/README.md) | Producer and analysis documentation from the recovered package. Its working-directory paths belong to that package; root wrappers perform the staging. |
| [Packet continuation guards](validation/ramp/packet_guard/README.md) | Historical continuation examples use that folder's relative layout. Regenerate complete raw histories first; use the root `resume-tests` wrapper for the supported replay. |
| [Packet spectrum](validation/ramp/packet_spectrum/README.md) | Historical stationary spectrum and packet-energy denominator. Its `reproduce.py` command assumes this folder; the root spectrum wrapper stages inputs and outputs safely. |
| [Pre-audit README](validation/ramp/history/before_full_audit/README.md) | Historical authority and paths. Do not use its old extraction-folder name or download instructions for this delivery; current corrections are in the ramp root. |
| [Sealed R4 capsule](validation/ramp/r4_certificate/README.md) | Original reference-taper certificate and candidate manuscript. Relative paths start at `validation/ramp/r4_certificate/`. Its older manuscript/PDF references belong to the retained capsule, not the current snapshot. Use the root `baseline` wrapper to replay safely. |
| [R4 export](validation/ramp/r4_certificate/analysis/full_domain_r4/export/README.md) | Detailed saved-state audit and strict continuation-runtime requirements. “Outer package root” means the R4 capsule root above, not the combined repository root. Use the public wrapper for a read-only replay in an external workspace. |
| [R4 integrator](validation/ramp/r4_certificate/analysis/full_domain_r4/integrator/README.md) | Low-level staged integration API. `run_stage.py` assumes the integrator folder. This is an advanced producer interface requiring valid settings and its pinned proof environment, not a quick validation command. |

## Third-party notices and preserved instructions

| File | Purpose |
|---|---|
| [REVTeX notice](manuscript/third_party/README-REVTeX.txt) | Unchanged upstream documentation and license information; not instructions for running this repository. |
| [R4 REVTeX notice](validation/ramp/r4_certificate/analysis/full_domain_r4/manuscript/dependencies/revtex-tds/doc/latex/revtex/README) | Byte-identical copy retained inside the sealed capsule. |

Original copies of the nine edited guidance files and the three refreshed hash inventories are in [docs/history/2026-11-13/](docs/history/2026-11-13/). These added archival copies are not new entry points. [DOCUMENTATION_UPDATE.json](DOCUMENTATION_UPDATE.json) records the original/current hashes and the permitted scope of this update.

For an equation, percentage or proof claim, consult the included manuscript, [EVIDENCE_BY_SECTION.md](EVIDENCE_BY_SECTION.md) and the specific scientific record. A README summary is not an additional proof or a new physical claim.
