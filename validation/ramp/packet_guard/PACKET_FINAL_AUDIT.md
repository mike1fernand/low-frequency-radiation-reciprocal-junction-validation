# Final packet audit — 8 November 2026

The bounded packet validation passes. This source/evidence review found no new defect in the stored physics results. It found a genuine **replay precondition** in the archived engine, addressed by an additive guard, and clarified what a compact reproducible repository contains. The archived physics engine, analyzer and scientific results are unchanged.

## Scientific scope checked

| Question | Finding |
|---|---|
| Was the entire selection preregistered? | No. The 56-box run failed the unchanged duration requirement, then the 80-box run followed those diagnostics. The selector is coefficient-blind, and its 57.4–64.2 interval was frozen before the half-timestep and larger-box controls. “Predeclared” must be read in that sequence; it does not mean the numerical-domain adaptation was globally preregistered. |
| Were control windows changed to improve agreement? | No. All three reported cases use the frozen interval; the existing independent reviewer recomputed all five predicates and all 35-row means. |
| Were the physical 48 contact springs retained? | Yes. There are eight springs per motif, six motifs, and actual per-link stiffnesses `(1,2,3,4,5,6)/7`. Independent full-force/current one-step controls checked the exact y/z reduction. There is no axial reflection assumption. |
| Does the direct packet measure monochromatic 73.6% transmission? | No. It measures finite-time stored energy in the two lower bulk branches of a periodic Fourier embedding, divided by initial physical open-line energy. The bandwidth-weighted stationary prediction is about 71.63%. |
| Are the lower branches exact transverse Maxwell fields? | This calculation makes no such claim. It does not establish a local qT flux, exact native-displacement transversality at all wavevectors, or an electromagnetic field. |
| Is the tiny modal allocation interval a total error bar? | No. It addresses degenerate-subspace allocation only. Timestep, domain, cadence and method differences remain separate sensitivities. |
| Is the secondary cadence endpoint correction current? | Yes. `packet/CURRENT_CADENCE_CHECKS.json`, `CADENCE_FIELD_CORRECTION.json` and the current analyzer supersede the old secondary common-mean fields. The primary cadence-0.2 means, frozen window and dynamics are unaffected. |

The direct main midpoint is `0.7164564184629316`; the larger-box control is `0.7162667759527709`. The stationary value with the matching open-energy denominator is **`0.7163064862075107`**. The old periodic-energy-normalized value, `0.7163064796290152`, remains legitimate when labeled with its own denominator. `packet_review/FROZEN_CONTROLS_CHECK.json` already uses the matched open normalization. The tail-only spectral interval is not a total numerical or dynamical uncertainty enclosure.

## What the reproducible repository must contain

The archived raw manifest lists **58 files totaling 3,162,169,048 bytes**. All 58 paths were present and their sizes matched in this pass. The existing manifest retains the individual SHA-256 values; this final pass did not redundantly hash all 3.162 GB again.

The approximately 9.2 MB main checkpoint at τ66 contains complete line and reduced-bulk displacement and velocity, plus the model and integrated-current metadata. Its SHA-256 is `9bc2744e79274ef57de60949416edce588b8c79e2d8ef24e57cd9046925f69ac`, rechecked here. Exact transverse reflection reconstructs its full physical state. **That state does not reconstruct its earlier trajectory.** The standard analyzer needs the history beginning at step zero.

A complete *reproduction repository* can include all source, configurations, preparation, mathematical inputs, results, manifests and one continuation checkpoint while omitting the multi-gigabyte regenerated traces. The download must describe this omission clearly. It must not describe the single checkpoint as all raw validation data. The main three production controls can be regenerated from the bundled initial state using `packet/README.md`; the first 56-box result remains a preserved failed-domain control.

Packet replay requires Python, NumPy with `trapezoid` (NumPy 2 or later), a C++17 compiler with OpenMP, and the POSIX `resource` module used by the wrapper. The stationary spectrum replay also requires SciPy and its bundled `inputs/` tree. Linux is the tested platform. Native binary files use little-endian float64/int32; floating-point bitwise identity on another platform is not promised. The 18 load-bearing small-file presence/hash records are in `PACKET_FINAL_AUDIT.json`.

## Additive correction to staged replay

The archived C++ engine appends to `line.bin`, `corner.bin`, `snap.bin`, `snap_times.bin` and `monitor.csv`. It verifies checkpoint parameters but does not verify that those histories end at that checkpoint. Consequently, an old checkpoint appended to later histories can duplicate data; a checkpoint by itself cannot supply missing earlier traces. This did not invalidate the saved runs: the preserved continuous-versus-staged small-box control already matched byte for byte.

`safe_resume.py` closes that avoidable replay problem without changing the engine:

- It checks little-endian layout, the pinned source, complete checkpoint size/finite state, six-motif profile, bond inventory, model constants, geometry and timestep.
- It requires all preceding histories, exact step-zero-inclusive byte counts, the complete saved-time grid, monitor sequence, and matching line/corner endpoints. It also matches a saved bulk frame to the checkpoint when their endpoints coincide.
- It rejects analyzer-ready final targets at or after τ48 that fall between the 0.2 snapshot times.
- Execution is explicit. It takes an exclusive wrapper lock, compiles the pinned source, rechecks before appending, preserves the exact old configuration, and updates the completed target fields and active binary hash only after successful history/state validation.
- It retains failure receipts and output. It neither attempts rollback after an interrupted append nor claims to detect arbitrary corruption inside earlier trace records. A failed or interrupted directory must be diagnosed before any retry.

The final **15 controls pass**: 12 expected rejection/failure controls, two positive structural checks, and one H=2 continuation step from 20 to 21. The synthetic snapshot fixtures test file layout only. They are not new physical trajectories. A second source review caught and corrected the configuration-target, final-cadence and failure-log issues before the final controls. See `SAFE_RESUME_VALIDATION_FINAL.json`.

## Repeating the independent reviewers

Saved reviewer receipts remain usable without raw snapshots. Re-executing `check_box80_analysis.py` requires regenerated raw histories and the expected `packet/candidate_dt004_box80` basename. Re-executing `check_physical_step.py` requires fresh `random_step_H2` and `random_step_H3` output directories. Work in a **fresh replay workspace**, copy only the needed source/configuration/preparation there, and regenerate the primary under its original basename. Do not overwrite or remove the frozen validation records to make the scripts run.

No new optimization, long dynamics, manuscript change or additional physics acceptance criterion was introduced by this audit. This pass is limited to packet validation and reproduction; it does not certify all Paper Three claims or empirical applicability.
