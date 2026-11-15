# Full validation audit and corrections

8 November 2026. This release audits the Paper One high-conversion candidate and makes the validation replay self-contained. It does not revise the manuscript, change the physical model, or start a new optimization campaign.

## Scientific outcome

No new load-bearing mathematical or packet-physics defect was found. The exact response derivations, original contour evidence, fresh independent numerical checks, packet spectrum, time evolution and controls retain their recorded results. The nominal candidate estimate is about 73.5874%. The rigorous containing interval is 73.2166%–73.9652%, displayed outward. The continuous six-motif ±5% lower bound is 0.6905008642… > 69%. These are three different statements: the interval does not certify one-decimal rounding to 73.6%.

Each motif contains eight physical springs. Its one coefficient changes all eight together. The robustness theorem covers six independent motif factors, not 48 independently perturbed springs. The contact has 48 springs and total physical stiffness 24 at the nominal point. Neither bandwidth, global optimality, experimental efficiency nor microscopic DST selection is established by these validations.

The packet is broadband. Its direct qT stored-energy result is compared with the energy-weighted stationary packet prediction, with the same open-line energy denominator. That prediction is approximately 0.7163064862, rather than the monochromatic 0.7358742884. The primary packet result, half-timestep control and larger-box control remain approximately 0.7164564185, 0.7164540867 and 0.7162667760. Their observed differences are sensitivity checks, not a combined rigorous error enclosure.

## Confirmed corrections

| Finding | Correction | Scientific effect |
|---|---|---|
| Earlier README supplied a relative ledger path with a duplicated package prefix. | The new baseline runner supplies an absolute staged-ledger path. | No change to the earlier successful absolute-path replay or bounds. |
| Earlier root file checker allowed missing manifest entries and extra files; Python `-O` disabled its assertions. | Explicit integrity errors, exact inventory comparison, malformed-path/symlink rejection and adverse fixture tests. | Stronger delivery verification; no changed physics. |
| Earlier compact core included only selected upstream R4 files. | All 683 original R4 files, including the inherited archive, are included here. | Removes the separate source-part dependency. |
| Historical replay scripts can write within their own source trees. | `replay.py` stages numerical work externally and verifies frozen evidence before and after. Child optimization is disabled. | Preserves evidence and failed-run logs. |
| Fresh quadrature drivers retained historical input paths. | Additive relocation adapter stages bundled inputs and changes exactly the I/O root. | Relocated pilot reproduced all nine responses and internal arrays exactly. |
| Proof and fresh-quadrature environments used different mpmath versions. | Separate pinned requirements and instructions. | Avoids claiming one exact environment for both. |
| A complete checkpoint was easily mistaken for complete analyzer history. | Additive safe-resume guard rejects missing, stale, truncated and misaligned histories, preserves configuration and failure records, and checks snapshot-grid compatibility. | Original engine unchanged; 15 guard controls pass, including a real small-box continuation. |
| Report generator depended on system font paths. | Bundled licensed fonts and an external report adapter. | Report content unchanged; current root documentation governs delivery. |

## Corrections and limits already present, preserved here

- The old secondary every-third-sample cadence field omitted the 64.2 endpoint. The primary mean and chosen interval were unaffected. `packet/CURRENT_CADENCE_CHECKS.json` and `packet/CADENCE_FIELD_CORRECTION.json` remain authoritative over that legacy field; the current analyzer includes the endpoint.
- The initial 56-box run failed to produce a valid interval. The 80-box adaptation was diagnostic-informed and coefficient-blind. The 57.4–64.2 interval was frozen before the two controls. This is not a globally preregistered experiment.
- The packet estimator measures global periodic-embedding lower-branch stored energy in a finite open evolution. It is not a local qT flux or an electromagnetic field.
- Baseline replay verifies 222 saved compact job states, fresh compact physical checks and the bridge. It does not reintegrate all 449,698 contour leaves. The raw ledger's pending flag and the successful combined physical/bridge flag refer to different stages.
- Kernel target registration is exact for the actual cached binary64 G/Q arrays. No recentering budget is silently added. The authenticated bounds are rounded upward to eG = 1/250 and eQ = 1/500 in the exact-rational response proof.

## Completeness and raw-data policy

This is the full replay validation repository: all contour certificate sources, candidate proofs, scientific inputs, numerical and packet source, archived results, independent reviews, one complete primary continuation checkpoint, current portability adapters and audit tests are included. The report is the same dated four-page scientific report; its old split-download paragraph is superseded by this release's README.

The raw inventory contains 58 generated packet files totaling 3,162,169,048 bytes. The repository includes the 9,218,168-byte primary checkpoint and omits the remaining 57 regenerable files totaling 3,152,950,880 bytes. Their sizes and hashes are in `packet/RAW_DATA_MANIFEST.json`. Exact preparation, configurations and engine source regenerate them. The included primary τ66 checkpoint is complete physical state, but its earlier analyzer histories must be regenerated before normal history-based continuation or analysis. The safe guard deliberately rejects the compact saved directory when those histories are absent. Cross-platform bitwise raw-file identity is not promised.

Historical source bytes are preserved. `SOURCE_PRESERVATION.json` maps original files to their unchanged locations or explicit historical copies. `FULL_RELEASE_CHECKS.json` identifies the checks actually executed for this release. The root `PACKAGE_SHA256SUMS.json` covers every regular file except itself; it detects accidental changes, and is not an external signature or a substitute for the mathematical proof.

## Improvements and stopping rule

The improvements concern replay reliability and evidence preservation. Independent reviews checked the additions; relocated and clean-extraction runs check that they work together. No evidence in this audit justifies changing the harmonic physics model or repeating the optimization. Future manuscript work should cite the containing interval, motif-level tolerance scope and packet normalization correctly. New physics would need its own stated claim and evidence; it is not a prerequisite for using these candidate validations.
