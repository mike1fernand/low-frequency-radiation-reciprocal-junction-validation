# Prediction and upstream validation

This directory restores the calculation paths named in Supplemental Material SM7. Each path in the supplement is relative to this directory. The upstream example uses a 24-motif contact with line speed 1.5 and frequency 2.9. Its reported upstream power of approximately 71% uses incident line power as the denominator. The 73.6% qT ramp calculation uses a different contact and frequency.

Run the following from the repository root. Each command requires a new output directory outside this evidence directory. It checks the original source hashes, copies the scientific inputs to the output directory, executes the producers there, compares key results with the supplied reference and writes logs and a numbered list of executed commands in `REPLAY.json`. These regression tolerances measure reproducibility and supply no physical integration-error bound. Python 3.12, NumPy 2.3.5 and SciPy 1.17.0 were used. The sampled high-precision root check additionally uses mpmath 1.3.0.

```sh
python -B validation/predictions/replay.py predictions ../paper1-predictions
python -B validation/predictions/replay.py upstream ../paper1-upstream
python -B validation/predictions/replay.py fresh-cartesian ../paper1-cartesian
python -B validation/predictions/replay.py fresh-gb ../paper1-known-regulator
```

| Profile | Scope |
| --- | --- |
| `predictions` | Independent low-frequency coefficient and sparse-contact checks, cache-rejection controls, 31 interference rows, 2,123 surface and 579 tetrahedron bandwidth rows, finite-frequency angular checks, and the low-frequency cone at frequency 0.003125 and order 24. The bandwidth and cone replays use the archived bulk matrices. |
| `upstream` | Combines 16 saved nonoverlapping Cartesian chunks covering N384, checks the three zero-width fits, applies the left/right reflection test, tests 80 transverse root samples at 70 decimal digits, and freshly integrates one N768 real-axis case at 2.9. Files for all 11 real-axis checks are under `upstream/realaxis/`. |
| `fresh-cartesian` | Clears the Cartesian output in the external copy, regenerates all 16 N384 chunks from the force law and computes the three extrapolations. The chunks span consecutive groups of 24 axial grid planes. |
| `fresh-gb` | Clears the output kernels in the copy, regenerates three regulators and the real-axis residue at each of N384, N768 and N1536, and performs both extrapolations. The specified regulator is `z=(omega+i*eta)^2`. This independent corroboration does not establish the generator or regulator used for the historical comparison row in SM4. |

The September 26 prediction archive was internally revised September 27. Its retained arrays, rows, mathematical derivation and source routines are in `baseline/Prediction_Calculations_2026-09-26/`. Use `lowfreq_audit/low_frequency_guarded.py` when making new low-frequency kernels. The baseline producer is preserved for the original calculation and deliberate stale-cache regression. Use a copy when running original scripts directly, since many write beside their inputs.

The current cone producers and records are under `angular_audit/` and `collection_extension/`. The latter includes four frequency/order cases and their source-matrix hashes. `controls_audit/` contains the finite-contact interference and bandwidth checks. A further set of independently implemented checks under `independent_comparisons/` reproduces kernel fits, source assembly and cone comparisons. In an external working copy its `portable_checks/rev/bandwidth_check.py 1.1125 1.2375 1.3 1.29` command checks 193 settings at each of the first three frequencies. It reports that no comparison rows exist for 1.29. Its `green_audit/check_green_evidence.py` tests a finite-width 3D resolvent, which cannot bound the zero-width extrapolation.

`upstream/prior_staged/` preserves the earlier 2.6 shell and 2.9 principal-value reconstruction. To repeat another real-axis case in a working copy, run `python -B upstream/reconstruct_realaxis.py 2.9 1536` or use a frequency and grid size recorded in `upstream/realaxis/`. To repeat a single Cartesian stage, run `python -B upstream/reconstruct_cartesian.py 384 0 24`. Resume the other disjoint 24-plane groups as needed and run `python -B upstream/reconstruct_cartesian.py 384 combine` once all 16 stage files are present. The combiner checks their producer identity and exact axial coverage. Keep assertions enabled for direct producer calls.

`SOURCE_FILES.json` records the original archive, exact member name and SHA-256 for each of the restored files. Historical absolute paths occurring in data are run provenance. Superseded editorial material is excluded. The upstream records reproduce the scientific checks described in SM7E. The 80 root samples have a specified seed and selection rule. General root completeness and a certified high-frequency integration error remain open.

The [control at frequency 2.6](upstream/control_26/README.md) adds saved-matrix and fresh-grid comparisons for the same 24-motif geometry, including the real-part and shell swaps. Its launcher stages the inherited sources outside this repository.
