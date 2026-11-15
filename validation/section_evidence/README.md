# Model, reference-shell and packet evidence

Paths in this guide are relative to this folder. The recovered scientific sources and data retain their archived bytes. `ORIGIN.json` identifies each source archive, nested path and SHA-256. Reports embedded in those sources refer to the editions and programs available at their dates.

| Location in the paper | Recovered evidence |
| --- | --- |
| SM1 central-link forces and moments | `source_mechanics/code/sectionIII_independent_audit.py` implements symbolic differentiation, 5,000 five-point gradient comparisons, work and torque identities, Fourier sums and long-wave checks. Its archived output in `source_mechanics/records/sectionIII_independent_checks_v1.json` gives the observed scaled maximum as 4.204130575684652e−14. The 4.21e−14 in SM1 is a conservative upward rounding. |
| Sec. IV and SM2 shell extrema and symmetry reduction | `bulk_modes/records/P04_SHELL_BRANCH_{0,1,2}.json` report each search with its coordinates, gaps and residuals. `P04_CORE_CHECKS.json` in that folder and `bulk_modes/derivations/EXACT_CORE_RESPONSE_REDUCTION.md` give the finite-core comparison and exact argument. `bulk_modes/WITNESS_CORRECTION_EVIDENCE.json` supplies the full-precision coordinates and $Q/a=1$ normalization of both printed source-overlap witnesses. The searches give numerical extrema estimates. The analytic proof in current SM2 separately bounds the upper/lower squared-frequency gap. |
| SM5C independent C++ scalar diagnostics | `independent_packet/v32/results/target_blind_global_td_120x56_dt004_v32.json` gives the original Python trajectory observations. `independent_packet/earlier_v4/evidence/sectionVII_scalar_closure_rows.csv` gives the separate OpenMP C++ observations. Six quantities can be compared at ten common integer times, 48 through 57. The maximum absolute difference across those 60 pairings is 2.825564047502013e−15 in the recorded units. Each comparison is of scalar quantities. |
| SM5C checkpoint statements | `independent_packet/v32/code/validate_segmented_checkpoint_equivalence.py` and the corresponding JSON in `v32/results/` give the 40 × 12 × 12 Python regression, $\Delta\tau=0.04$ to $\tau=10$. The later 120 × 56 × 56 result compares monolithic and restarted C++ trajectories through $\tau=57$. Both C++ implementations used in that restart check descend from the same C++ force kernel. |

The P07 record `P07_INDEPENDENT_REPRODUCTION.json` describes zero difference for a 46-row C++ replay. Its reference CSV is byte-identical to the earlier v4 C++ output, and a source comparison shows P07 preserved the v4 force kernel as it added configurable arguments and outputs. The independent comparison between languages uses the ten common times and six observables identified above. The v4 C++ file retains machine-specific output names. The equivalent P07 version in `independent_packet/run/p07_trajectory.cpp` accepts both output filenames as arguments.

From this directory, compare the two archived time series with:

```sh
python -B independent_packet/compare_diagnostics.py
```

To compare against a freshly generated reference-taper C++ CSV, append its path to that command. To use a new Python trajectory, also pass `--python-json PATH`. The threshold of 2.83e−15 checks reproduction of the historical discrepancy on this tested build. It is not a bound on physical model error or a platform-independent floating-point tolerance.

For fresh shell and force tests, first copy the entire `bulk_modes` and `source_mechanics` directories into a separate work directory. From that work directory run:

```sh
python -B bulk_modes/code/p04_run_all.py
python -B source_mechanics/code/sectionIII_independent_audit.py
```

These commands use the NumPy, SciPy and SymPy versions in the main reproduction guide. The shell driver writes its reports into the work copy's `bulk_modes/records/`. The force driver writes into `source_mechanics/results/`. Do not run these in the evidence copy.

To repeat the historical small Python checkpoint fixture, use a work copy of `independent_packet/`, install the packages in its `source/phase12_core/requirements.txt` including Numba, then run from that copy:

```sh
python -B v32/code/validate_segmented_checkpoint_equivalence.py --phase12-root source/phase12_core --output small-check.json
```

The original Python producer and the three required Phase 12 Python modules are included. The module hashes match between v32 and P07. With Numba 0.68.1, llvmlite 0.50.0, NumPy 2.3.5 and Python 3.12.14, the 40 × 12 × 12 uninterrupted and segmented results agreed exactly. The generated fixture summary was byte-identical to its archived counterpart. The repository's reference C++ smoke test is a separate replay.

For the full Python reference evolution on the 120 × 56 × 56 free-boundary lattice, use the same work copy and environment. The archived v32 snapshots occur every 1.0 from 48 to 57. Set `--trace-step 0.2` explicitly to recover its 285 trace times through 57.

```sh
python -B v32/code/run_target_blind_global_td.py --phase12-root source/phase12_core --trace-step 0.2 --output python-reference.json
python -B compare_diagnostics.py /path/to/trajectory.csv --python-json python-reference.json
```

Replace the CSV path with the current C++ reference-taper output. The fresh Python run reproduced all 285 archived scalar trace rows exactly. Comparing the fresh Python and C++ trajectories gave the same six-field maximum difference as the original cross-language comparison. The v32 Python driver uses a later analytic projector implementation to partition modal energies. Its modal-energy fields and metadata can differ from the pre-patch stored file. The six fields used for the cross-language comparison depend on the trajectory and local observables.
