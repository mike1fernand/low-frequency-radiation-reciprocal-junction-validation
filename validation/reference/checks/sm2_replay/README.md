# SM2 symmetry-reduction check

From the **validation repository root**, run the archived-record check with standard-library Python:

```sh
python -B validation/reference/checks/sm2_replay/verify.py --records-only
```

For the fresh P04 replay, use the reference environment's Python and omit `--records-only`:

```sh
../paper1-ref-env/bin/python -B validation/reference/checks/sm2_replay/verify.py
```

The second example uses a POSIX environment path. PowerShell users can substitute `..\paper1-ref-env\Scripts\python.exe`; see the root [run guide](../../../../RUNNING_VALIDATIONS.md). If you instead work from `validation/reference/`, the corresponding relative script path is `checks/sm2_replay/verify.py`. The old standalone-manuscript starting directory is not the root of this combined repository.

The default run requires NumPy and SciPy. It verifies the five archived files against their original SHA-256 digests and runs the original P04 producer in a temporary directory. `--records-only` checks the saved values using the Python standard library. `SOURCE.json` gives the SHA-256 of the original ZIP, its member prefix and member digests. The original producer, its two dependencies and both records are preserved byte for byte.

The passage concerns the six-segment reference taper at ω = 1.2, broadening η = 0.03 and a periodic 32 × 32 × 32 wavevector grid. In `records/P04_CORE_CHECKS.json`, read `metrics.reduction["32"]`. The lower-resolution grid has a separate entry. The producer uses seed 20260813 and 50 complex seven-component line inputs for the full and reduced response comparison on each grid.

- `V_leak`, `U_leak`, `G_leak` and `G_commutator` are the four relative Frobenius residuals defined in SM2.
- `full_reduced_field_max` and `effective_self_energy_max` are maxima of relative vector errors, with the Euclidean norm of the full core displacement or full self-energy action as the denominator (floored at 10⁻³⁰ in code).
- The separately recorded `orth` is the absolute Frobenius error in the 12-column basis orthonormality check.

All seven archived values are below 5 × 10⁻¹⁶. This is an observed result for the stated regression. The analytical justification for the exact reduction is the symmetry argument in SM2. Fresh numerical values may vary with the numerical libraries. The verifier reports those values and compares them to the original output. The original producer applies its own cross-grid regression tolerances.

A separate archived comparison is `reduction.G_match` in `ARCHIVED_SECTIONIV_INDEPENDENT_CHECKS_v1.json`. It compares the projection of the full core Green matrix with a separately constructed reduced matrix. Its relative Frobenius difference is 1.6311799907932728 × 10⁻¹⁵, with the projected full matrix norm as the denominator. This entry belongs to the independent historical record. The present numerical replay covers the P04 producer. The formula and conditions for the independent comparison are documented in §S4-10 of the July 19 Section IV audit.
