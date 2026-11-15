# Numerical inputs for the Section VIII sweep replay

The bundled inputs support new solutions of the loaded line and contact equations for the 85 one-factor cases and the reference taper’s frequency-and-tolerance screen. Run `OPENBLAS_NUM_THREADS=1 python tools/verify_robustness.py` from the cumulative package. The one-factor ledger is already included at `evidence/vi_results/ACTIVE85.json`. The second ledger, `REFERENCE_FREQUENCY_SCREEN.json`, is an unchanged copy of the prediction release’s `bandwidth/FINAL_SCREEN.json`.

`KERNELS.npz` contains 49 base-medium frequency arrays, five stiffness-ratio arrays, 11 surface-method screen arrays and three tetrahedron-method screen arrays. Every array stores the real part of the projected bulk Green matrix followed by the spectral matrices of branches 0, 1 and 2. The solver assembles the eight Cartesian link extensions of each complete motif, solves simultaneously for line and bulk contact motion, and computes qT power using the sum of the two lower-branch matrices. The one-factor replay uses the original lead regularizer of 10⁻¹². The joint screen uses the physical real-frequency lead.

The associated source archives and SHA-256 digests are:

| Archive | SHA-256 |
|---|---|
| PAPER1_VALIDATION_INTEGRATED_R7(1).zip | e2d3a219de57202001ccd1b20dd985079067c920959c8823368145f73d67137a |
| Prediction_Calculations_2026-09-26.zip | 5a0f9fd7259be6a899fd904253fda62b60ea973f9c79a5e7799316d801bb9629 |

The source member names and reconstruction equations are specified in `tools/prepare_robustness_kernels.py`. To regenerate the compact inputs from the two archives, run

```sh
OPENBLAS_NUM_THREADS=1 python tools/prepare_robustness_kernels.py \
  --r7 'PAPER1_VALIDATION_INTEGRATED_R7(1).zip' \
  --screen Prediction_Calculations_2026-09-26.zip \
  --output regenerated_kernels.npz
```

The script authenticates both archives. It reconstructs the real principal-value matrices of the common-medium N = 80 surface calculations from the two saved order-256 Jackson spectral histograms at each requested frequency. The changed-stiffness and three finer screen calculations use their own saved principal-value matrices. The tetrahedron inputs supply their own complex Green matrices, with their imaginary parts checked against the three stored branch contributions. Rebuilding the compact matrices requires NumPy. Tiny machine-dependent floating-point differences can arise across versions or BLAS backends.

The maximum newly solved channel discrepancy was 2.34 × 10⁻¹⁵ across 2,787 row comparisons on the tested installation. The 579 tetrahedron rows evaluate a subset of the same physical settings used by the 2,123 surface rows. The numerical replay uses retained bulk data and covers the sampled settings. The packet results and the continuous coupling-box certificates have separate checks.
