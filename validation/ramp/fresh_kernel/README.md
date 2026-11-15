# Fresh numerical quadrature checks for the 73.6% P1 candidate

8 November 2026. One fixed physical candidate: omega=1.2, actual per-link complete-motif stiffnesses g=(1,2,3,4,5,6)/7, same six-motif graph/48 springs, bare bulk a=m=K_axis=K_diag=1 and line K=1.69. No optimization, frequency search, manuscript or force-law change.

## What actually ran

1. Fresh marching-surface spectral matrices N=48,80,128, using the preserved basis-invariant nodal producer with strict guards against hidden marching failures or skipped small-gradient cells.
2. Fresh scrambled Sobol spectral histograms at2^18 and2^20 samples, two NEW seeds20261108/20261109, PV orders192/256/320, retaining individual seed and seed-mean results. These produce27 candidate response rows. The selected old N192/M256 candidate is rechecked only as a comparison; it is never called a fresh result.
3. The R8 archive's genuinely independent corrected linear-tetrahedron route was recovered. Five bound omega1.2 matrices (N96,128,128halfshift,160,192) were checked against kernel/config/source hashes. Each old reference12.6% row reproduces exactly before the new ramp is evaluated by that route's independent solver. No new fine-grid tetrahedron integration is claimed for these five kernels. Note that omega_000 is omega1.0 in several folders; the correct omega_001 was explicitly selected there.
4. ONE genuinely NEW corrected tetrahedron integration atN48 was executed after the producer became available. It completes the110,592cubes in two exact saved ranges [0,65536) and[65536,110592). No scalar factor was negative/nonfinite, branch matrices were PSD within roundoff, and all candidate state/source vectors are saved.

Fresh surface/PV results use the same archived algorithms as the selected method but new grids/samples. Tetrahedron integration is a different singular-cell integration route, with its own independently written bulk/contact/lead/scattering implementation. The fine-grid ladder is archived; the N48 kernel was generated during this validation.

## Main results

- Original selected cached candidate:73.58742884% qT.
- Fresh N128 surface with M256/two2^20-seed mean:73.58519712% qT.
- Fresh N48 corrected tetrahedron:73.59112487% qT.
- Independent archived tetrahedron N96:73.57884866%; N128:73.58198447%; N128halfshift:73.57406602%; N160:73.57315187%; N192:73.57847434%.
- All27fresh surface/PV order/seed/grid rows are in SUMMARY.json; they remain near73.6%. Their spread is observed sensitivity, not an error bar.
- Corrected tetra reference regressions have zero numerical difference from their saved rows. The NEW coarse N48 reference control gives12.52833481%; it is a new coarse-grid control, not an exact regression of the fine reference value.
- Fresh surface/PV computation:about21.7s total. Fresh tetrahedron:6.21s. The restored finer tetrahedron kernels were not regenerated.

These results strongly corroborate the stationary approximate73.6% point within the specified harmonic model. They do not prove the continuum/infinite-lattice decimal with a rigorous global error enclosure; do not assign a statistical95% confidence, bandwidth, packet or tolerance certificate from mesh/seed agreement. Full channel closure is internal consistency, not quadrature accuracy.

## Files

- run_fresh_kernel.py: executed staged fresh surface/PV script and exact predeclarations.
- inputs/: exact copied producer/API/core sources, selected old reference kernel and its metadata.
- results/: fresh spectra, full fresh spectral histograms,27response rows and every saved coupled state/system matrix.
- check_independent_tetra.py / independent_tetra/: bound original R8 kernels/metadata/configurations, independent source modules and new candidate response data.
- run_fresh_tetra_N48.py / fresh_tetra_N48/: exact original corrected tetrahedron sources, two restartable range checkpoints, full new kernel and candidate/reference results.
- SUMMARY.json: compact distinction between new quadrature and reused kernels.

Scripts retain their execution-time upstream source paths for provenance. Inputs required for audit are included. For a fresh rerun, stage into a new output tree and adapt only the input/output paths; do not overwrite preserved results. scikit-image0.26.0, lazy_loader0.6 and mpmath1.4.1 were installed task-locally. These dependencies are NOT part of the evidence and vendor/ must be excluded from the continuation archive. All package versions are recorded in ENVIRONMENT.json.
