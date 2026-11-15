# Independent packet corroboration of the new Paper One contact

8 November 2026. This is a new finite open-junction calculation for the six-motif ramp `g = (1,2,3,4,5,6)/7`. The manuscript, original packet producer, and earlier numerical records are unchanged. No further dynamics is selected by this record.

## Main result and the correct comparison

The stationary carrier coefficient is approximately **73.5874%** at `omega = 1.2`. The original Gaussian packet has a finite frequency spread, so its stationary packet prediction is approximately **71.63065%**, rather than 73.6%. The direct evolution must be compared with that energy-weighted prediction.

The predeclared main run uses a 120-site free scalar line and a `120 x 80 x 80` vector bulk, with timestep 0.04. The unchanged five-predicate, coefficient-blind selector chooses **57.4 <= tau <= 64.2**. On that interval, its invariant lower-pair stored-energy average is

`[0.7164564182492906, 0.7164564186765725]` of the initial open-line energy,

with midpoint **71.6456418463%**. The narrow allocation interval only covers the choice of the qL subspace at exact modal degeneracies; it excludes integration, timestep and finite-box errors.

The original `56 x 56` transverse box and `48 <= tau <= 57` search were tested first. They did not produce a passing interval of duration at least 2: only the isolated point 53.8 passed all predicates. Stronger-contact ringdown and later boundary exposure explain the failed predicates. That run and its failure are preserved, and no qT average was claimed from it. The explicitly authorized larger-box search retained every predicate and threshold.

## Controls on the frozen interval

The 57.4–64.2 interval was frozen in `FROZEN_WINDOW.json` before running the controls. It was not reselected to match a predicted coefficient.

| Run | Timestep | Transverse box | Mean qT / initial open-line energy | Relative change from main run |
|---|---:|---:|---:|---:|
| Main, predeclared | 0.04 | 80 x 80 | 71.6456418463% | — |
| Half timestep | 0.02 | 80 x 80 | 71.6454086714% | -0.00032546% |
| Larger box | 0.04 | 96 x 96 | 71.6266775953% | -0.02646951% |

All original five predicates pass throughout the same frozen interval in both controls. The main result remains the 80-box run; the 96-box value is a box sensitivity control, not a replacement chosen because it lies closer to the stationary prediction. Relative to the stationary in-band packet contribution, the main difference is +0.02093%; the larger-box difference is about -0.00554%. These are observed comparisons, not a combined error bar.

The total node-assigned energy in the unique union of all six two-cell boundary faces is independently checked after selection. Its maximum fraction is 0.18743% in the main box and 0.03758% in the larger box. The axial-only fraction in the larger box is 0.01077%. The selector's transverse kinetic-shell predicate is distinct from this full-energy recertification.

## Physics, preparation and normalization

The line stiffness is 1.69; the bulk has unit masses and unit axial/face-diagonal central springs, with all outer boundary bonds omitted. The six contact segments are 60–65. Every motif retains its eight physical springs; each `g_j` is its actual per-link stiffness. There are no wrapping bonds in the evolution and no external drive after preparation.

The initial line displacement is `0.002 exp[-(j-18)^2/98] cos[k0(j-18)]`, with `k0 = 2 asin(6/13)`. The 120 spectral velocities use `-i sign(k) omega_q(k)` with DC/Nyquist velocity zero. The bulk is initially exactly quiet. The prepared open-line energy is `3.5972978303311774e-5`; the new initial contact energy is `1.064211764392121e-22`. The periodic Parseval diagnostic adds a tiny closing-bond energy, producing a relative normalization difference of `9.184e-9`.

The stationary packet comparison retains every signed propagating bin in its numerator and the full-grid initial energy in its denominator. Nyquist energy remains in the denominator. The retained 14-frequency in-band contribution is `0.7163064796290152`; its tail-only interval is approximately `[0.7163064796290152, 0.716622986870909]`. This interval bounds omitted spectral weight only. It is not a rigorous quadrature, dynamics or finite-time uncertainty interval.

The direct diagnostic is global modal energy of a periodic Fourier embedding of each finite open-box snapshot. Both lower modes are counted once. Full upper/middle degenerate subspaces and triple points use basis-invariant energy-coherence allocation bounds. The calculation does not define a local qT flux or an electromagnetic field.

## Energy and numerical checks

The producer assembles every physical bulk and contact force before reducing only the exact transverse reflection symmetries. There is no axial-reflection assumption, so the asymmetric ramp is retained. Weighted reciprocity is exact in the assembled matrix. Independent full-grid spring-energy summation agrees with the matrix energy within the stated floating-point tolerances. The conservative `dt * omega_max` bound is 0.1711 at timestep 0.04, well inside the Verlet stability range.

The physical open energy balance uses line energy, open-bulk energy and the complete contact energy. In the main interval their mean fractions are approximately `0.2173146078 + 0.7826023803 + 0.0000830716 = 1.0000000598`. These must not be replaced by a sum mixing periodic branch energy with open-box energy. The periodic bulk mean is `0.7834437501`, a separately named diagnostic.

The main integrated interface-current/half-contact-energy balance residual is 0.03062094%. Halving the timestep reduces it to 0.00765353%, approximately the expected factor of four for this quadrature/integrator comparison. The corresponding larger-box residual is 0.03061924%.

`FINAL_SUMMARY.json` separates selected-window drift measured at cadence 0.2 from the largest drift at saved checkpoints across the whole run. The checkpoint-sampled whole-run drift is about `4.12e-5` at timestep 0.04, during the interaction. It is not a maximum over every integrator level. The later selector still uses its original `1e-5` drift predicate, which passes. No unrecorded all-step drift check is claimed.

## Minor cadence-field correction

The original independent reviewer found that the copied secondary `common_mean_low_high_total` field used every third sample and stopped at 64.0 instead of 64.2 for the new 35-sample window. The primary cadence-0.2 result, physical evolution and chosen interval are unaffected.

The raw main and half-step result JSONs remain frozen. **Their old secondary common-mean fields are historical and superseded by `CURRENT_CADENCE_CHECKS.json` and `CADENCE_FIELD_CORRECTION.json`.** The corrected grid explicitly includes the 64.2 endpoint; its main qT midpoint is approximately `0.7164548324`. The current analyzer is corrected for future replays. A separate regular cadence-0.4 comparison on the same complete interval is also recorded in `CURRENT_CADENCE_CHECKS.json`.

## Complete checkpoints and raw-file policy

Checkpoints contain the complete line and symmetry-reduced bulk displacement and velocity arrays, step, timestep, profile, bond inventory, model constants, cumulative current, last current and initial energy. The transverse reflection map reconstructs the complete physical bulk state. Restarting from a checkpoint recomputes acceleration from the stored displacement and does not infer a state from summed diagnostics. A nontrivial small-box control gives identical complete checkpoint and line/corner-history bytes for continuous versus staged evolution.

`RAW_DATA_MANIFEST.json` records every raw trace and checkpoint with size and SHA-256. The multi-gigabyte snapshots can be excluded from a compact core package because the exact source, preparation, configurations and diagnostics regenerate them in seconds to minutes. A useful optional continuation artifact is the main run's `candidate_dt004_box80/checkpoint_1650.bin` (about 9.2 MB, tau = 66). Generated binaries and private transfer helpers are excluded. Hashes identify these executions; cross-platform floating-point bitwise identity is not promised.

## Reproduction

Requires Python 3, NumPy with `trapezoid`, and a C++17 compiler with OpenMP. The core dynamics/analysis does not need the large R8 archive or its stationary kernels. Use a fresh output name so that results are not overwritten.

```bash
g++ -std=c++17 -O3 -fopenmp evolve_candidate.cpp -o evolve_candidate
python3 run_control.py --half-box 40 --dt .04 --end 66 --name replay_box80
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python3 analyze_candidate.py replay_box80 > replay_box80_result.json
python3 run_control.py --half-box 40 --dt .02 --end 64.2 --name replay_halfstep --frozen-window 57.4 64.2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python3 analyze_candidate.py replay_halfstep > replay_halfstep_result.json
python3 run_control.py --half-box 48 --dt .04 --end 64.2 --name replay_box96 --frozen-window 57.4 64.2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python3 analyze_candidate.py replay_box96 > replay_box96_result.json
```

The measured evolution wall times were about 12.7, 21.6 and 20.7 seconds, respectively, on two OpenMP threads in this environment. Peak producer memory remained below 1.3 GiB. Fourier postprocessing adds time and memory. `run_candidate.py pilot` separately preserves the first 56-box pilot and checkpoint-resume control.

## Provenance and scope

The unmodified original producer is in `recovered/p1_packet_reproduction/`, recovered from `Paper1_Packet_Reproduction_and_Erratum_2026-11-01.zip` (public copy). The original and public archive hashes are recorded in `ORIGINAL_PRODUCER_IDENTITY.json`. The C++ and Python producers and their inputs retain their original bytes. The public guide preserves the labels of four private source links. The two patch files document the new parameterization, checkpoint handling and selector-driven analysis changes.

This supports the new high-conversion contact through independent time evolution and appropriate spectral averaging. It does not turn the broadband mean into a monochromatic 73.6% measurement, establish a combined rigorous infinite-lattice error enclosure, or inherit the old reference tolerance certificate. Candidate-specific stationary and tolerance validation are separate records in the parent validation package.
