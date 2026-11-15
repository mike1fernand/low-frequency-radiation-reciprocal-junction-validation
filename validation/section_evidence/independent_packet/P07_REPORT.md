# P07 audit report — real-time packet corroboration

## Decision

**PASS_WITH_CORRECTIONS_FROZEN.**

The real-time calculation independently corroborates the two bounded scalar Paper-1 claims assigned to it:

1. energy leaves the prepared one-dimensional line and enters a bulk that starts from exactly zero displacement and velocity; and
2. the resulting finite-time global periodic branch content is predominantly in the invariant lower-two qT acoustic cluster.

It does not establish a unique local qT current, a remote detector-plane complex field, an asymptotic time-domain scattering coefficient, a continuum limit, polarization control, gauge dynamics, or electromagnetism. The authoritative v5 manuscript and supplement were not edited.

## Scope and authority

P07 groups Main Section VII with the scalar real-time controls in SM5. It uses the current v5 source as wording authority, fresh producer-backed calculations as numerical authority, and the older v32/v4 records only as historical corroboration and provenance inputs.

## Fresh execution summary

The audit performed:

- one full double-precision OpenMP C++ trajectory on the reference `120x56x56` geometry;
- invariant modal analysis of every cadence-0.2 snapshot;
- direct re-execution of the coefficient-blind selector;
- all 66 integration-step current, coupling, drift, and kinetic-shell checks over the selected interval;
- independent Parseval packet normalization and translated open-line controls through `tau=55.4`;
- invariant common-time runs at `dt=0.03`, `0.04`, and `0.05`;
- a fresh `120x64x64` transverse-box run on the same final interval;
- node-assigned six-face shell calculations for both boxes;
- independent current-versus-assigned-energy closure;
- a bounded checkpoint/resume code fixture; and
- clean-extraction replay through the package verifier.

## Packet preparation and initial purity

The exact spectral velocity construction gives a periodic-Parseval right-moving fraction

\[
0.9999997423126528,
\]

with stationary Nyquist fraction

\[
2.5768734723\times10^{-7}.
\]

Every bulk coordinate and velocity is exactly zero initially. The initial source-link energy is

\[
E_c(0)/E_0=2.05020\times10^{-18}.
\]

The full-DFT Parseval energy and the open-line initial energy differ by only

\[
9.18\times10^{-9}
\]

relatively. This is the periodic closing-bond convention and must be stated when the right-moving fraction is reported.

The computed frequency curve `0.8 <= omega <= 1.6` contains `99.9683493%` of the periodic spectral energy. The uncomputed fraction `3.16507e-4` remains a one-sided passivity bound. The independently recomputed discrete 1% and 99% frequency crossings are

\[
\omega_{1\%}=0.93175667,
\qquad
\omega_{99\%}=1.47265622.
\]

The historical interpolated median is not needed by the current manuscript and is not promoted by this audit.

## Open-line preparation control

Relative to the converged translated 160-site control, the 120-site production line changes cumulative incident flux by:

- `2.34619e-4` at `tau=52`;
- `2.22156e-4` at `tau=55`;
- `2.20256e-4` at `tau=55.4`.

The 160-, 200-, and 240-site controls agree within `6.8e-14` relatively at the final endpoint. The periodic spectral preparation therefore does not create a material open-line bias over the accepted interval.

## Coefficient-blind interval

The exact selector contract uses only direct-dynamics quantities on `48 <= tau <= 57`, cadence `0.2`, minimum duration `2.0`, and the five published thresholds. The qT result and both frequency-domain predictions are absent from every predicate.

Two passing runs are found:

- `51.8--52.0`, rejected because it is too short;
- `52.8--55.4`, accepted as the unique eligible longest run.

The selected interval is

\[
\boxed{52.8\le\tau\le55.4}.
\]

The path-independent scientific-content hash is

`9470a0c1d579621a3e8752c61506d1a2886201150c8e722e0787f4702789f1ff`.

At every one of the 66 integration steps on that interval:

- maximum coupling-energy fraction: `0.000849437679279`;
- maximum absolute instantaneous current: `0.000447848848174`;
- maximum kinetic-shell fraction: `0.00192617445806`;
- maximum energy drift: `3.54166927635e-07`.

All pass their frozen direct thresholds.

The total node-assigned six-face shell remains a post-selection recertification because the supplied chronology does not establish at 95% confidence that its numerical threshold was fixed before output inspection.

## Invariant qT result

The robust decomposition identifies 132 upper-middle doublets and 231 threefold points. It yields

\[
0.12762638643388616
\le
\overline{\mathcal E}_{\rm qT}^{\rm TD}
\le
0.12762638738500592.
\]

The invariant midpoint is

\[
0.12762638690944604.
\]

The production fallback value differs only below the displayed precision. The maximum instantaneous ambiguity is `1.851e-09` and the averaged ambiguity width is `9.511e-10`.

The defensible prose result remains

\[
\boxed{\overline{\mathcal E}_{\rm qT}^{\rm TD}\simeq0.1276}.
\]

The companion averages are

\[
\overline{\mathcal E}_{\rm qL}^{\rm TD}=0.0113020560,
\qquad
\overline{\mathcal E}_{\rm bulk}^{\rm per}=0.1389284430.
\]

Thus the lower cluster contains approximately `91.8648%` of the finite-time periodic modal bulk energy. This is a spatially nonlocal global branch-content diagnostic, not a unique local branch current.

## Packet-weighted comparison

The direct result is compared with incident-spectrum-weighted frequency-domain curves, not with the carrier coefficient. The central predictions are

\[
T_{\rm qT}^{\rm packet,GB}=0.1277094,
\qquad
T_{\rm qT}^{\rm packet,surf}=0.1277468.
\]

They differ from the direct result by `0.0650%` and `0.0943%`. All four retained interpolation prescriptions pass the declared 2% criterion. The outside-band energy is kept as a separate one-sided bound rather than silently extrapolated.

## Current and energy accounting

The exact bond-current identity is sound. The discretized trapezoidal current integral and independently evaluated half-bond assigned-energy gain differ by at most

\[
0.02968\%.
\]

This is a numerical quadrature and trajectory-accounting residual, not missing physical energy.

## Timestep and cadence controls

On exactly common times `52.8, 53.4, 54.0, 54.6, 55.2`, the invariant qT averages are:

- `dt=0.03`: `0.127615787172`;
- `dt=0.04`: `0.127616682197`;
- `dt=0.05`: `0.127617833511`.

Their relative span is `1.603503e-05`. The maximum degeneracy-bound width is `1.012e-09`. Cadence `0.2` versus `0.4` changes the selected-window qT average by `0.00749%`.

## Refined transverse-box control

On the fixed authoritative interval `52.8--55.4`, enlarging the transverse box from `56x56` to `64x64` changes:

- qT average by `0.04398%`;
- periodic total by `0.10498%`;
- qL by `0.7938%`;
- current by `2.064e-08` relatively;
- open bulk energy by `4.378e-08` relatively.

The qT quantity remains stable well inside the packet-comparison tolerance. The less stable qL diagnostic receives no stronger precision claim. The maximum node-assigned shell fraction falls from `0.001857839` to `0.000706344`.

## Independent implementation and checkpoint scope

A freshly compiled, separately written C++ trajectory reproduces every scalar field in the archived trajectory CSV exactly on this host. Complete state arrays were not cross-compared between independent implementations, so the result is correctly described as scalar-diagnostic reproduction.

The checkpoint control gives exact full-versus-segmented equality only on a `40x12x12`, `tau_max=10` fixture. A production-scale Python checkpoint replay did not complete within the execution window. This is an infrastructure noncompletion, not a scientific failure; the full production trajectory itself was independently rerun in C++.

## Confirmed corrections

1. Define the exact Parseval spectral-energy weight and normalization convention.
2. Qualify the right-moving fraction as periodic-Parseval energy and distinguish it from open-line `E0`.
3. Replace “physical-time average” by “continuous-`tau` (structure-time) average.”
4. Use invariant qT bounds and their midpoint as the canonical machine-readable central record.
5. Bind the exact selector algorithm and path-independent content hash.
6. Qualify the checkpoint result as a small-fixture code test.
7. Replace the historical earlier-window box number by the fresh refined-window comparison.
8. Extend the open-line preparation control through `tau=55.4`.
9. Add the actual invariant timestep values, not only their span.
10. Preserve the post-selection status of the total-shell recertification.
11. Narrow independent-reproduction wording to scalar trajectory diagnostics.
12. Bind the 2% packet criterion to its exact prospective plan or contract record.

## False positives rejected

- the bulk really is exactly unexcited initially;
- the source-link tail is negligible;
- the coefficient-blind selector does not use qT or frequency-domain targets;
- the current sign and half-bond assignment are consistent;
- the periodic qT quantity is already correctly described as nonlocal;
- the finite-time result is not called an asymptotic transmission coefficient;
- invariant degeneracy ambiguity is negligible;
- qT timestep, cadence, and box controls pass;
- no remote detector calculation is required for the present bounded scalar claim;
- no Paper-2 or Paper-3 gate belongs in P07.

## Remaining uncertainty

- no separate-person or cross-platform reproduction was performed;
- a full production-scale checkpoint replay was not completed;
- the exact pre-result chronology of the total-shell threshold remains below 95% confidence;
- the global periodic qT diagnostic is not a local branch-current measurement;
- the time-domain calculation does not establish a continuum limit.

## Gate state

- P00 — PASS
- P01–P07 — PASS_WITH_CORRECTIONS_FROZEN
- P08 — AUTHORIZED NEXT
- authoritative v5 source — unchanged
- consolidated manuscript edit — deferred

## Exact next task

P08 will audit the finite vector-guide/core extension and decide, from its scientific value and referee risk, whether it should remain in the main text, move to the supplement, be narrowed, or be split into a separate handoff paper. It will verify the physical incident channel, WFM normalization, channel tracking, spectral amplitude and covariance, packet comparison, angular evidence, source rank, growing-core status, and all stated open gates.
