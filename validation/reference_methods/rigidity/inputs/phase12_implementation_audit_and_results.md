# Paper 1 Phase 12 — implementation, audit, and results

**Date:** 2026-07-13  
**Model:** positive harmonic scalar spine coupled to a positive harmonic three-dimensional vector lattice.  
**Claim ceiling:** mechanical longitudinal-channel to transverse mechanical-wave energy transfer.

## Executive finding

The optimized Paper 1 program passes its stated numerical and independent-dynamics gates. At the reference point

\[
\omega=1.2,\qquad
(g_1,\ldots,g_6)=(0.1,0.2,0.2,0.2,0.2,0.1),
\]

the zero-broadening continuum calculation gives

\[
R_q=0.00357625,
\qquad
T_q=0.85979075,
\]

\[
T_L=0.01103646,
\qquad
T_T=0.12559653,
\]

\[
T_{\rm bulk}=0.13663299,
\qquad
\frac{T_T}{T_{\rm bulk}}=0.91922550,
\]

with

\[
|R_q+T_q+T_L+T_T-1|=1.30\times10^{-11}.
\]

A conservative envelope over principal-value samplers, reconstruction order, surface resolution, and the independent finite-broadening extrapolation is approximately

\[
\Delta T_{\rm bulk}=1.22\times10^{-4},
\qquad
\Delta T_T=1.35\times10^{-4}.
\]

These are deterministic method envelopes, not statistical confidence intervals.

## 1. Harmonic-mechanics correction

The word “void” is not a degree of freedom in this model. The exact local object is the signed segment strain \(\Delta q_i\). For one corner node,

\[
H_c=\frac g3\left[(u_x-\bar q)^2+\left(s_yu_y+s_zu_z+\frac{\Delta q}{2}\right)^2\right].
\]

The analytic force agrees with an independent finite-difference gradient to better than \(2\times10^{-10}\) in the test suite. At zero bulk displacement, the transverse force is exactly

\[
(F_y,F_z)=-\frac g3\Delta q(s_y,s_z).
\]

This is an inward or outward radial breathing force on four corner nodes, depending on the strain sign. No vacancy, missing site, density field, or externally imposed transverse drive is used.

## 2. Exact active-sector reduction

The original six-segment core contains 72 bulk coordinates. Symmetry restricts the source-reachable sector to 12 coordinates: one equal-axial pattern and one symmetric radial pattern per segment.

Fresh full-versus-reduced tests give:

- basis orthonormality residual: \(2.72\times10^{-16}\);
- source leakage outside the active sector: \(2.17\times10^{-16}\);
- local-stiffness leakage: \(2.00\times10^{-16}\);
- projected Green-operator leakage: \(1.25\times10^{-16}\);
- full/reduced Green-matrix relative difference: \(1.63\times10^{-15}\).

This is an exact symmetry reduction, not a truncation. It also makes explicit that Paper 1 accesses one radial cross-sectional source symmetry. It does not establish two independently controllable polarizations.

## 3. Principal-value reconstruction

Five independent spectral samplers were used: two differently shifted periodic-grid ensembles and three independently scrambled Sobol ensembles. At Jackson order 192:

- maximum pairwise matrix disagreement: \(6.62\times10^{-3}\);
- median pairwise disagreement: \(3.07\times10^{-3}\);
- maximum sampler deviation from the mean: \(4.47\times10^{-3}\).

Changing the Jackson order from 160 to 256 moves the mean matrix by about \(9\times10^{-4}\) to \(1.15\times10^{-3}\) relative to the order-192 result.

A separate direct calculation at

\[
(N,\eta)=(96,0.050),(112,0.040),(128,0.032),(144,0.026),(160,0.021)
\]

followed by a quadratic \(\eta\to0\) extrapolation differs from the Jackson result by \(9.82\times10^{-3}\) in matrix norm. Despite that conservative matrix-level difference, the resulting reference coefficients differ by only

- \(9.31\times10^{-5}\) in \(T_{\rm bulk}\);
- \(8.41\times10^{-5}\) in \(T_T\).

The direct binwise Hilbert transform was tested and rejected as load-bearing because modes landing near the target energy made it overly grid-sensitive. This rejection is part of the audit, not an omitted negative result.

## 4. Zero-broadening on-shell matrix

The imaginary spectral matrix was integrated directly over iso-frequency surfaces at grids from \(48^3\) through \(80^3\). Its norm changes from 1.41331 to 1.40947 over that sequence, with systematic stabilization. The final matrix is positive semidefinite; the smallest density eigenvalue is \(6.76\times10^{-12}\).

Before enforcing exact real reciprocity, the \(80^3\) surface quadrature has a relative forbidden-sector residual of \(1.11\times10^{-4}\). The final matrix is projected onto the exact real-symmetric harmonic sector. The projection is therefore controlled and smaller than the overall numerical error envelope.

At the reference dressed source, branch diagnostics are:

| Sorted bulk branch | Transmission | Longitudinal fraction | Transverse fraction |
|---|---:|---:|---:|
| 0 | 0.0621044 | 0.000291 | 0.999709 |
| 1 | 0.0635203 | 0.000167 | 0.999833 |
| 2 | 0.0110083 | 0.999950 | 0.000050 |

The two lower sorted branches are mechanically transverse-like at this frequency; the upper branch is longitudinal-like. The load-bearing observable remains their basis-invariant total transverse content.

## 5. Flux, passivity, and microscopic current

The scalar leads use an analytic retarded surface self-energy. The bulk broadening matrix is positive semidefinite and the reference closure is \(1.30\times10^{-11}\). The exact zero-coupling null gives

\[
T_{\rm bulk}=0,
\qquad
R_q\simeq1.1\times10^{-24},
\qquad
T_q=0.999999999988.
\]

For three generic wavevectors and every propagating bulk branch, the microscopic pair-current formula reproduces energy-per-cell times group velocity with a maximum relative error of \(4.40\times10^{-16}\).

## 6. Independent open-boundary dynamics

The time-domain implementation evolves the original real-space scalar/vector Hamiltonian with initially exact zero bulk coordinates and velocities. It does not call the Green-function, KPM, or surface-integration routines.

### Total converted energy

For a damped open calculation with a width-9 incident packet, the incident-spectrum-weighted frequency-domain prediction is

\[
T_{\rm bulk}^{\rm weighted}=0.13825438.
\]

The integrated microscopic coupling current gives

\[
T_{\rm bulk}^{\rm TD}=0.13822814,
\]

a relative difference of \(1.90\times10^{-4}\).

### Independent transverse share

An undamped causally isolated snapshot uses a basis-invariant kinematic Helmholtz projection. Its transverse fraction is

\[
0.91920772.
\]

For that packet, the snapshot values are

\[
T_{\rm bulk}^{\rm snap}=0.14019225,
\quad
T_L^{\rm snap}=0.01132645,
\quad
T_T^{\rm snap}=0.12886580.
\]

The corresponding incident-spectrum-weighted predictions are

\[
0.13894864,
\quad
0.01115289,
\quad
0.12779575.
\]

Relative differences are 0.895%, 1.56%, and 0.837%, respectively. The microscopic current integral and the bulk energy plus half of the residual interface-spring energy differ by only 0.0324%.

The snapshot L/T split is a kinematic mechanical diagnostic, not a claim about independently addressable photon polarizations.

### Numerical convergence

- timestep \(0.05,0.04,0.03\): coefficient range \(4.02\times10^{-5}\);
- transverse boxes 18, 22, 26: coefficient range \(6.15\times10^{-4}\);
- damping strengths 0.4, 0.6, 0.8: coefficient range \(1.58\times10^{-4}\);
- width-12 packet: spectrum-weighted discrepancy 0.173%;
- exact zero coupling: zero integrated bulk current.

## 7. Robustness

The exact/invariant calculation was swept over a nonzero parameter region:

- frequency \(0.6\le\omega\le1.8\): \(T_T\) remains positive, ranging approximately from 0.10491 to 0.21785;
- coupling scale: exact zero at zero coupling and nonzero conversion for every tested positive scale;
- weak-coupling log-log slope: 2.0315, consistent with quadratic leading behavior;
- slab lengths 2 through 10: \(T_T>0\);
- reference, uniform, soft-ramp, short-core, and edge-heavy profiles: \(T_T>0\);
- \(k_{\rm diag}/k_{\rm axis}\) from 0.75 through 1.25: \(T_T\) ranges from 0.09811 to 0.16908;
- lattice spacing \(a\) from 0.25 through 4 leaves the dimensionless transmission unchanged while dimensional long-wave speeds scale with \(a\).

The last point is a similarity/unit-scaling result, not a continuum-limit proof.

## 8. Audit corrections found during the loop

1. Removed an unnecessary “void dynamics” interpretation.
2. Corrected the geometry from four independent side oscillators to four corner nodes in one radial breathing symmetry.
3. Found and proved the exact 72-to-12 active-sector reduction.
4. Replaced a basis-sensitive two-transverse-label claim with invariant total \(T_T\).
5. Rejected an unstable direct binwise Hilbert transform.
6. Corrected negative-separation assembly: conjugate the spectral measure before applying the same retarded denominator; do not conjugate the denominator.
7. Enforced the exact real-symmetric reciprocal sector after recording quadrature residuals.
8. Replaced iterative lead self-energies with the analytic scalar-chain result.
9. Compared finite packets with spectrum-weighted coefficients rather than carrier-frequency values.
10. Derived symmetric interface-energy current accounting, including all vector force components.
11. Separated causal snapshot times from reflected-return times.
12. Removed an unsupported continuum-limit interpretation of lattice-spacing variation.

## 9. Claim that is now supportable

> In an explicit stable positive harmonic line-defect lattice, a prepared one-dimensional longitudinal incident wave transfers nonzero energy into outgoing three-dimensional mechanically transverse bulk waves.

This is a converged numerical existence result for the declared discrete model. It is not an analytic closed-form theorem, a homogeneous self-organization result, a two-polarization control result, or an electromagnetic identification.

## 10. Optimal next action

The next action is manuscript reconstruction from the locked claim, Hamiltonian, exact symmetry reduction, zero-broadening flux result, independent dynamics, and robustness data. No new conceptual mechanism is needed for Paper 1. Any additional work should improve presentation or tighten numerical envelopes rather than reopen the “void” interpretation.
