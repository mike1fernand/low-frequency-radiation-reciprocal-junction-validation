# P3 selected band versus the P1 frequency-1.2 qT candidate

8 November 2026. Independent analytic/source check, with bounded numerical controls. No scattering, Green-kernel, optimization or trajectory calculation was performed.

## Controlling sources and units

- Current P3 text: `P3_CURRENT.txt`, extracted from `PAPER3_WORKING_DRAFT2.pdf`, version 0, eleven pages, dated 7 November. The actual radial-spring Hamiltonian and nine undirected bond families are in section 2, Eq. (1), PDF pp. 2–3. The exact tangent symbol, registered selected band and field map are in section 4, Eq. (8) and (8a), PDF p. 7. The leading continuum comparison and Maxwell energy are section 6, Eqs. (11)–(13), PDF pp. 9–10.
- Actual P1 API: `numerical_probe/inputs/model_api_v1_1.py` in the previous coupling assessment. SHA256 `aa4128940dc7d39682af53c87cf04ae2532de4af880d25393ef4a4fe0b067980`. Its bond inventory and normalization are lines 70–108; `bulk_modes` explicitly sorts ascending eigenvalues. Current candidate uses a = m_bulk = K_axis = K_diag = 1, omega = 1.2. These are the matched bare-bulk units of the P3 calculation.
- Current P3 text SHA256: `7410acf62e2f99c23777dfebdfbcd777a1abd989be9e44a31b815b4fcfaed65d`.

## Exact deduction

Let B+ contain the three positive axial directions and the six face-diagonal families

    (1,0,0), (0,1,0), (0,0,1),
    (1,1,0), (1,-1,0), (1,0,1), (1,0,-1), (0,1,1), (0,1,-1).

The unit vector is bhat = b/|b|, so the face-diagonal dyads contain the factor 1/2. In the declared unit-mass/unit-stiffness units,

    D(k) = sum_b 2[1-cos(k·b)] bhat bhat^T.

Every dyad is positive semidefinite. The scalar inequality

    0 <= 2(1-cos theta) <= theta^2

therefore gives the *global* Loewner inequality

    0 <= D(k) <= C(k),
    C(k) = sum_b (k·b)^2 bhat bhat^T = |k|^2 I + 2 k k^T.

The final equality is an exact quadratic-polynomial identity for these nine families, not a small-k truncation of the inequality. The script checks all six monomial coefficient matrices with exact rational arithmetic. The continuum matrix C has ascending eigenvalues (rho^2, rho^2, 3rho^2), with rho = |k|. Eigenvalue monotonicity for real symmetric matrices gives

    lambda_0[D(k)] <= rho^2,
    lambda_1[D(k)] <= rho^2.

These are exactly the two lower branches used by P1's qT projector. Thus each lower-branch free frequency obeys

    omega_lower(k) <= |k|.

It follows that the lower-branch on-shell support at omega = 1.2 requires

    |k| >= 1.2.

It is disjoint from P3's registered selected band 0 < |k| <= 1. Consequently P3's Eq. (8a) and the particular band-restricted long-wave/Maxwell comparison cannot directly certify all of the new P1 candidate's outgoing qT radiation at omega = 1.2, nor manufacture an efficiency gain at that point.

This conclusion is only about the qT lower pair. The upper branch can have omega > |k|; the explicit adverse control demonstrates that extending the inequality to all three branches would be wrong. With common stiffness K and mass m, the dimensionless Bloch variable theta = a k_phys instead gives omega_lower <= sqrt(K/m)|theta| = c_T |k_phys|. The numerical cutoff 1 is therefore unit-dependent. Unequal axial/diagonal stiffness also requires a different continuum matrix.

## What the conclusion does not say

1. It does not invalidate P1's outgoing lower-two-branch qT definition or its calculated 73.6% coefficient. P1 directly evaluates its own bare on-shell branch matrices at omega = 1.2.
2. It does not imply that a harmonically forced *near field* has zero Fourier amplitude in |k| <= 1. Off-shell/reactive response can occur there. The support statement concerns the free bulk shell contributing to outgoing radiative power.
3. It does not say P3's exact coordinate/isometry or full energy/current identities cease to be algebraically meaningful for a driven low-k component. It says their particular selected band does not encompass the candidate's lower-branch radiative shell.
4. It does not establish a higher coefficient at another frequency. A different operating frequency/material/contact would need its own coupled response and appropriate wavevector support check.
5. It does not reopen microscopic DST or EM correspondence. Those questions are separate from this concrete model/band comparison.

## Maxwell energy: no additional reservoir or degeneracy factor

P3 section 4 explicitly states M†M = Q and retains the exact transformed stiffness Dy. It is an isometric coordinate change on the registered selected two-mode space, not two new physical modes. Its energy, source and complementary coordinates are retained.

P3's Eq. (13) is a leading-continuum transverse identity. For any Fourier v satisfying k·v = 0,

    E_m = -v_tau,
    B_m = i k cross v,
    |B_m|^2 = |k|^2 |v|^2,
    (|E_m|^2+|B_m|^2)/2
        = (|v_tau|^2+|k|^2|v|^2)/2.

The electric-like and magnetic-like terms are the kinetic and restoring portions of the same mechanical energy. They do not add a factor two to P1's power. P1 already sums both lower branches in Im_branch_0 + Im_branch_1. The exact finite-k Dy/current should not be silently replaced with the continuum E/B identity outside its approximation scope.

## Executed checks

`check_frequency_bound.py` and `FREQUENCY_BOUND_CHECK.json` record:

- Exact rational verification of all 54 entries of the six quadratic monomial coefficient matrices.
- 5,463 numerical wavevectors: 4,096 fixed-seed random BZ points, 343 structured grid points and 1,024 points in the selected rho <= 1 band.
- Independent bond-symbol versus actual API maximum difference 5.33e-15.
- The minimum eigenvalue of C-D is -2.52e-16; maximum lower-eigenvalue excess over rho^2 is 2.32e-16. These tiny signs are floating-point roundoff; the analytic inequality supplies the proof.
- Adverse controls reject using the same bound for the upper branch, omitting face-diagonal unit normalization, or quadrupling stiffness without rescaling the frequency relation.
- 128 complex transverse Fourier fixtures verify continuum mechanical E/B energy equality to 2.95e-16 relative discrepancy; doubling the energy fails.

No API file was modified. The matrix/energy fixtures verify the stated identities and conventions, not a new physical output efficiency or numerical quadrature enclosure.
