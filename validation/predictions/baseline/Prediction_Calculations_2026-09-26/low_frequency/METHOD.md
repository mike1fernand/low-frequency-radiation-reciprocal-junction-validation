# Exact-lattice low-frequency quadrature and checks

The executable `low_frequency.py` uses the reference lattice parameters and the manuscript's orthonormal X/R core basis. It evaluates the full projected 12-coordinate Green matrix and all three on-shell branch matrices, then solves the complete dressed seven-site problem. No fit supplies the asymptotic coefficients.

## BZ integral and poles

Parity reduces the cube to its positive octant. Split that octant into three face pyramids and write k=t p, where p is a permutation of π(1,u,v), with t,u,v in [0,1]. The factor eight from parity and the Jacobian π³t² cancel the Fourier denominator (2π)³, leaving t² dt du dv. The parity-summed projected blocks are implemented directly with their signed X/R off-diagonal separation dependence.

For each angular ray and sorted eigenbranch j, bisection locates t_j with λ_j(t_j)=ω². Set R_j=t_j² P_j(t_j)/λ'_j(t_j), with P_j the parity-projected eigenprojector. The on-shell imaginary part is π times the angular integral of R_j. For the real part subtract 2t_j R_j/(t²−t_j²) from the matrix integrand, integrate the smooth remainder, and add R_j log[(1−t_j)/(1+t_j)]. This symmetric-pole subtraction handles the nearby negative root without a frequency-dependent radial mesh. At ω=0 the t² Jacobian cancels the integrable acoustic singularity.

The exact D is evaluated in stable trigonometric form. An independent builder sums all three axial and six face-diagonal undirected spring bonds. Its values agree at 137 random wavevectors to 1.90e-16 relative, and the analytic ray derivative agrees with the independent finite-difference derivative to 9.35e-11.

## Why the shells are complete for 0 < ω ≤ 0.4

In reference units write D=2ssᵀ+diag(d_i), s_i=sin k_i, c_i=cos k_i, and

\[
 d_i=2[(1-c_i)(2-c_i)+c_i((1-c_j)+(1-c_l))].
\]

If all c_i≥0 and m=min c_i, then D≥2(1−m)I. For the minimum-index diagonal this follows immediately; for another index the residual bound is 2(1−c_i)(1−c_i+m)≥0. If some cosine is negative, a negative-index diagonal is ≥2(2+c_i+c_i²)≥3.5 and every nonnegative-index diagonal is ≥2, so D≥2I. Therefore every shell below ω=0.4 lies in the small acoustic neighborhood, with |k_i|≤2 asin(ω/2).

At radial distance r=2ω the lower bound gives λ_min≥2[1−cos(2ω/√3)]≥(4/3)ω²−(4/27)ω⁴>ω². All bond phases throughout 0<r≤2ω satisfy |k·d|≤2√2ω<π. Hence the bond-sum derivative D'(r) is positive definite (the active bond directions span three dimensions), each sorted eigenvalue increases strictly, and each branch has exactly one root in that interval. The analytic bound excludes any remote shell outside it. The selected bisection roots therefore exhaust the low-frequency propagating shells. This argument is specific to the tested reference lattice and range; the script is not a general arbitrary-frequency shell finder.

## Independent checks

1. **Spherical shell:** use positive-octant μ=n_x and azimuth φ, the explicit bond-built D, and roots in physical radius r. The measure is r²/(π²λ'_j) dμ dφ. Angular orders are twice each face order. This checks the parameterization, Jacobian, D implementation, derivative and shell power. The parity projection/assembly code is shared. Only the clustered lower doublet is claimed accurate near degeneracy.
2. **Complex contour:** independently integrate direct inverses of the bond-built D along a lower-half-plane radial detour for t∈[0,4ω/π], then the real remainder. Detour z=(4ω/π)[u−ih sin(πu)] uses h=0.12 and 0.24. No shell eigenvectors, root finding, residue formula, principal-value subtraction, or artificial finite broadening enters this check. It tests the total retarded Green matrix at ω=0.2,0.05,0.003125. The parity projection and face angular quadrature are shared. Agreement of two heights is a convergence control, not a complex-analysis proof of all possible singularities for arbitrary parameters.
3. **Augmented solve:** solve the complete 19-coordinate line/core block problem and compare power with the independently coded Schur elimination at two frequencies and two coupling scales.
4. **Analytic normalization:** the point-force spectral angular integrals, translation identities and unit-power force conversion independently fix the low-frequency coefficients. Flux closure alone is not taken as evidence of correct branch partition or power normalization.
5. **Held-out geometry and surrogate controls:** change the profile and number of segments using principal submatrices of the exact finite-core kernel; separately replace the full local Hessian by the consistent X-only or R-only Hessian before redressing.

See `verification.json` and the grid-specific JSON files for actual residuals and numerical differences. The root, spectral-positivity and closure residuals are roundoff diagnostics; none is a substitute for a rigorous quadrature error bound.
