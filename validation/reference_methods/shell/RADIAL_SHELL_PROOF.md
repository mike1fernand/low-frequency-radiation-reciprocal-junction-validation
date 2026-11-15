# A model-specific global shell reduction for Paper 1

This closes a genuine part of the exact-lattice feasibility problem: the complete real-frequency shell at the reference frequency has a regular radial description, and the highest branch is separated from the lower-two cluster. It does **not** enclose any Green-function integral or establish the existing conditional error budgets.

The result uses the published bare bulk only: lattice spacing and bulk mass one, axial spring constant one, all six positive face-diagonal bond representatives with spring constant one, and energy (E=\omega^2=36/25). It is not asserted for the stiffness sweep, a defect-loaded medium, or arbitrary frequencies. The finite loading remains in the subsequent exact finite response construction.

## 1. Exact decomposition and exterior exclusion

Set (c_i=\cos k_i), (s_i=\sin k_i). The actual central-force bond sum gives

\[
D_{ii}=6-2c_i(1+c_j+c_k),\qquad D_{ij}=2s_i s_j\quad(i\ne j),
\]

and hence

\[
\frac D2=\operatorname{diag}(d_x,d_y,d_z)+ss^T,
\qquad d_i=2+c_i^2-c_i(1+c_j+c_k).
\]

Fix (0\le q\le1). If at least one component, relabelled (a=c_x), satisfies (a\le q), then every (d_i\ge1-q). To see this without interval sampling, write (b=c_y,c=c_z\in[-1,1]).

For (a\le0),

\[
d_x-\tfrac74=(a+\tfrac12)^2-a(2+b+c)\ge0.
\]

For (0\le a\le q),

\[
d_x-(2+q^2-3q)=(q-a)(3-a-q)+a(2-b-c)\ge0,
\]

while (2+q^2-3q-(1-q)=(1-q)^2\ge0).

For another component (b\ge0),

\[
d_y-(1-q)=(1-b)^2+q(1-b)+b(q-a+1-c)\ge0.
\]

For (b\le0),

\[
d_y-\tfrac74=(b+\tfrac12)^2-b(2+a+c)\ge0.
\]

The same holds for (d_z). Since (ss^T\succeq0), this proves

\[
\min_i\cos k_i\le q\quad\Longrightarrow\quad D(k)\succeq2(1-q)I.
\]

Choose (q=1/4), (R=\arccos(1/4)). In the principal Brillouin zone, outside the open cube (C=(-R,R)^3),

\[
D(k)-EI\succeq(3/2-36/25)I=\frac3{50}I.
\]

Thus every branch's (E=36/25) shell is strictly inside (C). This is a global exclusion statement about the exact symbol; it does not rely on a mesh search. In particular, the exterior resolvent is uniformly bounded by (50/3) in operator norm. Its integral has not been enclosed.

## 2. Radial monotonicity without assuming smooth ordered branches

For any unit vector (n), let (D_n(r)=D(rn)) and (r_C(n)=R/\max_i|n_i|). For (0\le r\le r_C(n)), each axial or face-diagonal bond phase obeys (|\theta_b|=|rn\cdot b|\le2R<\pi).

For (0<|\theta|<\pi),

\[
\frac{\theta\sin\theta}{1-\cos\theta}=\theta\cot(\theta/2)
\]

is positive, even, and decreasing with (|\theta|); its derivative for positive \(\theta\) is

\[
\frac{\sin\theta-\theta}{2\sin^2(\theta/2)}<0.
\]

Its continuous value at zero is two. With

\[
\beta=2R\cot R=\frac{2R}{\sqrt{15}}>0,
\]

the bond representation therefore proves the matrix inequality

\[
rD_n'(r)=2\sum_b\theta_b\sin\theta_b\,n_bn_b^T
\succeq\beta D_n(r).
\]

Here (n_b=b/|b|), and the sum uses the nine actual positive bond representatives. For every (r>0) in this cube, (D_n(r)\succ0). Indeed, a null vector must be orthogonal to every active bond direction because every bond contribution is PSD. A nonzero coordinate (k_i) activates its axial bond and forces (v_i=0); for each other coordinate (j), at least one of (k_i+k_j,k_i-k_j) is nonzero and has magnitude below \(\pi\), so its face-diagonal bond forces (v_j=0).

Integrating \(\frac{d}{dr}[r^{-\beta}D_n(r)]\succeq0\) gives, for (0<r_1<r_2\le r_C(n)),

\[
D_n(r_2)\succeq(r_2/r_1)^\beta D_n(r_1).
\]

The min–max principle applies directly to the ordered eigenvalues, including at crossings:

\[
\lambda_j(r_2)\ge(r_2/r_1)^\beta\lambda_j(r_1)>\lambda_j(r_1).
\]

Each ordered eigenvalue is continuous, starts at zero, and is at least (3/2>E) at (r_C(n)). Consequently each of the three ordered branches crosses (E) exactly once along every ray, counting coincident lower-branch roots with their multiplicity. This is a complete radial census; it does not assume three distinct scalar roots or differentiable eigenvectors at a crossing.

At a simple root, Hellmann–Feynman gives

\[
\partial_r\lambda_j\ge\beta E/r
\ge\frac{\beta E}{\sqrt3R}=\frac{24}{25\sqrt5}
\approx0.429325.
\]

At a degenerate root the same lower bound holds for every eigenvalue of the restricted radial derivative on that eigenspace, and hence for the corresponding one-sided eigenvalue slopes. This excludes radial tangency; it does not claim differentiability of the individual ordered lower branches in angular coordinates.

## 3. A quantitative gap above the lower-two cluster

Use ordered eigenvalues \(\lambda_1\le\lambda_2\le\lambda_3\). In the cube, put \(u_i=1-c_i\in[0,3/4)\), \(S=\sum_i u_i\), and \(U_2=\sum_i u_i^2\). The exact decomposition above gives

\[
d_i=S+u_i(2u_i-S),\qquad
 d_i-d_j=(u_i-u_j)(u_i+u_j-u_k),\qquad
 \|s\|^2=2S-U_2.
\]

The range of the diagonal entries is at most \(U_2\). To prove this, order a chosen pair so \(u_i\ge u_j\). If \(u_i+u_j\ge u_k\), then

\[
|d_i-d_j|\le(u_i-u_j)(u_i+u_j)=u_i^2-u_j^2\le U_2.
\]

If \(u_i+u_j<u_k\), then

\[
|d_i-d_j|\le u_i u_k\le(u_i^2+u_k^2)/2\le U_2.
\]

Let \(A=\operatorname{diag}(d_i)\). The rank-one interlacing inequality gives \(\lambda_2(A+ss^T)\le d_{\max}\); a Rayleigh quotient in the \(s\) direction gives \(\lambda_3(A+ss^T)\ge d_{\min}+\|s\|^2\). Consequently

\[
\lambda_3(D)-\lambda_2(D)
\ge2\bigl(\|s\|^2-(d_{\max}-d_{\min})\bigr)
\ge4(S-U_2)\ge S,
\]

because \(U_2\le(3/4)S\). At the origin this simply gives zero; at every other point of the cube \(S>0\). Moreover,

\[
\operatorname{tr}D=10S+2U_2-2S^2\le10S.
\]

On any of the reference energy shells, \(E\le\operatorname{tr}D\) because \(D\succeq0\). Therefore the highest branch has the exact uniform gap

\[
\boxed{\lambda_3-\lambda_2\ge E/10=18/125=0.144.}
\]

This bound separates the entire lower-two projector from the highest projector on all three shells. It allows crossings within the lower pair. It is a global analytic inequality, not the minimum of the finite direction sample.

As an independent limiting check, on an axis the eigenvalues are \(2(1-c),2(1-c),6(1-c)\). With two nonzero components and \(c_z=1\), the isolated third diagonal is below an active diagonal by \((D_{xx}-D_{zz})/2=(1-c_x)(1+c_y)>0\). These limits agree with the uniform gap and preserve the lower-pair degeneracy. More generally, on \(k_z=0\), \(D_{zz}=4-2c_x-2c_y\) and \(\det(D_{xy\text{ block}}-D_{zz}I)=4(1-c_x^2)(1-c_y^2)-4s_x^2s_y^2=0\). The other active eigenvalue is strictly higher away from the origin, so the lower pair is tied throughout that coordinate plane. A future primitive check must include a generic plane-touching cell; the axis and body-diagonal cases alone do not test this distinct geometry.

## 4. Consequence for a future certified integrator

The previously open possibility of missed closed sheets or radial folds at the reference energy can now be removed from the plan. A sphere cover with one root for each ordered branch on each ray is complete by the lemma. The exterior cube has no on-shell contribution and has a certified nonsingular resolvent bound.

At a multiple radial root with eigenspace projector (P), the matrix pole residue is obtained from the inverse of the restricted derivative (PD_n'(r_*)P\) on the eigenspace. It must not be replaced by arbitrary equal-rank weights. This exact-residue requirement is distinct from the already accepted basis-invariant nodal approximation. A complete implementation must justify and enclose these residues, including their angular transition regions. The proof above supplies radial invertibility but no implemented residue integrator.

What still needs to be implemented and proved before the existing exact-lattice conclusion follows:

1. A complete angular cover and outward root/residue enclosures, using the known radial census and the proved quantitative cluster gap.
2. Enclosure of the 17 source-projected moments, with all signed cross moments retained. This includes the regular exterior volume integral, interior principal values, and lower-two spectral residues.
3. Principal-value subtraction with the actual radial Jacobian and endpoint logarithms; asymmetric cells do not have automatic residue cancellation.
4. Outward quadrature errors plus exact-to-selected center displacement, assembled as a full 12-by-12 operator-norm error with every Toeplitz occurrence retained.
5. Verification of \(\|1.69\Delta g\|_2\le0.020\) and \(\|1.69\Delta Q\|_2\le0.005\), followed by the already implemented conditional finite-response bridge. Until then the exact-lattice floor is unproved.

Paper 2's complete projective-root census and whole-cluster current construction (Sec. IV, equations `projective-pencil` and `cluster-current`; SM SIV) motivate retaining every root and treating tied eigenspaces jointly. Operational Duration's outward event certificate (SM equation `supp-e3-outward-interval`) separately illustrates that root brackets alone do not prove completeness. Paper 3A Sec. V.B motivates transforming source and resolvent together under a checked symmetry. These are methodological transfers; the inequalities here, including the top-branch gap, are derived from Paper 1's actual nine-bond symbol.

## 5. Checks performed and their limits

`check_radial_shell.py` verifies the ten polynomial identities exactly with SymPy, computes the three positive constants with 160-bit Arb arithmetic, and performs independent floating-point falsification checks against the explicit bond sum. Its direction set includes coordinate axes, face/body diagonals, near-axis and near-plane directions, and a fixed spherical sample. The root solves check all three ordered branches for every direction.

The exact inequalities and the min–max argument prove the global result. The finite sample is a check for an implementation or derivation mistake, not a substitute for that proof, not the justification of the global gap, and not an integration-error estimate. The local proof does not change the selected Green matrices or the manuscript's numerical values.

Run with the supplied NumPy/SciPy/SymPy/python-flint environment:

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python -B check_radial_shell.py --bulk-source /path/to/current_numerics/analysis/math_resolution/tetra/combined_tetra_diagnostic.py --output /path/to/fresh/RADIAL_SHELL_CHECKS.json
```

The output binds the checker and accepted bulk-source hashes. `RADIAL_SHELL_CHECKS.json` records the actual run.
