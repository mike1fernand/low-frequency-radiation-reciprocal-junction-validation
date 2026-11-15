# Low-frequency-first angular collection law

New deduction and discrete-lattice checks, 26 September 2026.

## Statement

Fix a finite connected line contact with a nonnegative complete paired-spring profile and at least one positive spring strength. Hold all remaining parameters fixed in the reference infinite three-dimensional bulk. First take frequency to zero at a fixed forward energy-flow cone half-angle \(0<\Theta\le\pi/2\). Write \(c=\cos\Theta\).

Let \(C_T,C_L\) be the full/X-only low-frequency coefficients already derived:

\[
C_T=\frac{2Z_q}{3\pi\rho c_T^3},\qquad C_L=\frac{Z_q}{3\pi\rho c_L^3}.
\]

Then the actual dressed-junction cone powers, normalized by incident line power, obey

\[
T_{T,\Theta}=C_T\omega^2\,\frac{2-3c+c^3}{4}+o(\omega^2),
\qquad
T_{L,\Theta}=C_L\omega^2\,\frac{1-c^3}{2}+o(\omega^2).
\]

The geometric factors are fractions of each channel's own full-sphere power; they are not qT shares of mixed collected output.

The complete/X-only result permits the nonnegative profiles stated above. For the separately redressed R-only surrogate, the result here uses a fixed strictly positive profile and its dipole coefficients \(D_T,D_L\) from the static derivation:

\[
D_T=\frac{|M|^2}{60\pi\rho c_T^5},\qquad
D_L=\frac{|M|^2}{15\pi\rho c_L^5}.
\]

The analogous cone powers are

\[
T^{R}_{T,\Theta}=D_T\omega^4\,\frac{2-5c^3+3c^5}{4}+o(\omega^4),
\]
\[
T^{R}_{L,\Theta}=D_L\omega^4\,\frac{8-15c+10c^3-3c^5}{16}+o(\omega^4).
\]

For the reference profile \(D_T=0.001766271088378\), \(D_L=0.000453225372594\). These R-only coefficients depend on the coupling profile. Identical leading powers hold in the corresponding backward cone because all four limiting angular weights are even under direction reversal. Finite-frequency forward/backward asymmetry is not excluded.

## Derivation and limit order

For the reference lattice, uniformly over direction,

\[
D(k)=|k|^2 I+2kk^T+O(|k|^4).
\]

The upper/lower acoustic cluster gap is order \(|k|^2\). The upper projector and the lower cluster projector converge uniformly to the longitudinal/transverse continuum projectors. Hellmann–Feynman differentiation gives a leading radial group velocity for every branch; the direction correction is \(O(|k|^2)\), including the limiting resolved-branch velocity at a lower-sheet degeneracy. An ordered scalar branch need not be differentiable exactly at the tie; the bounded limiting velocities give the same leading direction. The within-qT eigenbasis itself need not have a unique limit.

The complete effective source has a fixed leading axial force resultant after incident normalization. Its angular powers are proportional to \(\sin^2\theta/c_T^3\) and \(\cos^2\theta/c_L^3\). The R-only source has a leading moment proportional to diag(0,M,M); its weights are \(\sin^2\theta\cos^2\theta/c_T^5\) and \(\sin^4\theta/c_L^5\). Integrating these weights with \(2\pi\sin\theta\,d\theta\), and dividing by the corresponding full-sphere integrals, gives the four factors above.

For a fixed cone boundary, the exact group-direction indicator tends to the continuum one away from a set of zero solid-angle measure. The finite source, resolved acoustic shells and bounded cluster weights allow that limit inside the angular integral. This establishes the **frequency-first** cone law without treating a finite-frequency wave-normal cone as an energy-flow cone.

Only after this frequency limit may the following narrow-cone expansions be taken:

\[
\eta_{T,\Theta}^{\rm full}\sim\frac12
\left(\frac{c_L}{c_T}\right)^3\Theta^2\longrightarrow0,
\]
\[
1-\eta_{T,\Theta}^{R}\sim\frac23
\left(\frac{c_T}{c_L}\right)^5\Theta^2\longrightarrow0.
\]

Here \(\eta\) is the qT share of the cone's qT+qL output and angles are in radians. The absolute qT powers behave as \(3C_T\omega^2\Theta^4/16\) and \(15D_T\omega^4\Theta^4/16\) in this ordered asymptotic sense. R-only can have nearly pure axial qT collection while carrying much less absolute qT power.

This is not a uniform joint \((\omega,\Theta)\to(0,0)\) theorem and not a proved complete small-cone law at the reference frequency 1.2. It is a quantitative link between the audited point-force/dipole limits and the chosen collection observable. It uses ordinary elastic radiation, not a new universal source mechanism.

## Checks

An independent analyst derived all four polynomials and the small-angle coefficients without using this script. The new numerical script uses the bound refined low-frequency Green matrices, an explicit nine-bond dynamical matrix, literal Fourier sums over all 24 corner force sites, and group-velocity boundary solves. It integrates full and R-only powers directly. Whole-polar scans cover all three sheets and the opposite hemisphere at each tested azimuth; these remain numerical searches, not interval certificates.

At 10°, 20° and 30°, the maximum relative finite-frequency correction among all full/R-only qT/qL powers is:

| Frequency | Maximum relative difference from leading cone law |
|---:|---:|
| 0.05 | 0.003594 |
| 0.0125 | 0.0002241 |
| 0.003125 | 0.00001401 |

Changing angular order 24 to 40 at the last frequency changes the powers by at most \(2.48\times10^{-9}\) relative. That is numerical convergence evidence, not a rigorous error enclosure.

The limiting qT cone shares at 10°, 20° and 30° are respectively:

| Source | 10° | 20° | 30° |
|---|---:|---:|---:|
| Complete / X-only | 7.3857% | 24.6053% | 43.2682% |
| R-only | 99.8684% | 99.4573% | 98.7149% |

They differ from the existing omega=1.2 table; the table here is explicitly a low-frequency limit. At a fixed cone and coupling, the R-only/full qT-power ratio vanishes as omega squared despite the R-only purity advantage.

## Reproduction

Keep this folder beside baseline/Prediction_Calculations_2026-09-26. Python and NumPy suffice.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python cone_limits.py --omega .05 --order 16
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python cone_limits.py --omega .0125 --order 24
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python cone_limits.py --omega .003125 --order 24
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python cone_limits.py --omega .003125 --order 40

Every output records its exact Green-input hash, frequency, order, branch boundary ranges and raw powers. The role of stored low-frequency matrices is explicit; this script does not regenerate their principal-value integral.
