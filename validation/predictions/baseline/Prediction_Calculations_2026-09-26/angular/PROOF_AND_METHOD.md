# Dressed selection rules and angular-collection method

## Exact model and source

This calculation uses the latest Paper1 reference model: a = m_bulk = m_line = K_axis = K_diag = 1, K_line = 1.69, omega = 1.2, and six complete segment motifs with g = (0.1,0.2,0.2,0.2,0.2,0.1). Source Hessians are independently assembled from the 48 physical springs and checked against the X+R decomposition. The selected archived N192/M256 Green kernel determines the dressed incident-unit-power source. All core arithmetic is in `core.py`; no archive implementation module is imported.

The full effective source is eta = (I+U G0)^(-1) V chi, with chi = G_D sqrt(Gamma_L) e_L. The redressed X-only and R-only controls replace the complete source Hamiltonian consistently: U, V and K_self all change. The same-full-state terms instead use eta_X = (I+U G0)^(-1) V_X chi and eta_R analogously; eta = eta_X + eta_R. These are different controls.

Sources/locations in updated main TeX: X/R basis and exact symmetry identities 865–904; segment Hessians 958–981; Green/operator levels and shell form factors 983–1012; dressed core 1014 onward; source normalization and power 1050–1085; full-state versus redressed controls 1259–1291. The exact source hashes are saved in `inputs/MANUSCRIPT_BINDING.json`.

## Axial qT exclusion

For arbitrary complex coefficients in the centered A1 source space, suppressing an irrelevant common translation phase, the physical Cartesian Fourier source is

\[
f(k)=\sum_s e^{-ik_xs}\left[
2\eta_{X,s}\cos(k_y/2)\cos(k_z/2)\hat x
-i\sqrt2\eta_{R,s}\{\sin(k_y/2)\cos(k_z/2)\hat y+
\cos(k_y/2)\sin(k_z/2)\hat z\}\right].
\]

This identity is checked against a separate literal Fourier sum on all 72 core Cartesian coordinates (`SELECTION_CHECKS.json`, maximum difference 3.7e-15). The bare radial node itself was already stated in main 602–650 and SM1 169. The additional full-response conclusion follows because **every dressed effective source remains in this A1 space**, including its X part.

On k = (kx,0,0),

\[
D=2(1-\cos k_x)\operatorname{diag}(3,1,1),\qquad
P_{qT}=\operatorname{diag}(0,1,1),\qquad f\parallel\hat x.
\]

Thus f^dagger P_qT f = 0 exactly. The result uses the complete lower-cluster projector at the doublet; it is unchanged by rotating an eigenbasis in that subspace. At omega = 1.2 the qT roots are kx = +/-1.287002217586569. The qL roots are +/-0.707483211779343 and have nonzero full-state axial source weight (approximately 2.1423 forward, 0.008146 backward in the Caroli source convention), so the source does not vanish altogether on the axis. These weights are not themselves normalized powers.

Any longitudinal segment-profile change that retains complete centered motifs preserves this theorem, including the six-factor box. Breaking a transverse motif symmetry is outside both this reduction and the existing box certificate.

## Mirror-plane exclusion and its important degeneracy qualification

On kz = 0, the source has f_z = 0 and the bulk matrix separates z from its x/y block. Consequently its overlap with the pure z-polarized shear eigenmode vanishes. This is a definite polarization selection rule, not the statement that all qT radiation vanishes in the plane.

For this particular equal-stiffness lattice, the lower pair is actually tied **throughout each coordinate plane**, not just on the axes: on kz = 0 its repeated eigenvalue is 4-2 cos(kx)-2 cos(ky). Therefore an arbitrary eigensolver's two individual lower-branch weights are not invariant on that plane. Report the complete qT projector sum; use “z-polarization exclusion” only when that physically specified polarization is explicitly selected. Four regular on-shell mirror points are checked, with zero z overlap and invariant cluster weight under basis rotation. Off-plane group-velocity splitting cannot be inferred from a chosen degenerate eigenbasis at the plane.

## Complete radial shell coverage

`inputs/REFERENCE_RADIAL_CENSUS.md` preserves the previously audited exact reference-symbol proof from the R4 authority. Only its shell geometry statements (sections 1–3) are used here, not any later stage/status language in that historical note.

For E = 36/25, all three shells lie strictly inside |k_i| < arccos(1/4). Along every ray each ordered eigenvalue crosses E exactly once, counting coincident lower roots with multiplicity. The radial derivative is positive and bounded away from zero. Thus this spherical parametrization omits neither a disconnected sheet nor an outer-zone pocket. The upper/lower gap also has a global positive analytic bound. This census does not prove completeness of an angular **group-velocity** preimage search.

At simple points the radial coarea formula, in the manuscript's Caroli source normalization, is

\[
dT_j = \frac{2\pi}{(2\pi)^3}
\frac{r_j^2\,|e_j^\dagger f(r_j n)|^2}
{|\partial_r\lambda_j|}\,d\Omega_n.
\]

The denominator is strictly positive for this reference symbol, so code uses the positive radial derivative. No absolute derivative or outward channel is discarded. Direction is the **signed** group velocity grad(lambda_j)/(2 omega), obtained from analytic symbol derivatives and independently checked by finite differences. Group-direction cones select signed v_x/|v| >= cos(Theta); wave-normal cones instead select k_x/|k| >= cos(Theta). Both are computed to expose their difference.

Gauss-Legendre polar and azimuthal quadratures avoid exact coordinate-plane/body-diagonal degeneracies. Their measure-zero boundaries are treated by the separate projector tests, not an arbitrary half assignment of cluster power. Full qT **cluster** totals converge very rapidly; the individual ordered-sheet contributions converge more slowly across lower-sheet crossings and are retained only as diagnostics.

## Finite-cone preimage search and limit

Before integrating a cone, the entire polar hemisphere is scanned for every azimuth and every one of the three radial sheets. The code finds all sampled inside/outside transitions and also checks for opposite-kx points whose signed group velocity would enter the forward cone. Integration assumes a single interval only after that global scan reports exactly one exit, no re-entry, and no opposite-kx cone point. Each boundary is then solved directly in polar angle; no paraxial angular map or assumed Theta exponent is used.

The finest production scan uses 1025 polar nodes and 128 Gaussian azimuths in a C4v sector. A separate check uses 1537 polar nodes and 192 midpoint azimuths, explicitly including the opposite-kx test. All three sheets, all seven cone half-angles (0.5,1,2,5,10,20,30 degrees), both grid families give the same simple support topology. A 240-start numerical stationary-ray search finds only the axial qT doublets and axial qL points. These checks found **no additional components**; they do not constitute an interval proof excluding arbitrarily narrow unsampled angular islands. No exact real-space dark cone or rigorously exhaustive infinitesimal-ray theorem is claimed.

Transverse reflections and quarter turns allow integration over one eighth of azimuth with factor eight; the source, symbol, measure and cone are all invariant under these transformations. Both signs of kx are separately integrated. The independent midpoint scan confirms that reversed group velocity of the lowest branch does occur near transverse planes, but its magnitude never enters the tested opposite axial cones; therefore simply assuming sign(v_x)=sign(kx) would have been unjustified.

## What the collector means

The output is an **ideal asymptotic angular modal-power collector**: a partition of the radiated channel power by propagation direction. It has no oscillator, impedance, absorption, backreaction, or coherent receiving-mode overlap. A real finite-distance detector or a coherent Paper2-style receiver may collect a different quantity and needs its own definition and loading calculation. Far-field/incoherent-channel interpretation precedes the shrinking-aperture interpretation. In particular, an exact wave-normal source node does not establish a finite-width dark cone.

## Ordinary elastic comparator

For an axial point force in an isotropic continuum with cT = 1 and cL = sqrt(3), the angular weights are proportional to sin^2(theta)/cT^3 and cos^2(theta)/cL^3. Integrating one forward cone, and suppressing common 2pi and force factors, gives

\[
W_S=\frac{(1-c)^2(c+2)}{3},\qquad
W_P=\frac{1-c^3}{3(\sqrt3)^3},\qquad c=\cos\Theta.
\]

Only W_S/(W_S+W_P) is compared. No fictitious incident line power is assigned to a prescribed point force. The common axial node is ordinary elastic source radiation, not evidence unique to paired links. The finite lattice junction differs quantitatively through source extent, phase, dressing and lattice dispersion; the comparison does not isolate any one cause by itself.
