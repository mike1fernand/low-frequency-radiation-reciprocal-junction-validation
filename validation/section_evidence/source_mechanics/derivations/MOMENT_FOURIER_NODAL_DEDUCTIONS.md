# Exact deductions from the four-corner radial source

## 1. Moment-to-Fourier relation

For one active segment, let

\[
\mathbf f_c=-\frac{g\,\Delta q}{3}\mathbf t_c,
\qquad
\mathbf r_c=\frac a2\mathbf t_c,
\qquad
\mathbf t_c=(0,s_y,s_z)^T.
\]

The source is balanced, \(\sum_c\mathbf f_c=0\), and its first moment is

\[
M_{\alpha\beta}=\sum_c f_{c,\alpha}r_{c,\beta}
=-\frac{2ga\Delta q}{3}\operatorname{diag}(0,1,1).
\]

With the manuscript Fourier convention \(e^{-i\mathbf k\cdot\mathbf r}\),

\[
\widetilde{\mathbf F}(\mathbf k)
=\sum_c\mathbf f_c e^{-i\mathbf k\cdot\mathbf r_c}
=\sum_c\mathbf f_c
-i\sum_c\mathbf f_c(\mathbf k\cdot\mathbf r_c)+\cdots.
\]

The constant term vanishes, so for \(\mathbf k_\perp=(0,k_y,k_z)^T\),

\[
\boxed{
\widetilde{\mathbf F}^{\rm src}_{R,i}(\mathbf k_\perp)
=-i\mathbf M_i^R\mathbf k_\perp
+\mathcal O\!\left(g_i|\Delta q_i|a^3|\mathbf k_\perp|^3\right)
}.
\]

This exactly reproduces the leading term of the trigonometric source factor. It makes explicit that the vanishing zero-wave-number term and the leading linear transverse-wave-number behavior are the Fourier signature of the balanced dipole moment.

For a longitudinal Bloch component, \(\Delta q_i=\mathcal O(\kappa a)\); hence the complete local source scales as \(\mathcal O[(\kappa a)(k_\perp a)]\).

## 2. Exact source-factor norm and nodal set

Write

\[
\mathbf s_\perp(\mathbf k)=
\begin{pmatrix}
0\\
\sin(k_ya/2)\cos(k_za/2)\\
\cos(k_ya/2)\sin(k_za/2)
\end{pmatrix}.
\]

Then

\[
\boxed{
\|\mathbf s_\perp(\mathbf k)\|^2
=\frac{1-\cos(k_ya)\cos(k_za)}{2}
}.
\]

In the closed transverse first Brillouin zone, the exact zeros are:

- \((k_y,k_z)=(0,0)\); and
- the four zone corners \((\pm\pi/a,\pm\pi/a)\), which are periodically identified on the Brillouin torus.

Thus “symmetry-enforced nodal directions” can be replaced by an exact nodal-set statement. This is useful for the later shell-overlap audit but is not needed to prove local source generation.

## 3. Complete R-channel torque

For an arbitrary state,

\[
\mathbf F^{\rm tot}_{R,c}
=-\frac{2g}{3}B_c\mathbf t_c
\]

is collinear with \(\mathbf r_c=(a/2)\mathbf t_c\). Therefore

\[
\mathbf r_c\times\mathbf F^{\rm tot}_{R,c}=0
\]

for every corner separately. The complete R-channel torque about the segment midpoint vanishes even off the centered sector. Its resultant does not generally vanish off symmetry, and the complete X+R bulk torque can be nonzero off symmetry. The centered-sector qualification therefore remains essential for complete force closure and for zero transverse endpoint reaction.
