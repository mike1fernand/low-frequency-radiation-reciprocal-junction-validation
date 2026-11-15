# Basis-free projector and finite-aperture deductions

## 1. Upper spectral projector without a lower-doublet eigenbasis

Let the real-symmetric bulk dynamical matrix have ordered eigenvalues

\[
\lambda_0=\omega_0^2,
\qquad
\lambda_1=\omega_1^2,
\qquad
\lambda_2=\omega_2^2,
\]

and suppose the upper branch is isolated:

\[
\lambda_2-\lambda_1>0.
\]

The Lagrange interpolation polynomial that equals zero on `lambda_0` and `lambda_1` and one on `lambda_2` is

\[
p_2(x)=
\frac{(x-\lambda_0)(x-\lambda_1)}
{(\lambda_2-\lambda_0)(\lambda_2-\lambda_1)}.
\]

By the spectral theorem,

\[
P_{\rm qL}=p_2(D)
=
\frac{(D-\lambda_0 I)(D-\lambda_1 I)}
{(\lambda_2-\lambda_0)(\lambda_2-\lambda_1)}.
\]

This formula does not select an eigenbasis inside the lower pair. It therefore remains well defined when

\[
\lambda_0=\lambda_1
\]

provided the upper/lower gap is nonzero. The lower-cluster projector is

\[
P_{\rm qT}=I-P_{\rm qL}.
\]

Fresh validation over 20,000 random Brillouin-zone points gave maximum error below `4.70e-14` on the declared operating-like subset with upper/lower squared-frequency gap above `0.3`.

## 2. Exact six-segment aperture

For the symmetric profile

\[
(g_0,\ldots,g_5)=(0.1,0.2,0.2,0.2,0.2,0.1),
\]

define

\[
\mathcal A(q)=
\sum_{j=0}^{5}g_j e^{iq(j+1/2)},
\qquad
q=(\kappa-k_x)a.
\]

Pairing terms about the midpoint at `j+1/2=3` gives

\[
\mathcal A(q)
=e^{i3q}
\left[
0.2\cos\frac{5q}{2}
+0.4\cos\frac{3q}{2}
+0.4\cos\frac q2
\right].
\]

The bracket is real and even. Consequently the aperture power is even in `q`, and the phase factor does not affect its magnitude.

The first positive zero is exact:

\[
q_0=\frac{2\pi}{5}.
\]

The direct numerical root of the normalized power half-maximum gives

\[
\operatorname{FWHM}_{\rm power}
=1.07590112974806.
\]

The half-amplitude width is different:

\[
\operatorname{FWHM}_{\rm amplitude}
=1.47206933405548.
\]

This distinction must be retained in captions and tables.

At the reference line wavelength

\[
\lambda_q/a=6.548686917692805,
\]

the six-segment interaction length is

\[
\frac{6a}{\lambda_q}=0.916214208346,
\]

while the midpoint-to-midpoint span is

\[
\frac{5a}{\lambda_q}=0.763511840288.
\]

These are finite-profile compactness diagnostics. They are not a continuum limit and do not establish profile optimality.

## 3. Radial-source norm and exact transverse nodal set

The dimensionless transverse source factor is

\[
\mathbf s_\perp(k_y,k_z)=
\begin{pmatrix}
\sin(k_ya/2)\cos(k_za/2)\\
\cos(k_ya/2)\sin(k_za/2)
\end{pmatrix}.
\]

Using elementary trigonometric identities,

\[
\|\mathbf s_\perp\|^2
=
\frac{1-\cos(k_ya)\cos(k_za)}{2}.
\]

Therefore the source vanishes exactly when

\[
\cos(k_ya)\cos(k_za)=1.
\]

Inside the closed first transverse Brillouin zone this gives the origin and the periodically identified zone corners. Away from those sets, the radial form factor is nonzero.

## 4. Regular witness implies an open shell patch

Let `Sigma_lambda(omega)` be a smooth branch shell near `k_*` and assume

\[
\omega_\lambda(\mathbf k_*)=\omega,
\qquad
|\nabla_{\mathbf k}\omega_\lambda(\mathbf k_*)|>0,
\]

with a nonzero branch gap and

\[
W_\lambda(\mathbf k_*)
=
|\mathbf e_\lambda^\dagger
\widetilde{\mathbf F}_{X+R}^{\rm slab}|^2>0.
\]

The implicit-function theorem gives a smooth local shell. The isolated branch projector and the finite-source form factor are continuous on that neighborhood, so `W_lambda` remains positive on a sufficiently small relatively open shell patch.

This proves only an open set of nonzero **source overlap**. The physical outgoing-power statement still requires the branch-resolved spectral measure, density-of-states factor, and flux normalization introduced in Section V.
