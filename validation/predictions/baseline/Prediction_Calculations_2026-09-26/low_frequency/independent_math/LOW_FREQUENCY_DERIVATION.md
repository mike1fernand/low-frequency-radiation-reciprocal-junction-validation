# Independent low-frequency derivation

Analytic audit, 2026-09-26. Source: current Paper1 main-text projected response, paired-link Hessian, lead self-energies, and flux identities. No manuscript or archive input was changed. Only elementary coefficient arithmetic was evaluated; no spectral quadrature or dynamical simulation was run.

## Assumptions and notation

The physical complete paired coupling is retained, with a finite number of segments, positive spring constants, fixed lattice spacing, and identical scalar leads. The bulk is the infinite three-dimensional reference lattice with isotropic leading acoustic speeds cT and cL. The low-frequency bulk compliance uses the outgoing/decaying infinite-medium Green function, rather than the uniform global-translation solution.

Let G = g0^R(0) be the static projected bulk compliance; U, V, and Kself denote their actual values at a fixed coupling scale s > 0. Let e be the all-ones line vector, eL and eR the endpoint unit vectors, and p = (2,0,2,0,...,2,0)^T the reduced representation of a uniform axial core displacement. For any reduced core force f, p^T f is its total axial resultant. Define

\[
 Z_q=\sqrt{K_qm_q},\qquad \gamma_q=2Z_q\omega+o(\omega),
 \qquad L=K_D-K_q(e_Le_L^T+e_Re_R^T).
\]

L is the finite path Laplacian and Le = 0. All leading statements below take omega -> 0+ with geometry and material parameters fixed.

## 1. Static invertibility

G is finite and positive definite: the only bulk zero is the acoustic point, whose three-dimensional k^{-2} singularity is locally integrable, and a nonzero finite-core force cannot have an identically zero Fourier transform. Set

\[
 H=G^{-1}+U,\qquad
 A=L+K_{\rm self}-V^T H^{-1}V,\qquad
 D=(I+UG)^{-1}V.
\]

The quadratic form of the block matrix with diagonal blocks L+Kself and H and off-diagonal blocks V^T,V is the sum of the nonnegative line energy, the nonnegative paired-link coupling energy, and u^T G^{-1}u. A zero requires u=0; for the present positive complete segment profile, Kself is positive definite on the line coordinates, so q=0. The block matrix, and therefore its Schur complement A, are strictly positive definite.

This proof does not apply unchanged to the R-only/strain-only surrogate: its uniform line displacement can remain uncoupled, leaving a static zero mode.

## 2. Exact static force balance

Uniform common translation of the line and core costs no coupling energy. The explicit segment blocks therefore obey

\[
 Ve+Up=0,\qquad K_{\rm self}e+V^Tp=0.
\]

Using I-UH^{-1}=(I+UG)^{-1} and e^T L=0 gives

\[
 e^TA=-p^TD,\qquad
 p^TDA^{-1}e_L=-1.
\]

Thus a unit static device load transfers a unit net force to the bulk. This identity fixes the leading radiating resultant without evaluating G numerically.

## 3. Acoustic spectral measure and fixed-coupling transmission

At zero wave number the X patterns have resultant 2 x-hat and the R patterns have zero resultant. For rho_m = m_b/a^3, angular integration of the isotropic acoustic projectors gives

\[
 \operatorname{Im}_H g_{0,C}^R(\omega)
 =\omega b_C pp^T+o(\omega),
 \qquad
 b_{\rm qT}=\frac{1}{6\pi\rho_m c_T^3},\quad
 b_{\rm qL}=\frac{1}{12\pi\rho_m c_L^3}.
\]

Indeed the radial shell integral is omega/(2c_C^3), while the axial angular projector integrals are 8pi/3 for the two shear polarizations and 4pi/3 for compression.

Since G_D^R -> -A^{-1},

\[
 \chi_D^{(L)}=-\sqrt{\gamma_q}\,A^{-1}e_L+o(\sqrt\omega),
 \quad
 \eta_c=-\sqrt{\gamma_q}\,DA^{-1}e_L+o(\sqrt\omega),
 \quad
 p^T\eta_c=\sqrt{\gamma_q}+o(\sqrt\omega).
\]

Substitution into the exact Caroli power formula yields the universal leading coefficients

\[
 T_C(s,\omega)=4Z_q b_C\omega^2+o(\omega^2),\qquad s>0\ \text{fixed}.
\]

In particular,

\[
 T_{\rm qT}\sim\frac{2Z_q}{3\pi\rho_m c_T^3}\omega^2,
 \qquad
 T_{\rm qL}\sim\frac{Z_q}{3\pi\rho_m c_L^3}\omega^2,
 \qquad
 \frac{T_{\rm qT}}{T_{\rm bulk}}\longrightarrow
 \frac{2c_L^3}{c_T^3+2c_L^3}.
\]

No generic noncancellation assumption is needed for this leading bulk coefficient: static force balance forces a nonzero monopole resultant. The physical unit-incident-power source resultant has magnitude tending to 2sqrt(2Z_q); its squared magnitude 8Z_q independently reproduces these coefficients through the continuum point-force power formula.

For the remaining channels,

\[
 T_q=4Z_q^2\left|[A^{-1}]_{RL}\right|^2\omega^2+o(\omega^2),
\]
\[
 R_q=1-\left(4Z_q^2|[A^{-1}]_{RL}|^2
                  +4Z_q(b_{\rm qT}+b_{\rm qL})\right)\omega^2
       +o(\omega^2).
\]

The line coefficient retains coupling/profile dependence. Reflection tends to one. In the manuscript's reflection-amplitude convention, r_L -> 1; its overall phase convention does not affect this statement about power.

Reference units Z_q=1.3, rho_m=1, cT=1, cL=sqrt(3):

| Quantity | Leading value |
|---|---:|
| TqT / omega^2, fixed s > 0 | 0.275868568025952 |
| TqL / omega^2, fixed s > 0 | 0.0265454653351233 |
| Generated-bulk qT share | 0.912221450042867 |

## 4. Weak coupling first

Write U=sU1, V=sV1, Kself=sKself,1 and define

\[
 k_\Sigma=e^T K_{{\rm self},1}e=\frac83\sum_i g_i.
\]

For the uncoupled transparent line,

\[
 [\chi_{D,0}^{(L)}]_j=-\frac{i}{\sqrt{\gamma_q}}e^{i\kappa ja}.
\]

Therefore p^T V1 chi_D,0 = i kSigma/sqrt(gamma_q)+o(omega^{-1/2}), and the exact weak coefficients satisfy

\[
 \lim_{\omega\to0^+}\mathcal A_C(\omega)
 =\frac{b_C k_\Sigma^2}{Z_q}.
\]

For the reference profile sum(g_i)=1, these limits are 0.290197047233085 (qT) and 0.0279242238896761 (qL). They are limits of the weak coefficients, not the coefficients of omega^2 at fixed positive coupling.

The raw bulk transmission tends to zero in both iterated limits. The generated-bulk share also has the same point-force limit in the two orders, provided it is taken through positive s and omega. However, reflection and line transmission do not commute:

\[
 \lim_{s\to0^+}\lim_{\omega\to0^+}R_q=1,
 \qquad
 \lim_{\omega\to0^+}\lim_{s\to0^+}R_q=0,
\]

with reversed values for Tq. Thus it would be inaccurate to describe the raw bulk-transmission limits themselves as unequal.

## 5. Joint small-s, small-omega crossover

Let s=lambda omega with finite lambda >= 0. Projecting the device inverse onto the uniform zero mode of L gives the leading scalar denominator

\[
 -s k_\Sigma+2iZ_q\omega.
\]

Bulk static dressing enters at O(s^2), and nonuniform line modes and inertia change only higher orders in this joint limit. Define

\[
 F(\lambda)=\frac{(\lambda k_\Sigma)^2}
 {(\lambda k_\Sigma)^2+(2Z_q)^2}.
\]

Then

\[
 R_q\to F(\lambda),\qquad T_q\to1-F(\lambda),\qquad
 \frac{T_C}{\omega^2}\to4Z_q b_C F(\lambda).
\]

These leading Rq,Tq limits omit their O(omega^2) bulk-loss corrections; they must not be combined with the radiated term as an exact finite-frequency closure formula. The crossover scale is omega approximately s kSigma/(2Z_q), equal to 1.02564102564103 s in reference units. The omega^2 universal law is approached on the omega << s scale; the fixed-frequency weak law occupies s << omega.

## 6. Necessary-premise controls: X-only and R-only

The X-only surrogate retains both translation identities used above. Its static device matrix is also positive definite: a zero of the full static energy requires zero core displacement and a uniform line displacement (from L); the X coupling then forces that uniform displacement to zero. Consequently the X-only surrogate has the same leading fixed-positive-coupling qT/qL coefficients as the complete Hamiltonian, although its line coefficient and higher orders can differ.

The R-only surrogate has an exact uniform line zero mode and behaves differently. Let J be the N-by-(N+1) incidence matrix, (Jq)_i=q_i-q_{i-1}. Restrict to the radial core coordinates and define

\[
 C=G_{RR},\qquad S=\operatorname{diag}(4g_i/3),\qquad
 W=\tfrac12(S^{-1}+C)^{-1}>0.
\]

The actual positive coupling strengths, including any fixed scale s, enter S. Exact static elimination gives

\[
 A_R=J^T(K_qI+W)J,\qquad
 \eta_R=\sqrt2\,WJ\chi_D+o(\sqrt\omega).
\]

Here eta_R denotes the surviving radial-pattern coordinates in the R-only surrogate. Its axial-pattern coordinates vanish exactly: U_R has zero X rows and V_R has zero X rows, so the X rows of (I+U_R g0)eta=V_R chi give eta_X=0. This is distinct from the same-state eta_R decomposition inside the complete Hamiltonian discussed below.

The uniform part of chi is -i e/sqrt(gamma_q); balancing the endpoint forcing at the next order gives

\[
 J\chi_D=\frac{\sqrt{\gamma_q}}2(K_qI+W)^{-1}\mathbf1
          +o(\sqrt\omega),
 \quad
 \eta_R=\frac{\sqrt{\gamma_q}}{\sqrt2}
          W(K_qI+W)^{-1}\mathbf1+o(\sqrt\omega).
\]

The transverse dipole coefficient is guaranteed nonzero. Indeed

\[
 h=\mathbf1^T W(K_qI+W)^{-1}\mathbf1>0,
\]

because W(K_qI+W)^{-1} is positive definite. With the physical force normalization, the total moment tends to diag(0,M,M), where

\[
 M=-\sqrt{2Z_q}\,a h.
\]

Its continuum Fourier force is -i M k(0,n_y,n_z), to leading order. The P and S angular integrals are respectively 32pi/15 and 8pi/15. Hence

\[
 T_{\rm qT}^{R\text{-only}}
 =\frac{|M|^2}{60\pi\rho_m c_T^5}\omega^4+o(\omega^4),
 \quad
 T_{\rm qL}^{R\text{-only}}
 =\frac{|M|^2}{15\pi\rho_m c_L^5}\omega^4+o(\omega^4),
\]

and

\[
 \frac{T_{\rm qT}^{R\text{-only}}}{T_{\rm qL}^{R\text{-only}}}
 \longrightarrow\frac14\left(\frac{c_L}{c_T}\right)^5,
 \qquad
 \frac{T_{\rm qT}^{R\text{-only}}}{T_{\rm bulk}^{R\text{-only}}}
 \longrightarrow\frac{c_L^5}{4c_T^5+c_L^5}.
\]

At cL/cT=sqrt(3), these are 3.89711431702997 and 0.795798109812865. The same regular low-frequency expansion gives Rq -> 0 and Tq -> 1 (with Rq=O(omega^2)); the bulk loss is order omega^4. These are fixed-positive-coupling R-only limits, not the bare fixed-Q Bloch power scaling.

## Scope

The derivation concerns finite contacts in an infinite three-dimensional bulk with the physical complete coupling. It does not assert a uniform limit in contact length, physical radius, stiffness ratio, or spatial dimension. A numerically retained kernel must resolve its static part and the shrinking acoustic shells to test these results. The universal coefficient does not require a new numerical value of the static compliance; the static compliance is needed for the nonuniversal line-transmission coefficient and subleading terms.

The static dressed line displacement is not a prescribed Bloch wave, and its segment differences need not scale with omega. In addition, the manuscript's same-state source-sector vector eta_R=(I+Ug0)^{-1} V_R chi can acquire axial X-pattern components through core feedback. Thus a bare fixed-Q Bloch radial-source power law must not be assigned to that dressed eta_R contribution. The universal result above is for the total effective source and its complete qT/qL outputs.
