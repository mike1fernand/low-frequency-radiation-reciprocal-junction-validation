# Matching-plane and weak-shell deductions

## Matching-plane invariance

Adding pristine lead cells changes only propagation phases. If `n_inc` cells are
added on the incident side and `n_L,n_R` are the total left/right padding, the
active field, reflection, and transmission amplitudes are compared after

\[
\chi_{\rm active}\mapsto e^{-i\kappa n_{\rm inc}}\chi_{\rm active},
\quad
r\mapsto e^{-2i\kappa n_{\rm inc}}r,
\quad
t\mapsto e^{-i\kappa(n_L+n_R)}t.
\]

Physical powers are unchanged. A wrong self-energy attachment is not a matching-
plane change and fails strongly.

## Weak positive shell patch

At a simple shell point, `omega_lambda^2(k)-omega^2=0` with nonzero gradient.
The implicit-function theorem gives a local smooth shell. The finite-source form
factor and branch projector are continuous there. A strictly positive source
overlap therefore remains positive on a sufficiently small relatively open
shell patch.

The fresh 49-point exact-shell neighborhood verifies this condition directly.
Combined with the positive on-shell measure, it yields a positive first-order
coefficient and

\[
T_{\rm qT}(s)=s^2A_{\rm qT}+O(s^3)
\]

away from thresholds, stationary points, and dressed poles.
