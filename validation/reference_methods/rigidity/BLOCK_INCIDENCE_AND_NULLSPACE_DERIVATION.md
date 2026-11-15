
# Finite-truncation incidence form and nullspace deductions

## 1. Compatibility form

For any declared finite truncation, collect all scalar-line and vector-bulk coordinates in a real vector \(x\), and collect the canonical momenta in \(p\). Let \(B\) be the compatibility matrix whose rows are the signed linearized extensions of:

- line bonds;
- bulk axial bonds;
- bulk face-diagonal bonds; and
- endpoint-to-corner source links.

Let \(W\) be the diagonal matrix of the corresponding positive spring stiffnesses. Then the exact quadratic Hamiltonian is

\[
H=\frac12 p^{\mathsf T}M^{-1}p
 +\frac12 (Bx)^{\mathsf T}W(Bx)
 =\frac12 p^{\mathsf T}M^{-1}p+\frac12\|W^{1/2}Bx\|^2.
\]

Consequently,

\[
K=B^{\mathsf T}WB,\qquad K=K^{\mathsf T}\succeq0,\qquad \ker K=\ker B
\]

when all retained spring weights are positive. This is a compact proof of nonnegativity; it does not replace the physical bond inventory or the explicit source-link equations.

## 2. Conditional exact kernel theorem for the combined free truncation

Assume:

1. the finite bulk central-force framework is infinitesimally rigid, so its compatibility kernel consists of the six rigid motions
   \(u(r)=t+\Omega\times r\);
2. the finite scalar line is connected by positive axial springs, so zero line extension gives \(q_i=c\) at every line site;
3. one active line segment contains the complete four-corner, eight-link source motif, with positive source stiffness.

For a source corner at

\[
r=(x_m,s_y a/2,s_z a/2),\qquad s_y,s_z=\pm1,
\]

zero extension of the left and right source links gives

\[
u_x=c,\qquad s_yu_y+s_zu_z=0
\]

for every sign pair. Substituting the bulk rigid motion gives:

- from \(u_x=c\) for all four corners, \(\Omega_y=\Omega_z=0\) and \(t_x=c\);
- with those values, \(s_yu_y+s_zu_z=s_yt_y+s_zt_z=0\) for all signs, hence \(t_y=t_z=0\);
- \(\Omega_x\) remains free.

Therefore

\[
\ker K_{\rm combined}
=\operatorname{span}\{\text{common axial translation},\ \text{bulk rotation about the line}\}.
\]

The proof is exact under the stated hypotheses. It does not identify the unspecified historical finite fixture. The package's symbolic calculation gives source-constraint rank five on the seven rigid/line variables, and representative one- through four-segment fixtures reproduce nullity two. Setting the source stiffness to zero restores the disconnected nullity seven.
