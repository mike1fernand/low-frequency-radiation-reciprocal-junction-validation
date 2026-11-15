# Exact core-response reduction and the direct-rank distinction

Let `E_c` embed the 72 source-corner coordinates into the infinite bulk, and let

\[
P=BB^\dagger
\]

be the rank-12 `A1` projector on the finite core.

The full bulk stiffness commutes with the full transverse symmetry action but does not preserve finite spatial support. Therefore the exact statement is not

\[
K_b E_c\operatorname{Ran}B
\subset E_c\operatorname{Ran}B.
\]

Instead, after the exterior bulk is retained through

\[
G_{cc}^R
=E_c^\dagger
[K_b-m_b(\omega+i0^+)^2I]^{-1}
E_c,
\]

symmetry gives

\[
[P,G_{cc}^R]=0,
\qquad
PV=V,
\qquad
[P,U]=0.
\]

Consequently,

\[
G_{cc}^R B
=B(B^\dagger G_{cc}^R B),
\]

and every source-driven core response can be computed exactly in the 12-dimensional block.

Fresh `24^3` and `32^3` finite-eta regressions reproduce the full and reduced fields and their effective self-energies below `7e-16` relatively.

## Direct source rank versus response-sector dimension

The line device has seven retained scalar coordinates for a six-segment slab. Hence the direct reduced source map has rank at most seven; the fresh computation gives

\[
\operatorname{rank}V_{\rm red}=7.
\]

This does not contradict the 12-dimensional exact `A1` sector. Repeated action of the projected response/local-core operator produces a 12-dimensional Krylov span at the reference test point:

\[
\operatorname{rank}
\operatorname{span}
\{(G_{\rm red}U_{\rm red})^nV_{\rm red}:n\ge0\}
=12.
\]

Accordingly, “the 12-dimensional sector containing the scalar source” is universally safe. “The sector reached by the source” requires an explicit reachability definition and test point.
