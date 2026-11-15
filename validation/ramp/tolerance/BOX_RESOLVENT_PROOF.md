# Continuous six-motif ±5% certificate for the ramp candidate

8 November2026. This is a new candidate-specific proof. Each of the eight physical links in motif j has stiffness

    g_j=(j/7)(1+delta_j), j=1,...,6, independently |delta_j|<=1/20.

Frequency6/5, line stiffness169/100, bare media and geometry remain fixed. This does not cover48independent spring errors, changed frequencies, packet energy or graph changes.

## 1. Exact cached full-system notation

Use A0,y0,F,P_eta,R0 from POINT_RESOLVENT_PROOF.md, with y0=A0^-1b and R0=A0^-1. Since the complete-motif Hessian is affine in the six errors,

    A_delta=A0+sum_j delta_j A_j,
    B_delta=R0(A_delta-A0),
    y_c(delta)=(I+B_delta)^-1 y0.

The full uneliminated19coordinate system is used throughout; no intermediate core elimination is assumed in this proof. The cached qT power is ||F y_c(delta)||^2 with the exact same positive dyadic Q_c as in the point proof.

All the following auxiliary functions are convex in delta:

    ||B||_F, ||F B||_F, ||B y0||,
    ||P_eta B||_F, ||P_eta B y0||,
    ||P_eta R0 D_g(delta)||_F,
    ||F R0 D_g(delta)||_F,
    ||R0 D_g(delta)||_F,
    2 |Re[(F y0)^dagger F B y0]|.

Here D_g(delta)=[V(delta)^dagger;U(delta)]/k_star is affine too. Norm maxima on all64vertices consequently bound these auxiliary quantities on the entire box. Transmission itself is not assumed convex, and its extrema are not assumed to occur at vertices.

## 2. Sharper cached-energy lower bound

The exact resolvent identity gives

    F y_c=F y0-F B y0+R_out,
    R_out=F B(I+B)^-1 B y0.

Let a=||F y0||, r>=||B||_2, w>=||F B||_2, v>=||B y0|| and r<1. Then ||R_out||<=wv/(1-r). Define

    ell >= 2|Re[(F y0)^dagger F B y0]|.

Expanding the squared output and retaining the nonnegative remaining square gives

    T_c(delta)>=a^2-ell-2a wv/(1-r).

This keeps the signed first-order energy response instead of replacing it by the generally larger norm of the output perturbation. With outward caps a_low<=a<=a_high,

    T_c(delta)>=a_low^2-ell-2a_high wv/(1-r).

## 3. Uniform kernel-error propagation

The source amplitude obeys

    ||P_eta y_c(delta)||
      <= z0+d_eta+w_eta v/(1-r)=Z_c,

where z0>=||P_eta y0||, d_eta>=||P_eta B y0|| and w_eta>=||P_eta B||_2. Write R_delta=A_delta^-1. Using R_delta=(I+B)^-1R0,

    ||P_eta R_delta D_g(delta)||
      <= h0+w_eta d_g/(1-r),
    ||F R_delta D_g(delta)||
      <= j0+w d_g/(1-r),

where h0,j0,d_g bound the corresponding affine R0D_g quantities in section1. These bounds need no convexity of R_delta itself.

The point source/output theorem can now be applied uniformly for every delta, using these uniform H,J,Z bounds and a lower cached-output amplitude L_c. For scaled bare-kernel errors ||DeltaG||<=eG, ||DeltaQ||<=eQ, let

    Z=Z_c/(1-H_cap eG),
    d=J_cap eG Z.

If H_cap eG<1 and L_c>d,

    T_qT(delta)>=(L_c-d)^2-2eQ Z^2.

This accounts for contact and kernel uncertainty together. The same physical bare-kernel enclosure applies to all contacts in this fixed support and medium; the contact-dependent response bounds are recomputed for this ramp.

## 4. Outward verified rational constants

verify_box_rational.py uses exact Fraction arithmetic and the point verifier's exact inverse-residual enclosure at every64auxiliary vertex. The proposed inverse itself is not trusted. All sqrt2-dependent matrices and the physical line radicals are enclosed by rational intervals; cached entries are exact dyadics. No matrix square root is needed. The exact rational LDL proof establishes Q_c>0.

The following conservative caps pass for every vertex and therefore throughout the box:

    r<41/250=.164,
    w<17/250=.068,
    v<103/625=.1648,
    ell<173/10000=.0173,
    w_eta<649/10000=.0649,
    d_eta<63/1250=.0504,
    h0<497/500=.994,
    j0<6289/10000=.6289,
    d_g<218/125=1.744.

Together with the point proof's z0<.569 and .8578<a<.8579, these imply

    T_c(delta)>=3634096757/5225000000
                =.695520910430622...,
    ||F y_c(delta)||>.8339,
    Z_c<.633,
    H_cap<1.13,
    J_cap<.772.

For the authenticated rounded-up contour-error caps eG<=1/250 and eQ<=1/500, exact Fraction propagation yields

    T_qT(delta)>=42767052725496289/61936276900000000
               =.6905008642115569... >.69.

Thus incident-to-qT conversion exceeds69% throughout the entire continuous six-coupling ±5% box, at the specified real frequency. This lower floor is deliberately conservative; it is not a claim that the true minimum is69.05%.

With only the coarser SM5B caps eG<=1/50,eQ<=1/200, the same proof gives a still-positive lower floor .6746176877464624... . The sharper floor above depends on the authenticated tighter contour enclosures, not merely the coarse statement in the paper.

## Status and limits

The R4 saved contour certificate and physical/bridge checks have been freshly replayed by the source-audit task, and the current cached12x12G andQ entries match its targets exactly. The proof-chain record binds these upstream checks to the new point/box scripts and receipts. The box verifier itself does not repeat449,698contour integration leaves; it propagates the authenticated error bounds through the candidate's different response.

No bandwidth, nonlinear, experimental, packet or48-independent-link guarantee follows. The separate nominal interval establishes73.2166%<T_qT<73.9652%, which contains the approximate73.6% prediction but does not certify its one-decimal rounding.
