# Candidate-specific source/output bound for the 73.6% point

8 November 2026. This is a new bound for the specified ramp, not a transfer of the old reference-profile certificate. The theorem explicitly assumes bare-kernel error bounds. Their R4 saved contour certificate and physical/bridge checks have now been freshly replayed, and the current cached matrix entries match its exact targets. The proof-chain record binds that authenticated premise to this new candidate-specific propagation.

## Model and exact point

Use the same six available complete paired-spring motifs and twelve normalized X/R source coordinates. Set actual per-link stiffnesses g_j=j/7, j=1,...,6; a=m_q=m_b=K_axis=K_diag=1; K_q=k_star=169/100; omega=6/5. Both bulk lower branches are included once. The incident wave is from the left. There is no parameter optimization in this check.

At the real physical frequency, the retarded line self-energy is

    Sigma_q = -97/100 - i (3 sqrt(133))/25,
    gamma_q = 6 sqrt(133)/25.

This proof uses zero lead regulator. The earlier exploratory solver used eta_lead=1e-12. The small lead-regulator change is not silently interpreted as an exact real-axis identity.

Let G_c and Q_c be the exact dyadic entries of the selected N192/M256 cached Green matrix G0 and Im_qT, restricted to the leading12 coordinates of REFERENCE_KERNEL.npz. The latter is entrywise identical to the sum of the two lower-branch arrays here. Its exact real symmetry and positive definiteness are checked by a rational LDL decomposition.

Use the full uneliminated system from SM5B,

    A_c = [[h/k_star, V^dagger G_c],[-V/k_star,I+U G_c]],
    b = sqrt(gamma_q/k_star) [e_left;0],
    y_c = A_c^-1 b,
    P_eta = [0,I],  eta_c=P_eta y_c,
    F_c^dagger F_c=2 k_star P_eta^dagger Q_c P_eta.

Thus the cached normalized qT output is a_c^2=||F_c y_c||^2. The full block construction avoids requiring an intermediate eliminated core inverse in the uncertainty proof.

## Exact source/output resolvent

Define scaled true-kernel differences

    DeltaG = k_star (G-G_c), ||DeltaG||_2 <= e_G,
    DeltaQ = k_star (Q-Q_c), ||DeltaQ||_2 <= e_Q,
    D_g = (1/k_star) [V^dagger;U].

The true full matrix is A=A_c+D_g DeltaG P_eta. With R_c=A_c^-1, put

    H=P_eta R_c D_g,
    J=F_c R_c D_g.

Subtracting the two full systems gives

    eta=(I+H DeltaG)^-1 eta_c,
    F_c y=F_c y_c-J DeltaG eta.

If h e_G<1, where ||H||_2<=h, and ||eta_c||<=z, then

    ||eta|| <= Z=z/(1-h e_G),
    ||F_c(y-y_c)|| <= d=j e_G Z,

where ||J||_2<=j. No sign, direction or proportionality of DeltaG is assumed. The true output is

    T_qT = ||F_c y||^2 + 2 eta^dagger DeltaQ eta.

Consequently, if a_low<=a_c<=a_high and a_low>d,

    (a_low-d)^2-2 e_Q Z^2 <= T_qT
      <= (a_high+d)^2+2 e_Q Z^2.

This is a complete norm-ball kernel-uncertainty statement at this fixed physical contact, not merely a mesh comparison. Positive semidefiniteness of the exact spectral measure provides the physical nonnegative power interpretation; the displayed perturbation bound already permits either sign of DeltaQ.

## Conservative outward-verifiable caps

The exact-rational residual verifier establishes the convenient conservative bounds

    ||H||_2 <= ||H||_F < 947/1000,
    ||J||_2 <= ||J||_F < 3/5,
    ||eta_c|| < 569/1000,
    4289/5000 < a_c < 8579/10000.

No computed matrix square root is used in the proof. Weighted norms are verified from exact rational quadratic forms, in particular

    ||J||_F^2=2 k_star trace(H^dagger Q_c H),
    a_c^2=2 k_star eta_c^dagger Q_c eta_c.

The inverse is enclosed using an exact-dyadic numerical preconditioner and a rational residual bound. Rational intervals for sqrt(2) and sqrt(133) enclose every nonrational coefficient in the true physical system. The numerical preconditioner is a candidate, not a trusted exact inverse.

With the stated mathematical premise caps e_G<=1/50 and e_Q<=1/200, the rational constants above imply

    0.720565188645... <= T_qT <= 0.751346331479...,

and therefore T_qT>0.72.

The authenticated tighter archived enclosure establishes the rounded-up caps e_G<=1/250 and e_Q<=1/500 for these exact cached matrices. The same argument gives

    0.732166076068... <= T_qT <= 0.739651206235...,

and therefore 0.732<T_qT<0.740.

The decimal endpoints are explanatory displays. The script checks the strict rational conclusions by exact Fraction comparisons. This does not establish that the actual exact value rounds to73.6% at one decimal place. It establishes a rigorous interval containing the approximate prediction under the authenticated kernel premises.

## Scope and remaining dependency

This is a fixed-point certificate, not a new bandwidth, packet or plus/minus5% coupling-box certificate. The separate 64-vertex evaluation reuses the existing SM5B theorem and is labelled independently. The authenticated kernel enclosures match the same frequency, source geometry, normalization, actual cached entries and complete lower-branch metric. Exact target matching confirms that Q_c is the saved floating64 sum used as the R4 certificate center, so no extra recentering allowance is required. The upstream contour replay checks saved outward-rounded jobs; it does not re-integrate all449,698leaves. Neither a repeated cached solve nor this response proof alone could authenticate that contour chain.

No Hamiltonian, graph, force law, original source or accepted numerical record is edited. The new interval proof concerns only the actual ramp specified above.
