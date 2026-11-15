# Independent review of the continuous ±5% coupling-box extension

8 November 2026. **PASS**, after a complete source/algebra review and one
isolated exact-rational replay. No load-bearing issue was identified.

## Verified conclusion and domain

For

    g_j=(j/7)(1+delta_j),  j=1,...,6,
    delta_j independently in [-1/20,1/20],
    omega=6/5,

with the same bare bulk, line and source positions, the candidate-specific proof
gives the uniform infinite-lattice model bound

    T_qT >= 42767052725496289/61936276900000000
         = 0.6905008642115569... > 0.69.

It consumes the authenticated R4 kernel caps `eG<=1/250`, `eQ<=1/500` for the
same cached centers, already checked in `../KERNEL_CENTER_BINDING_CHECK.json`.
Every motif's eight equal links vary together. This is a continuous six-variable
box, **not 48 independently varied springs**. It is not a frequency-band,
packet-average, global-optimum or experimental claim.

## What was independently reviewed

- Full `verify_box_rational.py` source, SHA256
  `dce88db6b788ec1b0341c97debc1db24deb4be80f333f64dc037815638662d60`.
- Its reused point verifier, SHA256
  `8cc7821ddc6c2adb20f9ff3079cda17519086479e11b1ddeb833053dcdcf3359`,
  previously reviewed and independently replayed.
- Affine dependence of the physical matrices on all six errors, including the
  dependence of `D_g(delta)` on the changed couplings.
- Every auxiliary norm's convexity/vertex extension and the absolute linear
  energy scalar's extension. No claim that transmission extrema occur at
  vertices is used.
- Exact-Fraction radical and inverse-residual corrections for all 64 vertices.
- The nominal energy lower bound, output-amplitude lower bound and uniform
  source/output resolvent caps.
- Final propagation of the fixed bare-kernel uncertainty over the entire box.

## Why the energy lower bound is valid

Let `B=A0^-1(A_delta-A0)`, `y_delta=(I+B)^-1 y0`, `v0=F y0`, and

    e = F B(I+B)^-1 B y0,
    ||e|| <= w*v/(1-r) = rem.

Then `F y_delta = v0 - F B y0 + e`, so

    ||F y_delta||^2
      = ||v0||^2 - 2 Re<v0,F B y0> + 2 Re<v0,e>
        + ||-F B y0+e||^2
      >= a_low^2 - ell - 2*a_high*rem.

The **whole** final squared term is nonnegative and may be discarded. Thus the
proof has not omitted an uncontrolled cross term involving `F B y0` and `e`.
The exact scalar result is

    T_cached >= 3634096757/5225000000
              = 0.695520910430622...,

which implies the conservative amplitude floor `L_cached>0.8339`.

## Uniform kernel-perturbation extension

The exact source bound is

    ||P_eta y_delta,c||
      <= z0 + d_eta + w_eta*v/(1-r) < 0.633.

Writing `C_delta=R0 D_g(delta)` makes the contact dependence explicit. The
resolvent identity gives

    ||P_eta A_delta,c^-1 D_g(delta)||
      <= h0 + w_eta*||C_delta||/(1-r) < 1.13,

    ||F A_delta,c^-1 D_g(delta)||
      <= j0 + w*||C_delta||/(1-r) < 0.772.

The norms of `C_delta`, its projected/output versions, and the other affine
auxiliaries are bounded at all vertices and extended by convexity. `F` remains
constant because the bare cached bulk spectral metric, frequency and source
basis are fixed throughout the coupling box. These bounds can therefore be
inserted in the previously reviewed source/output point perturbation theorem
uniformly in every `delta`.

## Executed independent checks

1. Copied both verifiers into `replay/`, preserving the frozen originals.
2. Replayed all 64 exact auxiliary-vertex cases. The isolated run passed in
   about 29.4 seconds; its complete receipt and output are retained.
3. Independently reconstructed the scalar energy floor, the three uniform caps,
   and final rational bound from the declared auxiliary constants. The final
   fraction matches exactly, as recorded in `INDEPENDENT_SCALAR_CHECK.json`.

The displayed decimals are explanatory. The accepted inequalities are exact
Fraction comparisons with rational enclosures of the irrational physical
coefficients. Floating-point values serve only as proposals or displays.

No source model, original manuscript, accepted numerical result or physical
claim outside this stated domain was modified by this review.
