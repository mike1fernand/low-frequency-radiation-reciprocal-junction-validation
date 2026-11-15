# Independent review of the candidate-specific point certificate

8 November 2026. **PASS**, within the precise mathematical/model scope below.

## Conclusion

For the declared passive six-motif ramp with actual per-link stiffnesses
`g=(1,2,3,4,5,6)/7`, exact frequency `omega=6/5`, unit bulk parameters and
`Kq=169/100`, the newly derived point certificate composes with the replayed R4
bare-kernel enclosures to give

    0.7321660760687072... <= T_qT <= 0.7396512062343026...

The exact rational endpoints are retained in the receipts. In particular,
`0.732 < T_qT < 0.740`. This supports a rigorous lower bound above 73.2% at the
specified point. It does **not** prove that the exact value rounds to 73.6% at
one decimal place, nor certify a frequency band, packet average or coupling box.

## Actual independent work

1. Read the current SM5B definitions and the complete new rational verifier and
   source/output proof. Independently derived the perturbation identity, checked
   both full-block signs, `k_star` scaling, and unit-incident normalization.
2. Checked the real-axis lead explicitly:
   `Sigma=-97/100-i*3*sqrt(133)/25`, `gamma=6*sqrt(133)/25`. The new proof uses
   zero lead regulator. The exploratory 1e-12 regulator is not treated as exact.
3. Independently reconstructed an exact Fraction LDL of the cached dyadic Qc.
   All twelve pivots are strictly positive and reconstruction is exact. The
   cached `Im_qT` equals the stored floating-point lower-branch sum entrywise.
4. Independently evaluated the final scalar propagation with exact Fractions.
   Safe rational Green caps `eG=1/250`, `eQ=1/500` and the verified conservative
   matrix caps yield the displayed interval.
5. Independently inspected and replayed `verify_point_rational.py` in a new
   review directory. Source SHA256:
   `8cc7821ddc6c2adb20f9ff3079cda17519086479e11b1ddeb833053dcdcf3359`.
   Every acceptance comparison is an exact Fraction comparison. NumPy supplies
   only a proposed inverse; its residual is checked exactly. Irrational physical
   coefficients are enclosed by integer-isqrt rational intervals.
6. Independently compared every complex Gc entry and every real Qc entry with
   the R4 exact `TARGETS.json` entries. Both centers match exactly. The target Q
   is already the floating-sum dyadic center, so no additional centering shift is
   needed here. The diagnostic difference from the exact sum of branch dyads
   is not silently ignored or used as a different center.
7. Read the actual R4 saved-evidence replay receipt: all budget, audited ledger,
   physical-control and exact-lattice-bridge stages passed; the ledger was
   reproduced exactly. Verified its ledger and targets hashes and checked the
   exact rational operator bounds are below the safe caps above. This consumes
   the separately executed saved-evidence replay; I did not rerun global contour
   quadrature or independently reconstruct its entire physics proof.

## Core derivation checked

With `A=Ac+Dg*DeltaG*Peta`, `R=Ac^-1`, `eta=Peta*y`, and
`H=Peta*R*Dg`, subtraction gives

    eta = (I+H*DeltaG)^-1 eta_c.

With `J=Fc*R*Dg`, the cached output difference is

    Fc*y - Fc*yc = -J*DeltaG*eta.

Hence `Z=z/(1-h*eG)` bounds the source norm and `d=j*eG*Z` bounds its output
amplitude change. The true power differs from `||Fc*y||^2` by
`2*eta^dagger*DeltaQstar*eta`; its magnitude is at most `2*eQ*Z^2`.
The lower/upper squared-amplitude bounds therefore follow without assuming a
sign or direction for either kernel error.

The matrix-square-root issue is handled correctly: the verifier proves Qc SPD
and uses exact weighted quadratic forms, rather than trusting a rounded
factorization. Specifically `||J||_F^2=2*k*trace(H^dagger*Qc*H)`. The inverse
and source residual corrections in the verifier conservatively cover the
physical sqrt(2), sqrt(133) coefficients and the proposed numerical inverse.

## Skeptical-pass findings

No load-bearing error was confirmed. One wording erratum was sent to the
author: use "stated/archived" instead of "published" for the old premise caps,
because Paper One's publication status has not been established. This does not
alter the mathematics or require another computation cycle.

The rigorous interval is a specified-model result, not experimental validation,
a proof of microscopic DST ancestry, a global optimum, or a calibrated 95%
statistical confidence statement. It is separate from the new tolerance and
packet checks being performed by the other workers.
