# Final bounded mathematical audit

8 November 2026. **PASS: this additional skeptical pass found no new load-bearing mathematical or scientific error.** No scientific source or result was changed by this audit.

## What was checked

This pass read both complete point/continuous-box proofs, both complete exact-rational verifier sources, the kernel-centering record, the validated proof chain, the current-state record, and the saved R4 combined-replay receipt. It reviewed the model and normalization, resolvent identities, inverse-residual corrections, continuous-box coverage, real-frequency premises, and the scope of the conclusions.

The new `verify_audit_bindings.py` was actually executed successfully. It independently checked all 19 proof-chain SHA-256 bindings, all 144 complex G target entries, all 144 real Q target entries, the exact lead polynomial identity, upward rounding of both kernel caps, and the exact scalar propagation for both advertised certificates. Its machine-readable receipt is `MATH_AUDIT.json`.

This pass did not rerun the 64-vertex matrix certificate, the 222 saved-contour-job replay, or global contour integration. Those earlier checks have their own bound receipts. Reading their receipts and checking their hashes is distinct from rerunning them.

## Findings

1. **Physical point and normalization.** The certified point uses exact per-link motif stiffnesses 1/7 through 6/7, real frequency 6/5, line stiffness 169/100, and zero lead regulator. Each motif has eight links with its one common stiffness. The 19-coordinate uneliminated system and source normalization give the same incident-power transmission formula as the physical Green solve. The scaled expression is 2 k_star eta-dagger Q eta; the equivalent unscaled expression is 2 Gamma_q t-dagger Q t. No omitted factor of two or duplicate stiffness scaling was found.

2. **Real frequency and poles.** The exact lead identity gives Sigma_q = -97/100 - i 3 sqrt(133)/25 and Gamma_q = 6 sqrt(133)/25. Its polynomial radicand was independently verified with rational arithmetic. The certification uses the full uneliminated matrix; it does not assume an intermediate eliminated block has no pole. Its verified residual and perturbation norm inequalities establish invertibility of the matrix actually used. The archived exploratory small regulator is not silently substituted for the certified real-frequency limit.

3. **Kernel centering.** The cached G and Q used by the response verifiers agree exactly with the R4 target dyadics. Q also agrees with the saved binary64 branch sum. An exact rational sum of the two separately stored branch arrays need not equal their already-rounded binary64 sum; the certificate correctly binds the actual Q target instead. No unbudgeted recentering error was found. The tighter authenticated caps are safely rounded upward to eG = 1/250 and eQ = 1/500.

4. **Point proof.** The source/output resolvent identities, weighted Q norm, positive-definite exact LDL metric check, radical enclosure, and inverse-residual error corrections are consistent. The independent scalar recomputation reproduces exactly

   1135358893260054489/1550684920225000000 <= T_qT <= 4587863886935057969/6202739680900000000,

   approximately 0.7321660761 to 0.7396512062. This supports an approximate 73.6% candidate prediction. It does not prove the exact value rounds to 73.6% to one decimal place.

5. **Continuous-box proof.** The vertex argument applies to norms and an absolute affine energy scalar, which are convex auxiliary functions. It does not assume transmission extrema occur at vertices. The expansion retains the nonnegative square of the complete output change, so the cached-energy lower bound does not omit a negative cross term. The response correction includes the variable contact-dependent D_g. The exact scalar recomputation reproduces the uniform floor

   42767052725496289/61936276900000000 = 0.6905008642... > 0.69.

   This covers six independent relative motif changes of plus or minus 5%, with all eight links in a motif changing together. It does not cover 48 independent link errors, changed frequency, changed support, bandwidth, nonlinear pulses, or experiments.

6. **Upstream reuse.** Reuse of the bare-medium kernel enclosure is legitimate for the same fixed medium, contact support and real frequency; the new contact-dependent response is separately bounded. The old reference transmission floor is not transplanted to the ramp candidate. The raw ledger's pending floor flag and the successful combined physical/bridge replay describe different stages and are not contradictory.

## Consequence for the full validation package

The remaining completeness concern is documentary/reproducibility packaging: the compact delivered core contains selected R4 bindings rather than every upstream replay payload. The full deliverable should include the original R4 manifest payloads and inherited archive required by its replay, alongside the candidate-specific point/box certificates. This is not a newly discovered defect in the proven transmission bounds.

No new model modification, simulation, optimization, or repeated planning cycle is warranted by this pass. The manuscript should retain the existing distinction between an approximate 73.6% numerical estimate, its rigorous containing interval, and the continuous six-motif robustness floor.

## Reproduce this bounded check

Run from any directory, replacing CORE with the unpacked compact validation core and OUT with a writable receipt filename:

```bash
python verify_audit_bindings.py --core CORE --out OUT
```

Python 3 and NumPy are required. Acceptance for the scalar and target comparisons uses exact `Fraction` arithmetic; NumPy only reads the saved arrays in this audit script.
