# R4 complete-domain error budget

The acceptance thresholds remain
`||(169/100)(g_exact-g_selected)||₂ ≤ 1/50` and
`||(169/100)(Q_exact-Q_selected)||₂ ≤ 1/200`, with the exact lower-pair
spectral matrix Q positive semidefinite. The selected dyadic matrix identities
are copied without alteration into `TARGETS.json`. Failure of a sufficient
upper bound does not establish failure of the physical floor.

## Actual acceptance test

1. Validate complete angular coverage, radial paths including the Green
   exterior, all 17 moment identities, and the specified source/model hashes.
   Reusing chart 1 for chart 2 counts two physical contributions. A chart 0
   tile may be reused under u↔v only with its distinct reflected domain counted;
   diagonal square tiles are counted once. Boundaries have zero measure.
2. Sum exact rational component endpoints with each physical multiplicity.
   Do not sum overlapping old/new refinements. A valid intersection of two
   enclosures for exactly the same complete region can sharpen that region.
3. Evaluate `budget_ledger.assess` on the complete 17 complex rectangles and
   every actual selected matrix entry. The tiny non-Toeplitz and b0 residuals
   in the cache remain in the comparison. Real and imaginary Green errors
   are combined as a complex modulus, not tested independently against0.020.
4. Accept only if both scaled bounds pass unchanged thresholds. Coverage,
   integral validity and exact-Q positivity remain separate prerequisites.
5. Replay the unchanged bridge producer against its hash-bound selected
   inputs in a fresh output directory, then bind the passing global budget
   and PSD premise to its conditional conclusion. The existing producer
   deliberately labels itself as a finite cached-matrix proof; an external
   linked conclusion must discharge its premises explicitly.

## More efficient norm bound

Artifact auditing requires identical mathematical source hashes and the
reviewed python-flint version. It preserves the recorded producer Python
identity but permits a different Python build string in the reader's
environment. Export verification records the producer and auditor runtime
identities separately, without treating a cross-build audit as same-platform
replay or independent-platform integration. Computation resume retains its
original strict full-identity requirement. Scientific ledger values remain
runtime-neutral so exact exported-ledger comparison needs no ignored fields.

Let C be the 17 interval midpoints assembled into the physical matrix,
T the actual selected dyadic matrix and R the nonnegative matrix of the
17 component-radius moduli in every Toeplitz occurrence. Then

`||M_exact-T||₂ ≤ ||C-T||₂ + ||R||₂`.

The uncertain entrywise envelope R has operator bound
`min(sqrt(||R||₁||R||∞), ||R||F)`. The fixed complex matrix C−T receives
an exact rational LDL* check of `u²I−(C−T)*(C−T)>0`. Floating SVD only
proposes rational u; it never establishes the acceptance decision.
The result takes the minimum with R3's direct full-difference envelope.
Thus no valid R3 bound becomes worse. The improvement preserves useful
sign/phase cancellation in the known center displacement.

Applied to R3's diagnostic dyadic centers, exact rational checks give scaled
center upper bounds 0.001279250908 for g and 0.000818200487 for Q. These
are bounds on those specified centers versus the selected matrices, not
bounds on their errors versus the unknown exact lattice integrals.

## Initial campaign and selective refinement

Starting at Q component density1/6000 and g component density1/1500 is an
exploratory cost decision. Even with the diagnostic center displacement,
the corresponding uniform worst-case bounds are approximately0.010114
and0.053861 respectively, above the unchanged thresholds. Therefore a
successful local integration at those tolerances does not itself close
the global budget. Actual achieved interval radii may be much smaller.

After the complete initial campaign, sum the actual intervals and run the
exact global gate immediately. If it passes, do not refine further. If it
fails, use `rank_regions` to identify the largest reducible uncertainty
contributions, and inspect their moment-specific row contributions. This
is a priority measure for reducing a sufficient bound, not a physical-error
estimate. Refine selected regions, retain unaffected proven regions, and
recompute the actual center displacement after every replacement batch.

`recommend_allocations` supplies a conservative sufficient allocation when
needed. It permits different moment radii based on their maximum row
multiplicities and positive rational cost weights, and rescales them to
satisfy every one of the12 row constraints. Default weights permit larger
errors for rarer moments. Default10% residual allowance is a computational
reserve for possible midpoint changes; it is not a theorem that midpoint
changes remain within that reserve. Final achieved intervals always face
the actual gate again.

With R3's diagnostic center reserve, uniform modulus density ceilings before
any extra reserve are approximately3.35678e−4 for g and7.49830e−5 for Q.
For g, an equal real/imaginary component tolerance must be smaller by√2.
These are sufficient uniform allocations only; nonuniform allocation and
actual achieved radii may close the bounds without meeting them locally.

`allocate_regions` converts a global17moment radius allocation to different
regional component densities. Positive area-times-cost weights set the
relative allowances; the exact normalization counts every physical reuse.
For complex g, the real and imaginary component allowances divide the
modulus allowance by a rigorous upper bound on√2. Thus their Euclidean
combination stays within the allocated modulus budget. These remain
planning allocations; actual output intervals face the final gate.

## Staged checks and replay command

- Before a campaign: source/target identities; positivity and rejection
  controls for exact LDL; coverage/multiplicity guards in the controller.
- After each completed region: exact endpoint format; hash and geometry;
  correct component count; all physical copies; resumable result recorded.
- After a replacement batch: no overlap, no omitted regions, recomputed
  aggregate radii/centers, and both actual matrix budget decisions.
- Before a floor statement: complete proof prerequisites and fresh bridge
  replay. If a gate remains open, report it without promoting the floor.

The unchanged producer is
`validation/current_numerics/analysis/resolution_execution/certificate/basis_invariant_n192_sobol23_m256_v1/finite_kernel_certificate.py`
inside the recovered extension1 repository. Its `staged_repo` sibling holds
the two hash-bound selected input arrays; the historical filenames there do
not identify their actual numerical resolution. `BRIDGE_READINESS.json`
checks its prepared source/input identities against the preparation record
and TARGETS. The bridge has **not** been replayed in this budget subtask.

From that producer's directory, invoke a validated python-flint0.8.0 Python:

```bash
python -B finite_kernel_certificate.py --repo-root staged_repo --output-dir /absolute/fresh/r4_bridge_replay --precision 160
```

Use an absolute fresh output directory outside staged_repo and preserve its process exit and
stdout. The producer checks all64 coupling-box vertices and both inverse
conditions; re-evaluating only the final decimal expression is insufficient.
Twenty-six exact/algebraic/rejection checks for the new budget module were run
and passed in `BUDGET_CHECKS.json`. The campaign adapter passed31 further
checks in `CAMPAIGN_ADAPTER_CHECKS.json`, including all17 canonical moment
endpoint mappings and altered-source/settings/certificate rejection.
The guarded driver was run against the unfinished campaign: it rejected
the unclosed budgets before creating a proof output directory or launching
the original proof. This fail-closed check is in `BRIDGE_GUARD_CHECKS.json`.

The portable original producer and exact selected arrays are preserved in
`bridge_inputs/`. Once the campaign closes, the preferred command is:

```bash
python -B budget/replay_bridge.py --root . --output-dir /absolute/fresh/r4_bridge_replay
```

The driver recomputes the actual global gate from all source-bound campaign
checkpoints and inherited records, requires complete refined coverage,
then replays the unchanged original64vertex proof and binds its result to
the complete error bounds and the inherited exact-Q PSD argument.
