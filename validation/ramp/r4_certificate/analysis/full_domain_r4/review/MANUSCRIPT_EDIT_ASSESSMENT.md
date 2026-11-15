# Conditional manuscript edit assessment

No manuscript was edited in this assessment. The proposed replacements below
are authorized for consideration only **after** complete integral enclosures,
both original matrix budgets, exact-Q positivity and a fresh bridge replay
establish the stronger result.

## Verified current source

The authoritative recovered wave-mechanics manuscript is
`recovered_r3/paper1-validation-repository/manuscript/candidate/pass18-numerical-correction-rc4/`.
`manuscript/CURRENT_SOURCE.json` selects this candidate and the repository
`CURRENT_STATUS.md` agrees. I recomputed both selected file hashes:

- `paper1_main.tex`: `01acf82a2d5d3a01850ae25fdf5f155dc8ebecb3d6f20545f769ae41ab48b214`.
- `paper1_supplement.tex`: `66d8fe655b73f168110b3cb803504a8618def7d214ed4a0d50f891f2b1f7e3ac`.

Both match the authority record. Its target is Physical Review E. This is the
specific mechanical-junction paper, not the distinct observer-time manuscript
mentioned elsewhere in the user's research history.

## Existing wording that would become outdated

The main robustness paragraph (source line1424) presently certifies
`T_qT>0.108` only for the selected finite Green matrices, then says the
certificate does not enclose infinite-lattice integration error. The conclusion
(line1465) repeats the approximately10.8-percent computed-response floor.
The abstract (line55) describes only a perturbation bound for the computed
response.

The supplement subsection `sec:sm-coupling-box` begins at line430. Its opening
paragraph limits the whole-box result to fixed matrices; lines472–489 give the
finite-matrix certificate, numerical comparison table, and currently unproved
need for g0/Q integration-error bounds. The existing finite-matrix results
and comparison table remain correct and should be retained with their labels.
The final unresolved-dependency paragraph requires replacement if the new
bridge succeeds. The conditional bridge itself is not presently derived in
this subsection, so simply deleting the caveat would be insufficient.

## Minimal useful update after a successful bridge

Replace the main robustness paragraph with a paragraph of comparable or
shorter length. A suitable draft, contingent on the exact final bridge, is:

> The reduced response also gives a quantitative bound for simultaneous
> coupling changes. At the reference frequency, write
> $g_s=g_s^\star(1+\delta_s)$ with $|\delta_s|\leq0.05$ for all six segments;
> each segment's eight links vary together, with bare media and geometry
> fixed. An inverse-perturbation estimate in the qT output norm controls every
> point of this box, including reciprocal feedback. Validated integration of
> the infinite-lattice Green and lower-cluster spectral matrices, combined
> with the whole-box bound, certifies $T_{\qT}^{\rm modal}>0.1037$, compared
> with the nominal numerical value $0.126$. The Supplemental Material gives
> the bounds, and the validation archive supplies the executable certificate.

In the conclusion, replace only the sentence about the computed10.8-percent
floor with a statement that the infinite reference lattice retains more than
10.3-percent incident-qT conversion throughout the stated simultaneous
plus/minus5-percent coupling box. The abstract's current perturbation-bound
sentence can be replaced, not lengthened, to mention this certified robustness
result; this is worthwhile because it is a substantive new result.

In the supplement, retain the finite-matrix table and its original provenance.
Add the actual complete-integral operator bounds after comparison with the
selected cache, define the perturbations and1.69 scaling, and give the exact
bridge inequality and its evaluation. State exact Q positivity from its
positive lower-pair spectral measure, including tied modes. Detailed contour
proofs, adaptive subdivision logs and replay instructions belong in the
validation addition; they do not require another main-text safeguard section.

## Meaning and limits of the stronger result

The current12.6-percent headline is the nominal **monochromatic incident-power
coefficient** at omega1.2. The packet section separately studies a finite
frequency spectrum and a finite-time stored-energy estimator. Certifying the
reference-energy whole-coupling-box floor does not certify that packet
average, create an exact finite-window/asymptotic identity, cover an entire
frequency band, permit arbitrary bondwise errors, or identify the elastic
mode as electromagnetic. The existing definitions make these distinctions;
they need not be repeated as a long list in the paper.

The new floor is numerically lower than the stored-matrix floor because it
also accounts for integration uncertainty. Its scientific improvement is the
exact infinite-lattice guarantee under the stated assumptions, not a higher
conversion value. All numerical benchmark values can remain unchanged.
