# Calculations accompanying the November revision

These files recover the scientific code, machine-readable inputs and expected outputs from revision 18 of the November enhancement package. The source archive retained November 6 in its filename after its November 8 update. `SOURCE.json` gives the original archive hash, member prefix and member hashes. The paper under `../../manuscript/` gives the current equations and the scope of each conclusion.

With the reference Python environment installed, run this from the repository root. Choose a new output directory.

```sh
python -B validation/significance/replay.py --output ../paper1-significance-run
```

The driver copies the fixtures before running 18 checks. Each producer enforces its own recorded mathematical, numerical or physical controls. The driver also compares the two fixed-total profiles with the archived incident-power fractions and emitted-bulk shear shares. The output contains a log for each producer and the resulting JSON files.

| Calculation | Where to look | What the calculation tests |
| --- | --- | --- |
| Two positive-total profiles at the same total stiffness | `tradeoff_validation/TWO_PROFILE_RESULT.json` | Using the retained reference kernel, the endpoint-only profile raises incident-to-qT conversion from 12.566% to 15.858% while reducing the qT share of emitted bulk power from 91.94% to 86.62%. |
| Fixed-total static concavity and one-active-motif result | `uniform_integration/source_evidence/`, `post_revision_audit/` | Synthetic positive-compliance and exact static fixtures test the algebra used in SM7A. The manuscript supplies the general derivation. |
| Cluster projection and bulk units | `checks/`, `bulk_proposition/`, `recheck/check_prediction_corollaries.py` | Fourier-space projection, mass and stiffness scaling, and bound checks on finite fixtures. |
| Low-frequency coefficients and literature comparisons | `significance_audit/`, `comparison/` | Exact normalization, model-specific coefficient algebra, sphere averages, and sampled shell and curvature checks. |
| Dependence on source and profile | `recheck/` | Retained-kernel comparisons and controlled changes to the contact distinguish an emitted share from an incident conversion efficiency. The producers state the assumptions for each extension. |
| Conditional mechanical correspondence | `microscopic_correspondence/`, `correspondence_next/` | Checks of a specified mechanical parent and the homogeneous bulk. These fixtures do not derive the physical DST junction or identify an electromagnetic field. |
