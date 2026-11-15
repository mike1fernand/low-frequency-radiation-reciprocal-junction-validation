# Fixed-contact qT conversion relative to its long-wave prediction

The right panel of the frequency-partition figure uses two established series for the same six-motif taper. All frequencies are in lattice units. The incident-power fraction is summed over both lower branches and all directions. Its normalization is the leading fixed-contact prediction, with line impedance 1.3, bulk mass density 1 and transverse speed 1. The coefficient is 2(1.3)/(3π).

The four open circles select frequencies 0.05, 0.1, 0.2 and 0.4 from the refined low-frequency calculation in SM7A. The 49 filled circles are the nominal N=80, M=256 scan at 0.6 through 1.8 in steps of 0.025. Those recorded frequencies are plotted individually. No result is interpolated between the two series. All values use fixed coupling scale s=1. The circle around the scan value at 1.2 marks the benchmark. Its ratio is about 0.3164. The finer N=192 surface/PV reference gives 0.3163; both round to the manuscript's 0.316. The guide at one is the low-frequency limiting ratio. A strict frequency interval with a given error would require a separate calculation.

## Reproduce

With Python 3 and Matplotlib (3.10.8 used for the distributed artwork), from the validation repository root run:

```sh
python validation/artwork/low_frequency_ratio/render.py --check-sources --output /tmp/paper1_low_frequency_ratio.pdf
```

The optional source check verifies every frequency, incident fraction and derived ratio in `points.csv` against the controlling records and confirms the three-digit benchmark against the finer calculation. The source records are already part of the validation repository. To rerender using the self-contained figure data outside that repository, omit `--check-sources`. To verify a separate copy of the validation repository, pass `--repository PATH` with `--check-sources`.

| Record (relative to repository root) | Purpose | SHA-256 |
|---|---|---|
| `validation/predictions/baseline/Prediction_Calculations_2026-09-26/low_frequency/inputs.json` | Reference units, stiffness profile, power definition | `710cf5638687b6c1ed8e3d0b23b7e6fc6c6a2f42c8b6f0478bcc65efa5b4d4d3` |
| `validation/predictions/baseline/Prediction_Calculations_2026-09-26/low_frequency/refined/rows.json` | Filter s=1, nominal frequencies 0.05, 0.1, 0.2 and 0.4 | `7ead4aa47b536cbb39318702ca11abc28946468dff52ffa76448004103894ac3` |
| `validation/reference/evidence/vi_results/ACTIVE85.json` | All 49 entries of families.frequency with the original N=80 surface and M=256 PV settings | `93ae058159343720a54fb272b2d9f9dd5b0cce839d652bf1e0e66f69b8174ef5` |
| `validation/reference/evidence/vi_results/REFERENCE_METHODS.json` | N=192 check at 1.2 using partition.TT_modal | `956c9f0dd9d7ae3f0c39cd59f7dc681f43e197575e0cfadecb6f1e37dcdd04b9` |
