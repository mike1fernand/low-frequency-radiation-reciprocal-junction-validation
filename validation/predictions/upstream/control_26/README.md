# Directional control at frequency 2.6

SM7E uses the same 24 complete motifs, identical passive leads and line stiffness 2.25 at frequencies 2.6 and 2.9. Each calculation sums the powers of modes with negative axial group velocity to obtain the upstream power. The six-motif ramp at 1.2 is a separate configuration.

Use Python 3.12 with NumPy 2.3.5. From the repository root, choose a new directory for each run.

```sh
python -B validation/predictions/upstream/control_26/replay.py records ../paper1-directional-records
python -B validation/predictions/upstream/control_26/replay.py fresh ../paper1-directional-fresh --workers 3
```

The first command fits the saved Cartesian inputs, re-solves the loaded contact and compares method-specific outputs. Its real-axis summary reads the saved results. The second command repeats the bulk integrations and the same comparison. For each of the 256³, 384³ and 512³ grids it checks coverage of 16 disjoint groups of axial planes before solving the response. The five imaginary widths in squared frequency are 0.20, 0.15, 0.10, 0.075 and 0.05. Fits of degrees 2–4 extrapolate the finite-width directional powers to zero width. Direct real-axis integrations use four transverse resolutions at 2.6 and one at 2.9. Both integration routes use the same finite-junction solve.

| Calculation | Bulk power / incident | Upstream power / incident | Upstream power / bulk |
|---|---:|---:|---:|
| 2.6, archived real part and shell | 0.332496 | 0.171178 | 0.514828 |
| 2.6, Cartesian 512³ cubic fit to finite-width powers | 0.334653 | 0.172258 | 0.5147 |
| 2.9, direct real-axis at transverse resolution 768 | 0.720729 | 0.709743 | 0.984757 |

The old real part plus the new Cartesian shell at 2.6 gives a bulk fraction of 0.334655. The new real part plus the old shell gives 0.332492. These swaps identify the main numerical contribution to the difference between these two input sets. The span across fits and grids gives sensitivity information. The limiting values have no certified interval here.

Each run writes logs and `method-comparison.json` into the selected output directory. The real-axis script, shared contact code and stored kernel inputs are staged directly from the parent directory. Their project provenance and input digests are in `SOURCES.json`. The SHA-256 guard checks the historical trusted NumPy object file before loading it. The adapted Cartesian producer changes only the frequency in the predecessor's 2.9 script. All reference files for the comparison are included in the repository.
