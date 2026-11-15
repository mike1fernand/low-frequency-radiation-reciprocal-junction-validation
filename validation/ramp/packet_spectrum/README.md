# Stationary spectrum prediction for the new ramp packet

8 November 2026. New profile: actual per-link motif stiffnesses
`g = (1, 2, 3, 4, 5, 6) / 7`. The same six available motifs and all 48 springs
are retained. This calculation uses the original 120-site real Gaussian input
and already computed bulk kernels; it is not a time-domain simulation.

## Main result

The matched broadband packet prediction is approximately **71.63% qT**, rather
than the monochromatic 73.6% at the carrier frequency. The pulse has finite
bandwidth and its energy must weight the whole transmission curve.

| Existing kernel family | New ramp, in-band contribution divided by full packet energy | Original reference replay |
|---|---:|---:|
| Surface N128 / M256 PV | 0.7163304996129106 | 0.12768113293380567 |
| Surface N160 / M256 PV | 0.7163064796290152 | 0.12767301846182993 |
| Tetrahedron N96 | 0.7161938909753368 | 0.12756756192214755 |
| Tetrahedron N128 | 0.7162026653472909 | 0.12759466873725295 |
| Tetrahedron N160 | 0.7162269562660165 | 0.12760941962540592 |

The five central predictions span 0.0137 percentage points. This is observed
method/grid sensitivity, not a rigorous true-error bound. Tetrahedron matrices
retain the declared `1e-10` regulator; the source metadata remain included.

For surface N160, the full-energy in-band partition is:

- qT: 0.7163064796290152
- qL: 0.06633223302091823
- forward line: 0.08512617236127365
- reflected line: 0.13191860773507982

Its sum is the evaluated energy coverage, up to the recorded flux residual.

## Frozen weights and denominators

The fourteen distinct evaluated frequencies are exactly the positive-index
DFT bins 12 through 25, paired with their negative-index conjugates. Both
signed members describe one real right-moving packet. Every required node was
evaluated for the new physical contact profile; no interpolation or missing-row
substitution was used.

The numerator includes propagating bins only. The denominator retains all 120
DFT-bin kinetic-plus-potential energy weights, including the stationary Nyquist
component. The spectrum covers 0.9996834927581062 of that full energy. Thus the
unchanged conservative omitted-weight allowance is 0.0003165072418938042.

For surface N160, the resulting **tail-only** interval is

`[0.7163064796290152, 0.7166229868709090]`.

It bounds unevaluated spectral weight conditional on the input coefficient
values. It does not bound quadrature, finite-box, timestep, cadence or finite-time
stored-energy-estimator errors. The historical reference-centered tolerance
certificate also does not automatically apply to this profile or these other
frequencies.

The periodic Parseval energy is `3.597297863368444e-5`, whereas initial open-line
energy is `3.5972978303311774e-5`. The initial candidate contact energy is
`1.0642117643921177e-22`. The periodic-to-open-total normalization factor is
`1.000000009184...`; both versions are reported in the machine-readable result.

The contract preparation reproduces the frozen q/v hash exactly. The compact
time-domain preparation computes the algebraically identical carrier through
a slightly different floating-point expression: the observed q/v differences
are at most `1.82e-18` and `2.17e-18`, and normalized spectral-weight L1 change
is `1.46e-15`. This explicitly checked roundoff difference cannot account for a
meaningful packet-efficiency discrepancy.

## Executed checks

- 70 new ramp coupled responses: five kernel families times fourteen nodes.
- 70 original-reference regression responses: all four channels agree with
  their source rows within `1.00e-15`.
- No new Green kernel, adaptive frequency search or optimizer.
- The largest four-channel closure residual is `1.49e-11`.
- All physical state, source and coupling matrices are retained in
  `INTERNAL_STATES.npz`; the packet q/v and weights are also preserved.
- Each input kernel matches its source metadata hash; the tetrahedron
  producer sources also match the recovered configuration hashes. The actual
  source archive identity is recorded in `INPUT_BINDING.json`.
- The archived physical API, reduced solver, frequency contract and preparation
  match pinned known hashes in `SOURCE_AND_PREPARATION_CHECK.json`.

These checks verify model/input binding and calculation consistency, not a new
kernel error enclosure or agreement with any particular time-domain result.

## Files and replay

- `PACKET_SPECTRUM.json`: all five predictions, all 140 responses, weights,
  partitions, regression and closure results, tail-only intervals and limits.
- `INPUT_BINDING.json`, `inputs/`: exact fourteen-node matrices and metadata
  plus the frozen preparation and model modules.
- `INTERNAL_STATES.npz`: q, source t and coupling/self-energy matrices for all
  140 responses.
- `PACKET_PREPARATION_AND_WEIGHTS.npz`: both floating-point implementations of
  the initial q/v data, frequencies and full spectral weights.
- `run_packet_spectrum.py`: the executed relative-path producer, which refuses
  to overwrite its frozen output.
- `stage_inputs.py`: original recovery/staging procedure; its old workspace
  paths are a provenance record and are not needed for replay.
- `reproduce.py`: portable replay from bundled inputs into a new output folder.

Requires Python 3, NumPy and SciPy. Example:

```bash
python3 reproduce.py --output-dir ../PacketSpectrumReplay
```

The direct simulation must be compared with this weighted quantity, with its
own estimator and error controls, rather than forced to equal 73.6%.
