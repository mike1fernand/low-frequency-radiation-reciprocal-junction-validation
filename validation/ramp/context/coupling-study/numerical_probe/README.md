# Bounded actual-model coupling-scale probe

8 November 2026. This is a new, separately stored diagnostic requested for the
question of performance beyond Paper One's sampled 52.2% incident-to-qT value.
It does not alter the manuscript, force law, engine or accepted data.

A subsequent, separately predeclared twelve-profile comparison is complete in
`profile_comparison/`. Its best sampled incident-to-qT point is approximately
73.5874% using an increasing profile at total strength 3. The two batches together
use 40 parameter/cache states and zero new kernels; the scale-only results below
remain the original eight-point comparison rather than a final optimum.

## Result

At the same frequency `omega=1.2`, same six-motif graph, same line/bulk parameters
and same reference taper, scaling every physical source-link stiffness by `s=2.5`
gives approximately **62.9193% incident-to-qT conversion** with the selected
N192/M256 kernel. The actual per-link motif strengths are
`[0.25,0.50,0.50,0.50,0.50,0.25]`; all 48 physical contact springs remain present.
The reflected line, transmitted line, qT and qL fractions are respectively
`0.1406877982, 0.1680478443, 0.6291934933, 0.0620708642`.
The qT fraction of total emitted bulk power is `0.9102067632`.

The same refined input gives `0.5221205514` at `s=2`, reproducing the paper's
rounded 52.2% result. The increase from that point is **10.7073 percentage
points**, or **20.5073% relative**. The benchmark at `s=1` is reproduced to
`2.78e-17` maximum absolute channel difference.

| Coupling scale s | qT incident fraction, N192/M256 | Reflected line fraction |
|---:|---:|---:|
| 1 | 0.1256593 | 0.0035721 |
| 1.5 | 0.3087931 | 0.0056407 |
| 2 | 0.5221206 | 0.0444182 |
| 2.5 | 0.6291935 | 0.1406878 |
| 3 | 0.6265639 | 0.2426165 |
| 4 | 0.5519801 | 0.3729884 |
| 6 | 0.4639022 | 0.4784921 |
| 10 | 0.4031332 | 0.5473282 |

This table shows overcoupling rather than monotone improvement. `s=2.5` is the
best of eight predeclared sampled values, not a continuous or global optimum.

## Evidence and limits

- Exactly 16 parameter/cache states were evaluated: eight declared scales on
  each of two existing bare-bulk kernels. There was no new Green quadrature,
  trajectory, optimizer or adaptive expansion.
- The primary cache is the existing N192 surface-spectrum/M256 PV selected
  reference kernel. Its input hash, selection metadata and model conventions
  match the retained Paper One preparation.
- The sensitivity cache uses historical N80 branch matrices and **M192** PV.
  It is not the published N80/M256 coupling-family producer. At `s=2.5` it gives
  qT `0.6290751489`, differing by 0.0118344 percentage points. This agreement is
  observed sensitivity, not a rigorous error enclosure or independent full-method
  verification.
- Maximum four-channel closure error is `1.29551e-11`. Maximum relative
  full-block equation residual is `7.62321e-17`; maximum full-block 2-norm
  condition number is 15.511. These test internal consistency, not the accuracy
  of the approximate Green integrals.
- Every q/t state and its line/full-block system are retained in
  `INTERNAL_STATES.npz` (224 arrays). All arrays are finite. Inputs and exact
  original modules are copied under `inputs/` with hashes in the predeclaration
  and result file. All 468 payloads of the source manuscript package were
  checked byte-identical before/after execution.
- No new frequency-band, wavepacket, robustness, tolerance-certificate,
  microscopic-correspondence, novelty or publication-significance claim follows
  from this point calculation. Existing reference-point certificates do not
  automatically transfer to `s=2.5`.

## Files and continuation

`PREDECLARATION.json` fixes the eight scales and stop condition before the
response loop. `SCALE_PROBE.json` is the full result; `INTERNAL_STATES.npz`
contains reproducible residual data. `run_scale_probe.py` is the exact executed
script; it uses the preserved source-package path named near its top and refuses
to overwrite its result. For a rerun, export a fresh copy and point `PACKAGE`
to the unchanged recovered Paper One package. Input hashes must still match.
No rerun is required merely to report these results.
