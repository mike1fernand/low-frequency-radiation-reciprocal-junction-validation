# Bounded compact tile integrator

Use the pinned task interpreter. Write a settings JSON such as
`BODY_SAMPLE_SETTINGS.json`, then repeat:

```sh
python run_stage.py --settings BODY_SAMPLE_SETTINGS.json --output result.json --max-seconds 20
```

`integrate_stage(settings, destination, max_seconds=20)` is the corresponding
Python API. Stage budgets must be positive and at most 40 seconds. A stage
can stop inside a quadrature attempt; its unaccepted leaf is preserved and
retried. Returned state `status` distinguishes a partial checkpoint, blocked
checkpoint, and `rigorous_complete_tile_enclosure`.

Settings require `bounds` as `[u0,u1,v0,v1]`, `chart` as 0/1/2, and `mode`
as `qt` or `green`. Fractions are strings. Other settings are `height`,
three `orders` (3 or 5), `bits` (at least64), and either one positive
`density_tolerance` or 17 positive `density_tolerances` in moment order
`a0..a5,c0..c5,b1..b5`. Optional radius multipliers, endpoint candidates, and
maximum depth are included in the checkpoint identity. Unknown settings
are rejected.

Complete outputs include `moments` in the R3 serialization format and
`aggregate_endpoints`. The latter uses rows
`[real_lower,real_upper,imag_lower,imag_upper]`, with each endpoint represented
by `[integer_mantissa_string,exponent2]`. For qT, use the imaginary interval
of each row. `completed_moments(state)` audits the state and returns its
outward Arb values. It rejects pending coverage. An enclosing tile does not
by itself establish a full-domain enclosure or an exact-lattice floor.

`check_compact.py --reference <R3 completed body JSON> --output checks.json`
performs fresh whole-cell scalar/vector integrations, interruption/restart
comparison, exact endpoint arithmetic comparison with all saved R3 leaves,
and corruption/coverage rejection checks. Source files should be frozen
before starting a campaign; changed producers intentionally cannot resume
old checkpoints. Complete states resume without recomputation.

`COMPACTION_PROOF.md` states what the compact records prove and the difference
between requested truncation targets and actual enclosure widths.
