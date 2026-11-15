# Compact checkpoint arithmetic and coverage

This controller changes storage and adaptive error allocation. The R3 files
`arbitrary_tile_sharp.py`, `certified_tensor_sharp.py`, and `source_moments.py`
are copied unchanged. Their contour, exact-radical Gauss rules, Cauchy tail,
analytic-neighborhood inverse certification, source projection, and radial
separator arguments therefore remain the numerical primitives used here.

For accepted leaf `j`, let the primitive's resulting complex enclosure for
moment `m` have endpoint rectangle
`[l_R(j,m),u_R(j,m)] + i[l_I(j,m),u_I(j,m)]`. These endpoints already include
quadrature remainders and all interval arithmetic rounding. Arb exposes each
endpoint as a finite dyadic integer times a power of two. The controller
adds lower endpoints separately from upper endpoints, using integer shifts
and additions without floating-point rounding. Consequently its four saved
aggregate endpoints are **exactly** the Minkowski sum of the accepted leaf
rectangles. They are neither sampled values nor sums of displayed decimals.

Only the 17 aggregate rectangles for each bounded stage are retained. Each
stage records its accepted-leaf count and the hash of that consecutive slice
of accepted IDs. Its exact endpoint increment, prior and resulting aggregate
hashes, and preceding stage hash form a chain. The auditor independently
re-sums all stage increments and compares them to the saved aggregate. This
retains sufficient information to check the compact summation and stage
membership; rerunning the fixed source and settings verifies the individual
quadrature work. Hashes detect inconsistency and accidental changes; they
are not a mathematical substitute for the validated quadrature producer.

The root box for each contour leg is `[0,1] x [u0,u1] x [v0,v1]`. Every split
bisects exactly one rational coordinate interval. The saved split tree
therefore reconstructs every generated leaf and its exact geometry. The
auditor requires that the accepted and pending ID sets are disjoint and
their union equals precisely the reconstructed leaves. Three roots are
required for qT, four for the full Green function including its exterior
real segment. A complete-tile claim is possible only when no pending IDs
remain. Thus accepted-leaf aggregation alone is never a complete-tile result.

For vector density targets `tau_m>0`, the unchanged scalar primitive is
applied to the holomorphic functions `f_m/tau_m` with scalar target
`1/number_of_legs`. Multiplication of the returned enclosure by `tau_m`
recovers an enclosure of the original integral. The exact rational targets
are enclosed by Arb during arithmetic, so non-dyadic target representation
only widens the returned enclosure; it is never treated as an exact binary
number. This preserves the original scalar proof and permits the limiting
moment to change across regions. Requested targets control the analytic
remainder decision. Final certification must use the actual outward
aggregate endpoints, which also include arithmetic rounding, rather than
assume the requested targets were met.

The public `moments` field is an outward Arb reconstruction of the exact
aggregate. It may widen slightly. On each completed-state audit it must
match a fresh reconstruction at the stated precision. Radius diagnostics
are computed from exact endpoint widths and explicitly remain diagnostics;
the endpoints are authoritative. A global matrix-error calculation must
also account for its own summation/assembly error and displacement from the
selected cached matrices.

Every stage checks source, Python version, python-flint version, precision,
settings, recomputed endpoint certificate, coverage, sum chain, and complete
public values before reuse. The next checkpoint is written to a temporary
file, flushed, and atomically replaced. An interruption before replacement
loses the current bounded stage but preserves its predecessor. There is no
claim of an independent numerical backend or second-platform replication.
