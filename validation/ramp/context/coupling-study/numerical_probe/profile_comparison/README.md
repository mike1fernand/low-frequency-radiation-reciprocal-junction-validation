# Fixed finite profile comparison after the scale probe

8 November 2026. This second batch tests exactly four predeclared shapes at
three predeclared total strengths on the same six available motifs and at the
same frequency `omega=1.2`. It uses the same actual positive/zero paired-spring
contact blocks, line parameters and bare bulk. It changes no manuscript or
accepted result.

## Main result

The increasing profile at total `sum(g)=3`, with actual per-link motif
stiffnesses

```text
g = [1,2,3,4,5,6]/7
```

gives **73.5874% incident-to-qT conversion** on the selected N192/M256
kernel. All six motifs and 48 contact springs remain active. Its complete
incident-power partition is:

| Channel | Fraction |
|---|---:|
| Line reflection | 0.1195482344 |
| Line transmission | 0.0763146256 |
| Bulk qT | 0.7358742884 |
| Bulk qL | 0.0682628515 |

Total bulk transfer is 80.4137%; qT share of that bulk output is 91.5110%.
The historical N80/M192 cache gives 73.5948% qT, differing by 0.0073268
percentage points. This is approximate point evidence, not a new accuracy
certificate.

Relative to the refined reproduction of the earlier `s=2` point (52.2121%),
this is a gain of 21.3754 percentage points or 40.9395% relative. That comparison
changes total contact strength from 2 to 3. A separate equal-total comparison
shows the profile contribution without that change:

| Profile at total `sum(g)=3` | Incident-to-qT fraction |
|---|---:|
| Earlier reference taper, scaled to total 3 | 0.6265639 |
| Uniform | 0.6777331 |
| Increasing linear | 0.7358743 |
| Decreasing linear | 0.4836785 |
| Endpoint-only | 0.3784441 |

Endpoint-only has 16 active contact springs, whereas the other listed profiles
have 48. It retains the same available span but is not a fixed-active-count
comparison. Increasing/decreasing profiles have equal total, active count and
peak stiffness; reversing their order changes radiation/reflection for a
fixed left incident wave. Their line transmission agrees within `1.67e-16`,
consistent with reciprocal line transmission; this directional radiation
difference does not imply nonreciprocity.

## Checks and exact budget

`PREDECLARATION.json` was written before the response loop. It fixes totals
2.5, 3 and 4 and shapes uniform, endpoint-only, increasing-linear and
decreasing-linear. No values were added after seeing results. This second batch
is 12 profiles times two existing caches, or 24 parameter/cache states; together
with the first batch the total is **40 states**. No new Green kernel, optimizer,
frequency sample or packet simulation was evaluated.

All q/t vectors, line/full-block systems and loading data are in
`INTERNAL_STATES.npz` (336 finite arrays; 560 across both batches). Maximum
four-channel closure in this batch is `1.24e-11`, maximum relative full-block
residual is `8.30e-17`, and maximum full-block condition number is 17.566. All
468 source-package payloads remain byte-identical. The batch computation took
about 0.062 seconds after imports/input checks.

## Qualification and stop

The primary cache is the existing N192 surface/M256 PV selected reference;
the sensitivity cache is historical N80/M192, **not** the published N80/M256
family. Their differences are sensitivity, not rigorous enclosures. Shared-model
flux closure tests algebraic consistency, not the exactness of those numerical
Green integrals. The largest qT disagreement across all twelve profile cases is
0.0341 percentage points; no statistical probability or total numerical error
bar is inferred.

73.5874% is the best of these predeclared point samples, not an optimum or a
bandwidth, tolerance, real-time-packet or practical-fabrication result. Previous
reference-point error/tolerance certificates do not transfer automatically.
The fixed diagnostic is complete and no further search is selected here.
