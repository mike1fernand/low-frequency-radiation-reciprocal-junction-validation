# Current physical12×12 matrix ladders

equal mean of two separately scrambled Sobol sequences, seeds20260920 and20260921,2^23 samples per seed,16384 energybins on[0,12],margin1.02; table is first12×12 physicalblock.

| M | Mean PV Frobeniusnorm | Relative difference fromM256 | s00 seedspan/mean (%) |
|---:|---:|---:|---:|
|192|0.896396574|0.000512981|0.0241446|
|256|0.896608386|0|0.0600599|
|320|0.89653018|0.000154902|0.177431|

| N | Total Im norm | Branch2 Im norm | min eig(Im/pi) | Raw symmetryresidual | qT usingselectedM256 |
|---:|---:|---:|---:|---:|---:|
|128|1.04826449|0.113662137|4.12468e-08|2.12756e-09|0.1256705826|
|160|1.04805511|0.113675592|4.06329e-08|3.26015e-09|0.1256661413|
|192|1.04793293|0.113683564|4.03476e-08|5.40104e-10|0.1256592554|

Histograms have10 sampled separationblocks; first12PV equivalence checked against prior6 support. Spectralnseg10 versus6 top12 equivalence checked onlyatN192.

Observed sensitivity diagnostics; no exact-lattice errorbound.
