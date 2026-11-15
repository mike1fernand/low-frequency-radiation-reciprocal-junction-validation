# Independent packet source and primary-analysis review

8 November 2026. Primary source/analysis review passed. Same-window timestep and box controls were still running when this first receipt was written.

I read the changed C++ force/evolution/checkpoint code against the preserved original, and the preparation, runner, selector and Fourier analyzer. The profile is the actual complete eight-link motif ramp, not an imposed transverse force. All48 contact springs, the full axial line and the free outer boundaries remain. Only transverse reflections are reduced; their fourfold bulk mass/energy weight is correct.

Two arbitrary-state one-step controls independently assemble every physical bond/contact force in the full transverse box and compare the actual engine's resumed Verlet output. Position/velocity errors are at most4.34e-19 and current errors at most1.28e-21. No new long simulation was run by this reviewer.

The five-predicate coefficient-blind selector was independently recomputed from all stored rows and selects57.4–64.2. Three snapshots(first/middle/last) were independently transformed with the closed-form bulk symbol and a total-minus-upper allocation, including full coherence matrices at top ties/triples. They reproduce qT bounds/total energy within3.34e-16. The35-row primary trapezoidal average reproduces exactly. Neither this comparison nor selection uses the stationary73.6%or packet-weighted prediction as an acceptance target.

One confirmed secondary-field issue was corrected: B[::3] stops at64.0 rather than64.2, so it is not a same-window cadence control. The corrected coarse grid includes the final endpoint. The primary qT result is unchanged. The producer preserved the SHA-bound original JSON and supplied CADENCE_FIELD_CORRECTION.json; the independent checker recomputes the same corrected value.

The finite-window global qT energy remains distinct from asymptotic modal power and a local current. Initial failed56box results are retained; they are not silently replaced by a passing result. The80box run changes a numerical boundary/control domain after a measured predicate failure, with the same fixed thresholds and an independent coefficient-blind window. Subsequent controls must use the frozen passing window and retain any predicate failures explicitly.
