# Full validation package audit

8 November 2026. This audit examined the original compact validation and complete upstream R4 sources without editing them. Temporary altered copies were used to test the old integrity checker. The resulting findings concern preservation and replay; they do not refute the stationary conversion or tolerance results.

The full R4 source contains **683 regular files, 93,922,169 uncompressed bytes**, including the inherited 53,490,635-byte R3 ZIP. A new full validation ZIP should include that complete R4 tree intact. Merely copying its three selected proof-binding files is not a complete independent saved-job replay.

| Finding | Minimal correction |
|---|---|
| Old README repeats the R4 prefix in `--expected-ledger`. | Pass an absolute ledger path from the staged full R4 worktree. |
| Old integrity verifier accepts an omitted manifest entry, an extra file, and modified data under `python -O`. These failures were actually exercised. | Add a strict outer inventory verifier with explicit checks and negative controls. Preserve the old verifier for provenance. |
| R4 bootstrap extracts the nested inherited archive into its own package. | Replay in a fresh external copy. |
| Fresh drivers contain host-specific input roots. | Use the new hash-pinned staging adapter, reconstruct bundled inputs, and change input roots only in staged copies. |
| Packet runners create output and sometimes overwrite receipts beside the code. | Compile and run a staged packet copy outside the sealed repository. |
| R4 and fresh quadrature used distinct pinned mpmath versions. | Supply separate documented environments. |
| Report rebuild depends on installed DejaVu fonts. | Bundle fonts with license, or document this optional report-only dependency. |

No required fresh-quadrature scientific input was found missing: six surface/PV inputs, four fresh-tetrahedron source modules, and all five fine-tetrahedron kernel/configuration/metadata sets are bundled. Historical staging/audit scripts with host paths need not be advertised as portable replay entrypoints.

The multi-gigabyte raw packet recordings may be regenerated from the supplied source and exact preparation. If omitted, the final README must plainly distinguish a complete validation/proof repository from a complete raw-recording archive. Include the raw inventory, complete initial state, configurations, frozen selector/window, source/patch provenance, diagnostics and correction receipts. Cross-platform byte identity is not promised for floating-point dynamics.

The unchanged R4 verifier replays 222 saved jobs with physical checks and the bridge; it does not reintegrate 449,698 leaves. A single relocated fresh pilot tests adapter portability only. Keep those scopes separate from the exact rational point and whole-motif box certificates.

`PACKAGE_AUDIT.json` records the exact findings, required inclusions, raw-recording policy and executed negative fixtures. The underlying receipts are preserved separately in this folder.

The new replay wrapper should also force `PYTHONOPTIMIZE=0` for historical proof subprocesses and disable bytecode generation. This is a defensive consequence of their use of `assert`, not a change to the scientific calculation.

The completed fresh adapter is in `fresh_portability/delivery`. One relocated declared N48/Sobol 2^18 pilot passed: all nine responses and saved internal arrays exactly matched the archived pilot in this environment. Overwrite, tamper, missing-dependency and optimization guards also passed. This demonstrates portability of the declared pilot; it does not recertify the full scientific campaign. Include the ten selected delivery files, not the temporary relocated worktrees.
