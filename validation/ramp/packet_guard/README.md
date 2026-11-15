# Packet audit and safe staged continuation

Read `PACKET_FINAL_AUDIT.md` for scientific scope, the diagnostic-informed 56→80 domain adaptation, current cadence authority, the correct stationary normalization, and the raw-history omission policy. The machine-readable receipt is `PACKET_FINAL_AUDIT.json`.

The archived C++ engine and scientific results remain unchanged. This folder adds a strict history check around its existing checkpoint mechanism. **A complete state checkpoint alone cannot replace preceding histories for the standard analyzer.**

## Dry-run checks and guarded continuation

From this folder, the following checks a hypothetical stage of a freshly regenerated main run. The history must currently end at checkpoint step 1500 (τ60), and the supplied monitor cadence must match that used to produce it:

```bash
python3 safe_resume.py --packet-dir ../packet \
  --run-dir ../packet/replay_staged \
  --checkpoint ../packet/replay_staged/checkpoint_1500.bin \
  --target-step 1650 --every 300
```

The default is read-only validation. Append `--execute` to compile the pinned C++17/OpenMP source and continue. For production stops at/after τ48, choose a target on the 0.2 snapshot grid. The wrapper preserves the exact prior `configuration.json`, then updates the successful final target and active binary hash while retaining all other fields, including any frozen window.

The large τ66 checkpoint bundled with the compact repository is a complete continuation state, but a dry-run against its directory will correctly reject missing raw histories. Regenerate those histories from the initial state before analysis or guarded continuation. Never append an earlier checkpoint to histories already extending beyond it.

An interrupted engine is not a transaction. Preserve the failed receipt, all partial output, and the exact configuration backup; diagnose the directory before retrying. The wrapper does not repair histories, infer missing states, or verify every interior trace value. Its lock only coordinates invocations of this wrapper, so do not run the legacy engine concurrently in the same directory.

## Cheap reproducible guard tests

The preserved H=2 fixture is included here, so these tests do not need the excluded multi-gigabyte packet recordings:

```bash
python3 test_safe_resume.py --packet-dir ../packet \
  --output GuardReplayResults.json
```

Use a new output name. The final shipped receipt is `SAFE_RESUME_VALIDATION_FINAL.json`: 15 controls pass. Only one small-box integration step is executed. Other fixtures test acceptance, corrupt/missing/truncated histories, wrong metadata, final-time cadence, or failure logging. `fixtures/FIXTURE_PROVENANCE.json` identifies the preserved executed fixture and each file hash. Synthetic snapshot-layout fixtures are created temporarily and are not claimed as physical histories.

The guard and tests use the Python standard library plus `g++`/OpenMP for the one-step execution. Full packet postprocessing additionally requires NumPy 2 or later; stationary spectrum replay requires SciPy. The tested native binary layout is little-endian float64/int32 on Linux.
