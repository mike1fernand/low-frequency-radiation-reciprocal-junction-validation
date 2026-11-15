# Strong validation of the approximately 73.6% Paper One candidate

8 November 2026. This repository records completed computational validation. The manuscript was not edited.

## Results and current authority

The nominal contact is g_j = j/7 for j = 1,...,6. Each motif retains eight physical springs. Frequency is 6/5 in structure time; line stiffness is 169/100. Both masses, lattice spacing, and both bulk stiffnesses are one.

- Rigorous nominal infinite-lattice interval: 0.7321660760687072... <= T_qT <= 0.7396512062343026... . Exact rational endpoints are in tolerance/VALIDATED_PROOF_CHAIN.json. The safe outward percent display is 73.2166%–73.9652%.
- Rigorous continuous tolerance box: T_qT > 69% for all six independent motif factors in [0.95,1.05]. Each motif's eight springs vary together. The exact floor is 42767052725496289/61936276900000000 = 0.6905008642115569... . This is not a certificate for 48 independently varied springs or a frequency band.
- Stationary numerical estimate: 73.58742884% on the original selected kernel; fresh surface/PV gives 73.58519712%; fresh independent tetrahedron integration gives 73.59112487%; the finest recovered independent tetrahedron kernel gives 73.57847434%.
- Broadband packet: the coefficient-blind window is 57.4–64.2. The primary 80-box result is 71.64564185%; half timestep gives 71.64540867%; the larger 96-box result is 71.62667760%. The matching open-energy-normalized stationary packet prediction is 71.63064862%. These are distinct from the monochromatic carrier coefficient.

Start with Paper1_736_Strong_Validation_2026-11-08.pdf. Current authority is this README, CURRENT_STATE.json, tolerance/VALIDATED_PROOF_CHAIN.json, the two resolvent proofs and certificates, final independent reviews, packet/FINAL_SUMMARY.json, packet/CURRENT_CADENCE_CHECKS.json, and packet_review/FINAL_PACKET_REVIEW.json.

Early conditional and provisional records remain historical. Their pending-status wording does not reopen completed checks. The old secondary common-cadence fields in frozen packet JSONs are superseded by CURRENT_CADENCE_CHECKS.json and CADENCE_FIELD_CORRECTION.json. Primary means and the selected window were unchanged by that correction.

## Compact verification and response replay

Run these commands from the extracted paper1-736-validation directory. Choose new output directories so frozen results are preserved.

```bash
python3 verify_repository.py
python3 tolerance/reproduce_response_certificates.py --output-dir ../ResponseCertificateReplay
python3 packet_spectrum/reproduce.py --output-dir ../PacketSpectrumReplay
```

The first command checks saved bytes and proof-input bindings. It does not rerun the proof. The second re-executes the new exact-rational point and continuous-box certificates. Python 3 and NumPy suffice; an inverse proposal is accepted only after an exact residual check. The third replays the fourteen-node weighted packet prediction and its 140 response/reference checks using bundled kernels. See the individual README files for dependencies and scope.

Fresh quadrature code, inputs and results are under fresh_kernel/. Installed vendor packages are not distributed. Packet code, preparation, configurations and compact metrics are under packet/; its README gives C++17/OpenMP compile, evolution and analysis commands. Use a working copy and fresh output names.

The raw packet snapshots total about 3.162 GB and are hashed in packet/RAW_DATA_MANIFEST.json. They can be regenerated in seconds to minutes from the included source and preparation. They are not duplicated in this compact core. One complete primary checkpoint at tau = 66 is included: packet/candidate_dt004_box80/checkpoint_1650.bin. It contains full line and symmetry-reduced bulk displacement and velocity, model/step metadata and current state. The recorded reflection map reconstructs the full bulk. Diagnostic sums are not used as restart states. Cross-platform floating-point byte identity is not promised.

## Complete upstream certificate: two source parts

The old reference-contact floor was not transferred to the new profile. Only its bare-kernel enclosure was reused, after exact target matching and a fresh saved-evidence replay. The new response and tolerance bounds were derived separately.

This compact core contains three exact source-binding files under r4_certificate/ and the complete replay records under baseline_audit/. It does not contain every original R4 manifest payload. For an independent replay of the entire upstream certificate, obtain BOTH files:

- Paper1_736_Proof_Source_Part01.zip
- Paper1_736_Proof_Source_Part02.zip

Each is about 36 MB. Extract both into the same folder; their baseline_certificate subdirectories merge. Run baseline_certificate/assemble_baseline_certificate.py. It reconstructs the exact unchanged 71,926,283-byte PAPER1_COMPLETE_DOMAIN_VALIDATION_R4.zip and verifies SHA-256:

8068c897d9780a81670f0fe3d98b9a70d7ed04b828e3fff3b269de9775f257ce

The 463 MB R8 archive is unnecessary. The reconstructed R4 archive already bundles its inherited R3 archive. Extract its contents directly into this repository's r4_certificate directory. The three existing files are identical and may be replaced. Retain every original manifest payload. Then run:

```bash
python3 r4_certificate/analysis/full_domain_r4/export/verify_saved.py --output ../FullBaselineReplay --expected-ledger r4_certificate/analysis/full_domain_r4/budget/FINAL_CAMPAIGN_LEDGER.json --physical-check --bridge
```

The proof runtime used was Python 3.12.14, NumPy 2.3.5, python-flint 0.8.0 / FLINT 3.3.1, SymPy 1.14.0, and mpmath 1.3.0. Dependency details are in baseline_audit/R4_STANDALONE_REPLAY_DEPENDENCIES.json. The unchanged verifier audits 222 saved-job states, exact coverage/budget arithmetic, fresh compact physical/restart/corruption controls and the original bridge. It does not reintegrate all 449,698 leaves.

## Scope and preservation

The Paper One LaTeX masters, original assets and earlier accepted results remain unchanged. New proofs, analyses and parameterized simulation copies are isolated here. The initial 56-box failure and the secondary cadence correction are preserved. The accepted 80-box window was selected without a stationary coefficient and frozen before the timestep and larger-box controls. All five original reporting predicates pass in both controls.

The certificates establish deterministic mathematical bounds in the stated harmonic model. They do not establish exact one-decimal rounding to 73.6%, a universal physical realization, an experimental result, a global optimum, a bandwidth guarantee, or a 95% statistical probability. Packet comparisons remain finite-time corroboration with separate numerical sensitivities, rather than a combined rigorous error bar.
