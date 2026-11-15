# Running the Paper 1 validations

## Before running anything

- Extract **all** files. On Windows, use a short path such as `C:\p1\p1val`; see [Windows setup](WINDOWS_GUIDE.md).
- Start in the directory containing `verify.py`. Unless a section explicitly changes directory, paths are relative to that repository root.
- Use Python 3.12 for the recorded software environment. Do not use Python optimization (`-O`, `-OO` or a nonzero `PYTHONOPTIMIZE`) for direct scientific scripts.
- Keep generated outputs and virtual environments outside the repository, normally as sibling folders. Some submodules enforce an exact inventory and will reject added files.
- The included manuscript is the frozen **13 November 2026 snapshot**. [MANUSCRIPT_BINDING.json](MANUSCRIPT_BINDING.json) records its identity. The 15 November update changes documentation and packaging only.

## Check the package first

Windows PowerShell:

```powershell
py -3.12 -B verify.py
```

Linux/macOS with Python 3.12:

```sh
python3.12 -B verify.py
```

This needs only the Python standard library. It verifies release-file hashes and manuscript binding, compares all five Table 3 rows and fourteen SM5D frequency/percentage pairs with saved records, checks the reference and ramp inventories, and runs an independent rational covering of the reference-frequency shells. It does not freshly solve every physical problem. A successful run prints `PASS` messages and exits with code 0; a traceback or nonzero exit is a failure.

For Git upload and the separate staged-byte check, use [GITHUB_GUIDE.md](GITHUB_GUIDE.md). The `--git-index` flag is for a local Git checkout **after staging**, not for a newly extracted download. A passing working-folder check alone cannot prove that Git staged every input.

## Pick the right check

| Check level | Needs | Interpretation |
|---|---|---|
| Basic command above | Standard-library Python only | Saved-byte/source integrity, printed-record comparisons and rational shell localization |
| Compact suite below | `requirements-reference.txt` | Re-solves selected responses and analytic/bound checks using archived inputs; checks saved dynamics/spectra |
| Response and saved-certificate replay | Separate proof environment | Exact response/box proof replay and saved contour-chain audit; not a repeat of every quadrature leaf |
| Fresh dynamics or quadrature | Separate documented environments; compiler for dynamics | New computational runs, with larger time/disk/memory needs |
| Paper/artwork rendering | TeX or plotting dependencies | Optional presentation/reproduction task, not needed for the first validation check |

See [evidence by section](EVIDENCE_BY_SECTION.md) for earlier finite-box, three-method taper, prediction and enhancement calculations. A small sampled check, a saved-file comparison and a certified enclosure have different scopes; the module guides state their limits.

## Run the compact physics suite

The compact run recomputes Table 3 from the two archived kernels: the surface calculation at N = 192 and principal-value order M = 256, and the N = 192 tetrahedron calculation. It checks all five reflected fractions and five qT fractions against the percentages printed in the manuscript. This is a replay of the coupled scattering solves at the carrier frequency. The exact-lattice certificates for the reference and ramp have their own commands below.

Create the reference environment **outside** the repository. Internet/package-index access is needed to install dependencies; the scientific inputs themselves are already included. These commands do not activate the environment, so PowerShell activation-policy changes are unnecessary.

### Windows PowerShell

Run from the `p1val` folder:

```powershell
py -3.12 -m venv ..\paper1-ref-env
& "..\paper1-ref-env\Scripts\python.exe" -m pip install -r requirements-reference.txt
& "..\paper1-ref-env\Scripts\python.exe" -B verify.py --compact
```

### Linux or macOS

Run from the `p1val` folder with Python 3.12:


```sh
python3.12 -m venv ../paper1-ref-env
../paper1-ref-env/bin/python -m pip install -r requirements-reference.txt
../paper1-ref-env/bin/python -B verify.py --compact
```

The launcher stops on the first failing child, enables assertion guards for every child process and selects a single BLAS thread. Its commands are printed as they run. For direct historical Python commands in this guide, keep optimization disabled. The following settings apply to subprocesses. Run these in each new shell used for direct commands below.

```sh
export PYTHONOPTIMIZE=0 PYTHONDONTWRITEBYTECODE=1 PYTHONUTF8=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

In PowerShell, use `$env:PYTHONOPTIMIZE = '0'`, `$env:PYTHONDONTWRITEBYTECODE = '1'`, `$env:PYTHONUTF8 = '1'`, `$env:OPENBLAS_NUM_THREADS = '1'` and `$env:MKL_NUM_THREADS = '1'`.

| Compact checks | Evidence exercised |
| --- | --- |
| Incident power, scattering, low-frequency coefficients and profile checks | Contact and bulk symbol, static fixtures at two cubature resolutions, archived Green and spectral matrices. The uniform asymptotic theorem is proved in SM7A. |
| Reference taper, doubled taper, increasing ramp and angular controls | Modal partitions, geometric projections, flux closure, line-mass and stiffness controls, denominators and octant sums. |
| Fixed-frequency robustness | Reference vertex bounds, exact rational ramp response and all 222 saved contour-state/settings pairs. The original integrations are preserved separately from this replay. |
| Packet and spectrum comparisons | The stored reference and ramp records, packet preparation, window selection, degeneracy allocation and spectrum weighting. Fresh evolution uses the commands below. |
| SM7E directional control | Fits the saved Cartesian widths and recomputes mixed-input powers using the 2.6 shell and real-part records. The independent real-axis integrations are a separate fresh run. |
| Equal-total ramp response | Re-solves the scaled taper, uniform contact, ramp from each end and an interior permutation at 1.2. |
| Additional SM2 and SM7A checks | Exact symbolic bulk identity and gap constant, second rational mesh, 375 numerical shell roots and static force-moment fixtures. |

## Advanced reproduction

The remaining command blocks use a **POSIX shell on Linux**. On Windows, use WSL with Python 3.12, a C++17 compiler and OpenMP for the full workflows. The full ramp launcher imports the POSIX `resource` module. macOS may need a separate OpenMP-capable compiler; the recorded compiled runs were Linux runs. Do not paste the following POSIX blocks directly into PowerShell.

For unprefixed reference-environment `python` commands below, first run:

```sh
source ../paper1-ref-env/bin/activate
```

Other environments are invoked by their explicit executable paths. Keep proof, spectrum, fresh-quadrature and packet environments separate as shown; some pin different dependency versions. Use a new output directory for every run. The compact suite does not regenerate the full packet histories or the original complete-domain integrations.

The integrity and replay-error controls can be run separately without scientific dependencies. The second command uses subprocess fixtures and writes a new receipt outside the release.

```sh
python -B validation/ramp/test_integrity_verifier.py
python -B validation/ramp/test_replay_guards.py --output ../paper1-replay-control-tests.json
```

## Reproduce the higher-conversion ramp result

Use a fresh output directory for every run. The replay wrapper checks the ramp package before and after execution and writes logs and a task receipt outside the scientific inputs. The proof environment uses NumPy, python-flint, SymPy and mpmath with the pinned versions in its requirements file.

```sh
python -m venv ../paper1-proof-env
../paper1-proof-env/bin/python -m pip install -r validation/ramp/environments/proof-requirements.txt
../paper1-proof-env/bin/python -B validation/ramp/replay.py response --output-dir ../paper1-response-run
../paper1-proof-env/bin/python -B validation/ramp/replay.py baseline --output-dir ../paper1-baseline-run
```

`response` propagates the certified bare-kernel bounds at angular frequency 1.2 through the contact response. The nominal qT conversion estimate is approximately 73.5874% of incident power. Its outward containing interval is 73.2166%–73.9652%. The continuous lower bound exceeds 69% when the six motif strengths vary independently by ±5%, with all eight springs in each motif varying together.

`baseline` replays the 222 saved contour states, campaign ledger, physical controls and exact-lattice bridge. It does not reintegrate the 449,698 contour leaves. The November 11 runs of the two commands took roughly half a minute each on the recorded machine.

SM5D tabulates the fourteen nominal ramp values from the surface_N160 family ($N=160$, $M=256$). The saved frequencies, powers and weights are in [PACKET_SPECTRUM.json](validation/ramp/packet_spectrum/PACKET_SPECTRUM.json). To reproduce the coupled responses across all five saved kernel families and the fourteen-frequency packet weighting, use the separate spectrum environment.

```sh
python -m venv ../paper1-spectrum-env
../paper1-spectrum-env/bin/python -m pip install -r validation/ramp/environments/spectrum-requirements.txt
../paper1-spectrum-env/bin/python -B validation/ramp/replay.py spectrum --output-dir ../paper1-spectrum-run
```

## Regenerate and check the ramp dynamics

Create the packet environment and start the primary and both control cases. Use `--case base` for the 1,650-step, 120×80×80 primary case alone.

```sh
python -m venv ../paper1-packet-env
../paper1-packet-env/bin/python -m pip install -r validation/ramp/environments/packet-requirements.txt
../paper1-packet-env/bin/python -B validation/ramp/replay.py packet --case all --output-dir ../paper1-ramp-packets-all
../paper1-packet-env/bin/python -B validation/ramp/replay.py resume-tests --output-dir ../paper1-resume-tests
```

The three cases run sequentially and generate approximately 3 GB of files. Allow additional space for working data and several GB of memory. The finite-packet result compares stored lower-pair energy to initial line energy. In the primary window 57.4–64.2 its midpoint is approximately 0.71645641846. The two control cases cover halving the time step and widening the transverse box. The resume tests check checkpoint and history consistency. A late checkpoint alone lacks the earlier histories required for full-window analysis.

The independent primary replay can also regenerate the 91 selector rows and 35 qT allocations, project three snapshots with an independent implementation and compare six raw histories/checkpoints.

```sh
../paper1-packet-env/bin/python -B validation/reference/tools/replay_ramp_packet.py --output-dir ../paper1-ramp-raw
python -B validation/referee_checks/verify_packet_bytes.py ramp ../paper1-ramp-raw/packet/candidate_dt004_box80
```

The first command checks numerical agreement and reports the individual byte comparisons. The second command requires exact agreement for all six files and exits unsuccessfully on a missing file or changed byte. It accepts the same raw case directory under `../paper1-ramp-packets-all/packet/`. Compiler, platform or floating-point library differences can affect raw bytes. Report numerical tolerance and byte identity separately.

## Regenerate the reference taper in three stages

The reference producer evolves a 120-site line and 120×56×56 bulk for 1,425 steps. Start in the release root. The subshell returns to that directory when the stages finish. Each stage starts only after its predecessor succeeds.

```sh
mkdir ../paper1-reference-run
g++ -O3 -std=c++17 -fopenmp validation/reference/evidence/packet_reference/p07_trajectory_checkpoint.cpp -o ../paper1-reference-run/replay
(
  cd ../paper1-reference-run &&
  ./replay 1 - checkpoint30.bin - - 120 56 56 0.04 30.0 48.0 0.2 &&
  ./replay 1 checkpoint30.bin checkpoint47_8.bin - - 120 56 56 0.04 47.8 48.0 0.2 &&
  ./replay 1 checkpoint47_8.bin final.bin trajectory.csv snapshots.bin 120 56 56 0.04 57.0 48.0 0.2
)
../paper1-packet-env/bin/python -B validation/reference/tools/verify_reference_snapshots.py ../paper1-reference-run/snapshots.bin ../paper1-reference-run/trajectory.csv ../paper1-reference-run/final.bin
python -B validation/referee_checks/verify_packet_bytes.py reference ../paper1-reference-run
```

The 46 snapshots occupy about 397 MiB. The numerical analysis rebuilds the nine-bond bulk stiffness, projects every snapshot and handles full upper–middle degeneracy spaces at the fourteen reporting times. The strict byte check independently checks the trajectory rows, snapshots and final checkpoint against their three published digests. The auxiliary box, timestep and cadence checks in the compact suite use archived records.

## Directional control and ramp explanation

From the reference environment, the first command re-solves the three equal-total profiles, the opposite incidence and the interior-permuted profile using the stored reference-frequency kernel. The two remaining commands compare 24-motif directional calculations at a different frequency and contact length.

```sh
python -B validation/referee_checks/verify_profile_reversal.py
python -B validation/predictions/upstream/control_26/replay.py records ../paper1-directional-records
python -B validation/predictions/upstream/control_26/replay.py fresh ../paper1-directional-fresh --workers 3
```

At 2.6 the archived shell and Cartesian extrapolation give roughly 33.2% and 33.5% into the bulk; the corresponding upstream fractions are 17.1% and 17.2% of incident line power. Both give about 51.5% of bulk emission upstream. These differences inform reported precision. The separate 2.9 calculation gives about 71% of incident line power upstream. [Directional methods, staged grids and fit choices](validation/predictions/upstream/control_26/README.md)

## Fresh reciprocal-space quadrature

The fresh numerical environment and proof environment pin different mpmath versions. Install them separately. Each mode requires a new external output directory.

```sh
python -m venv ../paper1-fresh-env
../paper1-fresh-env/bin/python -m pip install -r validation/ramp/environments/fresh-requirements.txt
../paper1-fresh-env/bin/python -B validation/ramp/replay.py quadrature --mode pilot --output-dir ../paper1-fresh-pilot
```

`pilot` computes an N48 surface kernel and fixed Sobol histograms. `refine` computes N80 and N128 surface cases. `fresh-tetra-n48` computes a fresh corrected tetrahedron case, and `independent-tetra` reevaluates responses using saved fine-grid kernels. These checks corroborate the stationary calculation. [Detailed quadrature instructions](validation/ramp/fresh_kernel/REPRODUCE_FRESH.md)

See [Manuscript artwork](validation/artwork/README.md) for the figure generators, software versions and comparison of the frequency-scan vector curves.

## Optional build of the archived paper

```sh
python manuscript/build.py --out ../paper1-current-pdfs
```

The build script stages the bundled REVTeX components and runs three pdfLaTeX passes per document. Required TeX packages are listed in [the build instructions](manuscript/README.md). Each successful replay has the scope of the commands actually run.
