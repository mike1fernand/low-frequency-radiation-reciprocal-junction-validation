# Running the Paper 1 validations

Run the commands from the directory containing `verify.py`. Python 3.12 and little-endian Linux were used for the numerical replays. The main article and supplement are supplied in editable form under `manuscript/source/`. `MANUSCRIPT_BINDING.json` identifies the matching revised source-build ZIP and the SHA-256 digests of its scientific sources. Historical records and comparison inputs retain their original dates and scope.

See [Evidence by manuscript section](EVIDENCE_BY_SECTION.md) for the earlier source calculations, including the finite-box tests, the three-method taper comparison, the predictions and the November enhancement calculations.

## Verify the downloaded or staged files

```sh
python -B verify.py
```

This uses only the Python standard library. It checks the delivered file inventory and manuscript binding. For Table 3 it compares all five rows with the exact input stiffnesses and saved surface results. It also checks all fourteen printed SM5D ramp frequency–conversion pairs, runs the reference and ramp package verifiers, and performs an independent rational covering of the reference-frequency shells. Saved-file integrity on its own does not reproduce a numerical result.

To put the files in Git, copy the contents into an empty repository checkout or into a dedicated subdirectory of an existing checkout. Include both `.gitattributes` and `.gitignore` exactly as delivered. A subdirectory preserves any existing project files with the same names. Keep the distribution's reference and ramp validation subdirectories free of added files. Use Git or GitHub Desktop for the initial commit. Some scientific assets exceed the browser upload limit. Before committing, stage the intended release and check the bytes Git would publish.

```sh
git add .
python -B verify.py --git-index
```

The optional Git check verifies the staged blob for each listed file and for the release inventory itself. A local working copy can pass a file check even when a Git ignore rule has left scientific inputs unstaged. The same command can be run from a fresh clone after committing. Existing unrelated repository files need not appear in the release inventory. [GitHub file-size guidance](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)

With this package at the repository root on the default branch, the included `Paper 1 validation` workflow can also be started manually from GitHub Actions. Its default runs the integrity checks. Its compact option also installs the reference dependencies and runs the numerical suite. The local checkout procedure works for installations in a subdirectory.

## Run the compact physics suite

The compact run recomputes Table 3 from the two archived kernels: the surface calculation at N = 192 and principal-value order M = 256, and the N = 192 tetrahedron calculation. It checks all five reflected fractions and five qT fractions against the percentages printed in the manuscript. This is a replay of the coupled scattering solves at the carrier frequency. The exact-lattice certificates for the reference and ramp have their own commands below.

Create a reference environment in a sibling directory. The commands use POSIX environment paths. In Windows PowerShell, replace `../paper1-ref-env/bin/python` with `..\paper1-ref-env\Scripts\python.exe`, and use the corresponding executable for other environments. The compiled evolution commands below were tested on Linux and require C++17 and OpenMP. The full ramp driver also uses the POSIX `resource` module, so use a Linux environment, including WSL on Windows, for that workflow.

```sh
python -m venv ../paper1-ref-env
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

## Build the current paper

```sh
python manuscript/build.py --out ../paper1-current-pdfs
```

The build script stages the bundled REVTeX components and runs three pdfLaTeX passes per document. Required TeX packages are listed in [the build instructions](manuscript/README.md). Each successful replay has the scope of the commands actually run.
