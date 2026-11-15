# Paper 1 validation repository

This repository contains the code, saved inputs and evidence supporting Mike Fernandes, *Low-Frequency Radiation at a Reciprocal Junction between a Line and a Three-Dimensional Harmonic Lattice*. It studies radiation from a longitudinal line wave through a finite reciprocal contact into a three-dimensional harmonic lattice.

**First-time users: extract the complete package, then run the basic check below. You do not need to run the large simulations or build the paper first.**

## 1. Extract the ZIP

On Windows, use a short destination such as `C:\p1`. The extracted repository will be `C:\p1\p1val`. In that folder you should see `README.md`, `verify.py`, `RUNNING_VALIDATIONS.md`, `manuscript` and `validation`, together with `.github`, `.gitattributes` and `.gitignore`.

If extraction reports **0x80010135: Path too long**, cancel that incomplete extraction and extract this ZIP into the short location. See [Windows extraction and setup](WINDOWS_GUIDE.md). Do not skip individual files or rename internal folders: both can break the validation paths or file checks.

The ZIP's short outer folder is only a delivery wrapper. Keep its contents and relative paths together.

## 2. Check the downloaded files

Use Python 3.12, the version family used for the recorded replays. Open a terminal **inside `p1val`, next to `verify.py`**, then run the command for your system.

Windows PowerShell:

```powershell
py -3.12 --version
py -3.12 -B verify.py
```

Linux or macOS, with Python 3.12 installed:

```sh
python3.12 --version
python3.12 -B verify.py
```

This first check needs only Python's standard library. No NumPy, compiler, TeX installation or Git account is needed. Several `PASS` messages are normal; the final line points to `RUNNING_VALIDATIONS.md`. A successful process exits with code 0. If an error or traceback appears, follow [troubleshooting](WINDOWS_GUIDE.md#troubleshooting) before attempting larger runs.

The basic check verifies saved-file hashes and the included manuscript snapshot, compares printed Table 3 and SM5D entries with saved records, checks the two validation inventories and runs the rational reference-shell check. **It does not rerun all numerical calculations or certify every result by itself.**

## 3. Choose how much to reproduce

| Your goal | Next step | What it checks |
|---|---|---|
| Confirm a complete, unchanged download | Basic command above | File integrity, source bindings, saved-table transcription and the rational shell check |
| Recompute the smaller numerical checks | [Compact physics suite](RUNNING_VALIDATIONS.md#run-the-compact-physics-suite) | Coupled solves on saved kernels, selected analytic/bound checks and saved packet/spectrum comparisons |
| Re-execute the saved certificate proof chain | [Ramp response and certificate replay](RUNNING_VALIDATIONS.md#reproduce-the-higher-conversion-ramp-result) | Exact response/box calculations and saved-state auditing; separate from fresh bulk quadrature |
| Generate new trajectories or bulk integrals | [Advanced reproduction commands](RUNNING_VALIDATIONS.md#advanced-reproduction) | Fresh computations with their own dependencies, new output directories and larger resource requirements |
| Find the evidence for a paper section | [Evidence by section](EVIDENCE_BY_SECTION.md) | Which files and commands support each stated result |
| Upload the complete repository to GitHub | [GitHub upload guide](GITHUB_GUIDE.md) | Desktop/Git instructions, staged-byte verification and the large proof input |

For the compact suite, the run guide gives separate copyable PowerShell and POSIX commands. Full packet, proof-wrapper and quadrature workflows use the documented Linux environment; Windows users can use WSL. Generated output belongs **outside** this frozen repository. Never run old producers in place unless their guide explicitly says it is safe.

## Which paper version is included?

The scientific and manuscript snapshot here is dated **13 November 2026**, identified by [MANUSCRIPT_BINDING.json](MANUSCRIPT_BINDING.json) and [ASSEMBLY_CHECK.json](ASSEMBLY_CHECK.json). It matches the scientific sources in `Paper1_Wave_Motion_Source_Build_2026-11-13.zip`.

This **15 November documentation and packaging update** improves instructions and shortens the ZIP wrapper; it changes no calculation, scientific input, result, manuscript source or figure. Later editorial paper builds may have different wording. The verifiers intentionally check the included 13 November snapshot, not a substituted manuscript. [Documentation change and preservation record](DOCUMENTATION_UPDATE.json)

## Which result do you want to reproduce?

| Result | Evidence and command |
| --- | --- |
| Low-frequency radiation | The analytical coefficients and their assumptions are in the main text and SM7A. Run `validation/reference/tools/verify_lowfreq.py` and the independent shell and force-moment scripts in `validation/referee_checks/`. |
| Reference taper | At the reference frequency the taper converts approximately 12.6% of incident power to the two lower modes. The recovered surface, tetrahedron and generalized-Bloch calculations and their convergence records are listed in the section guide. |
| Doubled taper | Doubling every reference-taper spring stiffness gives about 52.2% incident-to-qT conversion at 1.2. This case is in the 85-configuration survey and in the current paper’s contact-comparison table. |
| Five contacts in Table 3 | The six-motif contacts are evaluated at 1.2. At equal total link stiffness, the scaled taper, uniform profile and ramp convert about 62.7%, 67.8% and 73.6% of incident power into qT. The uniform profile in this table uses 0.5 per link in units of the diagonal bulk stiffness. The compact check [recomputes all five rows](validation/referee_checks/verify_contact_table.py) with two archived kernels. |
| Increasing ramp | From its softer end, the six-motif ramp at dimensionless angular frequency 1.2 gives an estimated 73.5874% incident-to-qT power conversion. In the proof environment run `validation/ramp/replay.py response` and `baseline` with separate new output directories. The certified interval is 73.2166%–73.9652%, rounded outward. The six-motif ±5% guarantee is above 69% at this frequency. Each motif's eight springs share one stiffness variation. |
| Ramp frequency samples | [The archived spectrum](validation/ramp/packet_spectrum/PACKET_SPECTRUM.json) gives the nominal estimate for each of the 14 frequencies in SM5D. The highest sampled estimate is 77.6% near 1.3585. The spectrum-weighted prediction is 71.63065% of initial line energy. |
| Independent time-domain test | For the increasing ramp the primary packet deposits about 71.6456% of initial line energy in the lower bulk branches under the paper's averaging prescription. `validation/ramp/replay.py packet --case all` runs the primary, smaller-step and larger-box cases. The guide also gives the reference-taper packet commands. |
| Ramp incidence and profile order | At equal total stiffness, the scaled taper, uniform contact, soft-entry ramp and reversed ramp reflect about 24.3%, 16.7%, 12.0% and 38.7% of incident power. [Power-budget replay](validation/referee_checks/verify_profile_reversal.py) solves the four channels on the same finite-frequency kernel, including an interior permutation. |
| Two directional operating points | For 24 consecutive motifs on the stiffer line (line speed 1.5), the 2.6 control gives about 17% upstream incident power, with approximately 51.5% of bulk power traveling upstream. At 2.9 the upstream power is about 71% of incident power, carried by the upper branch. [Real-axis and Cartesian checks](validation/predictions/upstream/control_26/README.md) identify the method dependence of the control digits. |
| Two geometries at the same total coupling | `validation/significance/tradeoff_validation/TWO_PROFILE_RESULT.json` records the four power fractions for each profile. The accompanying replay uses the same bulk kernel for both profiles. |


## Reading the specialist guides

[README index](README_INDEX.md) identifies every README, its starting directory and its scope. Some deeper guides describe historical packages or sealed certificate capsules. Their dates, old equation/page numbers and source paths are provenance, not replacements for this repository's entry instructions. All needed scientific inputs for the documented public entry points are included; separately downloading the old proof-source parts is unnecessary.

The evidence distinguishes incident-to-qT conversion, emitted-bulk qT share and finite-packet stored energy. Keep those denominators distinct when comparing outputs.

## GitHub and preservation

Copy the **contents** of `p1val`, including its dotfiles, into the intended local Git checkout. Use GitHub Desktop or Git for the complete package. The 53,490,635-byte nested proof ZIP exceeds the browser uploader's limit but is below GitHub's regular-Git file limit. Keep that ZIP intact; no splitting or Git LFS is needed for this release's files. Detailed steps and the `--git-index` check are in [GITHUB_GUIDE.md](GITHUB_GUIDE.md).

The current inventory is `REPOSITORY_SHA256SUMS.json`. Original versions of edited instructions and refreshed inventories are preserved in `docs/history/2026-11-13/`. Numerical evidence, the complete R4 capsule, inherited archives, model/solver code and third-party notices retain their original bytes. Retained history is intentionally dated; the root guides control current navigation.
