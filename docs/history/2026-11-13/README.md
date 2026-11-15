# Paper 1 validation repository

This repository accompanies Mike Fernandes, *Low-Frequency Radiation at a Reciprocal Junction between a Line and a Three-Dimensional Harmonic Lattice*. The revised [manuscript and supplement sources](manuscript/source/) for Wave Motion describe energy leaving a one-dimensional line through a finite reciprocal contact and entering the modes of a three-dimensional harmonic lattice. The calculations distinguish the fraction of incident power converted to shear-like bulk modes from the shear-like share of the power that reaches the bulk.

This November 13, 2026 edition matches `Paper1_Wave_Motion_Source_Build_2026-11-13.zip`. `MANUSCRIPT_BINDING.json` records its exact identity, and `ASSEMBLY_CHECK.json` records the assembly checks.

## Start here

From the directory containing this README, run the following with Python 3.12.

```sh
python -B verify.py
```

The default check uses the standard library to verify the release files and manuscript binding, compare the printed Table 3 and SM5D entries with the saved data, and run a rational check on the reference-frequency shells. For the compact numerical checks and the separate proof, quadrature and time-evolution commands, see [Running the validations](RUNNING_VALIDATIONS.md). For calculations restored from earlier releases, see [Evidence by manuscript section](EVIDENCE_BY_SECTION.md).

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

## Putting this package on GitHub

Extract the ZIP and place the **contents of its inner folder**, including the dotfiles, in a local Git checkout. Use a Git client for the full package. The largest proof input exceeds the browser's single-file upload limit. If replacing the previous release, use a clean directory for this package so omitted private drafts and superseded files do not carry forward. Preserve unrelated files in an existing repository.

After staging the files, run `python -B verify.py --git-index` before committing. It checks that the Git index contains the release bytes a reader will receive. The [run guide](RUNNING_VALIDATIONS.md) describes software requirements and generated output locations. Original third-party notices are included. The author can add the public release URL, tag and DOI once those exist.
