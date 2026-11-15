# Reference taper and finite-contact validation

The reference module retains the computational evidence from the November 11 cumulative revision. The original evidence, simulation producers, figure inputs and numerical checks keep their original bytes. The editable 13 November manuscript sources are in `../../manuscript/`; generated manuscript PDFs are not included in this snapshot. The command guide is [Running the validations](../../RUNNING_VALIDATIONS.md).

`evidence/` contains the input matrices and simulation records. `tools/` contains the contact, low-frequency, angular, packet and robustness checks. `checks/` contains dated scientific outputs. For the six reference-grid symmetry/reduction residuals, see `checks/sm2_replay/`. The archived figure inputs in `context/figures/` use the historical panel layout. The current manuscript distributes some panels across the main paper and supplement.

To verify the public inventory, run `python -B tools/verify_package.py` from this directory. Instructions for regenerating the reference packet are in [REFERENCE_REPRODUCTION.md](REFERENCE_REPRODUCTION.md). The historical editorial material stays in the author’s earlier revision archive. File hashes and public-edition changes are recorded in `../../PUBLIC_SOURCE_PROVENANCE.json`.
