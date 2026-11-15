# Checks added in the November 11 referee review

The archived paper sources are in `../../manuscript/`. Start with the repository-root `verify.py`; [the root run guide](../../RUNNING_VALIDATIONS.md) sets up the reference environment and includes these checks in its compact suite. The two standard-library checks can also be invoked from the repository root as follows:

```sh
python -B validation/referee_checks/verify_contact_table.py
python -B validation/referee_checks/verify_reference_shell.py
```

On Windows, use `py -3.12` in place of `python` for those commands. For numerical scripts, use the reference environment's Python executable. The pure Python shell check is independent of all saved bulk matrices. The symbolic identity check uses SymPy. The second shell check and static force-moment calculation use NumPy and SciPy. These calculations check the printed SM2 and SM7A results.

`receipts/` records the completed review runs. The ramp response and saved-contour-chain receipts refer to the authenticated proof inputs now at `../ramp/`. The numerical reference tests are in `../reference/`. Dated source-preservation records identify the checked editions. The build receipt predates this outer repository and gives the SHA-256 of the archived inner LaTeX build. No private path or external data service is needed for these checks.

`verify_contact_table.py` compares the five Table 3 rows with saved surface results using the standard library. With `--replay` it also uses NumPy to solve all five contacts on each of the two archived kernels.

`verify_profile_reversal.py` uses the released N192/M256 Green and branch matrices to recompute the equal-total power budgets at frequency 1.2. With identical axial leads a mirrored profile maps to incidence from the other end. The script checks the equal line through-transmission and tests rearranging interior strengths while leaving the first two and last strengths fixed. These are stationary fluxes for one contact family.
