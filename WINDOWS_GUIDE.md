# Windows: extract, check and choose a replay

## Fix the extraction error first

The screenshot reports **0x80010135: Path too long**. Windows is objecting to a complete destination path, which includes the folders above the displayed PDF. That message does not establish that the PDF itself is damaged. The screenshot does not show the original destination, so its exact length cannot be reconstructed.

1. Download `P1_validation_2026-11-15.zip` and choose **Extract All**.
2. Set the destination to **`C:\p1`**, or another short writable location. Use a new, empty destination so files from an earlier incomplete extraction are not mixed in.
3. Open **`C:\p1\p1val`** after extraction. It should contain `verify.py`, `README.md`, `manuscript` and `validation`. That is the repository root.

This update shortens the outer folder to `p1val`. It keeps every scientific filename and internal folder unchanged. With the recommended destination, the longest packaged file path has 169 characters, below the traditional 260-character Windows limit. Copying the repository later into a deeply nested folder can bring the same problem back.

If File Explorer still reports a path error, check that it is extracting to the selected short destination. An alternative ZIP extractor, such as 7-Zip, can also be used with that destination. Do not press **Skip** and treat the resulting folder as complete. No registry change is needed for the short paths in this package.

Microsoft explains both the traditional limit and application-specific long-path behavior in [Maximum Path Length Limitation](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation).

## Run the first check

Use Python 3.12. If it is not installed, obtain it from [Python's official Windows releases](https://www.python.org/downloads/windows/). The command below requires the Windows Python launcher and an installed 3.12 interpreter; it does not require any scientific packages.

Open the extracted `p1val` folder in File Explorer. Right-click an empty area and choose **Open in Terminal**, if available. Select PowerShell, or open PowerShell and change to the actual repository folder:

```powershell
Set-Location C:\p1\p1val
py -3.12 --version
py -3.12 -B verify.py
```

For a different extraction location, change the first command. Quote a path that contains spaces. Wait for the command to finish. A successful run prints several `PASS` messages, then a line pointing to `RUNNING_VALIDATIONS.md`, and exits with code 0. Immediately after it completes, `$LASTEXITCODE` displays that exit code in PowerShell. A traceback or a nonzero exit means the check did not pass.

If `py` is unavailable but `python --version` reports Python 3.12, use `python -B verify.py` instead. Otherwise install or locate the correct interpreter before proceeding. The basic check verifies saved files and selected source/analytic consistency; it does not rerun the large simulations.

## Run the compact numerical suite

From the same repository root, follow the **Windows PowerShell** block in [Running validations](RUNNING_VALIDATIONS.md#run-the-compact-physics-suite). It creates a separate environment in `..\paper1-ref-env`, installs the pinned reference dependencies, and runs `verify.py --compact` with that environment's executable. Installation requires internet access to a package index. The scientific input files are included in the ZIP.

The explicit executable paths avoid PowerShell activation-policy changes. Do not paste the Linux `bin/python`, `source` or `export` commands into PowerShell.

The full ramp proof wrapper, fresh quadrature and compiled packet workflows are documented for Linux. For those workflows on Windows, use WSL and follow [Advanced reproduction](RUNNING_VALIDATIONS.md#advanced-reproduction). A fresh packet campaign can generate about 3 GB of data; it is not needed to check an extracted download.

## Troubleshooting

| Message or symptom | What to check |
|---|---|
| `0x80010135: Path too long` | Extract into a new short destination, such as `C:\p1`; keep the internal names unchanged. |
| Python cannot open `verify.py` | The terminal must be inside `p1val`, next to that file, rather than in the ZIP or its parent folder. |
| `py` is not recognized, or no 3.12 runtime is found | Check the installed interpreter/launcher. Use `python` only after confirming its version is 3.12. |
| `Missing packaged file` or `Changed packaged file` | Re-extract the complete download into a new folder and rerun the check. Do not remove an entry from a hash inventory to hide the failure. |
| `ModuleNotFoundError: numpy`, `scipy` or `sympy` during a numerical replay | Use the reference environment and its exact executable from the compact-suite instructions; installing into a different Python does not help that environment. |
| `ModuleNotFoundError: resource` in the full ramp launcher | That workflow uses a POSIX module. Run it in Linux/WSL, not native Windows Python. |
| An output directory already exists | Choose a new sibling output folder. Several replay wrappers intentionally reject existing output directories. |
| `--git-index` fails in a plain extracted folder | First create/clone a local Git repository and stage or commit all release files. See [GitHub guide](GITHUB_GUIDE.md). |
| GitHub's browser rejects the large proof ZIP | Use GitHub Desktop or Git for the complete repository; keep the proof ZIP intact. |

For an unresolved error, retain the full command, full terminal output, Python version and operating system. The single final error line may omit the actual failing check.

## What was checked on Windows?

This release was inspected for Windows-invalid names, case collisions and destination-path lengths. The verification commands were executed on Linux with Python 3.12 and the reference dependencies; native Windows extraction, PowerShell execution and WSL execution were not performed here. The instructions distinguish those facts from the expected Windows behavior.
