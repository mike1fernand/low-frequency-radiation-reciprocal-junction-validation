# Validation package documentation review — 15 November 2026

This is a documentation and packaging update to the 13 November 2026 scientific snapshot. Its scope is usability, preservation and reproducibility instructions. It does not revise the manuscript, equations, calculations, scientific claims or validation code.

## Problems addressed

| Finding | Assessment | Change |
|---|---|---|
| The screenshot reports Windows `0x80010135: Path too long`. | Confirmed error message. The full original destination is not visible, so its exact length cannot be established. | Short ZIP name, a single `p1val` outer folder, extraction to `C:\p1`, and a Windows guide. Internal scientific paths are unchanged. |
| A first-time user has to distinguish several environments and historical entry points. | Confirmed navigation problem; the scientific differences between checks are meaningful. | A basic standard-library check first, a purpose/requirements table, an index covering all 27 inherited READMEs, and links to specialist workflows. |
| The quick start depends on POSIX conventions and incomplete substitutions for Windows. | Confirmed presentation problem. Full ramp workflows also have an actual POSIX dependency. | Separate copyable PowerShell and Linux/macOS compact commands, explicit interpreter paths, and Linux/WSL guidance for advanced workflows. |
| Successful saved-file checks could be mistaken for freshly rerun simulations or proofs. | Clarity risk; those checks do not have identical scope. | Explain integrity/transcription checks, compact saved-kernel numerical replays, saved-certificate replay, fresh quadrature and fresh dynamics separately. |
| The SM2 README points to a standalone package starting directory. | Confirmed command-path error for this combined repository. | Correct root-relative paths, with the standard-library records-only check first and the numerical replay separately identified. |
| The reference README says manuscript PDFs are included. | Inaccurate for the current source snapshot. | Identify the included sources and optional build instructions. No new PDFs are generated for this update. |
| Some historical READMEs name earlier folders, pages, source builds or authority records. | These are useful provenance, rather than uniformly erroneous current scientific statements. | Preserve their bytes and explain their scope/current navigation in the index. Do not rewrite sealed proof-capsule or third-party documentation. |
| The nested proof ZIP is rejected by GitHub's browser uploader. | Confirmed size mismatch for that upload method; it is below the regular-Git blocking limit. | Keep the input intact and document Desktop/Git upload plus index verification. Explain that Desktop selection checkboxes alone are not proof that all bytes are staged. |

## Change boundaries and preservation

Nine original guidance files are edited: seven README files, `RUNNING_VALIDATIONS.md` and `EVIDENCE_BY_SECTION.md`. Three hash inventories are refreshed to cover those edits and the added guidance. Each changed original, including the old inventories, is retained byte for byte in `docs/history/2026-11-13/`.

All other original files are required to remain byte-identical at their original relative paths. This includes the entire sealed R4 capsule, its nested proof archive, code, scientific data, manuscript sources, figures, software requirements, Git attributes/ignore rules and workflow. The original root README's complete result table is retained unchanged. The shorter delivery wrapper changes no internal scientific filename.

The inherited public package already omits large regenerable raw histories, as its scientific guides explain. This update adds no omissions and does not silently replace those histories with claims that fresh dynamics were run. Full raw-output regeneration remains a separate optional workflow.

The authoritative per-file change record is [DOCUMENTATION_UPDATE.json](DOCUMENTATION_UPDATE.json). The current root hash inventory covers all packaged regular files except itself. Original inventories are historical records, not the inventory to use for the edited release.

## Verification status

The following checks were actually performed on Linux with Python 3.12.14. No physics code was edited to obtain a passing result.

| Check | Observed result |
|---|---|
| Original 13 November package's default verifier | Passed before editing. |
| Updated package's basic verifier in an isolated virtual environment created without pip or scientific packages | Passed, exit code 0. This tests the standard-library quick start. |
| Updated package's `verify.py --compact` | Passed, exit code 0, using the installed versions pinned in `requirements-reference.txt`. This was not a new package-installation test. |
| Corrected SM2 root-relative `--records-only` command | Passed, exit code 0. A new P04 numerical search was not run for this update. |
| Git index checks after local commits, at repository root and in a subfolder | Both passed. `core.autocrlf=true` was set in the test repositories; the supplied attributes preserved the required bytes. No remote was used. |
| Deliberately remove the nested proof ZIP from the test Git index while retaining its working file | Correctly failed with the missing archive named. Restaging the same archive restored a passing result. The delivered working files were not deleted. |
| Original-versus-updated file comparison | All 2,130 original paths retained. Exactly 2,118 originals byte-identical; only nine guidance files and three inventories changed. All 683 R4-capsule files remain byte-identical. Changed originals are also preserved in history. |
| Active guidance's relative Markdown links and anchors | All 103 checked links resolve. This is not a claim that every historical prose path or every external website was executed/available. |
| Windows filename, case-collision and full-path audit | Passed. At `C:\p1\p1val`, the longest packaged file path has 169 characters before its null terminator. No native Windows execution was performed. |
| Preflight ZIP CRC and fresh-extraction comparison | Passed for all 2,147 packaged files. Original file-mode attributes were preserved. The final delivery ZIP is checked again after this report and inventories are finalized. |

The 17 additions comprise five current guidance/change-record files and 12 preserved original files. There are no added scientific producers, datasets or manuscript figures. The reference-suite inventory changes only the two edited README digests; the ramp inventory is untouched. Basic verification and staged-byte checks are repeated after final documentation assembly. Full scientific numerical replays need not be repeated for changes confined to this report and its inventories.

## Limits of this review

The screenshot establishes the reported error type; it does not expose the exact failing full destination. No native Windows machine or GitHub Desktop session is available in this workspace. Static filename/path checks and Linux command execution must not be presented as observed Windows extraction or PowerShell execution.

This documentation update does not claim to re-establish every scientific result. Full contour integration, fresh packet campaigns, optional artwork rendering and manuscript PDF builds are outside its scope. Runtime-dependent numerical behavior remains subject to each scientific module's own tolerances and assumptions. No files were pushed to GitHub and no remote Actions run was performed.

## Source references for platform guidance

- [Microsoft: maximum path length](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation)
- [GitHub: large-file limits](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
- [GitHub Desktop: cloning](https://docs.github.com/en/desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop)
- [GitHub Desktop: first repository](https://docs.github.com/en/desktop/overview/creating-your-first-repository-using-github-desktop)

This is a focused usability and preservation review, not an exhaustive audit of all physics or every historical producer interface.
