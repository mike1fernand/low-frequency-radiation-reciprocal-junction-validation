# Archived Paper 1 LaTeX source snapshot

This folder contains the editable manuscript and supplement frozen for the 13 November 2026 validation snapshot. The 15 November documentation update does not replace these scientific sources with a later editorial manuscript. The abstract contains 249 words under the project’s whitespace count. The source, runtime figure inputs, standalone build script, bundled REVTeX components and verification records here match their counterparts in the author build.

See [the source binding](../MANUSCRIPT_BINDING.json) for the exact source-build ZIP and SHA-256, and [the validation guide](../RUNNING_VALIDATIONS.md) for the numerical analyses. The author build also contains the artwork exports and editable submission documents. The current revision patch and preservation check are included in `revision/`. The matching [panel-export recipe](artwork/submission/export_reference_panels.py) is included here. In the author ZIP, `artwork/submission/README.md` maps each submission export to the corresponding manuscript figure.

## Local build

Building the paper is optional and is not needed to run the validation checks. From the validation repository root, the following commands keep output outside the repository. Install pdfLaTeX and the dependencies listed below first.

**Windows**

```powershell
py -3.12 manuscript/build.py --out ..\paper1-paper-pdfs
```

**macOS or Linux**

```sh
python3.12 manuscript/build.py --out ../paper1-paper-pdfs
```

With these commands, the outputs are `../paper1-paper-pdfs/paper1_main.pdf` and `../paper1-paper-pdfs/paper1_supplement.pdf`. If instead you enter this `manuscript` folder and run `build.py` without `--out`, its default output is `manuscript/build/` relative to the repository root. The script uses Python 3.8 or later and an installed copy of pdfLaTeX on the command path. A full TeX Live or MiKTeX installation provides the standard LaTeX packages and fonts. For a smaller installation, include AMS math and fonts, mathtools, bm, booktabs, graphicx, hyperref, array, microtype, xcolor, PGF/TikZ, natbib, and their dependencies. Once these are installed, both documents build offline. Each document gets three pdfLaTeX passes.

## Locate the paper inputs

Keep the checked snapshot unchanged when validating this release. To revise the paper, work in a separate manuscript build; editing these inputs will intentionally fail the source-binding check.

| File or folder | Contents |
| --- | --- |
| `source/paper1_main.tex` | Manuscript, including the abstract, tables, TikZ figures and manuscript bibliography |
| `source/paper1_supplement.tex` | Supplement and its bibliography |
| `source/paper1_SM5D.tex` | Supplemental subsection included by the supplement |
| `source/paper1_methods_references.tex` | The three numerical-method references included by the supplementary document |
| `source/paper1_source_binding.tex` | Historical source identifiers retained by both documents |
| `source/figures/` | Vector and raster figure inputs. The additional right panel of Figure 5 is generated from the 53 records in `../validation/artwork/low_frequency_ratio/`; the original frequency-partition figure remains the left panel. Figure 4 and the SM4 diagnostic use separate panel files derived from the retained original two-panel asset. Figure 6 also has a supplementary diagnostic. |

Keep these filenames and the relative layout when moving the sources. References are typeset directly from the TeX files. REVTeX writes its own `*Notes.bib` control files during compilation.

## Check a build

The final logs, raw `.fls` recorder files and local-input lists are in the selected output folder's `logs/` subfolder. The release inventory in `verification/build-inputs.json` records which editable files, REVTeX components and distribution files pdfLaTeX read. For diagnostics, add `--work tex-work` to the build command to retain all intermediate files in a new directory. Use `--out OUTPUT_DIRECTORY` to choose the output location.

The current verification record identifies the compiler, inputs and observed page counts from the draft-mode checks. A local PDF build can exhibit different pagination depending on the TeX distribution and fonts.

The build script stages the four REVTeX files used by these sources from the included APS REVTeX 4.2f distribution. `third_party/revtex-4.2f-tds.zip` preserves the distribution files, including documentation and source. The upstream notices and license text are in `third_party/README-REVTeX.txt` and `third_party/LPPL-1.3c.tex`.
