# Preservation check for the 13 November source build

The reviewed v3 proposal has been applied to a fresh extraction of the 12 November source build. All 44 source substitutions matched their original text uniquely. The resulting source tree matches the reviewed v3 candidate exactly. Reversing the cumulative patch restores the complete predecessor source tree byte for byte.

- Every predecessor file remains represented. Earlier release records, the predecessor cover letter and highlights are retained under `revision/history/2026-11-12/`.
- All 254 displayed mathematical blocks were compared. Equation (135) ends in a period. Its expression and all other displayed formulas retain their original contents.
- Both documents passed three fresh draft-mode LaTeX passes through the supplied build script’s staging and compilation functions. The manuscript has 30 pages and the supplement 35. Equation (118) remains on page 17. All 252 label identifiers, numbers and page positions and the citation order match the predecessor.
- The final logs contain no unresolved references, duplicate destinations or overfull boxes. The same two underfull bibliography warnings and standard class/package notices remain.
- The abstract has 249 words under the project whitespace count. The cover-letter body has 309 words. The five highlights have 83, 83, 80, 85 and 84 characters.
- Only `word/document.xml` changes inside each submission DOCX. The original styles, relationships, properties, salutation and signature are preserved. The original author biography and author checklist are retained.
- Figure assets, their standalone exports and generating data, the historical TeX binding, numerical-method references, build script and bundled REVTeX distribution retain predecessor bytes.
- The matching validation repository contains identical current manuscript sources and current verification records. Its own manifests identify the author ZIP and protected scientific inputs.

`revision/preservation.json` records source identities and comparison results. `verification/build-inputs.json` records the actual compilation inputs and diagnostics. `revision/exact-replacements.json` and `revision/submission-changes.json` give the text changes.

The check uses draft-mode layout. No new manuscript or supplement PDF was generated. Author, deposit and companion-status confirmations remain in `submission/AUTHOR_ACTIONS.md`.
