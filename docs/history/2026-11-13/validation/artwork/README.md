# Manuscript artwork

Run from the repository root using the reference environment (NumPy 2.3.5, matplotlib 3.10.8 and PyMuPDF 1.26.6), with pypdf 6.21.4 for the Section IV vector labels.

```sh
../paper1-ref-env/bin/python -m pip install pypdf==6.21.4
../paper1-ref-env/bin/python -B validation/artwork/reproduce.py --output-dir ../paper1-artwork
```

The command renders the [normalized qT frequency panel](low_frequency_ratio/README.md) from all 53 selected or scanned reference-contact points. It also renders the reference power partition, packet window and selection diagnostics, radial-scale scan and continuum-share comparison. It rebuilds the Section IV panel annotations from the archived vector parent and embedded font. SHA-256 identifies all seven regenerated files against the assets in the manuscript source directory. For Figure 4 this check compares the original two-panel asset; the November 12 main paper and SM4 use separate single-panel derivatives. The [panel-export script](../../manuscript/artwork/submission/export_reference_panels.py) reproduces both derivatives from that original. Run it with `--output-dir ../paper1-artwork` from the repository root to keep outputs separate. Different plotting libraries or fonts can change output bytes. These are artwork replays using stored values. The physics checks and fresh-dynamics instructions are in [the main run guide](../../RUNNING_VALIDATIONS.md).

For the retained frequency-panel asset, the supplied frequency partition PDF, the verifier reads its vector coordinates and compares all four 49-point curves against the archived frequency scan. The original September 26 renderer for this particular PDF was not recovered in the compared releases. The command also checks the packet and selector vector paths against the trajectory and selection records. The TeX source supplies the diagrams drawn directly in LaTeX.

`PROVENANCE.json` identifies the restored plotting code, three numerical plot arrays, annotation script and its two input PDFs in their original archives. The R8 renderer contains historical routines for a superseded reference plot. The public entry point selects the three R8 renderings still used in this manuscript and uses the current reference renderer for the revised reference plot.
