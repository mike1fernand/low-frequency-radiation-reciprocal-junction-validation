# Portable fresh-quadrature reproduction

Place `reproduce.py` beside the original `fresh_kernel/SHA256SUMS.json`. It verifies all 64 frozen files and their original manifest, creates a new output tree, copies only bundled inputs into the expected upstream layouts, and adapts exactly one input-root assignment in each staged driver. Scientific calculations, original scripts, archived results, profile, frequency, grids, seeds, and tolerances remain unchanged. The original frozen scripts retain historical paths as provenance; use this adapter for portable replay.

Use Python 3.12 (the recorded environment used 3.12.14) and the packages in `requirements-fresh-pinned.txt`, in a separate environment from the baseline proof's `mpmath==1.3.0` environment. No 463 MB R8 archive, separate Paper One source bundle, vendored dependencies, or external scientific input files are required. Installing scientific Python packages is an environment prerequisite, not part of the evidence ZIP.

From the `fresh_kernel` directory:

```sh
python -m venv /path/to/fresh-venv
/path/to/fresh-venv/bin/python -m pip install -r requirements-fresh-pinned.txt
/path/to/fresh-venv/bin/python reproduce.py --mode pilot --output /path/to/new-pilot
```

Choose a new output directory outside this frozen `fresh_kernel` source tree. Existing output directories are rejected, so previous results cannot be overwritten. The pilot computes the declared N48 surface kernel and two fixed Sobol 2^18 histograms, then nine responses. It is a repeat of the already declared check and is not a profile search.

Other already declared modes, each in its own new output tree:

```sh
python reproduce.py --mode refine --output /path/to/new-refine
python reproduce.py --mode independent-tetra --output /path/to/new-independent-tetra
python reproduce.py --mode fresh-tetra-n48 --output /path/to/new-fresh-tetra
```

`refine` computes N80/N128 surface kernels and two Sobol 2^20 histograms. `independent-tetra` re-evaluates the candidate and reference on the five bundled fine-grid kernels without performing new quadrature. `fresh-tetra-n48` computes one N48 corrected tetrahedron kernel in its two archived integration ranges. These modes are independent and do not need to run in the same tree. To inspect staging without starting a numerical calculation, add `--stage-only`.

Each output contains `RELOCATION_MANIFEST.json`, original/staged driver hashes, exact I/O substitutions, copied input hashes, the command and runtime, and `EXECUTION.log`. Driver outputs retain their own predeclarations, provenance hashes, state matrices and results. Thread counts are set to one; `PYTHONOPTIMIZE=0` preserves the historical drivers' assertion guards and bytecode writing is disabled. The adapter verifies frozen originals again after execution, including unsuccessful numerical runs. Missing scientific dependencies are reported before an output directory is created.

The checks support numerical consistency, independent-method agreement and convergence/sensitivity at the frozen candidate. They do not establish a rigorous infinite-lattice error enclosure, measured physical efficiency, exact transversality, or equivalence to another parent Hamiltonian.
