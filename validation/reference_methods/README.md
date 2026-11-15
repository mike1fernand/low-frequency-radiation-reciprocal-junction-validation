# Earlier Paper 1 numerical methods

These files reproduce the **six-motif reference taper** and numerical controls cited in the current supplement. Its stiffness profile is (0.1, 0.2, 0.2, 0.2, 0.2, 0.1). The separate rising ramp reaches 73.6% at the same frequency. Its proofs and packet simulation use the repository's ramp directory.

Run the following commands from this directory with Python 3.12 and optimization disabled. Each output path is a new, writable directory outside the release. The commands below assume a POSIX shell and show example output locations under your home directory. Use the root reference environment for NumPy 2.3.5, SciPy 1.17.0 and SymPy 1.14.0. Historical programs retain the reviewed mathematical and physical assumptions of their dates. `ORIGIN.json` identifies the source archive, exact member and SHA-256 of each recovered file.

| Check | Command |
| --- | --- |
| SM4 high-resolution surface/PV and tetrahedron sequences at N=128, 160 and 192, three generalized-Bloch extrapolations | `python -B replay_reference.py ~/p1-reference-taper` |
| SM1 free-boundary and contact null spaces, independently assembled stiffnesses and perturbation controls | `python -B replay_fixtures.py rigidity ~/p1-rigidity` |
| SM3 finite-core scattering, matching-plane independence, power, degenerate-mode rotation, weak-coupling witness and nodal control | `python -B replay_fixtures.py scattering ~/p1-scattering` |
| SM4 continuum P and S weights for the finite source and axial point-force limit | `python -B replay_fixtures.py continuum ~/p1-continuum` |
| SM5C separate monolithic versus checkpointed C++ executions of the bounded fixture, comparing diagnostic bytes and snapshots (g++ with C++17/OpenMP) | `python -B replay_fixtures.py checkpoint ~/p1-checkpoint` |
| SM6 finite three-segment nonaffine coordinate fixture, including deliberately incorrect scaling controls | `python -B replay_fixtures.py clock ~/p1-clock` |
| Archived reference-shell radial derivative and gap checks with a hash-bound production symbol (also install python-flint 0.8.0) | `python -B replay_fixtures.py radial ~/p1-radial` |

`replay_reference.py` reconstructs each surface result by recomputing the Jackson order-256 principal-value matrices from the two archived Sobol spectral histograms and solving the loaded scattering problem with the corresponding N-grid shell matrices. It re-solves the tetrahedron scattering problem from each archived N-grid matrix and repeats all three Bloch extrapolations from their archived matrices and retained real-axis residue. The N=192 surface and tetrahedron results give qT fractions 0.1256592554 and 0.1256079633. The low-pair Bloch extrapolation gives 0.1256063100. The original generator and original regulator convention of those legacy Bloch matrices remain unidentified. Numerical agreement measures sensitivity across the three constructions.

`dynamics/PRODUCTION_CHECKPOINT_EQUIVALENCE.json` records the full reference run. The quick checkpoint command reruns the 80×12×12 bounded comparison. Restart identity compares executions of the original and checkpoint-capable programs. The claim about selected scalar diagnostics across separately implemented dynamics solvers requires its own evidence. The finite clock fixture tests a coordinate transformation on a 40-degree-of-freedom model. Its historical source hashes apply to the earlier version named in the fixture. Current manuscript identity is recorded at repository root.

## Reconstruct convergence tables and the 85 robustness cases

The producers write beside their inputs. **Make a working copy of the entire `analysis/` subtree outside the repository** to preserve its directory structure and the source evidence. Install scikit-image 0.26.0 and lazy-loader 0.6 for marching surfaces and the historical robustness producer.

```sh
mkdir ~/p1-replay-methods
cp -R analysis ~/p1-replay-methods/
cd ~/p1-replay-methods
python -B analysis/resolution_execution/surface/make_sm_ladders.py
python -B analysis/resolution_execution/robustness/postprocess.py --pv-contract analysis/resolution_execution/robustness/PV_CONTRACT.json
```

The ladder script calculates the PV-order and shell-resolution tables in SM4. The robustness producer verifies 54 input job records and their matrices, then reconstructs the 85 active configurations and the zero-coupling control in SM5A. Their qT values are fractions of incident power. The large shared histograms and numerical integrations are preserved inputs in this replay. `worker.py` and `run_checkpoints.py` in the robustness directory implement the original individual jobs and batch driver. A batch with the complete saved inputs will reuse accepted checkpoints.

## Produce new integrals

The commands below use the same external working copy. The first two compute a new surface and begin a fresh tetrahedral integration. The next invocation continues that tetrahedral campaign from the validated earlier cube ranges.

```sh
python -B analysis/resolution_execution/surface/surface_worker.py --omega 1.2 --grid 32
python -B analysis/resolution_execution/tetra/campaign.py --directory tetra12 --grid 12 --omega 1.2 --range-cubes 300 --max-new-parts 2
python -B analysis/resolution_execution/tetra/campaign.py --directory tetra12 --grid 12 --omega 1.2 --range-cubes 300
```

For the controlling surface run, use `--grid 192 --nseg 10`. The worker's historical fixed-PV diagnostic uses order 192. To reconstruct the paper's adopted result, combine the resulting shell matrices with the two-seed order-256 PV input. The `make_sm_ladders.py` and `replay_reference.py` examples use the six-segment shell array for this calculation. For the ten-segment result, use its leading 12×12 blocks and verify agreement with the six-segment arrays first. To regenerate the controlling tetrahedral run set `--grid 192 --omega 1.2 --range-cubes 65536`, use a new `--directory` and retain the same numerical settings for subsequent invocations. Each run saves its config hash, scope, completed ranges and result.

To regenerate the surface PV histograms, run `histogram_worker.py` inside the copied `analysis/resolution_execution/surface/` hierarchy with `--seed 20260920 --max-power 23 --bins 16384 --nseg 10`, then repeat for seed 20260921. The original reference used 8,388,608 samples per seed. Using a smaller grid or fewer samples changes the numerical problem. See the included scientific records for the separate roles of solver comparison, finite-resolution diagnostics and certified bounds.
