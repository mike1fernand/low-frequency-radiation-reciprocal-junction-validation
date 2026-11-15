# Reference taper packet replay

Run the compact check from the package root with Python and NumPy 2.x.

```sh
python tools/verify_packet_reference.py
```

The check reconstructs the Gaussian line displacement, its right-moving velocity, its open and periodic energies, and the signed DFT weights. It reselects the reporting interval using only five trajectory diagnostics. It integrates the branch-energy allocation bounds in that interval, compares the time-domain result with five independently weighted stationary grids, and checks the finite-box, cadence, time-step and interface-current diagnostics. The current check compares total assigned bulk gain against the recorded current integral from the start of the run. It also reintegrates every recorded step of the current within the reporting interval. The plotted reference trace and stationary overlays can be checked against the archived Figure 6(a) data. Full source-member names and digests are in `evidence/packet_reference/PROVENANCE.json`.

## Regenerate the complete finite-box reference trajectory

The three stages below use the original checkpoint-capable velocity-Verlet producer for the 120-site line and the 120 × 56 × 56 bulk. They retain the displacement and velocity of the complete system between stages. Run them from the package root. Choose a fresh sibling directory for the checkpoints and binary snapshots. Use a C++17 compiler with OpenMP, Python, and NumPy 2.x.

```sh
mkdir ../reference-run
g++ -O3 -std=c++17 -fopenmp evidence/packet_reference/p07_trajectory_checkpoint.cpp -o ../reference-run/replay
(
  cd ../reference-run
  ./replay 1 - checkpoint30.bin - - 120 56 56 0.04 30.0 48.0 0.2 &&
  ./replay 1 checkpoint30.bin checkpoint47_8.bin - - 120 56 56 0.04 47.8 48.0 0.2 &&
  ./replay 1 checkpoint47_8.bin final.bin trajectory.csv snapshots.bin 120 56 56 0.04 57.0 48.0 0.2
)
OPENBLAS_NUM_THREADS=1 python tools/verify_reference_snapshots.py ../reference-run/snapshots.bin ../reference-run/trajectory.csv ../reference-run/final.bin
```

The 1,425 integration steps produce 46 equally spaced snapshots, with a 32-bit snapshot handoff after double-precision evolution. The secondary analysis constructs the Bloch stiffness directly from the nine undirected bond families, diagonalizes its real symmetric matrix at every wavevector and transforms all 46 fresh snapshots. At the fourteen reporting times it treats each upper–middle degeneracy on its complete eigenspace and evaluates the extremal admissible qL allocations. This independently reproduces both qT bounds and their midpoint. The binary alone occupies about 397 MiB; the reproducible package contains the source, compact observables and replay receipt. Byte-identical file results refer to the tested build. The numerical comparison also accommodates small cross-build floating-point differences.

The production replay regenerates time evolution for the reference taper. The auxiliary box, timestep and cadence comparisons use retained numerical records, and the stationary coefficients are reweighted from the original spectral quadrature results.
