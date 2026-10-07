# Reference Benchmark Results

This directory stores reference performance data measured on the maintainer's
WSL Ubuntu environment (2026-10-07, Intel Ultra 9 185H, OpenMPI 5.0.10,
g++ 15.2, `-O3 -march=native`, `OMP_NUM_THREADS=1`), used for regression
comparison. See `docs/SCALING_REPORT.md` for the full report and methodology.

| File | Content |
|------|---------|
| `strong_scaling_wsl.json` | Strong scaling, 256×256 fixed, 1/2/4 MPI ranks, 10000 iterations |
| `weak_scaling_wsl.json` | Weak scaling, 128×128 per rank, 1/2/4 MPI ranks, 10000 iterations |
| `performance_256sq_4ranks.json` | Raw `performance_report.json` of the np=4 strong run |

Measurements are machine-specific; treat them as a baseline, not a promise.
