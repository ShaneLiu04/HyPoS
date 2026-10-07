# Reference Benchmark Results

This directory stores reference performance data measured on the maintainer's
WSL Ubuntu environment (2026, Intel Ultra 9 185H, OpenMPI 5.0.10, g++ 15.2,
`-O3 -march=native`), used for regression comparison. See
`docs/SCALING_REPORT.md` for the full report and methodology.

| File | Content |
|------|---------|
| `strong_scaling_wsl.json` | AR001 baseline, strong scaling, 256×256 fixed, 1/2/4 ranks |
| `weak_scaling_wsl.json` | AR001 baseline, weak scaling, 128×128 per rank, 1/2/4 ranks |
| `performance_256sq_4ranks.json` | AR001 baseline, raw np=4 report |
| `strong_scaling_optimized_wsl.json` | AR002 (fusion/persistent/memcpy/SIMD), strong scaling |
| `weak_scaling_optimized_wsl.json` | AR002, weak scaling |
| `performance_256sq_4ranks_optimized.json` | AR002, raw np=4 report |

Controlled head-to-head improvement (AR002 vs AR001, 3-run median iter_time):
np=1 -37.5%, np=4 -39.0%. Measurements are machine-specific; treat them as a
baseline, not a promise.
