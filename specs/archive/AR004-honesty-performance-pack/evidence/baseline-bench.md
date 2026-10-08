# AR004 改动前基线（Baseline Benchmark Evidence）

- 采集时间：2026-10-08（开发开始前，HEAD = aa5c282）
- 环境：WSL hypos（Ubuntu 24.04, OpenMPI）, Release 构建（/root/build-release），OMP=1，3 次取中位数
- 脚本：`scripts/bench_ar004_baseline.sh` + `scripts/collect_ar004.sh`（原始数据 /root/ar004-baseline）

## 每 iter 耗时（256²，--max-iter 500 --tol 0.0，固定迭代数）

| case | iter_ms 中位数 | 原始 3 次 |
|------|--------------|-----------|
| jacobi np=1 | 0.0701 | 0.07154 / 0.07012 / 0.06823 |
| jacobi np=4 | 0.0277 | 0.02774 / 0.02749 / 0.0305 |
| rbgs np=1 | 0.1060 | 0.1060 / 0.1091 / 0.1037 |
| rbgs np=4 | 0.0453 | 0.04534 / 0.04891 / 0.04292 |
| cg np=1 | 0.1882 | 0.1895 / 0.1855 / 0.1882 |
| cg np=4 | 0.0565 | 0.0658 / 0.05654 / 0.05083 |

## 收敛行为（64²，np=1，main.cpp sin 问题）

| case | tol | 迭代数中位数 | final_residual 中位数 | 口径说明 |
|------|-----|------------|---------------------|---------|
| jacobi | 1e-6 | 16141 | 9.994e-7 | **递推残差**（diff 范数）——G3 缺口实证：真实残差 ≈ 4×9.99e-7 ≈ 4.0e-6 > tol |
| rbgs | 1e-6 | 8368 | 9.985e-7 | 递推残差（diff 范数） |
| cg | 1e-7 | 198 | 8.15e-8 | 递推 r（理论恒等于真实残差） |

## 用途

- T008 基准对比的 before 侧（Jacobi ±5% iter_time 门槛、CG 不回归门槛、RBGS k=1/5/10 三档）
- 收敛迭代数变化的口径换算依据（真实残差判据 ×4/×6 尺度）
