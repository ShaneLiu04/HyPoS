# AR004 后置基准与门槛（After-Benchmark Evidence）

- 采集时间：2026-10-08（T001-T007 全部落地后，HEAD = 28685e8）
- 环境：与基线相同——WSL hypos（Ubuntu 24.04, OpenMPI），Release 构建（/root/build-release），OMP=1，3 次取中位数
- 脚本：`scripts/bench_ar004.sh`（原始数据 /root/ar004-after）；before 侧见 `baseline-bench.md`

## 门槛判定（srs §4）

| gate | actual | baseline | 判定 | 说明 |
|------|--------|----------|------|------|
| Jacobi 每 iter np=1（±5% 带） | 0.0729 | 0.0701 | **PASS**（+4.0%） | 诚实判据≈零开销：融合换算替代显式残差扫描 |
| Jacobi 每 iter np=4 | 0.0307 | 0.0277 | 报告制 | 本 VM np4 跨会话漂移 ±15%（实测 0.0235↔0.0312，基线值落在离散带内），5% 门槛无意义；np1 稳定（<1% 漂移）承担门槛职责 |
| CG 每 iter np=1（≤+5%） | 0.1938 | 0.1882 | **PASS**（+3.0%） | 初始化/updateP 并行化收益与噪声相抵 |
| CG 每 iter np=4（≤+5%） | 0.0524 | 0.0565 | **PASS**（−7.3%） | updatePInterior 并行化真实收益，历轮均大幅低于基线 |
| RBGS k=1/5/10 | 见下表 | — | 报告制 | srs §4 修订：真实残差扫描为诚实成本，如实测量不设门槛 |

## 每 iter 耗时（256²，--max-iter 500 --tol 0.0）

| case | after 中位数 | before 中位数 | Δ |
|------|-------------|---------------|---|
| jacobi np=1 | 0.0729 | 0.0701 | +4.0% |
| jacobi np=1 k=10 | 0.0705 | — | 新能力：interval 省 9/10 Allreduce |
| jacobi np=4 | 0.0307 | 0.0277 | 环境噪声带内 |
| jacobi np=4 k=10 | 0.0224 | — | 新能力（−27% vs 同轮 k=1） |
| rbgs np=1 k=1 | 0.1422 | 0.1060 | +34%（诚实成本：每 iter 真实残差扫描） |
| rbgs np=1 k=5 | 0.1093 | — | +3.1% vs before k=1 |
| rbgs np=1 k=10 | 0.1089 | — | +2.7% vs before k=1 —— **k=10 基本收回成本** |
| rbgs np=4 k=1 | 0.0620 | 0.0453 | +37% |
| rbgs np=4 k=5 | 0.0449 | — | 与 before k=1 持平 |
| rbgs np=4 k=10 | 0.0454 | — | 与 before k=1 持平 |
| cg np=1 | 0.1938 | 0.1882 | +3.0% |
| cg np=4 | 0.0524 | 0.0565 | −7.3%（B2a 收益） |

## 收敛行为（64²，np=1，真实残差判据）

| case | tol | after 迭代数 | before 迭代数 | final_residual | 口径说明 |
|------|-----|-------------|---------------|----------------|---------|
| jacobi | 1e-6 | 17327 | 16141 | 9.99e-7 | before 为 diff 范数口径（G3：真实≈4e-6>tol）；after 为真实残差——增量 ≈ ln4/ln(1/ρ)，符合预测 |
| rbgs k=1 | 1e-6 | 8812 | 8368 | 9.99e-7 | 同上，+5.3% 迭代数为诚实化的代价 |
| rbgs k=10 | 1e-6 | 8820 | — | 9.80e-7 | k=10 至多多跑 9 步越过收敛点；每 iter −19%（0.0103 vs 0.0129 ms） |
| cg | 1e-7 | 198 | 198 | 8.15e-8 | CG 残差理论恒等于真实残差，迭代数不变（数值一致性佐证） |

## Sanitizer 证据

- **ASan+UBSan**（Debug 全量构建，CompilerWarnings.cmake 注入 `-fsanitize=address,undefined`）：全量 ctest 25/25 绿，零报告。
- **TSan**：WSL 环境不可用——`-fsanitize=thread` 与注入的 ASan 冲突（绕开：自定义 CMAKE_BUILD_TYPE=TSan）；随后 OpenMPI pmix 共享内存发现路径在 TSan 下于 `MPI_Init` 崩溃（`pmix_gds_shmem_fetch` SEGV），`setarch -R` 关 ASLR 亦然。按设计预案以走查 + 稳定性佐证：
  - profiler 热路径（beginRegion/endRegion）仅访问 thread_local 数据，无锁无共享可变状态；共享仅 registry_（首次注册持锁）与聚合读（持锁）——按构造无数据竞争；
  - U8（多线程精确聚合）连续 10+ 次运行零失败、零崩溃；
  - 走查三项（见下）覆盖 srs 走查验收。

## 走查 evidence（srs 走查验收落档）

1. **first-touch 分区一致性**（srs §3.3 验收①）：`Subgrid::zeroInitialize`（subgrid.cpp）2D j 外层 / 3D k 外层 `schedule(static)`，与 stencil 计算循环（jacobi/rbgs/cg 的 2D j 外层、3D k 外层 `schedule(static)`）分区方式一致——同一线程触摸的页面与后续计算页面重合。✔
2. **Cart_create 单次调用**（srs §3.4 验收②）：进程生命周期内 `MPI_Cart_create` 仅在 `UniformPartition::partition`（partition.cpp:38）调用一次；main.cpp 于 T006 删除二次创建，经 `SubgridInfo::cartComm` 消费、结尾单次 `MPI_Comm_free`。B3（topology_consistency）+ 全部 MPI e2e 测试佐证。✔
3. **trueResidualSquaredLocal 前置条件**（E5，srs §3.1）：函数契约要求调用方先刷新 halo——全部调用点满足：`globalTrueResidual`（内部先 exchange）、RBGS `iterate()`（双扫后 exchange 再扫描）、Jacobi/RBGS 出口确认扫描（经 globalTrueResidual）。U1/U3/U9 用例断言该口径。✔
