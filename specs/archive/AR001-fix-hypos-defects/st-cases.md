# [AR001] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR001 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-07 |
| 执行环境 | WSL2 Ubuntu / Intel Ultra 9 185H / OpenMPI 5.0.10 + MPICH（ASan）/ gcc 15.2 / OMP_NUM_THREADS=4 |

## 测试用例列表

### ST-001：双模式构建（Release / Debug+ASan）

**关联需求：** srs §3.1
**测试类型：** 正常路径
**优先级：** High

**前置条件：** 干净源码树；WSL 已装 cmake/OpenMPI/MPICH/GTest

**测试步骤：**
1. `cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j`
2. `cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx.mpich && cmake --build build-debug -j`

**期望结果：** 两者退出码 0，产出 `hypos`/`test_hypos`；相对修复前基线（31 条警告）**零新增警告**

**实际结果：** 双构建成功；警告 0 条（`logs/T007-release-build.log`、`logs/T007-debug-build.log`、`logs/T-review-fix-*.log`）

**状态：** PASS

---

### ST-002：Halo 缓冲尺寸安全 + ASan 零报告

**关联需求：** srs §3.2
**测试类型：** 边界条件
**优先级：** High

**前置条件：** Debug+ASan 构建

**测试步骤：**
1. 运行 ctest Debug（unit/halo_2d_mpi/halo_3d_mpi/solver_mpi/overlap_consistency）
2. 自环 U1b（hw=2）与 U3/U4（hw=1）触发全部面打包/解包

**期望结果：** 无 ASan/UBSan 报告；`packFace/unpackFace` 断言 `idx==buf.size()` 通过

**实际结果：** Debug ctest 5/5（`logs/ST-ctest-debug.log`）；修复前 Red 证据显示 `std::vector` 越界断言失败（`logs/T002-red-evidence-asan.log`），修复后消失

**状态：** PASS

---

### ST-003：2D 非均匀场幽灵层精确匹配

**关联需求：** srs §3.3
**测试类型：** 正常路径
**优先级：** High

**前置条件：** 4 ranks，2×2 拓扑

**测试步骤：**
1. 单 rank 自环方向语义（U3，u=1000i+j 逐点断言）
2. np=4 哨兵场 M1：内点 `u=1000·gI+gJ`、halo 填哨兵，交换后断言

**期望结果：** 右 halo=邻居最左内点列，左 halo=最右列，上/下同理；PROC_NULL 侧保持哨兵

**实际结果：** U3/M1 通过（`logs/ST-ctest-release.log` 中 unit 与 halo_2d_mpi 条目）

**状态：** PASS

---

### ST-004：3D 六方向幽灵层精确匹配（含 back/front）

**关联需求：** srs §3.3/§3.4
**测试类型：** 正常路径
**优先级：** High

**前置条件：** 8 ranks，强制 2×2×2 拓扑

**测试步骤：**
1. U4/U2b 单 rank 自环 3D 六方向（hw=1/2）
2. np=8 M2：全局 4×4×4 哨兵场，6 方向断言

**期望结果：** 6 方向（含 k 方向 back/front）幽灵层逐点精确匹配

**实际结果：** 全部通过（`halo_3d_mpi` ctest 条目，`logs/ST-ctest-release.log`）

**状态：** PASS

---

### ST-005：3D z 不分解保持物理边界

**关联需求：** srs §3.4
**测试类型：** 边界条件
**优先级：** Medium

**测试步骤：**
1. np=4，3D 拓扑 2×2×1；断言 `neighborBack/Front == MPI_PROC_NULL`
2. 全局 64×64×8 跑 100 次迭代

**期望结果：** PROC_NULL 断言通过；解有限；无异常

**实际结果：** `SolverMpiTest.ThreeDNoZDecompositionProcNullAndRun` 通过（solver_mpi 条目含于 `logs/ST-ctest-release.log`）

**状态：** PASS

---

### ST-006：overlap on/off 数值一致性

**关联需求：** srs §3.5
**测试类型：** 正常路径
**优先级：** High

**测试步骤：**
1. O1 自环（np=1）与 O2（np=4）分别以 `JacobiSolver(true)`/`(false)` 求解至收敛
2. 断言迭代次数相同、解逐点差 ≤1e-12

**期望结果：** 迭代数一致；解 bit 级一致（独立残差 pass 保证）

**实际结果：** O1/O2 通过（`overlap_consistency` 条目含于 `logs/ST-ctest-release.log`、`logs/ST-ctest-debug.log`）

**状态：** PASS

---

### ST-007：Profiler 区域记录（on/off 模式）

**关联需求：** srs §3.5
**测试类型：** 正常路径
**优先级：** Medium

**测试步骤：**
1. P1（off）：断言 `jacobi_iteration`/`halo_exchange`/`stencil_interior`/`stencil_boundary`/`residual_allreduce` 均 >0 且 `halo_wait.callCount==0`
2. P2（on）+ I2 e2e：`halo_wait` 出现

**期望结果：** off 不记录 halo_wait；on 记录

**实际结果：** P1/P2 通过；`logs/ST-I2.log` 含 halo_wait（grep=1）

**状态：** PASS

---

### ST-008：comm_time_ms 真实且与 Profiler 对账（±5%）

**关联需求：** srs §3.6
**测试类型：** 正常路径
**优先级：** High

**测试步骤：**
1. np=1 与 np=4 运行 `--enable-profiling`，解析 stdout Profiler 与 JSON 报告
2. 比较 `comm_time_ms×iter/1000` 与 `halo_exchange(+halo_wait)` 总秒

**期望结果：** 偏差 ≤5%；comm_time_ms>0

**实际结果：** np=1 dev=0.12%、np=4 dev=0.00%；I1 comm_time_ms=0.001972>0（`logs/ST-I5-*.log`、`logs/ST-I1.log`）

**状态：** PASS

---

### ST-009：overlap_ratio 边界语义

**关联需求：** srs §3.6
**测试类型：** 边界条件
**优先级：** Medium

**测试步骤：**
1. off 运行（I1）→ 断言 ratio==0
2. on 运行（I2）→ 断言 ratio∈[0,1]

**实际结果：** I1 ratio=0；I2 ratio=0.09373（`logs/ST-I1.log`、`logs/ST-I2.log`）

**状态：** PASS

---

### ST-010：final_residual 回填（收敛/未收敛/边界）

**关联需求：** srs §3.7
**测试类型：** 正常路径 + 异常处理
**优先级：** High

**测试步骤：**
1. A2（未 solve → 0）、A3（maxIter=0 → 0）、A4（tol=∞ 首轮收敛 → iters=1）、A5（收敛后 0<residual≤tol）
2. E2 用例（maxIter=1,tol=1e-12 → residual>tol）
3. A6/U11（基类默认路径返回 0）
4. I1 报告 `final_residual>0`

**实际结果：** 全部通过；I1 final_residual=166.4195>0（`logs/ST-ctest-*.log`、`logs/ST-I1.log`）

**状态：** PASS

---

### ST-011：ctest 全绿（含多进程用例，无非法 skip）

**关联需求：** srs §3.8
**测试类型：** 回归测试
**优先级：** High

**测试步骤：**
1. Release：`ctest --test-dir build --output-on-failure`
2. Debug+ASan：`ctest --test-dir build-debug --output-on-failure`

**期望结果：** 5/5 通过（unit、halo_2d_mpi、halo_3d_mpi、solver_mpi、overlap_consistency）

**实际结果：** Release 5/5（1.87s）、Debug 5/5（3.31s）（`logs/ST-ctest-release.log`、`logs/ST-ctest-debug.log`）

> 说明：unit 在 np=1 下对多进程用例的 SKIP 是防误跑的守卫；ctest 已将多进程用例以 np=4/8 独立注册执行，属执行而非跳过。

**状态：** PASS

---

### ST-012：串行 vs 并行解一致性

**关联需求：** srs §3.8（串并行 ≤1e-10 追溯）
**测试类型：** 正常路径
**优先级：** High

**测试步骤：**
1. np=4 M3b：全局 128×128 固定 200 次迭代；每 rank 冗余运行全域名串行参考，按全局坐标逐点比较

**期望结果：** max diff ≤1e-10（设计目标 ≤1e-12）

**实际结果：** 逐点 bit 级一致（≤1e-12 断言通过，solver_mpi 条目）

**状态：** PASS

---

### ST-013：文档无失实宣称

**关联需求：** srs §3.9
**测试类型：** 正常路径
**优先级：** Medium

**测试步骤：** 核对 README/DESIGN/PERFORMANCE 对 collective/overlap/RMA/HDF5/PAPI 的表述

**期望结果：** README 存在"路线图（未实现）"节；DESIGN §2.4 标注 P2P 回退、§6 表注"路线图/未实现"；PERFORMANCE 无未实测数值

**实际结果：** D1 核查全过（T009 会话；`grep` 结果：README 路线图=1、失实 collective 声称=0、PERFORMANCE 假数据=0）

**状态：** PASS

---

### ST-014：SCALING_REPORT 实测数据回填

**关联需求：** srs §3.9
**测试类型：** 正常路径
**优先级：** Medium

**测试步骤：** 运行 1/2/4 进程强/弱扩展实测并回填报告与基准目录

**实际结果：** 强扩展 256² 数据（0.4507/0.2687/0.2043s）、弱扩展 128²/进程数据均已回填并含曲线图（`docs/SCALING_REPORT.md`、`benchmarks/reference_results/`、`scaling_plots/`）

**状态：** PASS

---

### ST-015：非法参数错误处理

**关联需求：** design §4.3.3 / §6.3 I6（srs §3.9 集成约束）
**测试类型：** 异常处理
**优先级：** Medium

**测试步骤：** `hypos --solver foo`、`hypos --comm-mode foo`

**期望结果：** 退出码 1 + 错误日志

**实际结果：** 两者退出码均 1（ST 运行记录）

**状态：** PASS

---

### ST-016：collective 回退与 p2p 结果一致

**关联需求：** design §6.3 I7
**测试类型：** 正常路径
**优先级：** Medium

**测试步骤：** np=4 分别以 `--comm-mode p2p/collective` 运行并输出 VTK，逐字节比较

**实际结果：** 完全一致（ST 运行 `I7 PASS (identical)`；回退告警按 rank 打印）

**状态：** PASS

---

### ST-017：小规模分解退化不崩溃

**关联需求：** design §6.4 E5
**测试类型：** 边界条件
**优先级：** Medium

**测试步骤：** `SolverMpiTest.TinyDecompositionSmoke` np=4：全局 4×4（rank 2×2）50 次迭代 + 全局 2×2（rank 1×1，内点盒退化为空）10 次迭代

**期望结果：** 迭代数正确、解有限、无越界（ASan）

**实际结果：** 通过（Release 与 Debug+ASan solver_mpi 条目）

**状态：** PASS

---

### ST-018：未收敛正常退出

**关联需求：** design §6.4 E2
**测试类型：** 异常处理
**优先级：** Medium

**测试步骤：** `SolverApiTest.NonConvergenceAfterMaxIterations`（maxIter=1，tol=1e-12）

**期望结果：** 返回 1；`lastResidual() > 1e-12`；不崩溃

**实际结果：** 通过

**状态：** PASS

---

### ST-019：maxIter=0 与首轮收敛语义

**关联需求：** design §4.3.3
**测试类型：** 边界条件
**优先级：** Medium

**测试步骤：** A3（maxIter=0 → 0 次迭代、残余 0）；A4（tol=∞ → 返回 1、残余有限）

**实际结果：** 通过

**状态：** PASS

---

### ST-020：lastResidual 基类默认路径

**关联需求：** srs §3.7 第 2 条
**测试类型：** 边界条件
**优先级：** Low

**测试步骤：** A6/U11：最小 `PoissonSolver` 子类不 override `lastResidual()`

**期望结果：** 返回 0，不崩溃

**实际结果：** 通过

**状态：** PASS

---

## 执行摘要

| 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|
| 20 | 20 | 0 | 0 |

---

## ST 执行报告

| 字段 | 内容 |
|------|------|
| 执行日期 | 2026-10-07 |
| 执行结果 | PASS |
| 执行轮次 | 第 1 轮 |

### 需求覆盖矩阵

| 需求 ID | 需求描述 | 测试用例 | 结果 |
|--------|---------|---------|------|
| §3.1 | F1 双模式构建 | ST-001 | PASS |
| §3.2 | F2 内存安全/尺寸 | ST-002 | PASS |
| §3.3 | F3 面-邻居语义 | ST-003, ST-004 | PASS |
| §3.4 | F4.1 3D back/front | ST-004, ST-005 | PASS |
| §3.5 | F4.2 overlap + 区域 | ST-006, ST-007 | PASS |
| §3.6 | F4.3 指标真实化 | ST-008, ST-009 | PASS |
| §3.7 | F4.4 收敛残差回传 | ST-010, ST-020 | PASS |
| §3.8 | F5 测试体系合规 | ST-011, ST-012 | PASS |
| §3.9 | F6 平台/文档对齐 | ST-013, ST-014（+ST-001、ST-011） | PASS |

**需求覆盖率：** 9 / 9（100%）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 10 | 10 | 0 | 0 |
| 边界条件 | 6 | 6 | 0 | 0 |
| 异常处理 | 3 | 3 | 0 | 0 |
| 回归测试 | 1 | 1 | 0 | 0 |
| **合计** | **20** | **20** | **0** | **0** |

### 覆盖率（gcov，修改模块）

| 模块 | 行覆盖 | 分支覆盖 |
|------|-------|---------|
| `src/comm/p2p_exchanger.cpp` | 100% | 100% |
| `src/solver/jacobi_solver.cpp` | 100% | 100% |
| `src/comm/collective_exchanger.cpp` | 100% | 100% |
| `src/main.cpp` | 86.9% | 73.8% |

### 遗留问题

| 严重性 | 描述 | 处理方式 |
|-------|------|---------|
| 信息级 | 退化子域（nxLocal ≤ 2·hw）时边界 slab 重叠双写同一单元；两次计算读取相同数据、结果 bit 级一致，无越界/竞态（E5 已覆盖验证） | 可选加固（内点盒端点归一化），延后处理 |
| 流程 | 仓库非 git 受控，未执行提交 | 用户如启用 git 后可一次性提交 |

### 结论

> **Go** — 全部 Go 条件满足：需求覆盖率 100%（9/9）；无 Critical/Major 未修复缺陷；单元/集成测试全通过且无回归（Release/Debug+ASan ctest 各 5/5、0 警告）；NFR 已验证（数值一致性、内存安全、指标对账 ±5%）。遗留仅 1 项信息级可选优化与 1 项流程备注，不阻塞归档。
