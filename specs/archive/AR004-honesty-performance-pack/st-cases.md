# [AR004] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR004 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-08 |

**执行环境：** WSL（OpenMPI），Release 构建 /root/build-release、Debug（ASan+UBSan）构建 /root/build-debug；基准 OMP=1、测试 OMP=4（srs §5 口径）。e2e 二进制 /root/build-release/hypos，输出目录 /root/st-ar004。

**执行方式标注：** E2E=本轮实际运行 hypos 二进制；AUTO=引用既有自动化测试（ctest/单测，本轮全量重跑通过）；WALK=代码/文档走查（evidence 落档）；BENCH=引用 T008 基准 evidence（≥3 次中位数）。

## 测试用例列表

### ST-001：Jacobi 真实残差收敛口径

**关联需求：** srs §3.1 R1 验收①
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 制造解（sin）问题，`--tol 1e-6`

**测试步骤：**
1. E2E 运行 `hypos --nx 64 --ny 64 --tol 1e-6 --max-iter 50000`
2. 读取 performance_report.json 的 iterations 与 final_residual

**期望结果：**
- 退出码 0；收敛（iterations < maxIter）；`final_residual ≤ 1e-6`（真实残差口径，非更新量范数）

**实际结果：** rc=0；iterations=17327（<50000）；final_residual=9.9874e-07 ≤ 1e-6。迭代数与 T008 bench evidence（17327）逐位一致。

**状态：** PASS

---

### ST-002：检查频率 k=10 的 Allreduce 节流

**关联需求：** srs §3.1 R1 验收②
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given `--residual-check-interval 10`，`--enable-profiling`

**测试步骤：**
1. AUTO：ctest `ar004_cli_verify`（interval10 模式，32² Jacobi）断言 `residual_allreduce.calls == ⌊iterations/10⌋` 且 `residual_confirm.calls == 1` 且 JSON 含 `residual_check_interval: 10`
2. E2E 交叉验证：`hypos --nx 64 --ny 64 --solver red_black_gs --tol 1e-6 --max-iter 50000 --residual-check-interval 10 --enable-profiling`，从 stdout profiler JSON 提取计数

**期望结果：**
- 调用数/迭代 ≈ 0.1（精确 ⌊N/10⌋）；确认扫描独立计数 == 1

**实际结果：** ctest 通过（Release 全量 25/25 内）；e2e rbgs：iterations=8820，residual_allreduce.calls=882==⌊8820/10⌋，residual_confirm.calls=1。

**状态：** PASS

---

### ST-003：默认 k=1 每步判定

**关联需求：** srs §3.1 R1 验收③
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 不传 `--residual-check-interval`

**测试步骤：**
1. AUTO：ctest `ar004_cli_verify`（default 模式）断言 JSON 默认 `residual_check_interval: 1` + CSV 17 列布局
2. AUTO：单测 I1s/R-I1s（solver 级精确计数：`residual_allreduce.calls == iterations`）

**期望结果：**
- 每步判定（语义等价改造前频率）；确认扫描独立计数

**实际结果：** ctest + 单测通过（Debug/Release 全量 25/25 内）。

**状态：** PASS

---

### ST-004：Jacobi overlap on/off 位级一致

**关联需求：** srs §3.1 R1 验收④、srs §4 NFR
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 同一问题 overlap 开关两档

**测试步骤：**
1. AUTO：ctest `overlap_consistency`（e2e 双跑逐位比较 u 场）
2. AUTO：单测 FusedResidualBitIdenticalAcrossOverlapModes

**期望结果：**
- u 场与残差逐位一致（真实残差语义下约束保持）

**实际结果：** ctest 通过（Debug/Release 全量 25/25 内）；期望迭代数已按真实残差口径更新（T001 记录）。

**状态：** PASS

---

### ST-005：非法 interval 值拒绝

**关联需求：** srs §3.1 R1 异常处理（前半）
**测试类型：** 异常处理
**优先级：** High

**前置条件：**
- Given `--residual-check-interval` 为 0 / -1 / abc / 缺失值

**测试步骤：**
1. AUTO：ctest `ar004_interval_zero/neg/nan/missing`（4 条 WILL_FAIL 条目：二进制须以退出码 1 失败，条目才 PASS）

**期望结果：**
- 退出码 1 + 可读错误（与既有 CLI 风格一致）

**实际结果：** 4/4 通过（Debug/Release 全量 25/25 内）。

**状态：** PASS

---

### ST-006：CG 传入 interval 告警忽略

**关联需求：** srs §3.1 R1 异常处理（后半）
**测试类型：** 异常处理
**优先级：** High

**前置条件：**
- Given `--solver cg --residual-check-interval 5`

**测试步骤：**
1. AUTO：ctest `ar004_cli_verify`（cg_warn 模式）断言退出码 0 + WARN「ignored」+ 迭代数与 k=1 一致
2. AUTO：main.cpp 静态确认不调用 setter（review S3 证据）

**期望结果：**
- WARN 提示 + 行为不回归（对齐 --overlap-comm 先例）

**实际结果：** ctest 通过（Release 全量 25/25 内）。

**状态：** PASS

---

### ST-007：CG 性能回归门槛

**关联需求：** srs §3.2 R2 验收①
**测试类型：** 正常路径（NFR）
**优先级：** High

**前置条件：**
- Given 256² np=1/4 基准，≥3 次中位数，OMP=1

**测试步骤：**
1. BENCH：scripts/bench_ar004.sh CG 档（T008 执行，evidence/ar004-bench.md）

**期望结果：**
- CG iter_time_ms 不劣于现状（单侧门槛：变快非回归）

**实际结果：** np1 +3.0%、np4 −7.3%，双档 PASS（+3.0% 在 ≤+5% 门槛内；−7.3% 为改进）；bench 脚本门槛判定 exit 0。

**状态：** PASS

---

### ST-008：CG 并行化数值不变

**关联需求：** srs §3.2 R2 验收②
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 同一问题并行化前后收敛行为

**测试步骤：**
1. AUTO：单测 U6 golden（tier=1 位级 hash 0x654fa5ddb3bd66b5 + 126 迭代；tier=4 迭代数+偏差<1e-11；Debug+ASan 绑定）
2. BENCH：cg64 收敛迭代数前后对比（198→198）

**期望结果：**
- 收敛迭代数与并行化前一致（逐元素运算不改数值）

**实际结果：** U6 通过（Debug 全量内）；bench cg64 198→198 不变。libgomp 归约到达序导致 OMP=4 跨运行低序位差异已按 tier 拆分记录（T004，design 前提修正）。

**状态：** PASS

---

### ST-009：全量 ctest 全绿

**关联需求：** srs §3.2 R2 验收③
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 双模式构建

**测试步骤：**
1. E2E：本轮 ST 就绪检查实际执行——`cd /root/build-debug && ctest`（串行）与 `cd /root/build-release && ctest`

**期望结果：**
- 全部测试通过，无回归

**实际结果：** Debug（ASan+UBSan）25/25 @ 91.67s；Release 25/25 @ 24.52s（2026-10-08 本轮执行）。

**状态：** PASS

---

### ST-010：first-touch 分区一致性走查

**关联需求：** srs §3.3 R3 验收①
**测试类型：** 正常路径（走查）
**优先级：** Medium

**前置条件：**
- Given zeroInitialize 并行实现

**测试步骤：**
1. WALK：evidence ①——zeroInitialize 循环分区（2D j 外层/3D k 外层，schedule(static)）与 stencil 计算循环分区对照
2. AUTO：单测 U7 全缓冲逐元素清零断言（2D halo=2 构造守护 z 平面覆盖）

**期望结果：**
- 初始化与计算同线程分区（AGENT_SPEC.md:54 first-touch 约定落地）

**实际结果：** 走查结论落 evidence/ar004-bench.md；U7 通过（Debug 全量内）。2D 分支含内层 k 循环（设计勘误已记录：2D 缓冲实含 ≥3 个 k 平面）。

**状态：** PASS

---

### ST-011：初始化语义不变（ctest）

**关联需求：** srs §3.3 R3 验收②
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given zeroInitialize 重构后全量测试

**测试步骤：**
1. E2E：同 ST-009 全量 ctest + B4 16³ 冒烟（单测内，901 迭代收敛）

**期望结果：**
- 全绿（初始化语义不变）

**实际结果：** Debug/Release 25/25；B4 通过。

**状态：** PASS

---

### ST-012：first-touch 文档宣称对齐

**关联需求：** srs §3.3 R3 验收③
**测试类型：** 正常路径（文档走查）
**优先级：** Medium

**前置条件：**
- Given README 与 AGENT_SPEC 的 first-touch 宣称

**测试步骤：**
1. WALK：README first-touch 段与实现对齐；OPTIMIZATION_GUIDE G1 勾选

**期望结果：**
- 无失实宣称；G1 勾选与完成状态一致

**实际结果：** T008 文档同步完成；review S1 R6 追溯确认（README/GUIDE 勾选含实测回填）。

**状态：** PASS

---

### ST-013：np=4 拓扑一致性

**关联需求：** srs §3.4 R4 验收①
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 单次 Cart_create 的拓扑，np=4

**测试步骤：**
1. AUTO：ctest `topology_consistency`（cartComm 非空/按行列分组尺寸和/offset 链连续/Cart_shift 四向互指/所有权释放）

**期望结果：**
- offsets/邻居/coords 与单次创建拓扑一致

**实际结果：** ctest 通过（Debug/Release 全量 25/25 内）。注：设计勘误——全局 Σ nxLocal==nx 在 2D tiling 下不变式错误，修正为按行/列分组求和（T006 记录）。

**状态：** PASS

---

### ST-014：Cart_create 单次调用走查

**关联需求：** srs §3.4 R4 验收②
**测试类型：** 正常路径（走查）
**优先级：** Medium

**前置条件：**
- Given 进程生命周期内 Cart_create 调用点

**测试步骤：**
1. WALK：evidence ②——partition.cpp 单次创建 + main.cpp 经 info.cartComm 消费 + 末尾单次 MPI_Comm_free

**期望结果：**
- `MPI_Cart_create` 仅调用一次；reorder=1 隐患消除（coords 取自 cart rank）

**实际结果：** 走查结论落 evidence；review D2 逐点核验（partition.cpp:38-49 cart rank 查 coords + main.cpp:398 单次 free）。

**状态：** PASS

---

### ST-015：拓扑重构无回归（ctest）

**关联需求：** srs §3.4 R4 验收③
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 拓扑单一来源重构后全量测试

**测试步骤：**
1. E2E：同 ST-009 全量 ctest（含 MPI 多进程条目）

**期望结果：**
- 全绿

**实际结果：** Debug/Release 25/25（含 mpiio_e2e、parse_error、mpibin 等 MPI 条目）。

**状态：** PASS

---

### ST-016：Profiler 并行区无数据竞争

**关联需求：** srs §3.5 R5 验收①
**测试类型：** 正常路径（NFR）
**优先级：** High

**前置条件：**
- Given OpenMP 并行区内调用 HYPOS_PROFILE

**测试步骤：**
1. AUTO：单测 U8（每线程独立区域名 + 全程锁定相位，10/10 稳定性）
2. AUTO：ASan+UBSan 全量（Debug 构建）零报告
3. WALK：profiler 热路径零锁走查（thread_local ThreadData；stats/report 持锁聚合）——TSan 实证不可用（WSL ASLR + OpenMPI pmix 于 MPI_Init SEGV，setarch -R 无效），已在 evidence 注记

**期望结果：**
- 无数据竞争（TSan 或走查口径记录 evidence）

**实际结果：** U8 通过；Debug 全量 25/25 ASan/UBSan 零报告；走查佐证落 evidence（srs §4「TSan 走查」口径满足）。

**状态：** PASS

---

### ST-017：剖面区名正确

**关联需求：** srs §3.5 R5 验收②
**测试类型：** 正常路径
**优先级：** Medium

**前置条件：**
- Given RBGS/CG 求解的 profile 报告

**测试步骤：**
1. AUTO：单测 U11（区名契约：rbgs_iteration/cg_iteration 与 jacobi_iteration 互斥）+ ProfilerAR004Test 全套
2. E2E：本轮 rbgs10 运行 stdout 的 profiler JSON 区名检查

**期望结果：**
- RBGS → `rbgs_iteration`，CG → `cg_iteration`，Jacobi → `jacobi_iteration`；G5 勾选

**实际结果：** 单测通过；e2e stdout 含 residual_allreduce/residual_confirm/stencil_interior 等正确区名；G5 已勾选。

**状态：** PASS

---

### ST-018：profile 报告字段兼容

**关联需求：** srs §3.5 R5 验收③
**测试类型：** 回归测试
**优先级：** Medium

**前置条件：**
- Given 既有报告消费方（performance_report.json、run_scaling_tests.py）

**测试步骤：**
1. AUTO：ctest `ar004_cli_verify`（default 模式）断言 JSON/CSV 输出结构（JSON additive 新字段、CSV 17 列）
2. AUTO：reporter 序列化单测

**期望结果：**
- 字段格式兼容（聚合不改变总时长语义；新字段 additive）

**实际结果：** ctest + 单测通过；review S3/C1 确认 JSON additive、CSV 无按位消费者。

**状态：** PASS

---

### ST-019：文档零失实宣称

**关联需求：** srs §3.6 R6 验收①
**测试类型：** 正常路径（文档走查）
**优先级：** Medium

**前置条件：**
- Given README/AGENT_SPEC/GUIDE 文档集

**测试步骤：**
1. WALK：GUIDE G1/G3/G5/G6、P1/P4/P6 勾选与 AR004 完成状态一致性；README CLI 帮助（`--help` 含 interval 行，ctest help 模式守护）

**期望结果：**
- 无失实宣称；勾选一致

**实际结果：** T008 文档同步 + review S1 R6 追溯全 YES；ctest help 模式通过。

**状态：** PASS

---

### ST-020：PERFORMANCE 基准数据落档

**关联需求：** srs §3.6 R6 验收②
**测试类型：** 正常路径（文档走查）
**优先级：** Medium

**前置条件：**
- Given PERFORMANCE §11

**测试步骤：**
1. WALK：检查 §11 含真实残差开销（RBGS k=1/5/10 三档）与 CG 并行化（前后对比）两组数据及口径说明

**期望结果：**
- 两组数据齐备，≥3 次中位数口径注明

**实际结果：** PERFORMANCE §11 落档（T008）：RBGS k1 +34%/k10 +2.7%；CG np1 +3.0%/np4 −7.3%；口径与门槛判定说明齐备。

**状态：** PASS

---

### ST-021：三求解器制造解收敛（B1）

**关联需求：** srs §3.1/§3.2 验收（design §6.3 B1）
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 64² 制造解，np=1，tol 1e-6

**测试步骤：**
1. E2E：本轮三求解器各一跑（jacobi/red_black_gs/cg），读 JSON
2. AUTO：np=4 档引用既有 ConvergedSolutionMatchesManufacturedSine（test_solver_mpi，ctest 内）

**期望结果：**
- 全部收敛；final_residual ≤ tol（真实口径）；CG 迭代数 < Jacobi（谱半径定性关系）

**实际结果：** E2E：jacobi 17327 iters/9.9874e-07、rbgs 8812/9.9886e-07、cg 183/8.8812e-07（183 < 17327 ✓）；np=4 既有测试通过。

**状态：** PASS

---

### ST-022：3D 冒烟（B4）

**关联需求：** srs §3.1 R1（design §6.3 B4）
**测试类型：** 边界条件
**优先级：** Medium

**前置条件：**
- Given 3D 16³ 问题（D=6 换算路径）

**测试步骤：**
1. AUTO：单测 B4（16³ np=1：901 迭代收敛 + 真实残差 ≤ tol + 正值检查）

**期望结果：**
- 收敛 + final_residual ≤ tol（3D 残差核正确）

**实际结果：** 通过（Debug/Release 全量内）。

**状态：** PASS

---

### ST-023：RBGS interval 三档收敛（B5）

**关联需求：** srs §3.1 R1 期望行为（design §6.3 B5）
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given RBGS k=1/5/10 三档

**测试步骤：**
1. E2E：本轮 k=10 档（8820 迭代 vs k=1 档 8812，差 8 < k=10）
2. AUTO：R-I1s/R-I5s/R-E4s 计数精确断言（solver 级 k=1/5 档）+ U4（k=1 收敛口径）
3. BENCH：evidence 64² 收敛对比（8812 vs 8820，差 8<10）

**期望结果：**
- 三档均收敛；k>1 允许多迭代 ≤k 步（差值 < k）

**实际结果：** 全部通过。注：B5 未作独立自动化用例，语义由上述分片覆盖（review D4 已补记替代口径，固化用例列后续 AR 候选）。

**状态：** PASS

---

### ST-024：CG breakdown 守卫（E2）

**关联需求：** srs §4 NFR 兼容性（design §6.4 E2）
**测试类型：** 异常处理
**优先级：** Medium

**前置条件：**
- Given CG 遇 pap ≤ 0（NaN 传播）

**测试步骤：**
1. AUTO：单测 BreakdownGuardFiresBeforeStateMutation（相邻双 inf rhs 确定性触发 pap=NaN，断言守卫先于 u 状态污染）

**期望结果：**
- WARN + 优雅退出（重构后既有行为保持）

**实际结果：** 通过（Debug/Release 全量内）。

**状态：** PASS

---

### ST-025：单 rank 退化子域（E3）

**关联需求：** srs §4 NFR（design §6.4 E3）
**测试类型：** 边界条件
**优先级：** Medium

**前置条件：**
- Given nLocal ≤ 2hw 退化分区

**测试步骤：**
1. AUTO：既有 test_solver 边界用例（TinyDecompositionSmoke 等）在全量 ctest 内

**期望结果：**
- 残差换算/扫描在退化分区下不越界

**实际结果：** 通过（Debug/Release 全量 25/25）。

**状态：** PASS

---

### ST-026：k ≥ maxIter 边界（E4）

**关联需求：** srs §3.1 期望行为（design §6.4 E4）
**测试类型：** 边界条件
**优先级：** Medium

**前置条件：**
- Given Jacobi maxIter=5 + interval=100

**测试步骤：**
1. AUTO：单测 E4s/R-E4s（循环内零判定：`residual_allreduce.calls == 0`；`residual_confirm.calls == 1`；lastResidual = 独立重算第 5 步后真实残差）

**期望结果：**
- 出口确认扫描兜底，lastResidual 恒为真实残差

**实际结果：** 通过（Debug/Release 全量内）。

**状态：** PASS

---

### ST-027：残差 helper 前置条件（E5）

**关联需求：** design §6.4 E5（走查口径）
**测试类型：** 异常处理（走查）
**优先级：** Low

**前置条件：**
- Given trueResidualSquaredLocal 的调用方契约

**测试步骤：**
1. WALK：evidence ③——全调用点核验「先 exchange + applyPhysicalBoundary」；residual.hpp:15-19 接口注释显式声明前置条件

**期望结果：**
- 运行时无法确定性断言的契约由走查 + 注释覆盖

**实际结果：** 走查结论落 evidence；review S5 核验通过。

**状态：** PASS

---

### ST-028：NFR 汇总验证

**关联需求：** srs §4 全表
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 双模式构建 + 基准脚本

**测试步骤：**
1. E2E/AUTO/BENCH：逐项核验 srs §4 表

**期望结果与实际结果（逐项）：**
- 残差语义：final_residual 独立验证一致 → ST-001/021 e2e + 确认扫描精确一致断言（U1-U5）PASS
- 位级一致：overlap on/off 逐位相同（ST-004）；CG 迭代数不变（ST-008）PASS
- 内存安全：ASan/UBSan 全量零报告（Debug 25/25）；Profiler 无新竞争（ST-016 走查口径）PASS
- 构建：双模式 0 新增警告（T008 Release 构建验证）PASS
- 兼容性：既有开关行为不变（parse_error 等既有 ctest 绿）；JSON 字段兼容（ST-018）PASS
- 性能门槛：Jacobi np1 +4.0%（±5% 带内）；RBGS 报告制如实（k1 +34%/k10 +2.7%）；CG 不劣于现状（ST-007）；bench_ar004.sh exit 0 PASS

**状态：** PASS

---

## 执行摘要

| 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|
| 28 | 28 | 0 | 0 |

## ST 执行报告

| 字段 | 内容 |
|------|------|
| 执行日期 | 2026-10-08 |
| 执行结果 | PASS |
| 执行轮次 | 第 1 轮 |

### 需求覆盖矩阵

| 需求 ID | 需求描述 | 测试用例 | 结果 |
|--------|---------|---------|------|
| §3.1 R1 | 真实残差收敛语义 + 检查频率 | ST-001, ST-002, ST-003, ST-004, ST-005, ST-006, ST-022, ST-023, ST-026, ST-027 | PASS |
| §3.2 R2 | CG 向量循环并行化 | ST-007, ST-008, ST-009, ST-024 | PASS |
| §3.3 R3 | First-touch 初始化 | ST-010, ST-011, ST-012 | PASS |
| §3.4 R4 | 拓扑一致性 | ST-013, ST-014, ST-015 | PASS |
| §3.5 R5 | 观测自愈（Profiler + 区名） | ST-016, ST-017, ST-018 | PASS |
| §3.6 R6 | 文档同步 | ST-019, ST-020 | PASS |
| §4 NFR | 非功能需求全表 | ST-004, ST-008, ST-016, ST-018, ST-024, ST-025, ST-028 | PASS |

**需求覆盖率：** 6 / 6 功能需求 + NFR 全表（100%）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 15 | 15 | 0 | 0 |
| 边界条件 | 4 | 4 | 0 | 0 |
| 异常处理 | 4 | 4 | 0 | 0 |
| 回归测试 | 5 | 5 | 0 | 0 |
| **合计** | **28** | **28** | **0** | **0** |

**关键执行记录（本轮实测）：**
- 全量 ctest：Debug（ASan+UBSan）25/25 @ 91.67s；Release 25/25 @ 24.52s（串行）
- E2E e2e（Release，OMP=1 默认）：jacobi64 17327 iters / 9.9874e-07；rbgs64 8812 / 9.9886e-07；cg64 183 / 8.8812e-07；rbgs10 8820 iters / 9.8037e-07 / residual_allreduce=882=⌊8820/10⌋ / residual_confirm=1
- 基准引用（T008 evidence，≥3 次中位数）：Jacobi np1 +4.0%（±5% 带内）；CG np1 +3.0% / np4 −7.3%；RBGS k1 +34% / k10 +2.7%；收敛口径 jacobi64 17327 / rbgs64 8812 / cg64 198

### 遗留问题

| 严重性 | 描述 | 处理方式 |
|-------|------|---------|
| Minor | B5「RBGS 三档迭代数差 < k」未作独立自动化固化用例（语义由 R-I* 计数 + U4 + e2e/bench evidence 分片覆盖） | 延后处理（review D4 已补记，列后续 AR 候选） |
| Minor | TSan 实证不可用（WSL ASLR + OpenMPI pmix 兼容性），Profiler 无竞争以 U8 稳定性 + ASan/UBSan + 走查佐证 | 延后处理（evidence 已注记，环境限制） |

### 结论

> **Go（建议）**：28/28 用例通过，需求覆盖率 100%（6/6 功能需求 + NFR 全表），无 Critical/Major 缺陷，无回归（双模式全量 25/25）。遗留 2 项 Minor 均为测试固化/环境限制类，已有记录与佐证，不影响核心功能。
