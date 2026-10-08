# [AR004] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR004 |
| AR 主题 | honesty-performance-pack（诚实性修复 + CG 并行化） |
| 关联 SR | 无（源自 docs/OPTIMIZATION_GUIDE.md §5 路线 ★1，条目 A1/A2/A3/A4/B2a） |
| 日期 | 2026-10-08 |
| 状态 | 已确认（用户授权自动化执行，依据 OPTIMIZATION_GUIDE.md §5 ★1） |

## 1. 背景与目标

OPTIMIZATION_GUIDE.md §3 诊断出 6 项「宣称 vs 实现」缺口（G 表）与 7 项性能级缺口（P 表）。本 AR 是路线 ★1「诚实性修复包」：消灭**可由本 AR 修复的**宣称缺口（G1/G3/G5/G6 + P6 拓扑隐患；G2 对应通信条目 C2、G4 为持续项，均在 Out of Scope），并修复最紧迫的性能缺口 P1（CG 向量循环单线程）。核心动机：教学样板的诚实性优先级最高，且 CG 是项目"快 1-2 个量级"的招牌，其串行更新循环是明显的性能漏洞。

## 2. 需求范围

**In Scope（对应指南条目）：**
- A1 First-touch 真实现（G1）
- A2 真实残差语义 + `--residual-check-interval`（G3、P4；覆盖 Jacobi 与 RBGS）
- A3 拓扑一致性：消除双重 `Cart_create`（P6）
- A4 观测自愈：Profiler 线程安全/per-thread 聚合 + 剖面区名修正（G5、G6）
- B2a CG 向量循环并行化（P1）
- R6 文档同步（README/AGENT_SPEC 宣称对齐、PERFORMANCE 新基准、GUIDE 勾选）

**Out of Scope：**
- B2b pipelined CG、B1 多重网格、B3 Chebyshev、B4 RBGS 通信（后续 AR）
- C/D/E/F 层全部条目（AR005+）
- CG 的残差语义变更（CG 的 r·r 本就是真实残差，仅需区名修复与并行化）

## 3. 功能需求

### 3.1 R1 真实残差收敛语义 + 检查频率（A2，Jacobi/RBGS）

**描述：** `--tol` 判据从递推残差（更新量范数）改为真实残差 ‖Au−f‖₂；新增 `--residual-check-interval k` 控制归约频率。

**触发条件：** 任意 Jacobi/RBGS 求解运行。

**期望行为：**
- 收敛判据与报告的 `final_residual` 均为真实残差口径；Jacobi 保留 AR002 的融合累积优化（更新量与残差有精确代数关系 r = D·diff，即融合累积 ×denom 即为**上一步迭代**的精确真实残差），判据通过后再做一次确认扫描报告当前迭代的精确 ‖Au−f‖；RBGS（就地更新无此代数关系）每 k 步做一次专用残差扫描
- `--residual-check-interval k`（默认 1）：每 k 步一次 Allreduce 判定；k>1 时可能多迭代 ≤k 步才停（文档注明）
- CG 不引入 interval（beta 计算依赖每步 dot(r,r)，无法跳过）；行为不回归

**异常处理：** `--residual-check-interval` 非法值（非正整数）→ 可读错误 + 退出码 1（与现有 CLI 风格一致）；k 对 CG 传入时告警忽略（对齐 `--overlap-comm` 先例）。

**验收标准：**
- Given 制造解（sin）`--tol 1e-6` 收敛，Then 报告的 `final_residual ≤ 1e-6`（真实残差口径，修复 OPTIMIZATION_GUIDE G3 记录的 204.6 失实）
- Given `--residual-check-interval 10`，When 查看 profiler 报告，Then `residual_allreduce` 调用数/迭代 ≈ 0.1
- Given 未传 interval（默认 1），Then 每步判定（语义等价现状频率）
- Given Jacobi overlap on/off，Then u 场与残差保持位级一致（延续既有约束）

### 3.2 R2 CG 向量循环并行化（B2a）

**描述：** CG 每迭代的向量更新循环（p/r/u 等 axpy/scale 类）全部 `omp parallel for + omp simd`。

**触发条件：** `--solver cg` 任意运行。

**期望行为：** 并行化范围：p 更新、r/u 更新、以及 iterate() 内对应循环；循环分区方式与既有 `axpyInterior`/`dotGlobal` 范式一致；数值结果位级不变（同分区同序，纯并行化不改浮点运算顺序的累加结构——向量内逐元素运算无归约）。

**异常处理：** 无新增异常路径。

**验收标准：**
- Given 256² np=1/np=4 基准（≥3 次中位数），Then CG `iter_time_ms` 不劣于现状（回归门槛），改进幅度如实记入 PERFORMANCE
- Given 制造解收敛测试，Then 收敛迭代数与并行化前一致（逐元素运算不改数值）
- Given 全量 ctest，Then 全绿

### 3.3 R3 First-touch 初始化（A1）

**描述：** `Subgrid::zeroInitialize()` 从单线程 `std::fill` 改为 OpenMP 并行，初始化分区与计算循环一致。

**触发条件：** 任意运行（构造子域后初始化）。

**期望行为：** 三个场（u/uNext/rhs）按与 stencil 计算相同的外层维度分区（j/k）并行初始化；全缓冲覆盖语义不变；`applyDirichletBC` 的 first-touch 一致性在 design 阶段评估（其填充为 setup 期单次执行，优先级低于逐迭代路径）。

**验收标准：**
- Given 代码走查，Then 初始化循环分区与 stencil 循环分区一致（AGENT_SPEC.md:54 first-touch 约定落地）
- Given 全量 ctest，Then 全绿（初始化语义不变）
- Given 文档走查，Then README:16 与 AGENT_SPEC 宣称有实现对齐，G1 勾选

### 3.4 R4 拓扑一致性（A3）

**描述：** 消除 `UniformPartition::partition()` 与 main.cpp 的双重 `MPI_Cart_create`，进程拓扑单一来源。

**触发条件：** 任意多/单进程运行。

**期望行为：** 笛卡尔拓扑只创建一次；分解 offsets、邻居关系、coords 来自同一 cart comm（消除 `reorder=1` 下两次创建可能映射不一致的隐患）。接口调整方式（partition 复用外部 comm 或返回自建 comm）由 design 决定并在 design.md 记录授权（AGENT_SPEC「接口变更须审批」）。

**验收标准：**
- Given 新增 np=4 一致性测试，Then offsets/邻居/coords 与单次创建的拓扑一致
- Given 代码走查，Then `MPI_Cart_create` 在进程生命周期内仅调用一次
- Given 全量 ctest，Then 全绿

### 3.5 R5 观测自愈（A4）

**描述：** Profiler 线程安全化（per-thread 统计 + report 聚合）；修正复制粘贴残留的剖面区名（RBGS/CG 误用 `jacobi_iteration`）；hpp 注释与实现对齐。

**触发条件：** `--enable-profiling` 运行；OpenMP 并行区内使用 profiler。

**期望行为：**
- Profiler 统计改为 per-thread 累积（无锁或低竞争），report 时聚合各线程；对既有报告消费方（performance_report.json、run_scaling_tests.py）**字段格式兼容**
- RBGS 区名改 `rbgs_iteration`，CG 区名改 `cg_iteration`（Jacobi 保持 `jacobi_iteration`）；依赖区名的测试/脚本同步
- hpp:13 "atomic counters" 注释更正为实际机制

**验收标准：**
- Given OpenMP 并行区调用 `HYPOS_PROFILE`（新增测试或走查），Then 无数据竞争（TSan 或走查口径记录 evidence）
- Given RBGS/CG profile 报告，Then 区名正确（G5 勾选）
- Given 并行区外使用场景，Then 报告数值与改造前口径一致（聚合不改变总时长语义）

### 3.6 R6 文档同步

**描述：** 文档与实现零缺口（延续 AR002/AR003 惯例）。

**期望行为：** README（first-touch 宣称对齐、`--residual-check-interval` 帮助文本）、AGENT_SPEC（如接口变更则同步）、PERFORMANCE（A2 开销对比 k=1/5/10、B2a CG 前后对比，≥3 次中位数）、OPTIMIZATION_GUIDE（G1/G3/G5/G6、P1/P4/P6 勾选 + 回填实测数据）。

**验收标准：**
- Given 文档走查，Then 无失实宣称；GUIDE 勾选与本 AR 完成状态一致
- Given PERFORMANCE，Then 含真实残差开销与 CG 并行化两组数据及口径说明

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 正确性 | 残差语义 | 报告的 final_residual 与独立计算 ‖Au−f‖₂ 一致（测试内独立验证；Jacobi/RBGS 为确认扫描的精确一致，CG 为递推残差——理论恒等 b−Ax、浮点下相对容差一致） |
| 数值稳定 | 位级一致 | Jacobi overlap on/off 位级一致约束保持；CG 并行化不改变收敛迭代数 |
| 内存安全 | ASan/UBSan | 全部测试零报告；Profiler 改造后无新竞争（TSan 走查） |
| 构建 | 双模式 | 0 新增警告 |
| 兼容性 | CLI/接口 | 既有开关行为不变；profile 报告 JSON 字段兼容 |
| 性能回归 | 基准门槛 | Jacobi iter_time_ms 不劣于现状 ±5%（换算零迭代内开销，仅终点一次确认扫描）；RBGS 真实残差存在必然开销（每 k 步一次扫描 + halo 交换），不设回归门槛——以 k=1/5/10 三档如实测量并记入 PERFORMANCE（诚实性的必要代价）；CG iter_time_ms 不劣于现状（中位数口径） |

## 5. 约束与假设

**约束：**
- 延续 AGENT_SPEC 全部约定；接口变更（R4）须在 design.md 记录授权
- 不引入外部依赖；C++17 + MPI + OpenMP
- 测试环境：WSL（OpenMPI 栈，AR003 结论：mpich 在 Ubuntu 24.04 有单进程退化缺陷）；基准 OMP=1、测试 OMP=4（延续既有口径）
- 单 AR 净新增代码 ≤ ~800 行（GUIDE §6 复杂度预算）

**假设：**
- Jacobi 保留融合累积的代数关系（r = D·diff，ω=1 纯 Jacobi 无加权）成立——design 阶段以代码核验
- 多节点 NUMA 收益在 WSL 单机不可测，A1 验收以代码正确性 + 宣称对齐为准（GUIDE A1 判据）
- 递推残差→真实残差的判据尺度变化会使既有固定迭代数测试的期望值偏移，属预期内测试更新

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| 递推残差 | 更新量范数 ‖uNew−u‖，当前实现误作收敛判据 |
| 真实残差 | ‖Au−f‖₂，A 为离散算子（2D 五点 / 3D 七点） |
| 检查频率 | 收敛判定的 Allreduce 归约间隔 k（`--residual-check-interval`） |
| first-touch | NUMA 页首次触摸策略：初始化与计算由同一线程分区执行（AGENT_SPEC.md:54） |
