# [AR009] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR009 |
| AR 主题 | pipelined CG（pcg）+ 残差历史输出（E2） |
| 关联 SR | GUIDE §4 B2b（P2）/ E2（路线图后续行） |
| 日期 | 2026-10-10 |
| 状态 | Draft |

## 1. 背景与目标

AR004（B2a）已并行化 CG 向量循环，但每迭代仍有 **2 次阻塞 `MPI_Allreduce`**（`cg_solver.cpp`：pap 与 rhoNew 两处 `cgDotGlobal`——GUIDE B2b 原文「3 归约/iter」为未核实口径，以实测 2 为准，本 AR 顺带修正文档）。阻塞归约在 np>1 时让全部 rank 空等，是 CG 迭代时间的同步成本主体。Gropp 管线化 CG（ICG）以非阻塞 `MPI_Iallreduce` + 延迟一步的标量更新把归约与向量运算重叠，理论上归约延迟被完全隐藏。

E2 残差历史是 GUIDE 明示「随任意算法 AR 顺带落地」的放大器基础设施：`--residual-history <file>` 把每迭代（cycle）残差写 CSV，配绘图脚本，任何算法改动的收敛效果即刻可视化（含 AR002 §9 对比叙事的复现基础）。

**目标**：① 新增 `--solver pcg`（管线化 CG，与裸 CG 同收敛迭代数 ±5%，每迭代阻塞归约 0 次）；② `--residual-history` CSV 全求解器可用 + 绘图脚本；③ pcg vs cg 对照入 PERFORMANCE §16（记录式）。

## 2. 需求范围

影响模块：`src/solver/`（新 pcg 求解器，复用 cg_kernels；CGSolver 零修改）、`src/main.cpp`（CLI 两项 + 分派 + 进度回调接线）、`src/utils/cmdline_parser`（不改）、`scripts/`（绘图脚本）、`tests/`（新增测试文件）、文档四处 + PERFORMANCE §16。

**In Scope（本 AR 要做的）：**
- Pipelined CG 求解器（`--solver pcg`）：Gropp ICG 结构（u/x/r/s/v/w/z/q 向量组，Iallreduce 重叠，alpha 延迟一步），复用 FP4 cg_kernels 内核
- `--residual-history <file>`：rank0 写 CSV（iter,residual），全求解器（jacobi/red_black_gs/cg/pcg/mg2/mgv/mgcg）统一出口
- 绘图脚本（`scripts/plot_residual_history.*`）出收敛曲线图
- pcg vs cg 基准入 PERFORMANCE §16（256²/512²、np=1/4、OMP=1、≥3 次中位、记录式）
- GUIDE B2b 勾选 + 「3→1 归约」口径修正为实测；E2 勾选

**Out of Scope（本 AR 不做的）：**
- 通信 Sensitivity/Chronopoulos-Gemignani 变体（仅 Gropp ICG 一种管线）
- CGSolver 本体行为任何变化（`--solver cg` 路径字节不动=复用契约，沿 AR008 FP4 惯例）

> **修订注记（2026-10-10，design D2 后）**：「cg_solver.cpp 零改动」措辞修正——W1 探针（阻塞归约 0 断言的非空洞自证）要求给 `cgDotGlobal` 的 `MPI_Allreduce` 加 profiler region，实现方式为拆分 `cgDotLocal`+`cgBlockingAllreduce`（净增 ~15 行，位于 cg_solver.cpp）。**`--solver cg` 路径浮点序位级不变**（循环体逐字搬移+同参 Allreduce，既有 cg golden 回归守护）——「字节不动」的实质契约（行为/数值零变化）保持，物理文件触碰属探针最小代价，沿 R1 修订注记机制如实记录。
- F3 CLI 类型校验统一（独立工程卫生项）
- 残差历史的实时可视化/在线分析（仅离线 CSV+图）
- Chebyshev/其他求解器（B3 独立 AR）

## 3. 功能需求

### 3.1 R1 Pipelined CG 求解器（pcg）

**描述：** 新增 `PipelinedCGSolver : PoissonSolver`（`--solver pcg`）：Gropp 管线化共轭梯度，每迭代 2 次 `MPI_Iallreduce` 与向量运算重叠，阻塞归约 0 次。

**触发条件：** `--solver pcg`（CLI 新值；既有 solver 值不受影响）。

**期望行为：**
- 数学等价 ICG（Gropp 2003）：延迟一步的 alpha/beta 更新，s/v/w/z/q 工作向量组；收敛判据为真残差口径（沿 G3 惯例——递推残差 s 范数仅作步进判据，退出前真残差复核）

> **修订注记（2026-10-10，设计阶段 D1 自审后）**：本条的「2 次 Iallreduce 与向量运算重叠」为文献口径的乐观转述——设计阶段推导证明 CG-1 单步结构下 p/q 更新依赖新鲜归约值算 β，归约窗口内无可用向量工作，**真实重叠需两步前瞻（CG-2），不在本 AR**。实际交付结构：**每迭代 1 次融合非阻塞 Iallreduce（[ρ, m] 双标量打包）**，收益=同步次数与消息数减半（2 阻塞→1 非阻塞）+ q=A p 递推免第二次 matvec；α 分母经 ν-递推（ν_m = m_m − β²ν_{m−1}）与 CG 精确等价。验收标准不受影响（阻塞归约 0 断言、迭代数 ±5%、np4 一致均照旧成立）；「pipelined」命名保留（标量更新的管线化重组），无失实宣称由 R4/R5 口径与 §16 诚实声明守护。
- 每迭代阻塞 `MPI_Allreduce` 调用数 = 0（profiler/计数断言）；2 次 `MPI_Iallreduce` 发出后与 axpy/matvec 内核重叠
- 数值注记（诚实记录）：延迟归约引入额外舍入漂移，长迭代链上残差递推偏差可能放大（Gropp 原文已知）；实测迭代数与裸 CG 对照 ±5% 判定，失配如实入档
- pap≤0 breakdown 防御（沿 cg/mgcg 惯例：WARN + 停止迭代）
- 初始猜测 0、r₀=−rhs、收敛后 applyPhysicalBoundary（沿 CGSolver 契约）
- `iterate()` 续态语义沿 PoissonSolver 契约（progress 回调、lastResidual）

**异常处理：** MPI 路径不抛异常；breakdown → HYPOS_WARN + 返回已完成迭代数，不产出失实收敛宣称。

**验收标准：**
- Given 制造解（256²/512²，多模式负载），When pcg 与 cg 各自求解至 tol，Then 迭代数一致 ±5%，且最终真残差 ≤ tol
- Given np=4，When pcg 求解，Then 与 np=1 解一致（L2 ≤ 1e-10·‖u‖ 口径，沿 U9 惯例）
- Given 任意既有负载，When `--solver cg`（及其他五值），Then 行为零变化（全量回归）
- Given profiler/计数探针，When pcg 每迭代执行，Then 阻塞 Allreduce 计数 = 0

### 3.2 R2 残差历史输出（E2）

**描述：** `--residual-history <file>`：求解过程每迭代（外层 cycle）记录残差到 CSV，rank0 单点写出。

**触发条件：** CLI 传入 `--residual-history <path>`（可选；缺省不写）。

**期望行为：**
- CSV 格式：首行 `iteration,residual`，其后每迭代一行（真实残差口径，与各求解器收敛判据同源）；mgv/mg2 为 cycle 计，jacobi/rbgs/cg/pcg/mgcg 为迭代计
- 全部六个 solver 值统一出口（进度回调链路或 driver 侧统一 hook，design 定；不改各求解器签名）
- np>1 时仅 rank0 写文件（残差本身是全局归约结果，各 rank 一致）
- 与 `--output-dir` 无耦合（独立路径；目录不存在 → HYPOS_ERROR 明示）
- 每 k 步采样沿 `--residual-check-interval` 既有语义（不新增 CLI）

**异常处理：** 文件不可写 → HYPOS_ERROR 明示退出（沿 I/O 惯例）；无该选项时零开销。

**验收标准：**
- Given `--solver cg --residual-history f.csv`，When 运行至收敛，Then CSV 行数=迭代数+1，~~残差列单调下降（容许数值噪声非严格）~~ 且末行 ≈ final_residual（**修订 2026-10-10/T004**：单调子句删除——cg 残差 2-范数实测非单调（45-47 行 0.153→0.204→0.246），理论上单调的仅误差 A-范数；判据以行数+末行打印精度内一致为准，见 design §6 E1）
- Given 六个 solver 值各跑一次，When 带 --residual-history，Then 六份 CSV 均合法且首末残差与各自收敛语义记录一致（mgcg=finish 行残差；其余=performance_report final_residual / Converged 行——**per-solver 口径与打印精度内比较见 design §6 E2**，finish 行大多不含残差值）
- Given np=4，When 输出残差历史，Then 仅 rank0 产出一份文件，~~内容与 np=1 轨迹一致~~（**修订 2026-10-10/T003**：np1 轨迹对照前提不成立——np1 fixture OMP=4 线程 vs np4 每 rank 1 线程，归约分块序不同使 cg 轨迹真实分岔（迭代数差 >±1）；判据改为 np4 工件内部一致性（header+行数==iterations+末行打印精度内==final_residual），np1 对照降级 evidence-only，见 design §6 E3）

### 3.3 R3 绘图脚本

**描述：** 消费 R2 CSV 出收敛曲线图（对比多求解器）。

**触发条件：** 手工/脚本调用（非运行时依赖）。

**期望行为：**
- `scripts/plot_residual_history.py`（matplotlib）：输入一或多份 CSV，输出 PNG（对数 y 轴收敛曲线，多文件叠加对比）
- 环境无 python3/matplotlib 时脚本不阻塞主流程（资产交付；本 WSL 环境验证口径见 §5 假设）

**验收标准：**
- Given 两份 CSV（如 cg vs pcg），When 运行脚本，Then 产出含两条曲线的 PNG（在可行环境验证并记录；不可行环境以脚本语法自检+CSV 契约测试代偿，如实注明）

### 3.4 R4 性能对照入档

**描述：** pcg vs cg 同 tol 对照，如实入 PERFORMANCE §16。

**期望行为：**
- 256²/512²、np=1/4、OMP=1、≥3 次中位；记录：迭代数、墙钟、（若可得）归约等待耗时；GUIDE B2b 判据②迭代数 ±5% 印证
- 记录式判定（沿 D6）：单机 WSL np4 下归约隐藏收益可能微小（共享内存 MPI 归约本就快）——如实记录，不预设倍数

**验收标准：**
- Given PERFORMANCE §16 新小节，Then 含 pcg vs cg 两组规模×np 数据与口径说明、迭代数一致性印证、无失实宣称

### 3.5 R5 测试矩阵与文档

**描述：** 单测+np 集成+e2e 全覆盖；文档零缺口。

**期望行为：**
- 单测：pcg 数学正确性（vs CG 迭代数/解）、breakdown 防御、阻塞归约 0 断言（计数探针）、CSV 契约（行数/格式/单调性）
- np 集成：pcg np=1/4 一致、残差历史 np4 单文件
- e2e：直跑 `--solver pcg` 与 `--residual-history`
- ctest 新增条目沿惯例（42 → 4x+；HYPOS_EXPECT_NP 守卫）
- 文档：README（solver 表+CLI 表+目录树）、GUIDE（B2b 勾选+口径修正、E2 勾选）、AGENT_SPEC（类层级 +PipelinedCGSolver）、PERFORMANCE §16、`--help`（solver 行两处提及）

**验收标准：**
- Given 全量 ctest（42 既有 + 新增），When 双模式运行，Then 全绿且 ASan/UBSan 零报告
- Given 文档走查 + `--help` 实跑，Then 无失实宣称

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 正确性 | 制造解收敛 | pcg 与 cg 迭代数 ±5%、解 L2 一致口径 |
| 兼容性 | 既有行为 | 六既有 solver 值与 42 条测试零变化；PoissonSolver/cg_kernels 接口签名不动（新增类不改旧文件语义；cg_solver.cpp 的 FP2 拆分触碰以 §2 修订注记为准——行为/数值零变化、物理触碰 ~15 行） |
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 构建 | 双模式 | 0 新增警告 |
| 性能 | 数据如实 | 记录式；单机收益微小时如实注明 |
| 回归 | 全量测试 | Debug+Release ctest 全绿 |
| 复杂度 | 预算 | 净新增 ≤ ~500 行（pcg ~220 + CSV 接线 ~120 + 脚本/测试另计沿惯例 ~600） |

## 5. 约束与假设

**约束：**
- 沿 AGENT_SPEC 全部约定（row-major、MPI 路径不抛异常、无外部依赖、C++17+MPI+OpenMP）
- CGSolver 零修改（`--solver cg` 字节级不变；pcg 为新增独立类，复用 cg_kernels 自由函数——沿 AR008 FP4 双消费惯例）【以 §2 修订注记为准：FP2 探针拆分致 cg_solver.cpp 物理触碰 ~15 行，浮点序位级不变由 cg golden 回归守护】
- CLI 新增恰 2 值（`pcg` solver 值 + `--residual-history` 选项）；不新增其他 CLI
- `MPI_Iallreduce`/`MPI_Wait` 错误码检查沿 MPI_ENV 惯例
- 绘图脚本为独立资产（无 python 运行时依赖进构建/测试主链路）

**假设：**
- OpenMPI 单机 np4 支持 Iallreduce 与计算真重叠（共享内存设备一般成立；收益大小实测入档）
- WSL 环境无 python3/matplotlib（AR008 会话实测）：绘图脚本以「可行环境验证」口径交付——若本环境确不可跑，PNG 生成以语法自检（py_compile 若可行）+ CSV 契约测试代偿并如实记录
- 延迟一步的 ICG 舍入漂移在 512² 制造解（迭代数 O(10³)）内不破坏 ±5% 判据（Gropp 原文与后续文献支持；实测验证）
- 残差历史开销可忽略（每迭代一次 rank0 fprintf；实测确认不劣化迭代时间）

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| Pipelined CG / ICG | Gropp（2003）管线化共轭梯度：dot 积用 MPI_Iallreduce 非阻塞归约，标量（alpha/beta）延迟一步使用，归约延迟与向量内核重叠 |
| 阻塞归约计数 | 单次迭代内 `MPI_Allreduce` 调用次数（本 AR 判据：pcg=0） |
| 残差历史 | 逐迭代（cycle）真实残差的 CSV 序列，收敛曲线数据源 |
| 延迟一步 | 第 k 迭代发出 Iallreduce、第 k+1 迭代 Wait 并消费——隐藏一整迭代的归约延迟 |
