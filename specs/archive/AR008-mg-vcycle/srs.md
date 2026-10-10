# [AR008] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR008 |
| AR 主题 | mg-vcycle（B1b：完整多层 V-cycle + MG-CG 预条件，复用 AR007 限制/延拓） |
| 关联 SR | 无（源自 docs/OPTIMIZATION_GUIDE.md §5 路线后续条目 B1b；依赖 AR007（B1a）的 R/P 算子与平滑入口；完成 §3「算法复杂度天花板」的完整 O(N) 叙事） |
| 日期 | 2026-10-10 |
| 状态 | 已确认（用户授权自动化执行，依据 OPTIMIZATION_GUIDE.md §5 路线依赖 B1a→B1b 与 B1 验收判据②） |

## 1. 背景与目标

AR007（B1a）落地了两层校正格式 `--solver mg2`：cycle 数规模无关（256² 制造解 14 cycles），但仅两层——最粗层 128² 仍用 CG 精解（每 cycle ~180 CG 迭代，占 mg2 耗时主体），且粗层绝对容差在极紧 tol 下饥饿（AR007 遗留 Minor，适用 tol ≳ 1e-8）。

B1b 把「两层」推广为「多层」：粗化链到 ≤8² 粗根，V-cycle 递归平滑-限制-粗解-延拓，渐近收敛率应与网格规模无关（经典 MG 理论）；再以 V-cycle 为预条件子构造 MG-CG（PCG），对标 DESIGN §6 路线图「CG 预条件子」——迭代数对规模与条件数双不敏感，是泊松 HPC 的标准答案形态。

GUIDE B1b 验收②：256²/512² 下 MG-CG 迭代数 vs 裸 CG 对比入 PERFORMANCE。

## 2. 需求范围

**In Scope（对应指南条目 B1b）：**
- 多层粗化链：2D 网格逐层减半（nx_l = nx/2^l）至粗根 ≤8² 量级（层数与粗根下界由 design 决定并记录）；逐层复用 AR007 的 restrictResidual/prolongateCorrection 纯内核与逐层 ×4 尺度补偿（(H_l/H_{l+1})² = 4）
- `--solver mgv`：完整 V-cycle 求解器（预平滑 → 限制 → 下层递归 → 延拓 → 后平滑；粗根精解）
- `--solver mgcg`：MG 预条件共轭梯度（PCG，预条件子 = V-cycle；对标 DESIGN 路线图）
- 粗层求解策略沿 AR007 复制式（逐层 Allgatherv，各 rank 冗余；粗层间无需新通信模式）——design 复核并记录
- AR007 遗留 Minor 顺带修复：粗层/粗根容差饥饿 → 相对容差（‖r‖/‖b‖ 口径，design 定具体形式）
- 顶层粗根求解器选择（CG 紧容差 vs 小规模直接法）：design 决定并记录
- 性能对照：256²/512²、np=1/4，mgcg vs 裸 cg、mgv vs mg2，迭代数与耗时入 PERFORMANCE §15（GUIDE B1b 验收②）
- 测试：V-cycle 多层正确性（层数链、逐层系数、制造解收敛、率 h 无关）、MG-CG vs CG 迭代数、np=1/4 布局一致性延续、既有 38 条全量回归（R/P 内核与 RBGS smooth 零变化）
- 文档同步：README/GUIDE/AGENT_SPEC/PERFORMANCE/--help

**Out of Scope：**
- 3D MG（延续 AR007 2D 限定，`--nz>1` → HYPOS_ERROR；3D 需 3D R/P 与七点粗算子，独立 AR）
- Neumann × MG（齐次 Neumann 校正边界 + 每层零空间投影是独立算法专题——AR007 srs 原文「留 B1b 及以后」，本 AR 拆出至后续 AR 以控规模（GUIDE §6 决策规则）；`--bc neumann` × mgv/mgcg → HYPOS_ERROR 明示）
- W-cycle / FMG / 自适应平滑策略（V(ν1,ν2) 固定调度；其他 cycle 类形后续 AR）

> **修订注记（2026-10-10，T002 实证后）**：本条与 R1 的「cycle 类形固定 V(ν1,ν2)」在实现阶段被 design D11 推翻——单 V-cycle 下逐层 ×4 补偿对插值杂散模态过补偿，深层链实证发散（齐次 [1,1] 稳态放大率 256²=2.01、512²=2.69/cycle；探针证据链归档 AR008 提交记录）；修复为 **W-cycle（γ=2，每层两次粗修正）**，其余调度（ν=2+2、×4 补偿、粗根 CG）不变。Out-of-Scope 的 W-cycle 条目由此改为「FMG / 自适应平滑策略」。R1/R3/R4 的功能验收标准不受影响（U3/U4/U5/§15 bench 全绿/实测达标）；性能语义反而更优（cycle 数 4× 细化反降）。srs 作为归档基线按实际交付结构再基线。
- 半粗化/代数 MG
- 更改 AR007 mg2 / 既有三求解器语义（mg2 保留为两层参照系，行为零变化）

## 3. 功能需求

### 3.1 R1 多层 V-cycle 求解器（mgv）

**描述：** 新增 `VCycleMGSolver : PoissonSolver`（`--solver mgv`）：粗化链上递归执行 V(ν1,ν2) 校正格式。**（修订：实现为 W-cycle γ=2，见 §2 Out-of-Scope 修订注记与 design D11；ν/补偿/判据不变）**

**触发条件：** `--solver mgv`（CLI 新值；既有 jacobi/red_black_gs/cg/mg2 不受影响）。

**期望行为：**
- 初始化：由全局 (nx, ny) 生成粗化链 [level 0(细) … level L(粗根)]，每层 2D Subgrid；粗根规模下界与层数上限由 design 定（≤8² 量级，≥4² 可解下限）；奇数维中途出现 → HYPOS_ERROR（粗化要求逐层可除 2——与 mg2 同防御口径）或由 design 定截断策略
- 每 cycle：ν1 预平滑 → 残差限制到 level 1 → 递归 V(level 1) → 延拓校正 → ν2 后平滑；level L 粗根精解（design 定 CG 紧相对容差或直接法）
- 逐层尺度补偿：level l+1 粗方程 rhs = 4 × R(残差)（沿 AR007 落地的 (H_l/H_{l+1})²=4 推导，见 PERFORMANCE §14「实现要点」与 mg_two_level_solver.cpp 粗 rhs 注释；逐层同因子，多层链累积 4^l）
- 收敛判据与既有求解器同口径（全局 L2 真残差 ≤ tolerance，外层 cycle 计）；iterate()/lastResidual()/进度回调/residualCheckInterval 沿 PoissonSolver 契约
- V-cycle 渐近率 h 无关（GUIDE B1b 理论）：512² 与 256² 同 tol cycle 数相近（验收②性能对照记录式印证）

**异常处理：** MPI/粗层路径不抛异常；不支持配置（3D、Neumann、奇数全局维、粗化链中途奇维、粗根过小）→ HYPOS_ERROR 明示退出，不产出失实结果。

**验收标准：**
- Given 制造解，When `--solver mgv` 运行至 tol（256² 与 512²），Then 两规模收敛且 cycle 数相近（h 无关定性判据，容差由 design 定）；解精度二阶（沿 AR007 二阶测试惯例）
- Given 256²/512² 制造解，When mgcg 与裸 cg 各自求解，Then mgcg 迭代数少于裸 cg（GUIDE B1b 验收②记录式对比，倍数如实入档）
- Given 任意既有负载，When `--solver jacobi/red_black_gs/cg/mg2`，Then 行为零变化（回归）
- Given 3D/Neumann/奇维，When `--solver mgv/mgcg` 启动，Then HYPOS_ERROR 且不产出结果

### 3.2 R2 粗化链与逐层算子复用

**描述：** 层级结构生成与逐层 R/P/A_l 组装，复用 AR007 纯内核。

**触发条件：** mgv/mgcg 初始化与每 cycle 的限制/延拓步。

**期望行为：**
- 粗化链生成：每层全局规模 (nx/2^l, ny/2^l)，层间 R/P 用 AR007 restrictResidual/prolongateCorrection（顶点重合约定不变），逐层布局沿用复制式（Allgatherv 各 rank 冗余粗层场）——粗层间无新通信模式（B1b 成本注记「粗层间通信用现有 exchanger 泛化」的复制式简化路线，design 复核论证）
- A_l：每层 rediscretized 五点模板 + 逐层 ×4 尺度补偿（与 AR007 mg2 的粗层方程完全同构，逐层传递）
- R/P 内核**零修改**：AR007 的已知系数场/P∘R 恒等测试不变绿转（复用契约）
- 粗根求解：规模极小（≤8²），design 定 CG 紧相对容差（‖r‖/‖b‖ ≤ 1e-12 量级）或直接消元，记录决策
- 粗层容差饥饿修复（AR007 遗留）：非粗根层的递归 V 不涉及 CG 容差（纯平滑+限制+延拓），仅粗根一处求解——相对容差口径消除饥饿（mgv/mgcg 适用 tol 范围较 mg2 放宽，实测入档）

**异常处理：** 纯内核无异常路径；链生成失败（奇维中途、层数上限）→ HYPOS_ERROR。

**验收标准：**
- Given AR007 全部 mg2 测试（38 条），When 本 AR 落地，Then 不变绿转（R/P 复用零回归）
- Given np=1 与 np=4 非均匀切分，When 各自执行多 cycle V-cycle，Then 校正场逐元素一致（跨布局一致性延续，沿 AR007 CrossLayout 模式）
- Given 粗化链（如 256→128→64→32→16→8），Then 每层规模/offset/尺度补偿逐层正确（单测可枚举验证）

### 3.3 R3 MG-CG 预条件求解器（mgcg）

**描述：** 新增 `MGPreconditionedCGSolver : PoissonSolver`（`--solver mgcg`）：PCG，预条件算子 = V-cycle（M ≈ A⁻¹ 的 MG 近似）。

**触发条件：** `--solver mgcg`（CLI 新值）。

**期望行为：**
- PCG 外层：沿既有 CGSolver 骨架（内积归约、方向更新、收敛判据），预条件步 z = M⁻¹r 以 V-cycle 执行（初始猜测 0 的校正格式近似 A⁻¹r）
- 收敛判据：真残差口径（沿 G3 修复惯例）≤ tolerance；迭代数对网格规模不敏感（MG 预条件把条件数从 O(h⁻²) 压到 O(1)，理论预期迭代数常数级）
- 既有 CGSolver 零变化（新增类组合/泛化，不修改 cg_solver.cpp 语义——沿 R3 复用契约）
- 预条件 V-cycle 的平滑步计数与 profiler 区名（沿 AR004 惯例，如 mgcg_iteration/mg_precycle）由 design 定

**异常处理：** 同 R1 防御口径；预条件 V-cycle 内部不抛异常。

**验收标准：**
- Given 制造解（256²/512²），When mgcg 求解至 tol，Then 收敛且迭代数少于裸 cg（记录式；GUIDE B1b 验收②，数据入 PERFORMANCE §15）
- Given np=4，When mgcg 求解，Then 与 np1 结果一致（真残差轨迹/解一致口径，沿 CrossLayout 惯例）
- Given `--solver cg` 任意负载，Then 行为零变化（回归）

### 3.4 R4 性能对照入档

**描述：** mgcg vs 裸 cg、mgv vs mg2 的同 tol 对比，如实入 PERFORMANCE §15。

**期望行为：**
- 256² 与 512²、np=1/4、OMP=1、≥3 次中位；记录：外层迭代数、总平滑步数、粗根求解次数/迭代、墙钟、真残差终值
- 512² 为 GUIDE B1b 验收②明示规模（srs §5 假设：单次运行超 60s 的项如实注明）
- 记录式判定（沿 D6 惯例）：不预设倍数，实测入档；理论定性（迭代数规模无关）文字印证

**验收标准：**
- Given PERFORMANCE §15 新小节，Then 含 mgcg vs cg（256²/512²）与 mgv vs mg2 两组数据及口径说明，无失实宣称

### 3.5 R5 测试矩阵

**描述：** 单测 + np 集成 + e2e 全覆盖（GUIDE B1b 验收②③④延续）。

**期望行为：**
- 单测（单进程）：粗化链生成正确性（层数/规模/offset/尺度补偿逐层）、mgv 制造解收敛 + 率 h 无关（256² vs 512² cycle 数相近断言，容差明确）、mgcg vs cg 迭代数、防御（3D/Neumann/奇维/粗化链失败）、二阶精度（沿 AR007 解析 rhs 惯例）
- np 集成：mgv/mgcg 制造解收敛 np=1/4、跨布局一致性（非均匀 {15,19} 沿 AR007 惯例）
- e2e：直跑 hypos `--solver mgv`/`--solver mgcg`
- ctest 条目沿 AR007 惯例（HYPOS_EXPECT_NP 探针守卫、专用 filter；38 → 新增条目）

**验收标准：**
- Given 全量 ctest（38 既有 + 新增），When 双模式运行，Then 全绿且 ASan/UBSan 零报告

### 3.6 R6 文档同步

**描述：** 文档与实现零缺口（沿 AR002-AR007 惯例）。

**期望行为：** README 求解器表加 mgv/mgcg + 路线图（B1b 勾选，后续条目仍列未实现）；GUIDE §3.4 天花板完整突破表述更新（O(N) 叙事完整兑现）+ ★路线表 B1b 勾选；AGENT_SPEC 类层级 +TwoLevelMG 两个新类；PERFORMANCE §15；`--help` solver 说明同步（solver 行 + residual-check-interval 行 + solver.hpp 接口注释）。

**验收标准：**
- Given 文档走查 + `--help` 实跑，Then 无失实宣称、帮助与行为一致

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 正确性 | 制造解收敛 | 二阶精度；np1/np4 布局一致；率 h 无关定性 |
| 兼容性 | 既有行为 | 既有四求解器（含 mg2）与全部 38 条测试零变化；PoissonSolver 接口签名不动 |
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 资源 | MPI 对象 | 逐层缓冲/通信对象生命周期完整 |
| 构建 | 双模式 | 0 新增警告 |
| 性能 | 数据如实 | 记录式判定，反直觉结果如实入档（沿 D6） |
| 回归 | 全量测试 | Debug+Release ctest 全绿 |
| 复杂度 | 预算 | 严控：目标净新增 ≤ ~800 行（GUIDE §6.4；mgv/mgcg 共享 V-cycle 核心，层级结构复用 mg2 的组装逻辑——design 论证复用路径） |

## 5. 约束与假设

**约束：**
- 延续 AGENT_SPEC 全部约定（row-major、MPI 路径不抛异常、snake_case、无外部依赖、C++17+MPI）
- 2D（nz=1）与偶数全局维延续 AR007 限定；粗化链中途奇维处理由 design 定（拒绝或截断，记录决策）
- 复制式粗层延续（逐层 Allgatherv；B1b 原文「现有 exchanger 泛化」的分布式粗层路线不采纳时须论证——设计阶段复核）
- ν1/ν2 沿 AR007 默认（2/2）；cycle 类形固定 V(ν1,ν2)
- CLI 新值最小化（mgv/mgcg 两值；层数/粗根阈值不新增 CLI 项，design 内定）
- R/P/平滑内核零修改（复用契约）；mg2 行为零变化（两层参照系保留）
- 测试环境：WSL 单机 OpenMPI；np>1 条目带 HYPOS_EXPECT_NP 探针守卫

**假设：**
- 制造解 rhs 由测试内构造（沿 AR007 解析 rhs 惯例）
- 粗根规模极小（≤8²），其求解成本可忽略（实测验证）
- 512² 单次运行可控（mgv/mgcg 预计秒级；裸 cg 512² 迭代数 ~O(n) 量级、秒-分钟级——超 60s 如实注明）
- 复制式逐层 Allgatherv 在 np=4 小规模下开销可接受（粗层场递减；实测验证，若反直觉如实入档）

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| V-cycle | 递归校正格式 V(l)：l 层预平滑→限制→V(l+1)→延拓→后平滑；粗根精解。V(ν1,ν2) 表示平滑步数 |
| 粗根 | 粗化链最粗层（≤8² 量级），直接精解不再递归 |
| MG-CG / PCG | 预条件共轭梯度，预条件子 M⁻¹ ≈ A⁻¹ 由 V-cycle 近似（z_k = V-cycle(r_k) 初始 0） |
| 逐层尺度补偿 | 粗方程 rhs = (H_l/H_{l+1})² × R(残差) = 4 × R(残差)，沿 AR007 偏差记录 6 推导逐层同因子 |
| h 无关收敛 | V-cycle 渐近收敛率不随网格加密劣化（多层 MG 理论核心；两层格式在最粗层仍 O(n) 求解时不完全成立） |
| 复制式粗层 | 每层粗问题在所有 rank 上冗余组装求解（逐层 Allgatherv），无粗层点对点通信 |
