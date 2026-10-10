# [AR007] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR007 |
| AR 主题 | mg-two-level（算法深水区第一步：B1a 两层几何多重网格校正格式） |
| 关联 SR | 无（源自 docs/OPTIMIZATION_GUIDE.md §5 路线 ★4，条目 B1a；修复 §3「算法复杂度天花板」：Jacobi/GS 类 O(N²) 迭代数，无 O(N) 级算法） |
| 日期 | 2026-10-09 |
| 状态 | 已确认（用户授权自动化执行，依据 OPTIMIZATION_GUIDE.md §5 ★4 与 B1 验收判据） |

## 1. 背景与目标

GUIDE §3 诊断：当前全部优化都在压「每步多快」（IO/通信/内存），而 Jacobi/RBGS/CG 对低频误差收敛慢（O(N²) 型迭代数）——多重网格改变「需要多少步」。B1a 是旗舰条目 B1 的第一步：**两层校正格式，以正确性验证为先**——细网格 RBGS 预平滑 → 残差限制到粗网格 → 粗网格 CG 精解 → 延拓校正 → 后平滑。收敛理论（2 层即消除全部可分辨高频）可直接对比实测迭代数（GUIDE B1a 验收①）。

现有积木：RBGS（现成 smoother，验收③要求复用且不改其位级一致测试）、CG（现成粗网格解法器）、均匀分解与 Subgrid offsets（现成限制/延拓的全局索引结构）、`trueResidualSquaredLocal`/`globalTrueResidual`（现成残差设施）——离 MG 只差粗网格生成、限制/延拓算子与两层驱动循环。

## 2. 需求范围

**In Scope（对应指南条目 B1a）：**
- 新求解器 `TwoLevelMGSolver`（实现 PoissonSolver 接口，`--solver mg2`）：预平滑 → 残差限制 → 粗网格 CG 精解 → 延拓校正 → 后平滑 → 收敛检查的外层循环
- 限制/延拓算子：full-weighting 限制 + 双线性延拓（2D；粗网格 = 全局 (nx/2, ny/2)）
- 粗层求解策略：粗网格数据布局与通信（gather/broadcast 复制式 或 分布式粗层——design 决定并记录，约束见 §5；GUIDE B1a 成本注记提示「粗层 gather/broadcast 与残差归约」）
- RBGS 平滑复用：新增轻量平滑入口（若干 sweep，不含残差扫描），**不改既有 solve/iterate 行为与位级一致测试**
- 收敛对照实验：制造解收敛 + 同 tol 下总迭代数（平滑步计）与两层理论定性一致，mg2 vs 纯 RBGS/Jacobi 迭代数对比入 PERFORMANCE（GUIDE B1a 验收①②口径）
- 新增粗层通信的 np=1/4 正确性测试（GUIDE B1a 验收④）
- 文档同步：README 求解器表/路线图、GUIDE B1a 勾选、AGENT_SPEC 类层级、PERFORMANCE 新小节

**Out of Scope：**
- B1b 完整 V-cycle / 多层粗化 / MG-CG 预条件子（后续 AR，复用本 AR 的限制/延拓）
- 3D MG（B1a 限 2D；`--nz > 1` 时 mg2 明示不支持——HYPOS_ERROR，不静默错误结果）
- Neumann 边界 × mg2（rediscretization 粗算子在 Neumann 下奇异——常数零空间，粗层 CG 停滞/发散风险；本 AR 仅支持 Dirichlet，`--bc neumann --solver mg2` → HYPOS_ERROR 明示，诚实优先于静默兼容；齐次 Neumann 校正边界 + 零空间投影留 B1b 及以后）
- 半粗化/代数 MG/自适应平滑策略
- 更改 RBGS/CG/Jacobi 既有求解语义

## 3. 功能需求

### 3.1 R1 两层校正求解器

**描述：** 新增 `TwoLevelMGSolver : PoissonSolver`，实现标准两层校正格式的外层循环。

**触发条件：** `--solver mg2`（CLI 新值；既有 jacobi/red_black_gs/cg 不受影响）。

**期望行为：**
- 外层循环：①ν1 次 RBGS 预平滑（细网格）②计算残差场 r_h = b - A·u ③限制 r_H = R·r_h ④粗网格上 CG 精解 A_H·e_H = r_H（紧 tol）⑤延拓校正 u += P·e_H ⑥ν2 次 RBGS 后平滑 ⑦全局真残差检查（复用 globalTrueResidual，含归约）——直至 tol 或 maxIter
- `iterate()` 语义：执行一个完整两层循环步并返回该步后全局真残差（满足 PoissonSolver 接口契约）；`lastResidual()`、进度回调、residualCheckInterval（外层步计）沿接口惯例
- ν1/ν2 默认值与可调性由 design 决定并记录
- 收敛判据与既有求解器同口径：全局 L2 真残差 ≤ tolerance

**异常处理：** MPI/粗层路径不抛异常（AGENT_SPEC）；不支持配置（3D、`--bc neumann`、奇数全局维、粗网格过小）→ HYPOS_ERROR 明示退出，不产出失实结果。

**验收标准：**
- Given 制造解（已知解析解的 rhs），When `--solver mg2` 运行至 tol，Then 收敛且解误差随网格加密二阶缩减（有限差分二阶精度佐证整体格式正确；测试落点见 R5 二阶收敛阶测试）
- Given 同 tol（如 1e-6）、256² 制造解负载，When mg2 与纯 RBGS 各自求解，Then mg2 总平滑步数严格少于纯 RBGS 迭代数（可判定底线），实际倍数如实记录不预设（记录式判定 + 两层理论定性结论文字为判据，GUIDE 验收①）
- Given 任意既有负载，When `--solver rbgs/jacobi/cg`，Then 行为零变化（回归）
- Given `--bc neumann --solver mg2`（及 3D、奇数维），When 启动，Then HYPOS_ERROR 且不产出结果

### 3.2 R2 限制/延拓算子与粗网格生成

**描述：** 2D full-weighting 限制 R 与双线性延拓 P，及粗网格 (nx/2, ny/2) 的生成。

**触发条件：** mg2 外层循环的②③⑤步。

**期望行为：**
- 粗网格全局规模 (nx/2, ny/2)（B1a 要求 nx, ny 偶数，奇数 → HYPOS_ERROR；故无取整歧义）
- full-weighting：r_H(I,J) = Σ w·r_h(2I+δi, 2J+δj)（9 点加权 1/16 系数），粗物理边界按 Dirichlet 0（校正量齐次边界）
- 双线性延拓：细点由相邻 4 粗点加权，细物理边界校正为 0（齐次）
- 子域边界处限制/延拓所需的邻居细网格数据经现有 HaloExchanger（任意 comm-mode）获得——粗层通信的正确性不依赖具体 exchanger 实现（验收④的 np=1/4 测试覆盖）
- R/P 为纯内核（无 MPI），可单测：与已知系数场逐点对照；P∘R 在粗点上的投影恒等性（r_H 经 P 回到细网格在粗点位置还原 r_H 值）——具体断言 design 细化
- 粗网格上的算子 A_H：与细网格同模板 rediscretization（2D 五点，D=4）——CG 直接在粗 Subgrid 上复用即得
- 粗网格规模下界：过小（<4² 量级）→ HYPOS_ERROR；具体阈值授权 design 决定并记录

**异常处理：** 纯内核无异常路径；越界/尺寸不符为调用方契约（HYPOS_ASSERT）。

**验收标准：**
- Given 已知系数场，When 限制/延拓，Then 逐点值与手算/解析期望一致（单测）
- Given np=1 与 np=4 非均匀切分同一全局残差场，When 各自执行限制→粗解→延拓，Then 细网格校正场在两种布局下逐元素一致（粗层通信正确性，GUIDE 验收④；串行参考实现对照）
- Given nx 或 ny 为奇数（或粗网格低于下界），When mg2 初始化，Then HYPOS_ERROR 且不产出结果

### 3.3 R3 RBGS 平滑复用（不改既有行为）

**描述：** 为 RedBlackGSSolver 增加轻量平滑入口（n 个 sweep，不含残差扫描），供 MG 平滑步调用。

**触发条件：** mg2 外层循环①⑥步。

**期望行为：**
- 新入口仅组合既有私有 sweep/iterateCore 逻辑（红黑各半 sweep + halo 交换 + applyPhysicalBoundary），**不复制粘改内核**
- 既有 `solve`/`iterate`/位级一致测试（OverlapTest 等）路径零变化——新增方法不触碰既有代码路径（验收③）

**异常处理：** 沿 RBGS 既有契约。

**验收标准：**
- Given 既有 RBGS 全部测试（含位级一致守护），When 本 AR 落地，Then 不变绿转（回归）
- Given 同一初始场，When 新平滑入口执行 n 步 vs 既有 iterate() n 次（扣除残差扫描），Then 更新结果位级一致（sweep 序完全相同）

### 3.4 R4 性能对照入档

**描述：** mg2 vs 纯 RBGS/Jacobi 的同 tol 迭代数与耗时对比，如实入 PERFORMANCE。

**期望行为：**
- 256² 制造解负载、np=1/4、OMP=1、≥3 次中位；记录：外层步数、总平滑步数（ν 加权计）、CG 粗层迭代数、墙钟时间、真残差终值
- 观测机制（步数/迭代数计数与 profiler 区名或 report 字段承载）由 design 决定，沿 AR004 区名惯例（如 `mg2_cycle`）记录
- 记录式判定（沿 D6 惯例）：不预设倍数门槛，实测如实；理论定性印证（步数不随网格加密而平方增长的可扩展性视角）可选加测 512² 佐证 O(N) 叙事——如加测则数据一并入档

**验收标准：**
- Given PERFORMANCE 新小节，Then 含 mg2 vs baseline 迭代数与耗时两组数据及口径说明，无失实宣称

### 3.5 R5 测试矩阵

**描述：** 单测 + np 集成 + e2e 全覆盖（GUIDE B1a 验收①③④）。

**期望行为：**
- 单测（单进程）：R/P 已知系数场对照、P∘R 投影恒等、粗网格生成尺寸/offset 正确性、奇数维防御、**二阶收敛阶测试**（R1 验收①落点：制造解在 ≥2 个网格尺寸如 32²/64² 下求解，解 L2 误差比 ≈4，容差明确——佐证整体格式正确性）
- np 集成：mg2 制造解收敛 np=1/4（含非均匀切分）、粗层通信跨布局一致性（R2 验收②）、RBGS 平滑入口位级一致
- e2e：直跑 hypos `--solver mg2`（沿 collective_mpi/datatype_e2e 模式）
- ctest 条目沿 AR006 惯例（np>1 带 HYPOS_EXPECT_NP 探针守卫、专用点名 filter）

**验收标准：**
- Given 全量 ctest（35 + 新增条目），When 双模式运行，Then 全绿且 ASan/UBSan 零报告

### 3.6 R6 文档同步

**描述：** 文档与实现零缺口（沿 AR002-AR006 惯例）。

**期望行为：** README 求解器表加 mg2 + 路线图无失实（B1b 仍列未实现）；GUIDE §5 ★4/B1a 勾选（附 AR 编号与证据）、「算法复杂度天花板」诊断更新；AGENT_SPEC 求解器类层级；`--help` solver 说明与实际一致，**且 `--residual-check-interval` help 文本与 solver.hpp:35-38 接口注释同步列 mg2**（stationary-iteration 语义——外层步计口径）；PERFORMANCE §14（R4）。

**验收标准：**
- Given 文档走查 + `--help` 实跑，Then 无失实宣称、帮助与行为一致

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 正确性 | 制造解收敛 | 二阶精度佐证；np1/np4 布局一致 |
| 兼容性 | 既有行为 | 既有三求解器与 overlap 位级守护零变化；PoissonSolver 接口签名不动（新增方法允许） |
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 资源 | MPI 对象 | 新增通信/缓冲生命周期完整（沿 exchanger 惯例） |
| 构建 | 双模式 | 0 新增警告 |
| 性能 | 数据如实 | 记录式判定，反直觉结果如实入档（沿 D6） |
| 回归 | 全量测试 | Debug+Release ctest 全绿 |
| 复杂度 | 预算 | **严控**：目标净新增 ≤ ~800 行（GUIDE §6.4；AR006 超支教训——测试模板化复用，断言循环抽 helper） |

## 5. 约束与假设

**约束：**
- 延续 AGENT_SPEC 全部约定（row-major、MPI 路径不抛异常、snake_case、无外部依赖、C++17+MPI）
- B1a 限 2D（nz=1）与偶数全局维；违反时 HYPOS_ERROR 明示（诚实优先于静默兼容）
- 粗层策略二选一由 design 决定并记录：A) 复制式粗层（残差 Allgather 到全体、各 rank 冗余串行 CG、无粗层 halo 交换——GUIDE 成本注记所指）；B) 分布式粗层（粗 Subgrid + 独立 exchanger，B1b 泛化更顺）。决策依据：B1a 正确性优先 + B1b 复用性 + 行数预算
- ν1/ν2、粗层 CG tol 等参数默认值由 design 决定；CLI 可调项最小化（避免 AR 膨胀）
- RBGS 平滑入口不得复制粘改既有内核（复用私有方法，必要时最小可见性调整——design 记录授权）
- 测试环境：WSL 单机 OpenMPI；np>1 条目带 HYPOS_EXPECT_NP 探针守卫

**假设：**
- 制造解 rhs 由测试内构造（沿既有测试 expectedValue 编码惯例），不新增 CLI 项
- 粗网格 CG 收敛到紧 tol 的迭代数在 128² 规模可忽略不计（相对外层步；实测验证）
- 512² 可扩展性佐证测试为可选（R4），如单次运行时间超 60s 则仅 256² 入档并注明

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| 两层校正格式 | 细网格平滑 + 粗网格校正的 MG 最小完整单元：u ← S^ν2(P·A_H⁻¹·R·(b-A·S^ν1(u)) + S^ν1(u)) |
| full-weighting 限制 | 9 点加权平均（角 1/16、边 2/16、心 4/16）把细残差压到粗网格 |
| 双线性延拓 | 细点由相邻 4 粗点按距离权重的插值 |
| rediscretization | 粗算子直接在粗网格上离散同模板（非 Galerkin），Poission 问题标准做法 |
| 平滑步 | 一次 RBGS 全迭代（红+黑各一 half-sweep）；「总平滑步数」= 外层步数 × (ν1+ν2)，MG 迭代数的标准计数口径 |
| 复制式粗层 | 粗问题在所有 rank 上冗余求解（Allgather 残差 + 串行 CG + 各自延拓），无粗层点对点通信 |
