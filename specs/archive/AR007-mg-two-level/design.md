# [AR007] 详细设计文档

| 字段 | 内容 |
|------|------|
| AR 编号 | AR007 |
| AR 主题 | mg-two-level（B1a 两层几何多重网格校正格式） |
| 关联 srs | specs/changes/AR007-mg-two-level/srs.md（R1-R6） |
| 关联指南 | docs/OPTIMIZATION_GUIDE.md §5 ★4 / B1a（:117-124，验收①③④；②为 B1b 遗留） |
| 日期 | 2026-10-09 |
| 状态 | 待门控（第 3 轮） |

## 1. 设计目标（从 srs 追溯）

| srs 需求 | 设计落点 | 验收锚 |
|------|------|------|
| R1 两层校正求解器 mg2 | FP4 外层循环 + FP5 CLI | 制造解收敛、步数严格少于纯 RBGS、既有求解器零变化 |
| R2 限制/延拓与粗网格生成 | FP1 算子内核 + FP2 粗层布局（含对角角交换） | 已知系数场对照、np1/np4 非均匀切分跨布局一致、奇数维防御 |
| R3 RBGS 平滑复用 | FP3 smooth 入口 + 残差场内核 | 位级一致、既有测试不变绿转 |
| R4 性能对照入档 | FP6 bench | PERFORMANCE §14 记录式 |
| R5 测试矩阵 | FP5 测试全表（§6） | 全绿 + ASan/UBSan 零报告 |
| R6 文档同步 | FP6 文档清单 | 零缺口走查 |

## 2. 功能点分解

| 序号 | 功能点 | 描述 |
|------|------|------|
| FP1 | 限制/延拓算子内核 | 2D 顶点重合 full-weighting 限制 R 与线性延拓 P（纯内核，无 MPI，布局无关）+ 粗胞所有权区间纯函数 |
| FP2 | 粗层布局与复制式粗解 | 所有权区间、对角角交换、Allgatherv 汇聚协议、MPI_COMM_SELF 复制式粗网格上 CGSolver 原样复用 |
| FP3 | 平滑入口与残差场 | RedBlackGSSolver::smooth（组合既有 iterateCore）+ residualFieldLocal（场版残差） |
| FP4 | TwoLevelMGSolver | 外层循环（预平滑→残差→限制→粗解→延拓→后平滑→收敛检查）、参数默认、防御四则、计数观测 |
| FP5 | CLI 接线与测试矩阵 | --solver mg2 + 前置校验 + ctest 条目 + e2e |
| FP6 | bench 与文档 | bench_ar007.sh（mg2 vs rbgs/jacobi 迭代数+耗时）→ PERFORMANCE §14 + 四文档同步 |

## 3. 影响范围（文件清单）

**修改：**
| 文件 | 变更 |
|------|------|
| src/solver/solver.hpp | +`TwoLevelMGSolver` 类声明；+`RedBlackGSSolver::smooth` 公开方法声明；residual-check-interval 注释列 mg2（R6/N4） |
| src/solver/red_black_gs_solver.cpp | +smooth 实现（n 次 iterateCore，~5 行；既有函数零变化） |
| src/solver/residual.hpp / residual.cpp | +`residualFieldLocal`（场版残差内核 + 契约注释） |
| src/main.cpp | --solver mg2 分支；前置校验（neumann/奇数维/3D → HYPOS_ERROR + return 1）；help 文本两处 |
| CMakeLists.txt | solver 源表 + mg_two_level_solver.cpp；**test_hypos 源表 + tests/test_mg_two_level.cpp**；3 个新测试条目（mg2_np1/mg2_np4/mg2_e2e） |

**新增：**
| 文件 | 内容 |
|------|------|
| src/solver/mg_operators.hpp | 限制/延拓/所有权区间纯函数声明 + exchangeCorners 声明（后者是 MPI 函数，**不适用 §4.1「无 MPI」内核纪律**，归此头仅为供求解器与测试共用） |
| src/solver/mg_two_level_solver.cpp | 算子实现 + 对角角交换 + TwoLevelMGSolver 实现 |
| tests/test_mg_two_level.cpp | §6 测试全表 |
| scripts/bench_ar007.sh | R4 bench |

**不动：** p2p/datatype/collective exchanger（角交换是 mg 模块自有通信，不改 exchanger）、jacobi_solver.cpp、cg_solver.cpp（原样复用）、RBGS 既有 solve/iterate/sweep、IO 栈。

## 4. 实现设计

### 4.0 方案对比：粗层策略（决策 D1）

**方案 A 复制式粗层（选定）：** 细残差场经现有 exchanger 面交换 + mg 自有对角角交换 → 各 rank 本地限制出自己负责的粗矩形片 → `MPI_Allgatherv` 汇聚成完整粗网格（每 rank 一份）→ 各 rank 在 **MPI_COMM_SELF** 粗 Subgrid 上**冗余串行** CG 精解（全 PROC_NULL 邻居 + PointToPointExchanger，CGSolver 原样复用）→ 各 rank 本地线性延拓（全员持全粗网格，**延拓零通信**）。
- 优点：新 MPI 原语仅 1× Allgatherv + 4× Sendrecv(1 Real)/cycle；无粗层 halo 生命周期管理；R/P 纯内核与布局解耦（B1b 复用保证）；粗解各 rank 严格串行同输入同代码 → **位级一致，无跨 rank allreduce 假设**；与 GUIDE B1a 成本注记「粗层 gather/broadcast 与残差归约」吻合。
- 缺点：粗解计算冗余（B1a 规模 128² 串行 CG ~毫秒级，可忽略）；粗内存每 rank O(nxH·nyH)（256² 粗层 128KB，可忽略）。

**方案 B 分布式粗层（否决）：** 粗 Subgrid 按同 cart 维度再分 + 独立 exchanger + 粗层 halo 交换 + 延拓前粗通信。
- 否决理由：额外 ~200+ 行，B1a 正确性优先目标下零收益；128² 粗网格切 4 份后每份 64² 反而放大通信占比；B1b 需要的是 R/P 内核复用而非粗层布局复用（已由方案 A 保证）。

### 4.1 FP1 限制/延拓算子内核（决策 D3：顶点重合约定）

**坐标系（顶点重合，全局 0 基）：** 粗值 r_H(I,J) 与细点 (2I, 2J) 位置重合；粗胞 I 关联细点对 {2I, 2I+1}。物理界外恒取 0（校正量齐次 Dirichlet）。

**full-weighting 限制（R，中心在细点 (2I,2J)）：**
```
r_H(I,J) = (1/16)·[ 4·r(2I,2J)
          + 2·(r(2I±1,2J) + r(2I,2J±1))
          + 1·(r(2I±1,2J±1)) ]     （物理界外项取 0）
```
- 内部粗胞权和 (4+8+4)/16 = 1（保常，**仅内部**——物理边界粗胞因界外取 0 而失常，与齐次校正自洽，断言范围见 §6）。
- 输入契约：细残差场 `r` 已完成**面** halo 交换 + **对角角**交换（§4.2）+ 物理面 applyPhysicalBoundary(r, 0.0)。
- 输出：写入调用方粗片缓冲（本 rank 所有权区间，见 FP2）。

**线性延拓（P，顶点重合）：**
```
细点 2I   （偶）:  e_f = e_C(I)                      （权 1）
细点 2I+1 （奇）:  e_f = (e_C(I) + e_C(I+1)) / 2     （I+1 物理界外取 0）
```
二维 = 行权 ⊗ 列权（1、1/2、1/4 组合）。延拓结果累加进 `u += P·e_H`。
- **P∘R 恒等（严格）**：P(R r)(2I, 2J) = r_H(I,J)——粗值经 R 落位、经 P 原位取回，单测锚。
- **保常性范围**：偶细点全部保常；奇细点在双粗父胞均在内时保常；右/上物理边界的奇细点（i = nx−1，I+1 = nxH 界外）只得 1/2 权——**记录的不对称**：界外取 0 与齐次 Dirichlet 校正自洽，对校正格式收敛性为边界局部效应（任何「合理」插值均保两层格式收敛），断言范围见 §6。
- 输入契约：粗解 e_H 为复制式完整场；内核直接读粗 Subgrid 的 u 内部区，零通信。

**粗胞所有权区间（纯函数，可单测）：**
```cpp
// mg_operators.hpp —— coarse cell I belongs to the rank owning fine point 2I
Index coarseRangeBegin(Index offset, Index nLocal) noexcept;   // = ceil(offset/2)
Index coarseRangeCount(Index offset, Index nLocal) noexcept;   // = ceil((offset+nLocal)/2) - ceil(offset/2)
```
**空粗片可达**：nLocal=1 且 offset 奇（如 x 向 8 rank 切 nx=8）→ count=0；Allgatherv count=0 合法、解包跳过——注记 + 测试覆盖（§6）。

**内核纪律：** 无 MPI、无异常路径；尺寸/越界为调用方契约（HYPOS_ASSERT）；OpenMP 沿 j 行 parallel for（行内串行，无竞争）。

### 4.2 FP2 粗层布局、对角角交换与复制式粗解

**所有权区间（coarseRange，零通信）：** rank 拥有细全局列 [ox, ox+nxL)，负责粗列 [ceil(ox/2), ceil((ox+nxL)/2))，y 同理。区间在 x/y 各自连续 → 粗片恒为**矩形**（任意 cart 布局）。验算（门控第 1 轮确认）：nx=10(3,3,4)→{0,1}/{2}/{3,4}；nx=8(2,3,3)→{0}/{1,2}/{3}——无缝不重叠。覆盖性：sum(局部粗片)=nx/2 ✓（奇数全局维在防御中拒绝）。

**halo 覆盖（含对角，门控 B1 修复）：** 粗胞 I 的 FW 模板需要细列/行 {2I−1 .. 2I+1} × {2J−1 .. 2J+1}（9 点）。
- **轴向**（面 halo）：区间下端 2I ≥ ox → 2I−1 ≥ ox−1（halo ✓）；上端 2I+1 ≤ ox+nxL（halo ✓）。haloWidth=1 面向充分——由现有 exchanger 面交换提供。
- **对角**（角 halo）：当 ox、oy 均偶（典型 2×2 偶切分）时，粗胞 (I_min, J_min) 的模板角点 (2I_min−1, 2J_min−1) = **对角邻居的内点**，落在本 rank 的**角 halo 单元**——三套 exchanger 均只交换面且横向范围为内部带，**角 halo 永不填充**；applyPhysicalBoundary 仅覆盖 PROC_NULL 面。**修复：mg 自有对角角交换**（不改 exchanger；守卫次序——**面邻居判定在 Cart 推导之前**，门控 R2-2）：
```
exchangeCorners(subgrid, data)：          // mg_operators.hpp 公开自由函数（供测试直调）
  对 4 个对角方向 d ∈ {(−x,−y),(+x,−y),(−x,+y),(+x,+y)}：
    先查面邻居：若 x 侧或 y 侧邻居为 MPI_PROC_NULL（物理边界）→ 跳过该角
    （该角 halo 已由 applyPhysicalBoundary(r, 0.0) 的 Dirichlet 全带填充置 0——物理界外=0，正确；
     且物理边界侧不存在对角 rank，Cart_rank 越界坐标非法，必须先守卫）
    否则：对角 rank = MPI_Cart_coords(本 rank) + d → MPI_Cart_rank（subgrid.comm() 是 cart comm）
    Sendrecv 1 Real：发本 rank 对应内角点值，收进对应角 halo 单元
```
每 cycle 新通信 = 4× Sendrecv(1 Real)（物理侧守卫跳过）+ 1× Allgatherv。tag 约定：沿 AR006 惯例独立 tag 段（避免与 exchanger tag 冲突；角交换用专用 tag 常量）。

**汇聚协议（每 cycle）：**
1. 首次 solve 时 `MPI_Allgather` 各 rank 粗矩形元数据 (cnxL, cnyL, coffX, coffY)，推导全局 nxH/nyH 与 Allgatherv counts/displs，缓存（契约：缓存键=细局部尺寸+粗全局尺寸，subgrid 尺寸变更时重推导——幂等 ensureInitialized）。
2. 限制内核写入本地打包缓冲（行优先 cnxL×cnyL；count 可为 0）。
3. `MPI_Allgatherv` → 全局粗缓冲 → 各 rank 按 (coffX, coffY) 把**全部**矩形解包进自己的复制式粗 Subgrid rhs 内部区（每 rank 持全量）。

**复制式粗 CG（决策 D-rev：MPI_COMM_SELF，门控 B3 修复）：** 粗 Subgrid 以 **MPI_COMM_SELF** 构造：(nxH, nyH, 1, halo=1, MPI_COMM_SELF)，全 PROC_NULL 邻居（默认），Dirichlet BC（默认）。CGSolver::dotGlobal 的 Allreduce 在 SELF 上**严格恒等**——无 ×np 缩放、无 √np 容差偏移；各 rank 对位级相同的输入执行相同串行代码 → e_H **位级一致**，不依赖任何跨 rank allreduce 一致性假设（D5 裕度由「冗余保险」降级为「纯保险」）。PointToPointExchanger 初始化于粗 Subgrid（全 PROC_NULL → 无消息）。

**符号链（门控验算通过）：** 本项目系统 `A u = b`，`A = D·I − S`，`b = −rhs`。残差场 `ρ = b − A·u = Σ(u_nb) − D·u − rhs = −(trueResidualSquaredLocal 内核式)`（范数等价）。粗层校正方程 `A_H·e_H = ρ_H`；CGSolver 内部取 `b = −rhs_buf` → 粗 Subgrid rhs 填 **−ρ_H**。解出后 `u += P·e_H`。

### 4.3 FP3 平滑入口与残差场

**RedBlackGSSolver::smooth（验收③核心）：**
```cpp
// solver.hpp（RedBlackGSSolver 公开区）
void smooth(Subgrid& subgrid, HaloExchanger& exchanger, Index sweeps) const;
// red_black_gs_solver.cpp
void RedBlackGSSolver::smooth(Subgrid& s, HaloExchanger& e, Index n) const {
    for (Index k = 0; k < n; ++k) iterateCore(s, e);   // 组合私有内核，零复制
}
```
既有 solve/iterate/sweep/iterateCore 逐行不动；新方法仅新增调用路径。sweep 序 = iterate()（iterateCore + 残差扫描，扫描不改 u）→ smooth(n) 与 iterate()×n 的 u 终态位级一致（单测锚；该性质与分解无关，np=1 验证即覆盖 srs R5 的 np 集成条目——注记闭合追溯）。

**residualFieldLocal（residual.hpp/cpp）：**
```cpp
void residualFieldLocal(const Subgrid& subgrid, Real* residual) noexcept;
// interior: residual[idx] = u[idx-1]+u[idx+1]+u[idx-nxT]+u[idx+nxT] - 4·u[idx] - rhs[idx]
// （2D；halo 单元不写；3D 走 D=6 分支沿既有内核惯例）
```
与 trueResidualSquaredLocal 同一算式取场不取范数。调用次序契约：u 的 halo 已同步（iterateCore 末尾 applyPhysicalBoundary 保证）→ 算 ρ → `exchanger.exchange(sg, ρ)` + `exchangeCorners(sg, ρ)` + `sg.applyPhysicalBoundary(ρ, 0.0)` → 限制可读 ρ 的面 halo、角 halo 与物理面。

### 4.4 FP4 TwoLevelMGSolver 外层循环

**类声明（solver.hpp，沿 CGSolver 成员惯例）：**
```cpp
class TwoLevelMGSolver : public PoissonSolver {
public:
    Index solve(Subgrid&, HaloExchanger&, Index maxIter, Real tolerance) override;
    Real  iterate(Subgrid&, HaloExchanger&) override;      // 一个完整 cycle + 真残差
    std::string name() const override { return "mg2"; }
    Real  lastResidual() const noexcept override;
    // 平滑步数：预/后平滑各 ν=2（决策 D2；不加 CLI 调节，防 AR 膨胀）
    static constexpr Index kPreSmoothSweeps = 2;
    static constexpr Index kPostSmoothSweeps = 2;
private:
    bool ensureInitialized(Subgrid&, HaloExchanger&);      // 防御 + 粗层元数据（幂等，尺寸变更重推导）
    void cycleOnce(Subgrid&, HaloExchanger&, Real tolerance);   // ①-⑧ 一个 cycle（不含检查）
    RedBlackGSSolver smoother_;
    CGSolver coarseSolver_;
    std::unique_ptr<Subgrid> coarse_;                     // 复制式粗网格（MPI_COMM_SELF）
    std::unique_ptr<PointToPointExchanger> coarseExchanger_;
    AlignedBuffer<Real> fineResidual_, coarsePack_, coarseAll_;
    std::vector<Index> coarseCounts_, coarseDispls_, rectMeta_;  // rectMeta: allgather 的 4-int 元数据
    Index nxH_ = 0, nyH_ = 0;
    Index coarseIters_ = 0;      // 累计粗层 CG 迭代数（日志/D6 观测）
    Real coarseTolerance_ = 0.0; // cycleOnce 粗层 tol 来源：solve() 传入 tolerance 缓存于此；
                                 // iterate() 独立调用（未经 solve）沿 CG stateReady_ 惯例用默认 1e-10
    bool valid_ = false;
    Real lastResidual_ = 0.0;
};
```

**solve() 伪码（计数口径：maxIter 与返回值均为外层 cycle 数；门控 B4 修复——先清场后填 rhs）：**
```
if (!ensureInitialized) return 0;            // HYPOS_ERROR 已记，u 不动
cycle = 0;
while (cycle < maxIter) {
    smoother_.smooth(sg, ex, kPreSmoothSweeps);      // ① 预平滑
    residualFieldLocal(sg, fineResidual_);            // ② ρ = b - A·u
    ex.exchange(sg, fineResidual_);                   // ③ halo：面 + 角 + 物理
    exchangeCorners(sg, fineResidual_);
    sg.applyPhysicalBoundary(fineResidual_, 0.0);
    restrictLocal(→ coarsePack_);                     // ④ 本地粗片
    MPI_Allgatherv(coarsePack_ → coarseAll_);
    coarse_->zeroInitialize();                        // 清 u/uNext/rhs 三场
    unpack(→ coarse_->rhs()) 后取负（rhs = -ρ_H）;     // ⑤ 填粗 rhs（在清场之后！）
    coarseIters += coarseSolver_.solve(*coarse_, *coarseEx_, nxH*nyH, tolerance);  // ⑥ 粗解
    prolongateCorrection(sg, *coarse_);               // ⑦ u += P·e_H（零通信）
    smoother_.smooth(sg, ex, kPostSmoothSweeps);      // ⑧ 后平滑
    ++cycle;                                          // 唯一自增点（门控 R2-4）
    if (cycle % residualCheckInterval_ == 0) {        // ⑨ 真残差（外层步计）
        lastResidual_ = globalTrueResidual(sg, ex);
        notifyProgress(cycle);
        if (lastResidual_ < tolerance) break;
    } else notifyProgress(cycle);
}
if (cycle > 0) lastResidual_ = globalTrueResidual(sg, ex);   // 出口确认扫描（沿 RBGS 契约）
HYPOS_INFO("TwoLevelMG finished in " << cycle << " cycles ("
           << cycle*4 << " smoothing sweeps, "
           << coarseIters << " coarse CG iterations), time = ...");
```
profiler 区名（沿 AR004 惯例）：`mg2_cycle`、`mg2_restrict`、`mg2_coarse_solve`、`mg2_prolongate`；平滑/交换/残差复用既有区名。

**防御四则（ensureInitialized，全局一致分支，MPI 路径不抛异常；门控 B5 修复）：**
1. **全局** nz > 1 → 3D 不支持（`MPI_Allreduce(MAX)` 推导 offsetZ+nzLocal 的最大值——`nzLocal>1` 单独判不完备：z 向被分解时全员 nzLocal==1 而全局 3D）；
2. 全局 nx/ny 奇数（同法 Allreduce(MAX) 推导）→ 粗化歧义；
3. `subgrid.boundaryCondition() == Neumann` → 粗算子奇异（srs B1）；
4. `nxH < 4 || nyH < 4`（粗网格 < 4²，srs N3 下界取 4）→ 粗解规模无意义。
处理契约（决策 D4）：HYPOS_ERROR 是**日志非终止**（logger.hpp:52）→ 各 rank 记 ERROR 后 `return 0` 迭代、**u 不被修改**、lastResidual 保持 0——可测契约，杜绝静默失实结果。四项条件均全局一致（全局维 by Allreduce、BC 全员相同）→ 无挂死风险。main.cpp 在求解器构造前用全局量直查另做同四项前置校验（HYPOS_ERROR + return 1，沿「Unknown solver」惯例；双保险）。

**iterate()：** 一个 cycle（①-⑧）+ `globalTrueResidual` 返回（沿 RBGS「单次调用契约」）。

### 4.5 FP5 CLI 接线与测试矩阵

**main.cpp：**
- 分支 `else if (solverName == "mg2") solver = std::make_unique<TwoLevelMGSolver>();`
- 前置校验（分支后立即）：`bcType=="neumann"` / `grid.nx%2||grid.ny%2` / `grid.nz>1` → HYPOS_ERROR + return 1（全局量直查，不依赖求解器内 allreduce 路径）。
- help 文本：`--solver` 行加 mg2；`--residual-check-interval` WARN/注释同步「jacobi、red_black_gs 与 mg2（外层步计）」（R6/N4，main.cpp:208 与 solver.hpp:35-38 两处）。
- `solverName != "cg"` 的 setResidualCheckInterval 分支已覆盖 mg2 ✓（不动）。
- **flopsPerCell 口径注记（门控建议 3）**：main.cpp:318-322 的 Est. FLOP/s 按「每迭代 8 flops/cell」折算，mg2 一个 cycle 含 4 平滑步+限制+粗 CG+延拓，实际算力高于折算值 → 该指标对 mg2 **系统性偏低**。处理：不改折算逻辑（超出本 AR 范围），在 PERFORMANCE §14 口径说明与 README 已知限制中注记，避免失实宣称。

**ctest 条目（沿 AR006 惯例）：**
```cmake
add_test(NAME mg2_np1  ... --gtest_filter=TwoLevelMGUnitTest.*:MpiEnvTest.*)          # 1 进程
add_test(NAME mg2_np4  ... --gtest_filter=TwoLevelMGMpiTest.*:MpiEnvTest.*)           # 4 进程
add_test(NAME mg2_e2e  ... $<TARGET_FILE:hypos> --solver mg2 --nx 32 --ny 32 --max-iter 20 ...)
# ENV：np4/e2e 条 "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4"，TIMEOUT 900（多迭代余量沿 solver_mpi）
```

### 4.6 FP6 bench 与文档

**bench_ar007.sh（沿 bench_ar006.sh 骨架）：** 256²、tol=1e-6、OMP=1、np∈{1,4}、solver∈{mg2, red_black_gs, jacobi}，各 ≥3 次取中位；从 stdout 日志提取：外层 cycle 数、总平滑步数（cycle×4）、粗层 CG 迭代数、墙钟、真残差终值。断言性检查（非门槛）：mg2 总平滑步数 < rbgs 迭代数（R1 验收②底线）。512² 可扩展性佐证为可选（单次 >60s 则跳过并注明）。

**文档清单：** PERFORMANCE §14（数据+口径+记录式结论+flopsPerCell 偏低注记）；README 求解器表+路线图（B1a 完成、B1b 仍列未实现）；GUIDE ★4/B1a 勾选（附 AR007+证据）+ §3「算法复杂度天花板」更新；AGENT_SPEC 求解器类层级（+TwoLevelMGSolver，标注平滑/粗解复用关系）；--help 与实际一致。

## 5. 接口描述

| 接口 | 签名 | 语义 |
|------|------|------|
| 新增 | `void RedBlackGSSolver::smooth(Subgrid&, HaloExchanger&, Index sweeps) const` | n 次 iterateCore（红黑半步+halo+物理边界），不含残差扫描 |
| 新增 | `void residualFieldLocal(const Subgrid&, Real* residual) noexcept` | interior 写入 ρ = b − A·u；halo 不写；无通信 |
| 新增 | `Index coarseRangeBegin(Index offset, Index nLocal) noexcept` / `Index coarseRangeCount(Index offset, Index nLocal) noexcept`（mg_operators.hpp） | 粗胞所有权区间纯函数；count 可为 0（空粗片） |
| 新增 | `void restrictResidual(const Subgrid& fine, const Real* fineResidual, Index coarseX0, Index coarseY0, Index cnx, Index cny, Real* coarsePacked)` | FW 限制到粗矩形 [coarseX0,coarseX0+cnx)×[coarseY0,coarseY0+cny)，行优先打包 |
| 新增 | `void prolongateCorrection(Subgrid& fine, const Subgrid& coarse)` | 复制式粗解 e_H 顶点重合线性延拓并**累加**进 fine.u() 内部区 |
| 新增 | `void exchangeCorners(Subgrid& sg, Real* data)`（mg_operators.hpp 公开自由函数，供求解器与测试共用） | 4 对角 Sendrecv(1 Real)，**先守卫面邻居 PROC_NULL 再 Cart 推导**；物理侧角由 Dirichlet 全带填充置 0 |
| 新增 | `class TwoLevelMGSolver : PoissonSolver` | solve/iterate/name/lastResidual；name()="mg2" |
| 修改 | main.cpp solver 分派 + 前置校验 + help | §4.5 |
| 不变 | PoissonSolver 纯虚接口、CGSolver、JacobiSolver、RBGS 既有方法、三 exchanger | 兼容性锚 |

## 6. 测试设计（全表）与交付物映射

tests/test_mg_two_level.cpp；`TwoLevelMGUnitTest`（np=1）+ `TwoLevelMGMpiTest`（np=4，GTEST_SKIP 守卫沿 SolverMpiTest 惯例）。

| 用例 | np | 断言要点 | srs 验收锚 |
|------|----|---------|-----------|
| CoarseRangePureFunction | 1 | coarseRangeBegin/Count 对照手算：偶/奇 offset、奇 nLocal、(1,1)→空片；(3,3)→1 胞；(0,4)→2 胞；无缝性 Σcount=nx/2 | R2② |
| RestrictionKnownCoefficients | 1 | 已知 r 场（常数、线性）逐点对照解析期望；**内部粗胞保常**；物理边界行/列/角粗胞按界外取 0 权手算（无常性断言） | R2① |
| ProlongationKnownCoefficients | 1 | 常数粗场：偶细点与内部奇细点保常、右/上边界奇细点 1/2 权（手算锚）；单位脉冲粗场 → 1/1/2/1/4 权形 | R2① |
| ProlongRestrictProjectionIdentity | 1 | 任意 r 场：P(R r) 在细点 (2I,2J) 逐点 == r_H(I,J)（位级） | R2① |
| ResidualFieldMatchesTrueResidualScan | 1 | ‖residualFieldLocal‖₂ == sqrt(trueResidualSquaredLocal)（同口径恒等，1e-15 相对） | R3/R1 |
| SmoothMatchesIterateBitExact | 1 | 同初始场：smooth(n) 终态 u == iterate()×n 终态 u（位级；性质与分解无关，np1 覆盖 R5 np 集成条目） | R3③ |
| DefenseRejectsUnsupportedConfigs | 1 | 3D 子域 / 奇数全局维（offset 布局模拟）/ Neumann BC / 粗网格 <4²：solve 返回 0、u 逐位不变、lastResidual()==0 | R1④ |
| ManufacturedSineConvergence | 1+4 | 64² sine 制造解（沿 SolverMpiTest:42 模式）：solve 返回 < 上限、lastResidual ≤ tol、l2Error < 1e-3 | R1① |
| SecondOrderConvergence | 1 | 32² 与 64² 制造解、tol=1e-8（≪ 离散误差 ~1e-3，3 量级隔离）：l2Error(32)/l2Error(64) ∈ [3.5, 4.5] | R1①（阶验证） |
| CrossLayoutCorrectionConsistency | 4 | **真非均匀切分**（srs R2 验收②）：2×2、全局 nx=ny=34、**手动构造 x/y 向 {15,19} 分片**（offsets {0,15}；两片均奇——奇端所有权+全四角触发；注意 34=17+17 的默认均匀切分不是非均匀，故测试手工 setOffsets/setNeighbors 构造，沿 test_solver_mpi 手工布局惯例）；固定 rhs+u=0；np4 单 cycle 后 u vs 测试内串行参考（34² 全场冗余复算，沿 FixedIterationsMatchesSerialReference:154 模式）逐元素 \|diff\| ≤ 1e-12·max\|u\|（D5：COMM_SELF 复制式粗解理论位级，留浮点裕度） | R2②/验收④ |
| CoarseOwnershipPartitionsGlobalRange | 4 | {15,19} 非均匀布局下各 rank 粗矩形（[0,8) 与 [8,17)，手算锚）allgather 后无缝不重叠覆盖 [0,17)×[0,17)；矩形性；**含空片退化断言**（构造元数据 (1,1)→count 0） | R2② |
| DiagonalCornerExchangeFilled | 4 | 测试**自建整场缓冲**（角点填已知 rank 相关值）→ 先 `applyPhysicalBoundary(buf, 0.0)`（前置调用，沿 §4.2 调用次序契约）→ 直调公开 exchangeCorners 后：内边界角单元 == 对角邻居内角点值；物理侧角 == 0（Dirichlet 全带填充） | R2②（B1 锚） |
| mg2_e2e（ctest） | 4 | hypos --solver mg2 --nx 32 --ny 32 --max-iter 20 正常运行出文件 | R5 |
| 既有全量回归 | — | 35 条不变绿转（RBGS 位级守护/overlap 零变化为关键锚） | R1③/R3③ |

**NFR 映射：** ASan/UBSan 全绿（双模式 ctest；B1 角 halo 未初始化读对 ASan 不敏感——由 DiagonalCornerExchangeFilled + CrossLayout 直接断言兜底）；0 新增警告；行数预算见 §7 D8。

## 7. 决策记录

| # | 决策 | 理由 |
|---|------|------|
| D1 | 粗层策略=方案 A 复制式 | §4.0：GUIDE 成本注记吻合、行数预算、B1a 正确性优先；B1b 复用性由 R/P 布局无关内核保证 |
| D2 | ν1=ν2=2；粗层 CG tol=外层 tolerance、maxIter=nxH·nyH | 文献标准平滑数；粗层 tol 由外层真残差检查兜底（粗解不完美只降效率不失正确性）；maxIter 取 CG 理论上界（实际 ~O(百)） |
| D3 | **顶点重合约定**：R 中心=细 (2I,2J)，P 偶细点权 1、奇细点 1/2+1/2、界外 0；右/上边界奇细点 1/2 权为记录的不对称 | 门控 B2 修复（原胞心 P 与顶点 R 错位半胞、恒等式不可通过）：顶点重合使 P∘R 恒等**严格成立**（srs R2 验收锚）；界外取 0 与齐次 Dirichlet 校正自洽；边界局部效应对两层格式收敛性无碍 |
| D4 | 防御契约=HYPOS_ERROR 记日志 + 返回 0 迭代 + u 不动（可测）；main 前置校验双保险；3D 判定用全局 nz（Allreduce） | HYPOS_ERROR 非终止（logger.hpp:52）；杜绝静默失实；nzLocal>1 判定不完备（z 向分解时全员 nzLocal==1，门控 B5） |
| D5 | 跨布局一致性断言 1e-12 相对而非位级 | COMM_SELF 复制式粗解 + 冗余串行同输入同代码使理论位级成立；保留浮点裕度防脆弱测试（保险性断言） |
| D6 | 观测=INFO 日志（cycles/sweeps/coarse iters/time）+ profiler 区名；report 不加新字段 | 沿 RBGS/CG stdout 取数惯例；bench 从日志提取 |
| D7 | mg2 计入 residual-check-interval 支持列表，步计=外层 cycle | stationary-iteration 语义（外层步）；main:208 WARN 文本与 solver.hpp:35-38 注释同步（N4） |
| D8 | 行数预算预估 **~1100**（实现 ~520 含角交换 ~50、测试 ~530、脚本/cmake ~50）；超 §6.4 ~800 建议走论证路径 | 门控第 1 轮修正预估（原 850-950 偏乐观）：测试占比 ~48% 沿既有模板；B1a 验收④硬性要求非均匀切分跨布局+串行参考复算（~100+ 行）、阶验证与防御四则各 60-120 行；对角角交换为 B1 修复必需（~50 行）。论证随 T006 实测数入 tasks.md |
| D9 | 对角角交换为 mg 模块自有通信，不改三 exchanger | exchanger 面语义是既有位级守护的契约（验收③）；角需求是 MG 特有（1 Real × 4），独立实现面最小 |
| D10 | 粗 Subgrid 用 MPI_COMM_SELF（非 cart comm） | 门控 B3 修复：dotGlobal Allreduce 在 cart comm 上是 ×np 缩放（√np 有效容差偏移、np1/np4 粗 CG 停点不同）；SELF 上严格恒等，粗层 tol 语义（D2）精确成立 |

## 8. 风险与缓解

| 风险 | 缓解 |
|------|------|
| 对角角交换实现错（cart 坐标推导/角单元索引错位） | DiagonalCornerExchangeFilled np4 直测角值来源；CrossLayout 34² 非均匀（奇端+全四角）逐元素对照串行参考；专用 tag 隔离避免与 exchanger 冲突 |
| 32²/64² 收敛阶测试被迭代误差污染 | tol=1e-8 ≪ 离散误差 ~1e-3（3 个量级隔离）；比值区间 [3.5,4.5] 容差 |
| 粗片非矩形/覆盖漏洞/空粗片（布局边角） | CoarseRangePureFunction 手算表 + CoarseOwnershipPartitionsGlobalRange np4 直测；空片（count=0）在 Allgatherv 合法、解包跳过，测试含退化断言 |
| 粗层 CG 未收敛到 tol（效率损失） | 外层真残差检查兜底正确性；coarseIters 计数入日志，异常如实记录 |
| RBGS smooth 触碰既有路径（验收③回归） | smooth 仅新增 5 行调用 iterateCore；SmoothMatchesIterateBitExact + 既有 35 条回归双锚 |
| 奇数维/3D 检测在 solver 内依赖 allreduce（初始化早期通信） | ensureInitialized 首次调用完成全部 allgather/allreduce（全局一致分支，无挂死）；main 前置校验先行拦截正常路径 |
| flopsPerCell 折算对 mg2 系统性偏低（Est. FLOP/s 失真） | 不改折算逻辑（超范围）；PERFORMANCE §14 口径说明 + README 已知限制注记（无失实宣称底线） |
| ensureInitialized 缓存与 subgrid 尺寸不符（同实例换网格） | 缓存键=细局部尺寸+粗全局尺寸，不匹配即重推导（幂等） |
