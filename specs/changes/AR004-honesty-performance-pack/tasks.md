# [AR004] 任务跟踪

| 字段 | 内容 |
|------|------|
| AR 编号 | AR004 |
| 关联 srs.md | ./srs.md |
| 关联 design.md | ./design.md（待生成）|
| 创建日期 | 2026-10-08 |

## 任务列表

| ID | 任务描述 | 依赖 | 状态 | 备注 |
|----|---------|------|------|------|
| T001 | A2-Jacobi 真实残差 + interval：新增 src/solver/residual.hpp/cpp（trueResidualSquaredLocal）；iterate() 返回 denom×‖diff‖；solve() 判据换算 + 出口统一确认扫描（区名 residual_confirm）。测试先行：U1-U3、U10、B2 期望值更新 | - | passing | 2026-10-08 完成；新增 solver 级 interval 用例 I1s/I5s/E4s 一并落地；maxIter=0 出口跳过确认扫描（保基类契约） |
| T002 | A2-RBGS 真实残差 + interval：iterateCore/iterate 拆分（iterate 每调用必扫描）；solve() 每 k 步判定。测试先行：U4/U5、B5 | T001 | passing | 2026-10-08 完成；新增 R-I1s/R-I5s/R-E4s；rbgs_iteration 区名随 solve 重写落地（原 T007 项）；新增共享 globalTrueResidual（design §4.3.1 已补录） |
| T003 | A2-CLI `--residual-check-interval`：基类 setResidualCheckInterval（0 钳 1）；main 解析/校验（<1 退出码 1；CG 传入 WARN 忽略）；RunConfig +residualCheckInterval 透出。测试先行：U9、I1-I6 | T001 | passing | 2026-10-08 完成；I1-I6 以 ctest 条目 + cmake/VerifyAr004Cli.cmake 四模式（interval10/default/cg_warn/help）落地；CSV 列插入 overlap_comm 后（17 列，位置由 default 模式守护）；拒绝类 4 项 WILL_FAIL；CG 不调用 setter（WARN 语义一致）；Red 期修正：verify 条目需 HYPOS_TEST_ENV（LSan 噪声） |
| T004 | B2a+A4-CG：updatePInterior 抽取（消除 solve/iterate 重复）；初始化循环并行化；区名改 cg_iteration。测试先行：U6（golden 迭代数，先跑基线固化）、E2 | - | passing | 2026-10-08 完成；U6 实测修正设计前提：libgomp 归约按线程到达序合并 → OMP=4 跨运行非位级确定，tier=1 保位级 hash、tier=4 改迭代数+相对差 1e-11（噪声带宽实测 ≤1.4e-12）；golden 绑定 Debug+ASan（NDEBUG 跳过）；E2 以相邻双 inf rhs 确定性触发 pap=NaN 守卫 |
| T005 | A1 First-touch：zeroInitialize() 外层维度并行（全缓冲，三场单遍历）；applyDirichletBC 保持串行并文档记录。测试先行：U7、B4 | - | passing | 2026-10-08 完成（微型任务，主代理直做）；设计勘误：2D 子域缓冲实含 nzTotal≥3 个 k 平面，设计原文"2D 仅 j/i 全范围"会漏 2/3 缓冲——2D 分支补内层 k 循环、外层仍并行 j（保线程利用率），U7 全缓冲断言（halo=2 构造）守护该语义 |
| T006 | A3 拓扑一致性：SubgridInfo +cartComm（所有权移交）；partition 用 cart rank 查 coords；main 删除二次 Cart_create。测试先行：B3 | - | passing | 2026-10-08 完成；Red=编译失败（info.cartComm 不存在）；B3 设计勘误：原"Σ nxLocal==nx"全局求和在 2D tiling 下不变式错误（Σ=nx×dims[1]），修正为按行/列分组求和 + 链端点检查（走查 evidence ② 单次 Cart_create 随 Green 落实） |
| T007 | A4 Profiler：per-thread ThreadData 注册表（热路径零锁）+ 聚合；hpp 注释对齐；test_performance.cpp:83 区名断言复核（RBGS 区名已随 T002 落地）。测试先行：U8、U11 | T004 | pending | design D7；JSON 格式兼容 |
| T008 | 全量回归（Debug/Release + ASan；TSan 不可用则走查佐证）+ 基准（scripts/bench_ar004.sh：Jacobi ±5% 门槛、RBGS k=1/5/10 三档、CG 前后，≥3 次中位数）+ 文档同步（README/AGENT_SPEC/PERFORMANCE §11/GUIDE 勾选 G1/G3/G5/G6、P1/P4/P6） | T001-T007 | pending | SCALING_REPORT 口径注记 |

## 状态说明

- `pending`：待开始
- `in_progress`：进行中（当前会话）
- `passing`：开发完成，测试通过
- `failed`：测试失败，需修复

## 进度记录

> 每个开发会话结束后追加，记录完成情况。

### 2026-10-08 会话记录（T006）

- 完成任务：T006 A3 拓扑单一来源
- TDD：Red（子代理：PartitionTopologyTest.Np4TopologyIsConsistent + GridTest.PartitionUniform 补所有权释放 + ctest 条目 topology_consistency；编译失败即合法 Red，错误仅限 test_grid.cpp）→ Green（主代理：partition.hpp SubgridInfo+cartComm（默认 MPI_COMM_NULL、所有权注释）；partition.cpp cart rank 查 coords（修 reorder=1 潜在错位）+ 删除内部 Comm_free；main.cpp 删二次 Cart_create，下游经局部 cartComm=info.cartComm 无缝衔接，末尾单次 free 即释放分区器通信器）
- B3 实现期勘误（设计 §6 B3 原文"Σ nxLocal==nx（各维）"有误）：2D tiling 下各 rank x 向全宽参与求和 → Σ=2048=2×1024；正确不变式 = 按行（固定 coords[1]）Σ nxLocal==nx、按列 Σ nyLocal==ny、nz 维不分解逐 rank 相等——已修正测试并在 tasks.md 记录
- B3 断言集：cartComm 非空、分组尺寸和、行/列 offset 链连续（链首 0、链尾恰为全局尺寸）、Cart_shift 四向邻居互指（Sendrecv 令牌核对）、coords 由 cartComm 查询、结尾 MPI_Comm_free（所有权契约演示）
- 修改文件：src/grid/partition.hpp/cpp、src/main.cpp、tests/test_grid.cpp、CMakeLists.txt（+topology_consistency）
- 测试：Debug 全量 25/25（61s）、Release 全量 25/25（14s）、Release 0 警告

### 2026-10-08 会话记录（T005）

- 完成任务：T005 A1 first-touch zeroInitialize 并行化（微型任务，主代理直做 Red+Green）
- Red：U7（test_grid.cpp：2D halo=2 与 3D halo=1 构造，三场投毒后断言全缓冲逐元素==0.0，回归固化语义基线即绿）+ B4（test_solver.cpp：16³ jacobi 901 迭代收敛 + lastResidual≤tol 真实口径 + 最大值原理正值检查）
- Green：subgrid.cpp zeroInitialize 由三次 fill 改为单遍历三场并行——2D 分支 omp parallel for over j∈[0,nyTotal) + 内层 k∈[0,nzTotal)（设计勘误：2D 缓冲含 z-halo 平面，原设计 2D 循环不含 k 会漏清 2/3 缓冲）；3D 分支 k 外层；schedule(static) 与 stencil 循环分区一致（走查 evidence ① 落实）；noexcept 保持
- 修改文件：src/grid/subgrid.cpp、tests/test_grid.cpp（+U7）、tests/test_solver.cpp（+B4）
- 测试：unit 全绿；Debug 全量 24/24（61s）；Release 0 警告 + unit 绿

### 2026-10-08 会话记录（T004）

- 完成任务：T004 B2a+A4-CG 并行化 + 去重 + 区名修复
- TDD：Red（子代理写 U6/E2；U6 初版正弦问题恰为离散算符特征向量 → CG 1 步收敛、不覆盖重构尾部，主代理改用 uniform rhs=-1 多模态问题 → 126 迭代全路径覆盖；golden 采集于基线代码，按回归固化语义 Red 期即绿）→ Green（主代理：updatePInterior 抽取 + omp parallel for/simd、三场置零与初始化循环并行化、区名 cg_iteration；一次编译修正 const 限定）
- 重要实测发现（修正设计 §6 U6 前提）：libgomp 归约按线程到达序合并，OMP=4 档跨运行 hash 不稳定（1.0~1.4e-12 波动）→ 测试口径改为 tier=1 位级 hash golden（0x654fa5ddb3bd66b5，两轮 10+ 次运行稳定）+ tier=4 迭代数 golden（126，收敛裕度 ~20% 不受 ulp 噪声翻转）与 tier=1 解最大偏差 <1e-11
- golden 绑定 Debug+ASan 构建：Release（O3/FMA）低序位合法不同，NDEBUG 下 GTEST_SKIP
- E2：相邻双 inf rhs → r=p=-inf → matvec 得 inf-inf=NaN → pap=NaN 确定性触发 breakdown 守卫（与线程数无关）；断言 0 迭代、u 全零（守卫先于 axpy）、lastResidual=inf、事后 iterate() 安全
- 修改文件：src/solver/cg_solver.cpp（updatePInterior + 三处并行化 + 区名）、src/solver/solver.hpp（私有声明）、tests/test_alt_solvers.cpp（CgAR004Test 两用例 + fnv1a64Interior 助手）
- 测试：CgAR004Test 5+10 次重复运行零失败；Debug 全量 24/24（79s）；Release 0 警告 + unit 绿；tier-1 hash 精确匹配证明并行化未改变单线程数值

### 2026-10-08 会话记录（T003）

- 完成任务：T003 A2-CLI + RunConfig 透出
- TDD：Red（子代理：cmake/VerifyAr004Cli.cmake 四模式 + CMakeLists 8 条目；8 项全失败、旧测试不受影响；两处偏差均合理——verify 条目补 HYPOS_TEST_ENV 抑制 LSan 噪声、拒绝类 Red 期以 Timeout 失败属预期）→ Green（主代理：main.cpp 解析/校验/WARN/setter/RunConfig、reporter.hpp/cpp JSON+CSV 序列化、printUsage；一次通过 8/8）
- 修改文件：src/main.cpp（:137 解析、:202-213 校验+WARN+setter、:332 RunConfig、printUsage+1 行）、src/perf/reporter.hpp（RunConfig +residualCheckInterval=1）、src/perf/reporter.cpp（JSON config 段 +1 字段、CSV +1 列）、CMakeLists.txt（+8 ctest 条目）、cmake/VerifyAr004Cli.cmake（新增）
- 测试：ar004_* 8/8 绿；Debug 全量 ctest 24/24 绿（340s）；Release 构建 0 警告 + unit 绿
- 设计决策落地：CG 与 k≠1 → WARN 含 "ignored"、不调用 setter（与 overlap-comm 先例 :198-200 同构）；CSV 新列插 overlap_comm 后（无按位消费者，经 grep 证实）

### 2026-10-08 会话记录（T002）

- 完成任务：T002 A2-RBGS 真实残差 + interval
- TDD：Red（子代理写 5 用例：U4/U5/R-I1s/R-I5s/R-E4s，运行期断言失败确认，24 既有测试保持绿）→ Green（主代理，一次通过）→ Refactor（新增共享 globalTrueResidual 并回改 Jacobi 确认扫描消除内联重复；design §4.3.1 补录该接口）
- 修改文件：src/solver/residual.hpp/cpp（+globalTrueResidual）、src/solver/solver.hpp（RBGS 文档 + iterateCore 声明）、src/solver/red_black_gs_solver.cpp（iterateCore/iterate/solve 重写 + rbgs_iteration 区名）、src/solver/jacobi_solver.cpp（确认扫描改用共享 helper）、tests/test_alt_solvers.cpp（+5 用例）
- 测试：unit 全绿（Debug 9.92s / Release 2.38s）；Debug 全量 ctest 16/16 绿（144s，较 T001 前增加——RBGS 扫描的诚实成本，Debug+ASan 口径）；Release 0 警告

### 2026-10-08 会话记录（T001）

- 完成任务：T001 A2-Jacobi 真实残差 + interval
- TDD：Red（子代理写 9 个用例：U1 2D/3D、U2、U3、U9、U10、I1s、I5s、E4s，编译失败确认）→ Green（主代理实现，一次通过）→ Refactor（无重复可消除；确认扫描 helper 待 T002 出现第二消费方时评估提取）
- 修改文件：src/solver/residual.hpp/cpp（新增）、src/solver/solver.hpp（基类 setter+成员+Jacobi 文档）、src/solver/jacobi_solver.cpp（iterate 换算+solve 重构）、CMakeLists.txt（+residual.cpp）、tests/test_solver.cpp（+9 用例）
- 测试：unit 全绿（Debug 9.16s / Release 1.76s）；全量 ctest 16/16 双构建绿；0 新增警告
- 实现注记：①maxIter=0 出口跳过确认扫描（ZeroMaxIterations 契约「未迭代返回 0」）；②setter 对 (Index)(-5) 回绕值经 >2^30 上界检测钳为 1；③判据 Allreduce 复用 iterate() 已换算返回值（平方回代 ~1 ulp，位级一致约束不受影响——overlap 两模式同路径）
- 基线采集：evidence/baseline-bench.md（256² 每 iter + 64² 收敛口径，3 次中位数）

### 2026-10-08（design 阶段）：srs §4 性能回归门槛修订——RBGS 真实残差 k=1 存在必然开销（每 k 步一次扫描 + halo 交换），±5% 门槛与诚实性目标物理上不可兼得，改为「不设门槛、k=1/5/10 三档如实测量」（design.md 决策 D10）。属 design 阶段对 req 的定向反馈，其余 srs 内容未动。

## 阶段门控记录

> 由 sdd-phase-gate skill 在阶段门控审查后追加，记录每轮审查结果（PASS/FAIL + 轮次）。格式见 sdd-phase-gate SKILL.md Step 6。

### 2026-10-08 design 门控记录（第 2 轮）

- 门控结果：PASS（第 2 轮）
- 审查项数：15 项（G10 SKIP-目录不存在；0 NO；2 项 Minor WARN）
- 修复的问题：第 1 轮 G12/G14（Important，区名用例 U11 + 走查 evidence 落 §6.5）及 5 项 Minor 全部修复到位；本轮新发现 W1（ThreadData 生命周期：堆分配永不释放 + reset 保 registry）、W2（residual_confirm「建议」改规范性定名）已补入 design.md
- 审查代理：sdd-gate-reviewer
- 亮点：行号引用 15 处与源码相符；代数换算经逐行核验；I1/I5 精确计数口径自洽

### 2026-10-08 design 门控记录

- 门控结果：FAIL（第 1 轮）
- 失败项：G12（功能点 10 区名修复无 §6 用例）、G14（srs §3.5 验收② 无追溯用例；§3.3①/§3.4② 走查验收未落 §6.5 evidence）——Important，同根因
- 严重性：Important；另 G3/G9/G11/G15 + CLI -1 可读性 5 项 Minor
- 审查代理：sdd-gate-reviewer
- 修复：§6.1 新增 U11 区名断言；§6.5-1 补三项走查 evidence；§4.3.2 补 noexcept 声明；§4.4 补 reporter.cpp；U6 按档固化 golden；I1/I5 改精确计数（确认扫描独立区名 residual_confirm）；B1 判据收紧；E4 改 k≥maxIter 边界 + E5 前置条件走查；CLI -1 空格形式口径注明

### 2026-10-08 req 门控记录

- 门控结果：PASS（第 1 轮）
- 审查项数：9 项（3 项 Minor WARN）
- 修复的问题：G3 §4 残差一致性补 CG 递推/容差口径；G4 §1「全部缺口」措辞收敛为实际范围；G4 §3.3 补 applyDirichletBC 评估点。特别核验 5/5 与代码相符（Jacobi ω=1 代数关系、CG dot(r,r) 不可跳、RBGS 无换算关系均实证）
- 审查代理：sdd-gate-reviewer
