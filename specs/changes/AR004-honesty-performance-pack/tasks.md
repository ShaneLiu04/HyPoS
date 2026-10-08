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
| T003 | A2-CLI `--residual-check-interval`：基类 setResidualCheckInterval（0 钳 1）；main 解析/校验（<1 退出码 1；CG 传入 WARN 忽略）；RunConfig +residualCheckInterval 透出。测试先行：U9、I1-I6 | T001 | pending | design D4/D9；U9 已随 T001 落地 |
| T004 | B2a+A4-CG：updatePInterior 抽取（消除 solve/iterate 重复）；初始化循环并行化；区名改 cg_iteration。测试先行：U6（golden 迭代数，先跑基线固化）、E2 | - | pending | design §4.2.2-4；位级不变 |
| T005 | A1 First-touch：zeroInitialize() 外层维度并行（全缓冲，三场单遍历）；applyDirichletBC 保持串行并文档记录。测试先行：U7、B4 | - | pending | design D8 |
| T006 | A3 拓扑一致性：SubgridInfo +cartComm（所有权移交）；partition 用 cart rank 查 coords；main 删除二次 Cart_create。测试先行：B3 | - | pending | 接口变更授权见 design §4.3.1；design D5/D6 |
| T007 | A4 Profiler：per-thread ThreadData 注册表（热路径零锁）+ 聚合；hpp 注释对齐；test_performance.cpp:83 区名断言复核（RBGS 区名已随 T002 落地）。测试先行：U8、U11 | T004 | pending | design D7；JSON 格式兼容 |
| T008 | 全量回归（Debug/Release + ASan；TSan 不可用则走查佐证）+ 基准（scripts/bench_ar004.sh：Jacobi ±5% 门槛、RBGS k=1/5/10 三档、CG 前后，≥3 次中位数）+ 文档同步（README/AGENT_SPEC/PERFORMANCE §11/GUIDE 勾选 G1/G3/G5/G6、P1/P4/P6） | T001-T007 | pending | SCALING_REPORT 口径注记 |

## 状态说明

- `pending`：待开始
- `in_progress`：进行中（当前会话）
- `passing`：开发完成，测试通过
- `failed`：测试失败，需修复

## 进度记录

> 每个开发会话结束后追加，记录完成情况。

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
