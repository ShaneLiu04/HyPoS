# [AR007] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR007 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-10 |

## 测试用例列表

### ST-001：制造解收敛与二阶精度（np1）

**关联需求：** srs.md §3.1 R1 验收①
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 制造解（已知解析解 rhs），np=1

**测试步骤：**
1. ctest -R unit 中 TwoLevelMGUnitTest.ManufacturedSineConvergenceNp1（64²，tol 1e-7）
2. TwoLevelMGUnitTest.SecondOrderConvergence（解析 rhs，32²/64² L2 误差比）

**期望结果：**
- Then 收敛且解误差随网格加密二阶缩减（比率 ∈ [3.5, 4.5]）

**实际结果：** np1 组 9/9 通过；SecondOrderConvergence 比率 ≈4（Release/Debug 双模式）
**状态：** PASS

---

### ST-002：mg2 平滑步数严格少于 RBGS 迭代数（256² 同 tol）

**关联需求：** srs.md §3.1 R1 验收②（GUIDE B1a 验收①）
**测试类型：** 正常路径（记录式）
**优先级：** High

**前置条件：**
- Given 256² 制造解、同 tol（1e-2，绝对真残差口径）、OMP=1、np=1/4、3 次中位

**测试步骤：**
1. scripts/bench_ar007.sh 逐组运行 mg2 / red_black_gs / jacobi

**期望结果：**
- Then mg2 总平滑步数（cycles×(ν1+ν2)）严格少于纯 RBGS 迭代数；倍数如实记录

**实际结果：** mg2 14 cycles / 56 平滑步 / 0.10s；red_black_gs 85540 迭代 / 11.6s；jacobi 166440 / 10.0s（np1）；np4 同量级（mg2 14/56/0.10s vs rbgs 83021）。cycle 数规模无关（tol 1e-6：64²/128²/256²=20/22/23）。PERFORMANCE §14 入档（含审查第 3 轮独立复测逐位复现）
**状态：** PASS

---

### ST-003：既有求解器行为零变化（回归）

**关联需求：** srs.md §3.1 R1 验收③
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 任意既有负载（--solver rbgs/jacobi/cg）

**测试步骤：**
1. 全量 ctest（35 既有条目，含 RBGS 位级一致守护、overlap 位级守护）

**期望结果：**
- Then 行为零变化、不变绿转

**实际结果：** 38/38 全绿（Release 15.4s / Debug-ASan+UBSan 35.7s）；既有 10 个测试文件零语义修改（diff --stat 仅新增 test_mg_two_level.cpp）；审查第 1 轮 S3/D3 确认三 exchanger/jacobi/cg/IO 栈零触碰
**状态：** PASS

---

### ST-004：不支持配置防御（CLI 层实跑）

**关联需求：** srs.md §3.1 R1 验收④
**测试类型：** 异常处理
**优先级：** High

**前置条件：**
- Given --bc neumann --solver mg2（及 3D、奇数维）

**测试步骤：**
1. 实跑 `hypos --solver mg2 --bc neumann`（32²）
2. 实跑 `hypos --solver mg2 --nz 4`（3D）
3. 实跑 `hypos --solver mg2 --nx 33`（奇维）
4. solver 层：TwoLevelMGUnitTest.DefenseRejectsUnsupportedConfigs（5 子场景）

**期望结果：**
- Then HYPOS_ERROR 明示且不产出结果

**实际结果：** 三条 CLI 实跑均输出 [ERROR]（neumann/2D/偶维前置校验）且输出目录无结果文件；DefenseRejectsUnsupportedConfigs 5 子场景全过（solve==0、u 逐位不动、lastResidual==0）
**状态：** PASS

---

### ST-005：限制/延拓已知系数场对照（R/P 纯内核）

**关联需求：** srs.md §3.2 R2 验收①
**测试类型：** 正常路径 + 边界条件
**优先级：** High

**前置条件：**
- Given 已知系数场（含边界行/列/角点）

**测试步骤：**
1. RestrictionKnownCoefficients / ProlongationKnownCoefficients / ProlongRestrictProjectionIdentity

**期望结果：**
- Then 逐点值与手算期望一致；P∘R 在粗点投影恒等

**实际结果：** 三个单测全过（含空粗片、粗边界 9/16 与 12/16 权、右/上奇细点 1/2 权、P∘R 位级恒等）
**状态：** PASS

---

### ST-006：np=1/4 非均匀布局校正一致性

**关联需求：** srs.md §3.2 R2 验收②（GUIDE B1a 验收④）
**测试类型：** 边界条件
**优先级：** High

**前置条件：**
- Given np=1 与 np=4 非均匀切分（{15,19} 手工布局）同一全局残差场

**测试步骤：**
1. CrossLayoutCorrectionConsistency：分布式 1 cycle vs 串行参考逐元素对照
2. DiagonalCornerExchangeFilled / CoarseOwnershipPartitionsGlobalRange

**期望结果：**
- Then 校正场两种布局逐元素一致（maxDiff ≤ 1e-12·尺度）

**实际结果：** np4 组 4/4 通过（CrossLayout 修复 pre-smooth 后 halo 陈旧 bug 后转绿——bug 由本测试捕获，过程记录于 tasks.md 偏差 4）
**状态：** PASS

---

### ST-007：奇数维/粗层过小/布局不铺满防御

**关联需求：** srs.md §3.2 R2 验收③
**测试类型：** 异常处理
**优先级：** Medium

**前置条件：**
- Given nx 或 ny 奇数（或粗网格 <4²、布局不铺满）

**测试步骤：**
1. DefenseRejectsUnsupportedConfigs 奇维/粗 3×3/不铺满子场景 + CLI 奇维实跑

**期望结果：**
- Then HYPOS_ERROR 且不产出结果

**实际结果：** 单测子场景 + CLI 实跑（ST-004 步骤 3）均按预期拒绝
**状态：** PASS

---

### ST-008：RBGS 既有位级一致守护不变

**关联需求：** srs.md §3.3 R3 验收①
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 既有 RBGS 全部测试（含位级一致守护）

**测试步骤：**
1. 全量 ctest 中 RBGS 相关条目

**期望结果：**
- Then 不变绿转

**实际结果：** 38/38 全绿；red_black_gs_solver.cpp 仅 +8 行 smooth（既有 solve/iterate/sweep/iterateCore 逐行不动，审查第 1 轮 C1 逐行核对）
**状态：** PASS

---

### ST-009：平滑入口位级一致（smooth(n) ≡ iterate()×n）

**关联需求：** srs.md §3.3 R3 验收②
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 同一初始场

**测试步骤：**
1. SmoothMatchesIterateBitExact
2. ResidualFieldMatchesTrueResidualScan（残差场内核口径）

**期望结果：**
- Then 更新结果位级一致（sweep 序完全相同）

**实际结果：** 两用例通过（Release/Debug 双模式）
**状态：** PASS

---

### ST-010：PERFORMANCE §14 入档无失实

**关联需求：** srs.md §3.4 R4 验收
**测试类型：** 文档走查（记录式）
**优先级：** Medium

**前置条件：**
- Given PERFORMANCE 新小节

**测试步骤：**
1. 文档数字与 bench 日志核对；口径说明完整性检查

**期望结果：**
- Then 含迭代数与耗时两组数据及口径说明，无失实宣称

**实际结果：** §14 含全量数据表（np1/np4 × 三 solver）、口径注（平滑步计法）、规模无关性数据、适用 tol 范围注记（粗层容差饥饿边界，诚实边界）、实现要点（×4 尺度补偿）。审查第 3 轮独立复测全部数字逐位复现并勘误 3 处口径滑误（已修复提交 5dfc46b）
**状态：** PASS

---

### ST-011：全量 ctest 双模式全绿 + ASan/UBSan 零报告

**关联需求：** srs.md §3.5 R5 验收
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 全量 ctest（35 既有 + 3 新增 = 38 条）

**测试步骤：**
1. Release 与 Debug(ASan+UBSan) 双模式 ctest

**期望结果：**
- Then 全绿且 ASan/UBSan 零报告

**实际结果：** Release 38/38（15.4s）、Debug 38/38（35.7s，detect_leaks=0 口径沿既有惯例）
**状态：** PASS

---

### ST-012：e2e 直跑 mg2（np4）

**关联需求：** srs.md §3.5 R5（e2e：直跑 hypos --solver mg2）
**测试类型：** 集成
**优先级：** High

**前置条件：**
- np=4 mpirun，32²，max-iter 20

**测试步骤：**
1. `mpirun -np 4 ./hypos --nx 32 --ny 32 --solver mg2 --max-iter 20 --output-dir ...`
2. 检查产出与 performance_report.json

**期望结果：**
- Then 正常运行、日志含 cycle 汇总、产出完整

**实际结果：** 17 cycles 完成（日志「TwoLevelMG finished」+ 逐 rank 粗层 CG），performance_report.json 产出（iterations=17）；ctest mg2_e2e 条目 Passed
**状态：** PASS

---

### ST-013：文档与 --help 一致性走查

**关联需求：** srs.md §3.6 R6 验收
**测试类型：** 文档走查
**优先级：** Medium

**前置条件：**
- Given 文档走查 + `--help` 实跑

**测试步骤：**
1. `hypos --help` 实跑（solver 行、residual-check-interval 行）
2. README/GUIDE/AGENT_SPEC/PERFORMANCE 走查

**期望结果：**
- Then 无失实宣称、帮助与行为一致

**实际结果：** --help 实跑：`--solver ... jacobi, red_black_gs, cg, mg2`、`--residual-check-interval ... (jacobi/red_black_gs/mg2 outer cycles; ...)` 与 solver.hpp 接口注释一致；README 特性表/--solver/--residual-check-interval/目录树、GUIDE B1a 完成标注/★4 勾选/§3.4 天花板部分突破（B1b 仍列未实现）、AGENT_SPEC 类层级+算子+residualFieldLocal 契约均更新且经审查第 3 轮无失实复核
**状态：** PASS

---

## ST 执行报告

| 字段 | 内容 |
|------|------|
| 执行日期 | 2026-10-10 |
| 执行结果 | PASS |
| 执行轮次 | 第 1 轮 |

### 需求覆盖矩阵

| 需求 ID | 需求描述 | 测试用例 | 结果 |
|--------|---------|---------|------|
| §3.1 R1 | 两层校正求解器 | ST-001, ST-002, ST-003, ST-004 | PASS |
| §3.2 R2 | 限制/延拓算子与粗网格生成 | ST-005, ST-006, ST-007 | PASS |
| §3.3 R3 | RBGS 平滑复用（不改既有行为） | ST-008, ST-009 | PASS |
| §3.4 R4 | 性能对照入档 | ST-010 | PASS |
| §3.5 R5 | 测试矩阵 | ST-011, ST-012 | PASS |
| §3.6 R6 | 文档同步 | ST-013 | PASS |

**需求覆盖率：** 6 / 6（100%）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 6 | 6 | 0 | 0 |
| 边界条件 | 3 | 3 | 0 | 0 |
| 异常处理 | 2 | 2 | 0 | 0 |
| 回归测试 | 3 | 3 | 0 | 0 |
| 集成/e2e | 1 | 1 | 0 | 0 |
| 文档走查 | 2 | 2 | 0 | 0 |
| **合计** | **17** | **17** | **0** | **0** |

（自动化承载：全量 ctest 38 条双模式 + bench 3×6 组；手工项：CLI 防御实跑×3、--help 实跑、e2e np4 实跑——本会话现场执行，记录于上）

### 遗留问题

| 严重性 | 描述 | 处理方式 |
|-------|------|---------|
| Minor | mg2 粗层绝对容差在极紧 tol 下饥饿（256² tol<~1e-9 不收敛）——§14/README 已注记适用范围 | 延后至 B1b（后续 AR 考虑粗层相对容差）；已文档化不阻塞 |
| Minor | 512² 可扩展性佐证为 R4 可选项，未加测（256² 已充分支撑规模无关结论） | 按 srs §5 假设条款「可选」处理，不阻塞 |
| Info | 预修复周期的 ∝n² 尾段为「欠校正×容差饥饿」耦合效应（tasks.md 偏差 6 已机制说明） | 已记录 |

### 结论

> **Go**。全部 Go 条件满足：需求覆盖 100%（6/6）、无 Critical/Major 缺陷、单元/集成/e2e/回归全绿（Release 38/38 + Debug-ASan/UBSan 38/38）、NFR 达标（行数 633 ≤ ~800、双模式 0 新增警告、MPI 路径无异常）、记录式性能判定如实入档且经独立复测。2 项 Minor 均已文档化/延后，不阻塞归档。
