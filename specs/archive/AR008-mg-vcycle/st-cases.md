# [AR008] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR008 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-10 |

## 测试用例列表

### ST-001：mgv 制造解收敛 + cycle 数 h 无关

**关联需求：** srs.md §3.1 R1（验收①）
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 制造解负载（解析 rhs），256² 与 512²，tol 收敛口径

**测试步骤：**
1. ctest 运行 U3（ManufacturedSineConvergence 64²/256²）与 U4（RateHIndependence 256²/512²）
2. bench 复证（scripts/bench_ar008.sh，tol 1e-6）

**期望结果：**
- Then 两规模收敛且 cycle 数相近（design D8 判据 |cyc512−cyc256|≤4 且 ≤40）；W-cycle 修订后判据不变

**实际结果：** U3/U4 绿：cyc256=8、cyc512=7（D8 实测）；bench §15：15/14 cycles（tol 1e-6）——4× 细化反降 1，h 无关成立
**状态：** PASS

---

### ST-002：mgv 解二阶精度

**关联需求：** srs.md §3.1 R1（验收①后半）
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 64²/128² 制造解收敛解

**测试步骤：**
1. ctest SecondOrderConvergence（沿 AR007 惯例 ratio∈(3.5,4.5)）

**期望结果：**
- Then 收敛阶 ratio 接近 4（L2 范数二阶）

**实际结果：** 断言绿（W-cycle 下解精度不受 cycle 类形影响）
**状态：** PASS

---

### ST-003：mgcg 迭代数少于裸 cg

**关联需求：** srs.md §3.1 R1（验收②）/ §3.3 R3
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 256²/512² 制造解，mgcg 与裸 cg 各自求解同 tol

**测试步骤：**
1. ctest U5（multiMode rhs：cg=45 vs mgcg 严格更少）
2. bench §15 记录式对比（256²/512²、np=1/4）

**期望结果：**
- Then mgcg 迭代数少于裸 cg；倍数如实入档（记录式，不预设门槛）

**实际结果：** U5 绿（多模式负载 cg=45，mgcg 更少）；§15：mgcg 8 迭代恒定 vs cg 700/1378（≈87×/172×）；np1 512² 耗时 5.1×
**状态：** PASS

---

### ST-004：既有求解器行为零变化（回归）

**关联需求：** srs.md §3.1 R1（验收③）
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 任意既有负载，--solver jacobi/red_black_gs/cg/mg2

**测试步骤：**
1. 全量 ctest 双模式（42 条含既有 38 条：mg2 冻结契约/cg golden/AR004-006）
2. git diff 冻结文件核对（solver.hpp/mg_operators.hpp/mg_two_level_solver.cpp/red_black_gs_solver.cpp 零改动，实现审查已核）

**期望结果：**
- Then 行为零变化

**实际结果：** release 42/42（18.88s）、debug ASan/UBSan 42/42（66.77s）零回归；冻结契约文件 diff 为空
**状态：** PASS

---

### ST-005：不支持配置防御（3D/Neumann/奇维/链失败）

**关联需求：** srs.md §3.1 R1（验收④）
**测试类型：** 异常处理
**优先级：** High

**前置条件：**
- Given 3D/Neumann/奇维负载，--solver mgv/mgcg 启动

**测试步骤：**
1. CLI 黑盒实跑：`--nx 64 --ny 64 --nz 2 --solver mgv` 与 `--solver mgcg`
2. ctest U7（Defense 四配置 × 双驱动，0 迭代+u 不动断言）+ U1 链生成拒绝（8×20 狭长）

**期望结果：**
- Then HYPOS_ERROR 明示退出，不产出失实结果

**实际结果：** 实跑日志 `[ERROR] --solver mgv requires 2D grids (--nz 1)` 与 mgcg 同款（2026-10-10 04:34 UTC）；U7/U1 断言绿
**状态：** PASS

---

### ST-006：R/P/内核复用零回归

**关联需求：** srs.md §3.2 R2（验收①）
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given AR007 全部 mg2 测试（38 条基线）

**测试步骤：**
1. 全量 ctest 双模式（既有 38 条含 mg2 系列）

**期望结果：**
- Then 不变绿转（R/P 复用契约，mg_operators.hpp 零修改）

**实际结果：** 38 条既有全绿（含 mg2 冻结契约与 P∘R 恒等）；mg_operators.hpp diff 为空
**状态：** PASS

---

### ST-007：跨布局一致性（np=1/4 非均匀）

**关联需求：** srs.md §3.2 R2（验收②）
**测试类型：** 边界条件
**优先级：** High

**前置条件：**
- Given np=1 与 np=4 非均匀切分 {15,19}

**测试步骤：**
1. ctest U8（CrossLayout，校正场逐元素一致）+ U9（mgcg np4 解 L2≤1e-10·‖u‖+轨迹采样）

**期望结果：**
- Then 校正场/解逐元素（或约定容差）一致

**实际结果：** U8/U9 双绿；U9 实测 np4 与 np1 迭代轨迹一致（§15：8 迭代均为 8）
**状态：** PASS

---

### ST-008：粗化链逐层结构正确

**关联需求：** srs.md §3.2 R2（验收③）
**测试类型：** 正常路径
**优先级：** Medium

**前置条件：**
- Given 粗化链 256→128→…→8（及截断链 [34,17]、12² 边界）

**测试步骤：**
1. ctest U1（链枚举：256/64/34/12 全对照 D2 表）+ U2（三层机器精度锚）

**期望结果：**
- Then 每层规模/offset/尺度补偿逐层正确

**实际结果：** U1 链枚举断言绿（含 [34,17] 截断与 [12,6] 边界）；U2 参考对比 ≤1e-12/1e-13（D12 机器精度级）
**状态：** PASS

---

### ST-009：mgcg 收敛 + np4 一致（R3 汇总）

**关联需求：** srs.md §3.3 R3（验收①②）
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 制造解 256²/512²，np=1/4

**测试步骤：**
1. ctest mgcg_e2e（直跑二进制）+ U5/U9；bench §15 np1/np4 对照

**期望结果：**
- Then 收敛、迭代数规模不敏感、np4 与 np1 一致

**实际结果：** mgcg_e2e 绿；§15：8 迭代在 256²/512²×np1/np4 四组恒定（规模与进程双无关=预条件后条件数 O(1) 直证）
**状态：** PASS

---

### ST-010：PERFORMANCE §15 数据完整如实

**关联需求：** srs.md §3.4 R4
**测试类型：** 文档走查
**优先级：** Medium

**前置条件：**
- Given PERFORMANCE §15 新小节

**测试步骤：**
1. 走查 §15：mgcg vs cg 与 mgv vs mg2 两组数据、口径说明（OMP=1/3 次中位/tol 语义）、无失实宣称
2. 核对 bench 脚本可复现（scripts/bench_ar008.sh）

**期望结果：**
- Then 数据、口径、结论注记齐备，记录式判定

**实际结果：** §15 含 16 行数据表（2 规模×2np×4solver）+6 条结论注记（含 W-cycle 结构注记、mg2 np4 反慢如实记录、四种迭代数口径不可互换说明）；脚本入仓可复现
**状态：** PASS

---

### ST-011：全量测试矩阵双模式

**关联需求：** srs.md §3.5 R5
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 全量 ctest（38 既有+4 新增=42）

**测试步骤：**
1. release ctest + debug（ASan+UBSan）ctest

**期望结果：**
- Then 全绿且 ASan/UBSan 零报告

**实际结果：** release 42/42（18.88s）、debug 42/42（66.77s，detect_leaks=0 口径），零报告
**状态：** PASS

---

### ST-012：文档与 --help 同步

**关联需求：** srs.md §3.6 R6
**测试类型：** 文档走查
**优先级：** Medium

**前置条件：**
- Given README/GUIDE/AGENT_SPEC/PERFORMANCE/--help

**测试步骤：**
1. `hypos --help` 实跑核 solver 行与 interval 行
2. 走查五处文档（README 特性表+CLI 表+目录树、GUIDE §3.4/★5/§5、AGENT_SPEC 类层级、§15）

**期望结果：**
- Then 无失实宣称、帮助与行为一致

**实际结果：** `--help` 实跑：`--solver <string> Solver type: jacobi, red_black_gs, cg, mg2, mgv, mgcg` 与 interval 行含 mgv、排除 mgcg（与 WARN 行为一致）；五处文档均已更新（T005 提交 8160fee）；B1b「≤8³」已澄清「每维 ≤8」
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
| §3.1 R1 | mgv 多层 MG 求解器（收敛/二阶/对照/回归/防御） | ST-001, ST-002, ST-003, ST-004, ST-005 | PASS |
| §3.2 R2 | 粗化链与算子复用 | ST-006, ST-007, ST-008 | PASS |
| §3.3 R3 | mgcg MG-CG 预条件 | ST-003, ST-009 | PASS |
| §3.4 R4 | 性能对照入档 | ST-010 | PASS |
| §3.5 R5 | 测试矩阵 | ST-011 | PASS |
| §3.6 R6 | 文档同步 | ST-012 | PASS |

**需求覆盖率：** 6 / 6（100%）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 6 | 6 | 0 | 0 |
| 边界条件 | 2 | 2 | 0 | 0 |
| 异常处理 | 1 | 1 | 0 | 0 |
| 回归测试 | 3 | 3 | 0 | 0 |
| **合计（ST 用例）** | **12** | **12** | **0** | **0** |
| 底层自动化证据 | ctest 42×2 模式 + bench 16 组 + CLI 黑盒 3 次 | 全绿 | 0 | 0 |

### 遗留问题

| 严重性 | 描述 | 处理方式 |
|-------|------|---------|
| Minor | 行数预算净增 ~908 超 ~800 约 13% | 归档记录背书（T005 论证在案：W-cycle 为正确性必需，10 组探针证据链证明单 V 不可交付） |
| Minor | GUIDE §5 已补登 Neumann-MG/3D-MG 后续条目 | 留待后续 AR 排期（req V1/V4 门控已落实） |

### 结论

> **Go**：需求覆盖 6/6（100%）；12 个 ST 用例全 PASS（含 3D 防御 CLI 黑盒实跑与 --help 实跑证据）；全量 ctest 42 条×release/debug 双模式零回归、ASan/UBSan 零报告；PERFORMANCE §15 记录式数据完整；无 Critical/Major 缺陷；2 个 Minor 已记录并确认延后。
