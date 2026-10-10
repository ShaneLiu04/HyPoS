# [AR009] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR009 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-10 |

> 执行口径：ST 为黑盒验证（直跑 `hypos` 二进制 + 工件断言），自动化承载于 ctest（fixture 链含直跑 e2e 条目与 verify 断言条目）；用例标注承载条目，执行记录见文末报告。srs 开发期修订注记（R2① 单调删除、R2③ np4 判据、R1 Iallreduce 融合口径）已纳入判据。

## 测试用例列表

### ST-001：pcg 与 cg 迭代数一致性（±5%）+真残差 ≤ tol

**关联需求：** srs §3.1 R1 验收①
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- 制造解负载（256²/512²，tol 1e-6，OMP=1，WSL Release 构建）

**测试步骤：**
1. `bench_ar009.sh` 逐配置跑 pcg/cg（3 次中位）
2. 对照迭代数与 finish 行真残差口径

**期望结果：**
- 四配置（256²/512²×np1/4）迭代数一致 ±5%（含 512² 半项代偿 R1①）；最终真残差 ≤ tol

**实际结果：** PASS——四配置迭代数**精确相同**（700/834/1378/1697，零差异远优于 ±5%）；pcg_e2e 真残差复核 ≤ tol(1+5e-4)（PERFORMANCE §16+ctest pcg_e2e/verify 绿）

**状态：** PASS

---

### ST-002：pcg np=4 与 np=1 解一致

**关联需求：** srs §3.1 R1 验收②
**测试类型：** 正常路径（np 集成）
**优先级：** High

**前置条件：**
- np4 mpirun（HYPOS_EXPECT_NP=4 守卫）

**测试步骤：**
1. np1/np4 各跑 `--solver pcg` 同负载
2. 解向量 L2 差 ≤ 1e-10·‖u‖（U9 惯例）

**期望结果：**
- L2 ≤ 1e-10·‖u‖

**实际结果：** PASS——ctest pcg_np4（HYPOS_EXPECT_NP=4 守卫）绿；P4 断言 L2 差 0（同轨迹同迭代数）

**状态：** PASS

---

### ST-003：既有六 solver 值行为零变化（全量回归）

**关联需求：** srs §3.1 R1 验收③
**测试类型：** 回归测试
**优先级：** High

**测试步骤：**
1. 全量 ctest（既有 42 条+AR009 前 60 条口径中的既有部分，含 cg golden U6 位级守护）

**期望结果：**
- 既有条目全绿；FP2 拆分后 `--solver cg` 浮点序位级不变（U6 golden 断言）

**实际结果：** PASS——Release 全量 60/60（含 AR009 前 42 条既有全绿）；U6 cg golden（FNV-1a 位级+迭代数）过

**状态：** PASS

---

### ST-004：pcg 每迭代阻塞 Allreduce 计数 = 0

**关联需求：** srs §3.1 R1 验收④
**测试类型：** 正常路径（探针）
**优先级：** High

**测试步骤：**
1. Profiler 计数：pcg 运行全程 `cg_blocking_allreduce` region 计数
2. 阳性对照：cg 同负载计数 == 2·iters+1

**期望结果：**
- pcg 阻塞归约计数 == 0（探针非空洞自证：cg 阳性对照成立）；`pcg_iallreduce` == iters+1

**实际结果：** PASS——P3 双向断言绿：pcg `cg_blocking_allreduce`==0+`pcg_iallreduce`==iters+1；cg 阳性对照==2·iters+1（探针通道非空洞自证）

**状态：** PASS

---

### ST-005：pcg 边界与防御（零迭代/截断/续态/breakdown 记录式）

**关联需求：** srs §3.1 R1（breakdown 防御+iterate 契约）
**测试类型：** 边界条件
**优先级：** Medium

**测试步骤：**
1. rhs=0 零迭代（P8）；maxIter=3 截断（P9）；iterate 续态/前置违反（P7）；ν≤0 防御代码走查（P5 记录式——正定负载不可稳定构造，沿 AR008 breakdown 先例）

**期望结果：**
- 零迭代计数边界==1；截断轨迹正确；前置违反 WARN+返回 0；防御分支在案（①② 顺序防误报）

**实际结果：** PASS——P8（rhs=0 零迭代+计数==1）/P9（maxIter=3 截断轨迹）/P7（前置违反 WARN+0）绿；ν≤0 防御代码在案（review D2 轮核 ①② 顺序防误报，P5 记录式沿 AR008 breakdown 先例）

**状态：** PASS

---

### ST-006：cg 残差历史 CSV 契约

**关联需求：** srs §3.2 R2 验收①（含 2026-10-10/T004 修订：单调子句删除）
**测试类型：** 正常路径
**优先级：** High

**测试步骤：**
1. 直跑 `--solver cg --residual-history f.csv` 至收敛
2. 断言 header==`iteration,residual`、数据行数==迭代数、末行 %.4e 打印精度内==JSON final_residual

**期望结果：**
- 三断言全过；无单调性断言（修订后判据）

**实际结果：** PASS——E1（ctest residual_history_cg+verify）绿：header==`iteration,residual`、数据行数==迭代数、末行 %.4e==JSON final_residual（reformat4e 打印精度内一致）

**状态：** PASS

---

### ST-007：六 solver 值 CSV 均合法且首末残差与收敛语义一致

**关联需求：** srs §3.2 R2 验收②
**测试类型：** 正常路径
**优先级：** High

**测试步骤：**
1. jacobi/rbgs/mg2/mgv/cg/pcg/mgcg 七值各带 --residual-history 直跑（fixture 链）
2. per-solver 口径断言（jacobi 首行‖b‖+末行 Converged 行 %g；其余末行 %.4e==JSON；mgcg ≥5e-6 相对）

**期望结果：**
- 六分支全过（mgcg 走 finish 行残差口径）

**实际结果：** PASS——E2 六分支 ctest 全绿（jacobi 首行‖b‖ 1e-12 相对+末行 %g==Converged 行；rbgs/mg2/mgv/cg/pcg 末行 %.4e==JSON；mgcg ≥5e-6 相对对 finish 行）

**状态：** PASS

---

### ST-008：np4 残差历史单 rank0 文件+内部一致性

**关联需求：** srs §3.2 R2 验收③（含 2026-10-10/T003 修订：np1 对照降级 evidence-only）
**测试类型：** 正常路径（np 集成）
**优先级：** High

**测试步骤：**
1. np4 直跑带 --residual-history（目录输出）
2. 断言目录恰一份 CSV；header+行数==iterations+末行 %.4e==final_residual（内部一致性）

**期望结果：**
- 单文件+三断言过；np1/np4 行数对照 RecordProperty 记录（evidence-only）

**实际结果：** PASS——E3（ctest residual_history_np4+np4_verify）绿：目录恰一份 CSV+header/行数==iterations/末行 %.4e==final_residual；np1 对照经 RecordProperty 记录（evidence-only）

**状态：** PASS

---

### ST-009：残差历史异常处理（目录不存在/不可写）

**关联需求：** srs §3.2 R2 异常处理
**测试类型：** 异常处理
**优先级：** Medium

**测试步骤：**
1. `--residual-history /nonexistent_dir/f.csv` 直跑（两变体：带/不带 --output-dir）

**期望结果：**
- HYPOS_ERROR 明示（[ERROR] 行）+非零退出码

**实际结果：** PASS——E6 两变体（带/不带 --output-dir）ctest 绿：[ERROR] 输出+非零退出码

**状态：** PASS

---

### ST-010：interval 采样语义

**关联需求：** srs §3.2 R2（每 k 步采样沿 --residual-check-interval）
**测试类型：** 正常路径
**优先级：** Medium

**测试步骤：**
1. `--solver jacobi --residual-check-interval 10 --residual-history f.csv` 直跑
2. 断言数据行数==iterations/10（整除前提）

**期望结果：**
- 行数符合采样公式

**实际结果：** PASS——E4（ctest residual_history_interval）绿：数据行数==iterations/10（整除前提）

**状态：** PASS

---

### ST-011：与 --output-dir 组合（双通道写出）

**关联需求：** srs §3.2 R2（无耦合+组合）
**测试类型：** 正常路径
**优先级：** Medium

**测试步骤：**
1. 同时传 --output-dir 与 --residual-history（jacobi 收敛负载）
2. 断言两产物并存（解文件+CSV）且 CSV 末行对齐 Converged 行

**期望结果：**
- 双通道并存，CSV 契约成立；jacobi 单调性成立（E7，定常收缩）

**实际结果：** PASS——E7（ctest residual_history_save_combo）绿：解文件+CSV 双产物并存；CSV 末行 %g 对齐 Converged 行；jacobi 单调断言过（checkMonotone 1e-12 相对）

**状态：** PASS

---

### ST-012：绘图脚本（S1 代偿口径）

**关联需求：** srs §3.3 R3
**测试类型：** 正常路径（资产交付，环境受限代偿）
**优先级：** Low

**测试步骤：**
1. `scripts/plot_residual_history.py` 语法自检（python ast/compile——若环境可行）
2. CSV 契约测试为验收主体（ST-006~008）；PNG 运行演示如实注明不可行（WSL 无 python3）

**期望结果：**
- 语法自检 OK（或同等只读自检）；代偿口径记录在案（srs §5 假设/W3）

**实际结果：** PASS（代偿口径）——Windows 侧 python `ast.parse` 只读自检=SYNTAX OK；WSL 无 python3（复测确认），PNG 运行演示不可行如实注明；CSV 契约测试（ST-006~008）为验收主体（W3 口径）

**状态：** PASS

---

### ST-013：PERFORMANCE §16 数据与口径

**关联需求：** srs §3.4 R4
**测试类型：** 文档验收（记录式）
**优先级：** Medium

**测试步骤：**
1. 走查 §16：两组规模×np 数据+口径说明+迭代数一致性印证+无失实宣称（单机无收益诚实记录）

**期望结果：**
- 四配置数据齐；口径说明含「单机下界」声明；无倍数宣称

**实际结果：** PASS——§16 含 256²/512²×np1/4 四组数据（迭代数精确相同+wall 对照）+口径注（pcg/cg 每 iteration 组成）+「单机下界」诚实记录（+0.7%~+19% 无收益，不外推跨节点）+归约计数佐证；无倍数宣称

**状态：** PASS

---

### ST-014：全量 ctest 双模式全绿（含 ASan/UBSan）

**关联需求：** srs §3.5 R5 验收①
**测试类型：** 回归测试（全量）
**优先级：** High

**测试步骤：**
1. Release 全量 ctest（60 条）
2. Debug（ASan+UBSan，detect_leaks=0）全量 ctest

**期望结果：**
- 双模式全绿；ASan/UBSan 零报告；0 新增警告（构建输出核对）

**实际结果：** PASS——Release 60/60（26.97 s）+Debug（ASan+UBSan，detect_leaks=0）60/60（94.04 s）零报告；构建输出 4 个 warning 均既有（mg_two_level_solver.cpp unused iB/jB×2+exception.hpp -Wterminate×3，皆非 AR009 触碰文件）

**状态：** PASS

---

### ST-015：文档走查 + --help 实跑

**关联需求：** srs §3.5 R5 验收②
**测试类型：** 文档验收
**优先级：** Medium

**测试步骤：**
1. `hypos --help` 直跑：solver 行含 pcg（pipelined 措辞）+--residual-history 行
2. 走查 README/GUIDE/AGENT_SPEC/PERFORMANCE 四文档新增段落

**期望结果：**
- help 两行在；文档无失实宣称（B2b「2 阻塞→1 非阻塞」实测口径、D1 无真实重叠声明在案）

**实际结果：** PASS——`hypos --help` 实跑：`--solver` 行含 `pcg (pipelined …)`、`--residual-history <path>` 行均在；四文档（README 特性/CLI 表、GUIDE B2b/E2/P2/路线图、AGENT_SPEC 类层次、PERFORMANCE §16）review 轮已逐项核，无失实宣称

**状态：** PASS

## 执行摘要

| 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|
| 15 | 15 | 0 | 0 |

## ST 执行报告

| 字段 | 内容 |
|------|------|
| 执行日期 | 2026-10-10 |
| 执行结果 | PASS |
| 执行轮次 | 第 1 轮 |

### 需求覆盖矩阵

| 需求 ID | 需求描述 | 测试用例 | 结果 |
|--------|---------|---------|------|
| §3.1 R1① | pcg/cg 迭代数 ±5%+真残差 ≤ tol | ST-001 | PASS |
| §3.1 R1② | np4 解一致 | ST-002 | PASS |
| §3.1 R1③ | 既有行为零变化 | ST-003 | PASS |
| §3.1 R1④ | 阻塞归约计数 0 | ST-004 | PASS |
| §3.1 R1 边界/防御 | breakdown/零迭代/截断/续态 | ST-005 | PASS |
| §3.2 R2① | cg CSV 契约 | ST-006 | PASS |
| §3.2 R2② | 六 solver CSV 语义一致 | ST-007 | PASS |
| §3.2 R2③ | np4 单文件+内部一致性 | ST-008 | PASS |
| §3.2 R2 异常 | badpath HYPOS_ERROR | ST-009 | PASS |
| §3.2 R2 interval | k 步采样 | ST-010 | PASS |
| §3.2 R2 组合 | 双通道并存 | ST-011 | PASS |
| §3.3 R3 | 绘图脚本（代偿口径） | ST-012 | PASS |
| §3.4 R4 | §16 数据与口径 | ST-013 | PASS |
| §3.5 R5① | 双模式全量绿 | ST-014 | PASS |
| §3.5 R5② | 文档+help 走查 | ST-015 | PASS |

**需求覆盖率：** 5 / 5（100%，R1-R5 全验收句对位）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 8 | 8 | 0 | 0 |
| 边界条件 | 1 | 1 | 0 | 0 |
| 异常处理 | 1 | 1 | 0 | 0 |
| 回归测试 | 3 | 3 | 0 | 0 |
| np 集成 | 2 | 2 | 0 | 0 |
| 文档验收 | 2 | 2 | 0 | 0 |
| **合计** | **15**（承载：Release ctest 60/60+Debug ASan/UBSan 60/60+bench 16 runs+help 直跑+ast 自检） | **15** | **0** | **0** |

### 遗留问题

| 严重性 | 描述 | 处理方式 |
|-------|------|---------|
| Minor | PNG 绘图运行演示不可行（WSL 无 python3/matplotlib）——S1/W3 代偿口径（脚本资产+ast 语法自检+CSV 契约测试为验收主体） | 已记录（srs §5 假设+tasks.md T004），有 python 环境时可复现 |
| Minor | 单机 wall 无 pcg 收益（+0.7%~+19%）——D1 诚实声明预期内（无真实重叠），跨节点未测 | 已如实入档（PERFORMANCE §16「单机下界」+GUIDE B2b），CG-2 两步前瞻列后续 AR |

### 结论

> **Go**——需求覆盖 100%（R1-R5 全验收句对位）；15/15 用例 PASS（双模式 ctest 120/120 全绿含 ASan/UBSan 零报告+回归零破坏）；2 项 Minor 均已如实记录不阻塞（环境受限代偿口径在案+诚实性能记录）；核心指标亮点：pcg/cg 迭代数四配置精确相同（零差异）+阻塞归约 0 断言带阳性对照自证。
