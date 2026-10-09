# [AR005] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR005 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-08 |

## 测试用例列表

### ST-001：已知 u 场 .vti 结构规范 + 独立解析器还原

**关联需求：** srs §3.1 R1（验收①）
**测试类型：** 正常路径
**优先级：** High

**前置条件：** 4×4 已知 u 场（毒化 halo 后 interior 编码赋值）

**测试步骤：**
1. 运行 `AltFeatureTest.VtkAppendedBinaryParsesBack`（U3，独立解析器：定位 AppendedData/`_`/解码 UInt64 长度头/逐位 memcmp/逐 double 还原）

**期望结果：** XML 结构含 header_type="UInt64"/format="appended"/offset="0"/AppendedData raw；长度头==128B；16 单元位级还原一致

**实际结果：** PASS（Release 构建运行 1 test PASSED；断言含结构①-⑤全项）

**状态：** PASS

---

### ST-002：np=4 拼装 .pvti 索引与 piece extents 正确

**关联需求：** srs §3.1 R1（验收②）
**测试类型：** 正常路径（端到端）
**优先级：** High

**前置条件：** WSL OpenMPI，np=4

**测试步骤：**
1. `mpirun -n 4 hypos --nx 32 --ny 32 --max-iter 5 --output-format vtk --output-dir /root/st-vtk-np4`
2. 检查产出文件集、pvti Source 引用、Extent 拼装、每 piece appended 标记

**期望结果：** 退出码 0；solution_5_r{0..3}.vti + solution_5.pvti；索引引用 4 piece；4 个 extent 无缝拼满 0-31²

**实际结果：** PASS（evidence/st-vtk-np4.log：exit 0、Source 计数 4、WholeExtent="0 31 0 31 0 0" + 4 象限 extent 无重叠拼满、每 piece 2 个 appended 标记；另 VtkPieceAndParallelIndexFiles 单测断言结构守护 GREEN）

**状态：** PASS

---

### ST-003：256²/512² vtk binary vs ASCII 体积耗时对比入 PERFORMANCE

**关联需求：** srs §3.1 R1（验收③）
**测试类型：** 正常路径（性能数据）
**优先级：** High

**测试步骤：**
1. `scripts/bench_ar005.sh baseline|after`（3 次中位数）
2. 核对 PERFORMANCE §12 与 evidence/（baseline|after）-io-bench.md

**期望结果：** 两组数据 + 口径说明入 PERFORMANCE

**实际结果：** PASS（PERFORMANCE.md §12：vtk 256²=36.2×/512²=48.7×、体积 -57.9%；evidence 双文件含 D6 门槛判定表）

**状态：** PASS

---

### ST-004：块写与逐元素版输出逐位一致

**关联需求：** srs §3.2 R2（验收①）
**测试类型：** 正常路径（正确性）
**优先级：** High

**测试步骤：**
1. 运行 `AltFeatureTest.BinaryOutputBitIdenticalToExpectedBytes`（U2，从零构造期望字节——比旧实现对比更硬）

**期望结果：** 56B 头 7 字段正确 + 16×8B 数据逐位一致；halo sentinel 不混入

**实际结果：** PASS（Release/Debug 双模式 GREEN；块写落地后字节布局逐位不变——after 基准字节数 524,344/2,097,208 与 baseline 逐位相同佐证）

**状态：** PASS

---

### ST-005：512² 块写前后耗时对比入 PERFORMANCE

**关联需求：** srs §3.2 R2（验收②）
**测试类型：** 正常路径（性能数据）
**优先级：** Medium

**期望结果/结果：** PASS（PERFORMANCE §12：binary 256²=3.35×/512²=3.79×，字节数不变；数据含 baseline 与 after 双侧）

**状态：** PASS

---

### ST-006：WSL OpenMPI np=4 全部 MPI ctest 探针全绿

**关联需求：** srs §3.3 R3（验收①）
**测试类型：** 回归
**优先级：** High

**测试步骤：**
1. `ctest`（Release + Debug 双模式全量）

**期望结果：** 全部条目 PASS（含 11 个 np>1 条目，其中 7 gtest 条目 filter 含 MpiEnvTest.* 消费探针）

**实际结果：** PASS（Release 30/30 + Debug(ASan+UBSan) 30/30；evidence/st-ctest-{release,debug}.log；探针不引入误报——四态中「未设」为 SKIP 不 FAIL）

**状态：** PASS

---

### ST-007：探针咬住单进程退化（Red 语义验证）

**关联需求：** srs §3.3 R3（验收②）
**测试类型：** 异常处理
**优先级：** High

**测试步骤：**
1. 单进程运行探针二进制（HYPOS_EXPECT_NP=4）——模拟 mpich 4.2.0 静默退化语义
2. 记录退出码

**期望结果：** FAIL 且退出码非 0（宁可红不可假绿）

**实际结果：** PASS（evidence/T001-probe-red-rc.log：Red-1 exit code 1，诊断含实际 size 与退化提示；Red-2 abc 非法值 fail-safe exit 1；Green np=1 exit 0——退出码传播链路闭合）

**状态：** PASS

---

### ST-008：CI MPICH 路径在修复后平台全绿（本地等价口径）

**关联需求：** srs §3.3 R3（验收③）+ §3.4 R4（验收①）
**测试类型：** 集成（环境受限替代口径）
**优先级：** High

**前置条件：** gitee 主仓库不跑 Actions（srs §5 约束已声明）

**测试步骤：**
1. 核对 ci.yml runs-on=ubuntu-22.04（脱离 24.04 mpich 4.2.0 缺陷平台）
2. 本地等价验证：Debug 全量 ctest（ASan+UBSan，CI Debug 段同款调试语义）+ 探针绿

**期望结果：** 配置正确 + 等价命令证据 + 口径记录

**实际结果：** PASS（ci.yml:14-17 ubuntu-22.04 + D4 注释；本地 Debug 30/30 含探针；口径记录于 srs §5 与 tasks.md T005 备注；无回退必要）

**状态：** PASS

---

### ST-009：文档与实际状态一致（无失实宣称）

**关联需求：** srs §3.4 R4（验收②）+ §3.5 R5（验收①）
**测试类型：** 走查
**优先级：** Medium

**测试步骤：**
1. 走查 README:28（mpich 缺陷记录 + 探针机制注记）、--output-format 描述、GUIDE P5/D1/F1 勾选、PERFORMANCE §12

**期望结果：** 宣称与实现/测试零缺口

**实际结果：** PASS（README mpich 段补探针注记属实（ctest ENV 注入+FAIL 语义）；--output-format vtk 描述改为 appended raw 与实现一致；GUIDE P5 ✅/D1 ✅/F1 探针+runner 完成标注与实际对应；审查第 2 轮独立核验文档项）

**状态：** PASS

---

### ST-010：PERFORMANCE 含数据及口径说明

**关联需求：** srs §3.5 R5（验收②）
**测试类型：** 走查
**优先级：** Medium

**实际结果：** PASS（§12 含 vtk/binary 四组合 baseline/after 数据表、io_write 剖面口径、save-interval 说明、D6 门槛判定、复测脚本指引入 scripts/bench_ar005.sh）

**状态：** PASS

---

### ST-011：ParaView 兼容口径（替代验收）

**关联需求：** srs §4 NFR 兼容性（读取方）
**测试类型：** 走查 + 解析器测试（MANUAL 替代）
**优先级：** Medium

**测试步骤：**
1. U3 位级解析器测试（独立最小解析器）
2. VTK XML appended raw 规范要点走查（header_type/byte_order/offset/`_` 标记/长度头宽度/行主序点序/pvti 结构）

**期望结果：** 规范合规 + 局限性声明记录

**实际结果：** PASS（U3 GREEN；走查覆盖已知易错点；局限性声明——未做真实 ParaView GUI 加载演示——如实记录于 evidence/after-io-bench.md「ParaView 可视化验证口径」节）

**状态：** PASS

---

### ST-012：NFR 验证——WARN 语义 / ASan+UBSan / 警告 / 行数预算

**关联需求：** srs §4 NFR（兼容性/内存安全/构建/复杂度）
**测试类型：** 回归 + 走查
**优先级：** High

**测试步骤：**
1. WARN 路径：`AltFeatureTest.{Vtk,Binary}WriteToInvalidPathWarnsAndSkips`（含恢复写断言）
2. ASan/UBSan：Debug 全量 ctest
3. 构建警告：双模式构建日志
4. 行数预算：`git diff 5cf1ed7..HEAD --shortstat`

**期望结果：** 不抛不崩 + 零报告 + 0 新增警告 + 净新增 ≤~800 行

**实际结果：** PASS（两用例 GREEN；Debug 30/30 零报告；双模式构建 0 error 0 warning；净新增 793-49=744 行（含 specs/evidence 文档，纯代码更少）≤~800）

**状态：** PASS

---

## ST 执行报告

| 字段 | 内容 |
|------|------|
| 执行日期 | 2026-10-08 |
| 执行结果 | PASS |
| 执行轮次 | 第 1 轮 |

### 需求覆盖矩阵

| 需求 ID | 需求描述 | 测试用例 | 结果 |
|--------|---------|---------|------|
| §3.1 | R1 VTK appended binary .vti | ST-001, ST-002, ST-003 | PASS |
| §3.2 | R2 BinaryIOBackend 块写 | ST-004, ST-005 | PASS |
| §3.3 | R3 CI mpich 假通过探针 | ST-006, ST-007, ST-008 | PASS |
| §3.4 | R4 CI 平台修复 | ST-008, ST-009 | PASS |
| §3.5 | R5 文档同步 | ST-009, ST-010 | PASS |
| §4 | NFR（正确性/兼容/内存安全/构建/性能/回归/复杂度） | ST-004, ST-011, ST-012, ST-003, ST-006 | PASS |

**需求覆盖率：** 6 / 6（100%，R1-R5 + NFR 全表）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 7 | 7 | 0 | 0 |
| 边界条件 | 1（ST-002 np=4 拼装） | 1 | 0 | 0 |
| 异常处理 | 2（ST-007 退化/非法值） | 2 | 0 | 0 |
| 回归测试 | 2（ST-006/ST-012 双模式全量） | 2 | 0 | 0 |
| **合计** | **12** | **12** | **0** | **0** |

自动化佐证：ctest Release 30/30 + Debug(ASan+UBSan) 30/30（evidence/st-ctest-{release,debug}.log）；np=4 vtk 拼装 e2e（evidence/st-vtk-np4.log）；探针退出码四态（evidence/T001-probe-red-rc.log）。

### 遗留问题

| 严重性 | 描述 | 处理方式 |
|-------|------|---------|
| Minor | test_alt_solvers.cpp 职责膨胀（~1140 行），AR005 IO 用例宜拆分至 test_io_layout.cpp | 已记入 GUIDE §5 AR008+ 工程卫生名单，延后处理 |
| Minor | ParaView 真实 GUI 加载演示未做（环境无 GUI） | srs §4 预授权替代口径已采用，局限性声明入 evidence |

### 结论

> **建议 Go**：12/12 用例 PASS，需求覆盖 100%，无 Critical/Major 缺陷；2 个 Minor 均已记录并确认延后口径（其一为环境限制的预授权替代，其二为跨 AR 工程卫生项）。全量回归双模式 30/30 零警告零 Sanitizer 报告。
