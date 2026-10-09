# [AR006] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR006（halo-comm-paradigms：C1 派生数据类型直传 + C2 真集合 halo） |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-09 |
| 执行环境 | WSL 单机 OpenMPI 4.1.6；Release=/root/build-release，Debug(ASan+UBSan)=/root/build-debug；OMP=1；证据见 ./evidence/ |

## 测试用例列表

### ST-001：datatype 路径复跑 halo 正确性矩阵全绿

**关联需求：** srs §3.1 R1 验收①（GUIDE C1 验收①「不变绿转」）
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 既有 halo 正确性测试矩阵（2D 方向语义/多 halo/PROC_NULL/np4 非均匀、3D 方向语义/np8）与新增 datatype 用例

**测试步骤：**
1. ctest release + debug 双模式全量（含 datatype_np4/np8 点名条目与 unit 全量条目）

**期望结果：**
- Then 全部通过（U1-U4/B1/B2 与既有矩阵不变绿转）

**实际结果：** release 35/35、debug 35/35（evidence/st-ctest-release.log、st-ctest-debug.log）；datatype_np4/np8 各自 PASS

**状态：** PASS

---

### ST-002：pack 与 datatype 路径 halo 带等价

**关联需求：** srs §3.1 R1 验收②
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 同一 Subgrid 与数据场（3D hw=2）

**测试步骤：**
1. 运行 `HaloExchangeTest.DatatypeEquivalentToPack`（U4：P2P vs Datatype 双 exchanger 同场 halo 带逐元素对比）

**期望结果：**
- Then halo 带内容一致

**实际结果：** PASS（ctest unit 条目内，两模式均过）

**状态：** PASS

---

### ST-003：np=1 全 PROC_NULL 安全（三 exchanger）

**关联需求：** srs §3.1 R1 验收③ + §3.2 R2 期望（np=1 空图）
**测试类型：** 边界条件
**优先级：** High

**前置条件：**
- Given np=1（全 PROC_NULL 边界）

**测试步骤：**
1. 运行 `HaloExchangeTest.NoNeighborNp1Safe`（E1：P2P/Datatype/Collective 参数化，exchange + begin/end 两场次快照对比）

**期望结果：**
- Then 无 MPI 通信调用且行为安全（halo 不变）；collective 空图路径（nEdges=0）被触达且安全

**实际结果：** PASS（release/debug 单跑 + ctest 全量均过；实现审查第 2 轮确认空图全链路真实执行）

**状态：** PASS

---

### ST-004：collective 模式 np=1/4/8 全绿且无委托 WARN

**关联需求：** srs §3.2 R2 验收①（GUIDE C2 验收①）
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given collective 模式，np=1/4/8（含非均匀切分、3D）

**测试步骤：**
1. ctest collective_np4/collective_np8（np4 非均匀 2D、np8 3D 点名条目）
2. np=1 直跑 collective（E1 覆盖空图）
3. np=4 直跑 `hypos --comm-mode collective` e2e 并 grep WARN/委托字样

**期望结果：**
- Then 全绿且无「falls back to PointToPointExchanger」WARN

**实际结果：** collective_np4/np8 PASS；e2e 输出零 WARN 匹配（evidence/st-collective-nowarn.log，NO-WARN）；委托 WARN 代码已物理删除（collective_exchanger.cpp 重写，diff 可见）

**状态：** PASS

---

### ST-005：委托测试改写非静默删除

**关联需求：** srs §3.2 R2 验收②
**测试类型：** 回归测试
**优先级：** Medium

**前置条件：**
- Given 旧 `CollectiveExchangerDelegatesToP2P` 测试

**测试步骤：**
1. 核对新用例 `CollectiveExchangerIsTrueCollective`（B5）存在且运行
2. 核对 git 历史（提交 2b92982）记录改写而非删除

**期望结果：**
- Then 覆盖升级为真集合语义断言（name 契约 + 自环 p2p 等价 + begin/end 拆分，零日志字符串断言）

**实际结果：** B5 PASS（两模式）；提交 2b92982 diff 可见旧用例→新用例改写

**状态：** PASS

---

### ST-006：`--comm-mode p2p` 既有行为零变化

**关联需求：** srs §3.2 R2 验收③
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given 基线提交 8d40c21（AR006 起点）

**测试步骤：**
1. `git diff --stat 8d40c21..HEAD -- src/comm/p2p_exchanger.cpp`（实现零改动）
2. 既有 p2p 路径测试全绿（halo_2d_mpi/halo_3d_mpi/solver/overlap 位级守护）

**期望结果：**
- Then p2p 实现零 diff；既有测试不变绿转

**实际结果：** p2p_exchanger.cpp 零 diff（evidence/st-p2p-zerodiff.txt，PowerShell git 核验空输出）；ctest 35/35 含全部既有条目（overlap 位级一致守护过）

**状态：** PASS

---

### ST-007：两路径同一 ctest 矩阵全绿

**关联需求：** srs §3.3 R3 验收①
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given 切换机制（--comm-mode 三值）与 ctest 矩阵

**测试步骤：**
1. ctest 全量（35 条目：既有 30 + datatype_np4/np8、collective_np4/np8、datatype_e2e）

**期望结果：**
- Then 双模式全绿

**实际结果：** release 35/35、debug 35/35（evidence/st-ctest-*.log）；datatype_e2e（np4 直跑 hypos --comm-mode datatype）rc=0（evidence/st-datatype-e2e.log）

**状态：** PASS

---

### ST-008：PERFORMANCE 两组对比数据 + 口径说明（记录式）

**关联需求：** srs §3.3 R3 验收②
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given PERFORMANCE §13

**测试步骤：**
1. 复测 `scripts/bench_ar006.sh`（200 迭代口径）
2. 稳定性口径：2000 迭代 ×3 中位（np4 三模式）
3. 文档走查：数据、口径、结论一致性

**期望结果：**
- Then 含 pack vs datatype 与 p2p vs collective 两组对比数据及口径说明；不预设结论、无失实宣称

**实际结果：** 三轮实测（首测/ST 复测/2000 迭代中位）排序不稳定——范式间差异与单机漂移同量级；**§13 已按实测修订为「同量级、无稳定优劣、不下结论」**，撤回首测「2.34× 劣化」强结论（判定为噪声），collective 空图 np1 固定开销两轮一致稳定可测（evidence/st-bench-rerun.log 含三轮原始数据）

**状态：** PASS

---

### ST-009：文档走查无失实宣称 + `--help` 一致

**关联需求：** srs §3.4 R4 验收①②
**测试类型：** 正常路径
**优先级：** Medium

**前置条件：**
- Given 文档集（README/GUIDE/AGENT_SPEC/PERFORMANCE/halo_exchanger.hpp 注释）

**测试步骤：**
1. grep「委托/roadmap item/falls back」残留
2. 实跑 `hypos --help` 核对 comm-mode 说明
3. 核对 README 路线图条目移除、GUIDE P7/G2 勾选、AGENT_SPEC 类层级

**期望结果：**
- Then 委托旧表述全部更新；--help 与实际行为一致；无失实宣称

**实际结果：** `--help` 输出三值 + datatype 注记（evidence/st-help.txt）；README 路线图真集合条目已移除（含委托注记）、:69 comm-mode 表更新；GUIDE P7/G2 ✅ 带 AR006 证据、§5 ★3 已完成；AGENT_SPEC 三子类层级；hpp 过时注释已随 T004 改写；实现审查第 1 轮已核文档清单全落地

**状态：** PASS

---

### ST-010：NFR——ASan/UBSan 零报告 + 资源释放 + 双模式零新增警告

**关联需求：** srs §4 非功能需求（内存安全/资源/构建/回归）
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given Debug 构建（ASan+UBSan）全量测试

**测试步骤：**
1. ctest debug 全量（35 条目，含三 exchanger 生命周期场景：重复 initialize/析构在 MPI_Finalized 前/在途请求防御均有用例或代码路径触达）
2. 构建日志核新增警告

**期望结果：**
- Then 零 ASan/UBSan 报告；零新增编译警告

**实际结果：** debug 35/35 零报告（evidence/st-ctest-debug.log）；release/debug 构建均零警告输出

**状态：** PASS

---

## ST 执行报告

| 字段 | 内容 |
|------|------|
| 执行日期 | 2026-10-09 |
| 执行结果 | PASS |
| 执行轮次 | 第 1 轮 |

### 需求覆盖矩阵

| 需求 ID | 需求描述 | 测试用例 | 结果 |
|--------|---------|---------|------|
| §3.1 R1① | datatype 矩阵不变绿转 | ST-001 | PASS |
| §3.1 R1② | pack vs datatype 等价 | ST-002 | PASS |
| §3.1 R1③ | np=1 安全 | ST-003 | PASS |
| §3.2 R2① | collective np=1/4/8 全绿无 WARN | ST-003, ST-004 | PASS |
| §3.2 R2② | 委托测试改写非静默删除 | ST-005 | PASS |
| §3.2 R2③ | p2p 零变化 | ST-006 | PASS |
| §3.3 R3① | 两路径同矩阵全绿 | ST-007 | PASS |
| §3.3 R3② | PERFORMANCE 数据+口径 | ST-008 | PASS |
| §3.4 R4①② | 文档无失实 + help 一致 | ST-009 | PASS |
| §4 NFR | ASan/UBSan/资源/警告/回归 | ST-010 | PASS |

**需求覆盖率：** 10 / 10（100%）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 6 | 6 | 0 | 0 |
| 边界条件 | 1 | 1 | 0 | 0 |
| 回归测试 | 3 | 3 | 0 | 0 |
| **合计** | **10** | **10** | **0** | **0** |

### 测量发现与处置（ST-008 专项记录）

- **发现**：bench 首测（T006）collective 2.34× 劣化触发 D6 评审条款；ST 复测同口径排序反转（0.87×）；2000 迭代 ×3 中位为 1.14×，rep 间漂移 ±30%。
- **处置**：按「无失实宣称」底线修订 PERFORMANCE §13——三轮数据全部入档、撤回强结论、判定 2× 差异为单机噪声、维持「同量级无稳定优劣」记录式结论；默认 p2p 未动，无回退动作。srs R3「不预设结论」的记录式判定设计经实测证明必要性。

### 遗留问题

| 严重性 | 描述 | 处理方式 |
|-------|------|---------|
| Minor | PERFORMANCE §13 结论为「无稳定优劣」——更大面/更大 halo 宽度/多节点收益空间未测 | 后续 AR（如 C3 RMA 或性能专项）可扩展测量口径 |
| Minor | 行数预算 1245 超 GUIDE §6.4 ~800 建议 | 已按「论证」路径记录于 design §8（测试占 59% + 共享头抽取 + 验收覆盖要求） |

### 结论

> **Go**。10/10 用例 PASS、需求覆盖 100%、无 Critical/Major 缺陷、双模式 35/35 零报告、p2p 基线零改动。测量类发现（bench 噪声）已按诚实优先原则处置并入档。用户已授权自动化 ST（无需人工确认），执行归档。
