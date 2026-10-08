# [AR003] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR003 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-08 |

> 执行环境：WSL（hypos 发行版，OpenMPI 4.1.x），构建目录 `/root/build-debug`、`/root/build-release`。
> 复用策略：开发阶段已固化的 ctest（16 项）作为自动化验收主体；srs 字面参数场景（如 §3.2-1 的 50/120）与跨格式等价性由 ST 专项手工执行补齐（`scripts/st_ar003.sh`），全部留 evidence。

## 测试用例列表

### ST-001：np=4 单文件产物与内容正确性

**关联需求：** srs.md §3.1-1
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given np=4 运行 `MPIIOBinaryBackend::write` 完成（4×4 子域拼接，制造解参考场）

**测试步骤：**
1. 运行 gtest `MpiIoOutputTest.FourRanksGlobalAssemblyMatchesReference`（mpiio_np4）

**期望结果：**
- Then 恰好 1 个 `solution_<step>.bin`；头部字段与运行参数一致；数据区与全局参考场逐点一致（≤1e-15）

**实际结果：** PASS——4 rank 实跑（OpenMPI，非跳过），全 rank 断言通过（evidence/T-st-ctest-debug.log #11；FourRanks 实跑证据 evidence/T-fourranks-np4.log）

**状态：** PASS

---

### ST-002：np=1 逐位一致

**关联需求：** srs.md §3.1-2
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given np=1、16×16 制造解

**测试步骤：**
1. 运行 gtest `MpiIoOutputTest.SingleRankHeaderAndDataMatchReference`（unit / mpiio_np1）

**期望结果：**
- Then 单文件数据区与串行参考场逐位一致（256 double）

**实际结果：** PASS——头部字段 + 256 double 逐位断言通过（T-st-ctest-debug.log #1/#10）

**状态：** PASS

---

### ST-003：np=1 mpibin 与 binary 分片数据区语义等价

**关联需求：** srs.md §3.1-2（后半句）
**测试类型：** 正常路径（跨格式等价）
**优先级：** Medium

**前置条件：**
- Given 同配置 np=1 分别以 `--output-format mpibin` 与 `--output-format binary` 运行

**测试步骤：**
1. WSL 手工执行两种格式各一次（32²，max-iter 5，`scripts/st_ar003.sh`）
2. `cmp` 对比 mpibin 数据区（跳过 72B 头）与 binary 分片数据区（跳过 56B 头）

**期望结果：**
- Then 两数据区字节级相同（同一全局行主序场）

**实际结果：** PASS——`ST-003: data areas byte-identical (32*32*8 = 8192 bytes)`（evidence/T-st-manual.log）

**状态：** PASS

---

### ST-004：MPI-IO 打开失败告警跳过

**关联需求：** srs.md §3.1-3
**测试类型：** 异常处理
**优先级：** High

**前置条件：**
- Given 输出目录不存在（`/nonexistent_dir_xyz/`）

**测试步骤：**
1. 运行 gtest `MpiIoOutputTest.WriteToInvalidPathWarnsAndSkips`（unit/mpiio_np1/mpiio_np4）

**期望结果：**
- Then 不抛异常不崩溃；恢复路径写正常产出文件且内容正确

**实际结果：** PASS——三种 np 下均通过，恢复写读回断言全过（T-st-ctest-debug.log #1/#10/#11）

**状态：** PASS

---

### ST-005：save-interval 中间单文件（srs 字面场景 50/120）

**关联需求：** srs.md §3.2-1
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given `--output-format mpibin --save-interval 50 --max-iter 120` np=4

**测试步骤：**
1. WSL 手工运行 hypos（32²，默认 tol 1e-6 下 120 迭代不提前收敛，`scripts/st_ar003.sh`）
2. 列出输出目录文件集

**期望结果：**
- Then 存在 `solution_50.bin`、`solution_100.bin` 与最终解 `solution_120.bin`，全部为单文件（无 `_r` 分片、无索引）

**实际结果：** PASS——exit=0，文件集恰为 `solution_50.bin; solution_100.bin; solution_120.bin`（evidence/T-st-manual.log）

**状态：** PASS

---

### ST-006：未传 save-interval 仅最终解

**关联需求：** srs.md §3.2-2
**测试类型：** 边界条件
**优先级：** High

**前置条件：**
- Given `--output-format mpibin` 且不传 `--save-interval`，np=4

**测试步骤：**
1. 运行 ctest `mpiio_e2e_final` + `mpiio_e2e_final_verify`（VerifyMpiioOutput.cmake 精确文件集断言）

**期望结果：**
- Then 输出目录仅有 `solution_5.bin` 一个文件

**实际结果：** PASS——`mpiio_e2e_final` 与 `mpiio_e2e_final_verify` 均通过（T-st-ctest-debug.log #15/#16）

**状态：** PASS

---

### ST-007：既有格式回归不变

**关联需求：** srs.md §3.2-3
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given AR002 基线（binary/vtk/json/csv 分片路径）未被破坏

**测试步骤：**
1. 全量 ctest（含 collective_mpi、ranks3_smoke、solver_mpi 等 binary/默认格式路径）

**期望结果：**
- Then 全部通过，与 AR002 基线一致

**实际结果：** PASS——Debug 16/16（49.66s）+ Release 16/16（12.96s），AR002 既有 12 项全部通过（T-st-ctest-debug.log、T-st-ctest-release.log）

**状态：** PASS

---

### ST-008：非法格式值退出码 1

**关联需求：** srs.md §3.2（期望行为末条）
**测试类型：** 异常处理
**优先级：** High

**前置条件：**
- Given `--output-format mpibinn`（非法值）

**测试步骤：**
1. 运行 ctest `mpibin_bad_format`（WILL_FAIL）

**期望结果：**
- Then 进程以非零退出码结束

**实际结果：** PASS——ctest #3 mpibin_bad_format Passed（WILL_FAIL 语义：非零退出被判定为预期）

**状态：** PASS

---

### ST-009：np=4 失败路径无挂死、告警有界

**关联需求：** srs.md §3.3-1
**测试类型：** 异常处理
**优先级：** High

**前置条件：**
- Given np=4、输出目录不可写

**测试步骤：**
1. 运行 ctest `mpiio_np4`（X-1 在 4 rank 下执行，COMM_SELF 独立失败面）
2. 走查 `warned_` 闩锁与 collective 一致性（evidence/T004-leak-walkthrough.md）

**期望结果：**
- Then 运行完成不挂死不崩溃；告警限流为代码走查结论（X-1 无日志捕获能力，已在 T004 以走查口径记录）

**实际结果：** PASS——mpiio_np4 0.4s 完成（无挂死）；走查确认 `warned_` 闩锁每 rank 仅首次告警、open 失败全 rank 一致提前返回（T004-leak-walkthrough.md）

**状态：** PASS

---

### ST-010：无 MPI_File 泄漏

**关联需求：** srs.md §3.3-2
**测试类型：** 代码走查（MANUAL）
**优先级：** Medium

**前置条件：**
- Given ASan 启用但 WSL root 下 LeakSanitizer 关闭（detect_leaks=0）

**测试步骤：**
1. 走查 evidence/T004-leak-walkthrough.md（close 无条件、Type_free 无条件、open 失败提前点无资源）

**期望结果：**
- Then 无泄漏路径；全 rank collective 对齐

**实际结果：** PASS——唯一提前返回点在 open 失败（无资源存在）；其余路径 `MPI_Type_free` ×2 + `MPI_File_close` 无条件执行；审查代理独立复核（tasks.md 审查记录第 2 轮）

**状态：** PASS

---

### ST-011：新测试 np=1/np=4 执行通过

**关联需求：** srs.md §3.4-1
**测试类型：** 正常路径
**优先级：** High

**前置条件：**
- Given test_mpiio_output.cpp 7 用例

**测试步骤：**
1. 运行 ctest `unit`、`mpiio_np1`、`mpiio_np4`

**期望结果：**
- Then 全部断言通过（np=1 运行 5 用例，np=4 跑 4-rank 用例 + COMM_SELF 用例）

**实际结果：** PASS——unit 7.46s、mpiio_np1 0.36s、mpiio_np4 0.42s 全部通过（T-st-ctest-debug.log）

**状态：** PASS

---

### ST-012：全量回归 Debug/Release + ASan 零报告

**关联需求：** srs.md §3.4-2
**测试类型：** 回归测试
**优先级：** High

**前置条件：**
- Given Debug 构建（-fsanitize=address,undefined）与 Release 构建

**测试步骤：**
1. `ctest --test-dir /root/build-debug`（16 项）
2. `ctest --test-dir /root/build-release`（16 项）

**期望结果：**
- Then 全绿、零 sanitizer 报告、0 新增警告

**实际结果：** PASS——Debug 16/16（49.66s）、Release 16/16（12.96s）；构建零警告零错误（会话记录）；ASan/UBSan 零报告（全部测试正常退出）

**状态：** PASS

---

### ST-013：3D 布局冒烟（64×64×8 np=4）

**关联需求：** srs.md §3.4-3
**测试类型：** 边界条件
**优先级：** Medium

**前置条件：**
- Given 64×64×8 全局网格、np=4

**测试步骤：**
1. 运行 gtest `MpiIoOutputTest.FourRanks3DLayout`（mpiio_np4）

**期望结果：**
- Then 单文件数据区尺寸 = 全局内点数 × 8 字节（3D 行主序布局正确）

**实际结果：** PASS——4 rank 实跑通过（T-fourranks-np4.log）

**状态：** PASS

---

### ST-014：基准参考数据入档

**关联需求：** srs.md §3.5
**测试类型：** 文档/数据走查
**优先级：** Medium

**前置条件：**
- Given T007 基准已执行（256²，np=1/4，OMP=1，3 次中位数）

**测试步骤：**
1. 核对 docs/PERFORMANCE.md §10 含 mpibin vs binary 对比表、环境与方法说明
2. 核对 evidence/T007-bench.log 原始数据（12 条无 MISSING）

**期望结果：**
- Then 表格存在、口径注明（无失实宣称）

**实际结果：** PASS——PERFORMANCE.md §10（:253 起）含对比表 + 取数口径 + 小文件固定开销主导的诚实说明；T007-bench.log 12 条原始数据完整

**状态：** PASS

---

### ST-015：文档全套同步无失实宣称

**关联需求：** srs.md §3.6
**测试类型：** 文档走查
**优先级：** Medium

**前置条件：**
- Given README/DESIGN/PERFORMANCE 已更新

**测试步骤：**
1. README：格式表含 `mpibin`（:72）；路线图节（:259-265）无「MPI-IO 单文件」残留（仅列 RMA/真集合/HDF5/PAPI）
2. DESIGN：IO 层含 `MPIIOBinaryBackend` 描述（:75，含布局/collective 语义/与分片取舍）
3. 全文走查无失实宣称（审查代理已独立走查通过）

**期望结果：**
- Then 三份文档与实现一致

**实际结果：** PASS——三处核查全部命中，无失实宣称（review 第 2 轮 PASS 佐证）

**状态：** PASS

---

## 执行摘要

| 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|
| 15 | 15 | 0 | 0 |

## ST 执行报告

| 字段 | 内容 |
|------|------|
| 执行日期 | 2026-10-08 |
| 执行结果 | PASS |
| 执行轮次 | 第 1 轮 |

### 需求覆盖矩阵

| 需求 ID | 需求描述 | 测试用例 | 结果 |
|--------|---------|---------|------|
| §3.1 | MPI-IO 单文件二进制输出后端 | ST-001, ST-002, ST-003, ST-004 | PASS |
| §3.2 | CLI 集成与中间保存 | ST-005, ST-006, ST-007, ST-008 | PASS |
| §3.3 | 异常处理与防御性 | ST-009, ST-010 | PASS |
| §3.4 | 验证与测试 | ST-011, ST-012, ST-013 | PASS |
| §3.5 | 基准参考数据 | ST-014 | PASS |
| §3.6 | 文档全套同步 | ST-015 | PASS |

**需求覆盖率：** 6 / 6（100%）

### 测试执行汇总

| 类型 | 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|------|
| 正常路径 | 7 | 7 | 0 | 0 |
| 边界条件 | 2 | 2 | 0 | 0 |
| 异常处理 | 3 | 3 | 0 | 0 |
| 回归测试 | 2 | 2 | 0 | 0 |
| 走查（MANUAL/文档） | 3 | 3 | 0 | 0 |
| **合计** | **15**（去重映射至 16 项 ctest + 2 项手工） | **15** | **0** | **0** |

### 覆盖率说明

项目文档（AGENT_SPEC.md / design.md）未声明数值型代码覆盖率阈值；以需求覆盖率（6/6 = 100%）与 design.md §6.6 追溯矩阵 100% 映射（review 已核验）作为覆盖口径，延续 AR001/AR002 既有 ST 实践。ASan/UBSan 全量零报告提供内存安全佐证。

### 遗留问题

无。

### 结论

> **Go 建议**：15/15 用例通过，需求覆盖 100%，无 Critical/Major/Minor 缺陷，无回归（AR002 既有 12 项测试全绿），NFR（ASan 零报告、零警告、双模式、接口兼容）全部达标。建议归档 AR003。
