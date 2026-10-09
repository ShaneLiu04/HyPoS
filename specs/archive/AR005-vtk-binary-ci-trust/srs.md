# [AR005] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR005 |
| AR 主题 | vtk-binary-ci-trust（I/O 快赢 + CI 可信度） |
| 关联 SR | 无（源自 docs/OPTIMIZATION_GUIDE.md §5 路线 ★2，条目 D1 + F1 紧急项；修复 §3 P5） |
| 日期 | 2026-10-08 |
| 状态 | 已确认（用户授权自动化执行，依据 OPTIMIZATION_GUIDE.md §5 ★2） |

## 1. 背景与目标

OPTIMIZATION_GUIDE.md §3 诊断 P5：VTK 输出为 ASCII Float64（512² 一次输出数百 MB 文本），binary 后端逐元素 `ofstream.write` 非块写。本 AR 落地路线 ★2「I/O 快赢 + CI 可信度」：D1 将 `.vti` 改 appended binary（raw offset 模式）+ `BinaryIOBackend` 改块写，输出体积/耗时预计降 1-2 个数量级；同时完成 F1 **紧急项**——CI 的 Debug(MPICH) job 跑在 ubuntu-latest（24.04），该平台 mpich 4.2.0 存在 PMI/PMIx 不匹配缺陷（AR003 实测确认，README:28 已记录），np>1 测试存在「多 rank 假通过」风险，需以进程数探针修复 CI 数据可信度。

## 2. 需求范围

**In Scope（对应指南条目）：**
- D1-a：`VTKIOBackend::write` 的 `.vti` 从 ASCII 改 appended binary（raw offset 模式：XML 头 `format="appended"` + `offset`，`<AppendedData encoding="raw">` 后跟 UInt64 长度头 + 数据镜像）
- D1-b：`.pvti` 并行索引适配（PDataArray 声明与 piece 一致性；pvti 本身无数据负载，仍为 ASCII）
- D1-c：`BinaryIOBackend::write` 数据区改块写（连续段一次 `ofs.write(ptr, bytes)`，消除逐元素调用）
- F1-紧急：ctest 进程数探针——独立探针测试 + 既有 np>1 MPI 测试激活 `MPI_Comm_size == 预期 np` 断言（假通过即 FAIL）
- F1-修复：CI ci.yml 的 Debug(MPICH) job 脱离缺陷平台（换发行版或 pin 修复版本，design 决定并记录）
- R4 文档同步：README（mpich 记录更新）、PERFORMANCE（ASCII vs binary 体积/耗时对比）、GUIDE 勾选（P5/D1/F1 紧急项）

**Out of Scope：**
- F1 完整矩阵（GCC/Clang × OpenMPI/MPICH、ccache、覆盖率 job）——规模超本 AR，列后续 AR
- D2 VTI 单文件 MPI-IO（依赖 D1，路线 ★后续）
- D3 异步输出双缓冲（指南定位为仅设计注记）
- E3 性能回归门禁（CI 报警比对基线，列后续 AR）

## 3. 功能需求

### 3.1 R1 VTK appended binary `.vti`（D1-a/D1-b）

**描述：** `.vti` piece 的数据负载从 ASCII 文本改为 VTK appended raw binary。

**触发条件：** `--output-format vtk` 任意运行。

**期望行为：**
- XML 头：`<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt64">`；`<DataArray type="Float64" Name="u" format="appended" offset="N">`（N 为 appended 数据流中的字节偏移）
- 数据：`<AppendedData encoding="raw">` 后 `_` 起始，每数组前置 UInt64 字节长度头，随后为行主序（x 最快）内存镜像（little-endian）
- interior 元素遍历序与旧 ASCII 版一致（k→j→i），保证同一 u 场字节序不变
- `.pvti` 索引仍为 ASCII XML；`<PDataArray>` 声明与 piece 的 DataArray 类型/名称一致
- 错误处理沿用现状：打不开文件 WARN + 跳过（不抛异常）

**异常处理：** 文件打开失败 → `HYPOS_WARN` + return（既有行为保持）。

**验收标准：**
- Given 已知 u 场（可复现编码），When 写 `.vti`，Then XML 结构符合上述规范且 appended 字节可被独立解析器还原为相同数值序列（测试内实现最小解析）
- Given np=4 拼装场景，When `.pvti` + 4 个 `.vti` 生成，Then 索引引用与 piece extents 正确（更新既有 VtkPieceAndParallelIndexFiles 断言）
- Given 256²/512² 输出，Then binary vs 旧 ASCII 的文件体积与写耗时对比数据入 PERFORMANCE（≥3 次中位数）

### 3.2 R2 BinaryIOBackend 块写（D1-c）

**描述：** `.bin` 数据区从逐元素 `ofs.write` 改为连续段块写。

**触发条件：** `--output-format binary` 任意运行。

**期望行为：**
- 利用 row-major 布局：每个 i 行（iBegin..iEnd）为连续内存段，一次 `ofs.write(ptr, rowBytes)`；k/j 平面循环保序
- 文件字节布局不变：头部 7×Index 字段 + interior 数据镜像（旧读码/测试兼容）
- 错误处理沿用 WARN 模式

**验收标准：**
- Given 同一 u 场，When 块写版输出，Then 文件与逐元素版逐位一致（Red 期固化 fixture）
- Given 512² 规模，Then 块写前后耗时对比入 PERFORMANCE

### 3.3 R3 CI mpich 假通过探针（F1 紧急项）

**描述：** ctest 增加进程数探针，任何「多 rank 测试实为单进程」的情形变为确定性 FAIL。

**触发条件：** 任意 np>1 的 MPI 测试运行（本地 WSL OpenMPI 与 CI MPICH 双路径）。

**期望行为：**
- 新增独立探针测试（ctest 条目，np=4）：`MPI_Comm_size(MPI_COMM_WORLD) != 4` 即 FAIL 并打印实际值
- np>1 的 gtest MPI 条目（halo_2d_mpi、solver_mpi、topology_consistency 等 8 条）在测试体内激活 size 断言（filter 追加探针用例，消费注入的预期 np）
- hypos 直跑条目（collective_mpi 等 4 条，无 gtest 测试体）由环境级探针守卫：mpirun 退化为 launch 环境的全有/全无属性，同 suite 探针运行在相同 mpirun 下，任一暴露退化即全 suite 不可信（design D7）
- 探针在正常环境（WSL OpenMPI）全绿；在 mpich 4.2.0 缺陷环境应 FAIL（宁可红不可假绿）

**异常处理：** size 不符 → gtest FAIL（携带「实际 size/预期 np」诊断信息）。

**验收标准：**
- Given WSL OpenMPI np=4，When 跑全部 MPI ctest，Then 探针全绿（含既有条目）
- Given 探针测试，When `mpirun -np 4` 静默退化为单进程的缺陷环境语义（模拟：直接单进程跑探针二进制），Then 探针 FAIL（Red 期验证探针真的会咬缺陷）
- Given CI ci.yml，Then Debug(MPICH) job 在修复后的平台上跑全量测试且探针绿（真 np4 证据）

### 3.4 R4 CI 平台修复（F1 修复）

**描述：** CI 的 Debug(MPICH) job 脱离 ubuntu-latest（24.04）mpich 4.2.0 缺陷平台。

**触发条件：** CI 运行。

**期望行为：** 换发行版（如 ubuntu-22.04）或 pin 已修复 mpich 版本——具体方案 design 阶段决定并在 design.md 记录依据；修复后 MPICH 路径 np>1 测试为真多进程（探针背书）。

**异常处理：** 若评估后无可用修复平台，则回退为 Debug(OpenMPI)（记录决策依据），MPICH 路径暂出 CI 矩阵——诚实优先，不留假绿。

**验收标准：**
- Given 修复后的 CI 配置，When CI 运行，Then MPICH（或替代）路径全绿且探针在日志中可见
- Given README:28 与 GUIDE F1 条目，Then 文档与实际 CI 状态一致

### 3.5 R5 文档同步

**描述：** 文档与实现零缺口（延续 AR002-AR004 惯例）。

**期望行为：** README（mpich 缺陷记录更新为 CI 实际状态、P5 宣称对齐）、PERFORMANCE 新增 I/O 对比小节（D1 数据 + 口径说明）、OPTIMIZATION_GUIDE 勾选 P5 + D1 + F1 紧急项（附 AR 编号与证据）、AGENT_SPEC 如有接口变更则同步。

**验收标准：**
- Given 文档走查，Then 无失实宣称；GUIDE 勾选与本 AR 完成状态一致
- Given PERFORMANCE，Then 含 binary vs ASCII 体积/耗时两组数据及口径说明

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 正确性 | 输出逐位一致 | `.bin` 块写版与逐元素版文件逐位一致；`.vti` appended 数据字节 = interior u 场行主序内存镜像 |
| 兼容性 | 读取方 | `.vti` 符合 VTK XML appended raw 规范（ParaView 可加载为目标；本环境无 ParaView 时以规范合规断言 + 独立解析器测试替代，evidence 记录口径） |
| 兼容性 | 既有行为 | 打不开文件 WARN+跳过语义不变；CLI 接口不变（无新 flag）；`.pvti` 仍是 ASCII |
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 构建 | 双模式 | 0 新增警告 |
| 性能 | I/O 提速 | 256²/512² VTK 输出体积与耗时较 ASCII 显著下降（数据如实记录，不预设倍数门槛）；binary 块写耗时如实记录 |
| 回归 | 全量测试 | Debug+Release ctest 全绿；探针不引入误报 |
| 复杂度 | 预算 | 单 AR 净新增代码 ≤ ~800 行（GUIDE §6） |

## 5. 约束与假设

**约束：**
- 延续 AGENT_SPEC 全部约定（row-major、MPI 路径不抛异常、snake_case 文件名）
- 不引入外部依赖；C++17 + MPI + OpenMP
- VTK appended raw 模式依赖 little-endian 宿主（与 mpibin 后端既有假设一致，README 已有口径）
- 测试环境：WSL（OpenMPI，AR003 结论：mpich 在 Ubuntu 24.04 有单进程退化缺陷——本 AR 探针即为其修复）；CI：GitHub Actions（gitee 镜像仓库不跑 Actions，CI 修复以配置文件变更为交付物，运行证据以本地等价命令 + 探针代替，evidence 记录口径）

**假设：**
- VTK appended raw + header_type="UInt64" 为 ParaView 5.x 支持的标准格式（VTK File Formats 文档规范）
- i 行连续（row-major + 无 padding）成立——AlignedBuffer 无填充，design 阶段以代码核验
- CI Actions 在本仓库的执行可用性受限（gitee 主仓库），修复以「配置正确 + 本地等价验证」验收

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| appended raw | VTK XML 格式变体：数据负载追加在 XML 之后（`<AppendedData encoding="raw">`），二进制原样存储，XML 内以 offset 引用 |
| header_type | VTKFile 属性，声明 appended 长度头的整型宽度（本 AR 用 UInt64） |
| 块写 | 以连续内存段为单位调用 `ostream::write`，区别于逐元素调用 |
| 进程数探针 | 断言 `MPI_Comm_size == ctest 启动 np` 的测试，捕获 MPI 实现静默退化单进程的缺陷 |
| 假通过 | 测试进程组实为单进程但仍报告 PASS 的失绿情形 |
