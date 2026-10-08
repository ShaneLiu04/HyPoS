# [AR003] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR003 |
| AR 主题 | mpiio-single-file-output（MPI-IO 单文件输出） |
| 关联 SR | 无（源自 AR002 design §4.1.1 方案 C 路线图） |
| 日期 | 2026-10-08 |
| 状态 | 已确认（用户授权自动化执行） |

> 注：`specs/component-detail-design/` 组件详设仍缺失；以 AR001/AR002 归档文档（`specs/archive/`）、`AGENT_SPEC.md` 与当前代码为领域基线。

## 1. 背景与目标

AR002 实现了全局解输出，采用方案 A（每 rank 分片 + rank0 索引：VTI/PVTI、per-rank `.bin`）。其 design §4.1.1 方案 C（MPI-IO 单文件自定义布局）当时以"自研格式生态差、工作量大"为由列入 AR003 路线图。

本 AR 目标：新增基于 MPI-IO 的**单文件**二进制输出能力，展示并行聚合 I/O 的 HPC 工程深度，同时保持现有分片行为完全不变（作为并行可扩展的默认路径）。

## 2. 需求范围

**In Scope：**
- R1 MPI-IO 单文件二进制输出后端（新格式值 `--output-format mpibin`、新类 `MPIIOBinaryBackend : IOBackend`、完整自描述头 + 全局行主序数据区）
- R2 CLI 集成与 `--save-interval` 中间保存（每 step 一个单文件）
- R3 异常处理（MPI-IO 失败告警跳过，不崩溃）
- R4 验证与测试（测试内读回：np=1 等价 + np=4 全局拼接一致；ctest 回归）
- R5 基准参考数据（np=1/np=4 写入耗时 mpibin vs binary 分片，不设硬指标）
- R6 文档全套同步（README/DESIGN/PERFORMANCE）

**Out of Scope：**
- VTK 单文件 MPI-IO 合并（`.vti` XML 头部需预留空间方案，复杂度独立评估）
- HDF5、移除/淘汰现有分片输出行为
- 独立读回工具脚本（Python 等）
- MSVC/Windows 支持

## 3. 功能需求

### 3.1 R1 MPI-IO 单文件二进制输出后端

**描述：** 全 rank 通过 MPI-IO collective 调用将各自子域写入单一二进制文件。

**触发条件：** `--output-format mpibin` 任意 rank 数运行（最终输出与 `--save-interval` 中间保存）。

**期望行为：**
- 新增 `MPIIOBinaryBackend : IOBackend`（src/io），`write()` 为 **collective 语义**：全 rank 以相同文件名调用（命名 `solution_<step>.bin`，无 `_r<rank>` 后缀）；`writeParallelIndex()` 为 no-op（单文件无需索引）
- 文件布局：**完整自描述头** + **全局行主序数据区**
  - 头部字段：magic（4 字节）、version、全局 nx/ny/nz、dx/dy/dz、边界类型（dirichlet/neumann）、数据区起始 offset
  - 数据区：按全局网格行主序（x 维最快，与 AGENT_SPEC 内存布局一致）连续存储全局内点场；各 rank 使用 MPI subarray filetype 将本子域直写到自身全局位置，无块拼接缝
- 读回视角：文件可被任意支持 MPI-IO/POSIX 读的程序解析——头部自描述，数据区即完整全局场
- 单 rank（np=1）时行为与全局等价（同一单文件，覆盖全局域）

**异常处理：** `MPI_File_open`/写入失败时 `HYPOS_WARN` 告警并跳过本次输出（对齐 AR002 §3.1 约定：不崩溃、不阻断求解）。

**验收标准：**
- Given np=4 运行完成，When 检查输出目录，Then 存在**恰好 1 个** `solution_<step>.bin`；C++ 测试读回后头部字段（magic/version/全局尺寸/间距/边界类型）与运行参数一致，数据区与全局参考场逐点一致（np=1 逐位一致；np=4 ≤1e-15）
- Given np=1 运行，Then 单文件数据区与同问题串行参考场逐位一致；与现有 binary 分片（np=1 时单分片）数据区语义等价
- Given 任一 rank 的 MPI-IO 调用失败（如目录不可写），Then 告警且求解正常完成

### 3.2 R2 CLI 集成与中间保存

**描述：** `mpibin` 作为新格式值接入 main 流程与 `--save-interval`。

**触发条件：** `--output-format mpibin [--save-interval N]`。

**期望行为：**
- main.cpp 按格式值实例化 `MPIIOBinaryBackend`；`writeSolution` lambda 对单文件后端使用无 rank 后缀的共享文件名（全 rank 一致）
- `--save-interval N` 中间快照同走 MPI-IO 单文件（每 step 一个 `solution_<step>.bin`）；N=0（默认）零额外开销
- 其他格式值（json/csv/vtk/binary）行为与 AR002 完全一致；非法格式值仍退出码 1

**异常处理：** 单次保存失败仅告警，不中断求解。

**验收标准：**
- Given `--save-interval 50 --output-format mpibin --max-iter 120`，When 运行完成，Then 存在 `solution_50.bin`、`solution_100.bin` 与最终解 `solution_<实际迭代数>.bin`，均为单文件
- Given 未传 `--save-interval`，Then 输出目录仅有最终解单文件
- Given 同配置 `--output-format binary`/`vtk` 回归运行，Then 行为与 AR002 基线一致（ctest 通过）

### 3.3 R3 异常处理与防御性

**描述：** MPI-IO 路径的错误处理对齐项目约定。

**期望行为：**
- 打开/写失败：告警跳过本次输出，不抛异常、不改变迭代流程（"Never throw in MPI communication paths"）
- collective 一致性：告警后所有 rank 走同一后续路径（不因单 rank 失败造成 collective 挂死）
- 文件句柄 RAII 管理或保证所有路径 `MPI_File_close`

**验收标准：**
- Given 输出目录不可写（模拟失败），Then np=4 运行完成不挂死、不崩溃，告警有界（不每迭代刷屏，限流或仅首次）
- Given 代码检视，Then 无 `MPI_File` 泄漏路径（异常/告警分支均 close）

### 3.4 R4 验证与测试

**描述：** 单文件正确性由 C++ 测试内读回验证。

**期望行为：**
- 新增 MPI 测试：np=1 与 np=4 下运行 `MPIIOBinaryBackend::write` → 以 POSIX/MPI-IO 独立读回 → 断言头部字段与数据区
- np=4 拼接正确性：读回的全局场与制造解参考场逐点一致（子域边界处无错位/缝隙/覆盖）
- np=1 等价性：与串行构造的全局参考场逐位一致
- ctest 注册新用例；全量回归（Debug + Release，含 ASan）零报告

**验收标准：**
- Given 新测试在 np=1 与 np=4 下执行，Then 全部断言通过
- Given 全量 ctest（Debug/Release + ASan），Then 全绿、零 sanitizer 报告
- Given 3D 冒烟（如 64×64×8 np=4），Then 单文件读回数据区尺寸与全局一致（3D 布局正确）

### 3.5 R5 基准参考数据

**描述：** 提供 mpibin vs binary 分片的写入耗时参考数据（不设硬性指标）。

**期望行为：**
- 同环境同配置（256²，np=1/np=4，OMP=1，≥3 次取中位数）对比两种格式的输出耗时
- 数据记入 PERFORMANCE 文档与 AR evidence（含环境说明）

**验收标准：**
- Given 基准运行，Then PERFORMANCE 含 mpibin vs binary 对比表，注明环境与方法（无失实宣称）

### 3.6 R6 文档全套同步

**描述：** 文档与实现一致（AR002 D1 文档规范延续）。

**期望行为：**
- README：输出格式表新增 `mpibin`（单文件、MPI-IO）；「路线图」节移除/更新 MPI-IO 单文件条目（已实现项不得留在"未实现"清单）
- DESIGN：IO 层新增 `MPIIOBinaryBackend` 描述（布局、collective 语义、与分片方案的取舍）
- PERFORMANCE：R5 基准数据与方法

**验收标准：**
- Given 文档走查，Then 无失实宣称；路线图节不含已实现能力

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 正确性 | 数据区 vs 参考场 | np=1 逐位一致；np=4 逐点一致（≤1e-15，浮点求和顺序不变应为逐位） |
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 构建 | 双模式 | 0 新增警告 |
| 兼容性 | CLI/接口 | 现有格式值与默认行为不变；`IOBackend` 接口不破坏（新增派生类） |
| 集体一致性 | collective | 失败路径全 rank 一致，无挂死风险 |

## 5. 约束与假设

**约束：**
- C++17 + MPI + OpenMP；POSIX/Linux（MSVC 显式拒绝）
- 保持 AGENT_SPEC.md 约定与既有分层/风格；不引入外部依赖（仅 MPI 标准 IO 接口）
- 头部字段使用固定宽度类型（`Index`=size_t 8B、`Real`=double 8B）；小端假设须在头注释/文档注明
- 测试环境：WSL；测试 OMP=4、性能基准 OMP=1（延续既有环境结论）

**假设：**
- 所用 MPI 实现（MPICH/OpenMPI）支持 MPI-IO 且对本地文件系统可用（WSL ext4 下 MPI_File_open 可用，AR002 已验证 MPI 栈）
- subarray filetype + `MPI_File_write_at_all` 在 np≤4 规模下正确性可依赖标准语义
- 256² 规模的输出耗时受调度噪声影响，参考数据以中位数给出（与 AR002 scaling 快照口径一致）

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| MPI-IO | MPI 标准第 9/13 章并行文件 I/O 接口（MPI_File_open、MPI_File_write_at_all 等） |
| subarray filetype | MPI_Type_create_subarray 描述的全局数组局部视图，配合 collective 写实现各 rank 直写全局位置 |
| 单文件布局 | 自描述头 + 全局行主序数据区；区别于 per-rank 分片（每 rank 一文件）+ 索引 |
| collective write | 全 rank 共同参与的同步写调用（write_at_all），文件名与调用次序须全 rank 一致 |
