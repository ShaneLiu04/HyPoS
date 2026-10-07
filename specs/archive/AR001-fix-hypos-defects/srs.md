# [AR001] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR001 |
| AR 主题 | fix-hypos-defects |
| 关联 SR | 无（源自代码审查缺陷清单） |
| 日期 | 2026-10-07 |
| 状态 | 已确认（用户授权自动化执行） |

> 注：`specs/component-detail-design/` 组件详设缺失，本 AR 以仓库实际代码、`AGENT_SPEC.md`（共享接口契约）与 `docs/DESIGN.md` 作为领域知识基线。

## 1. 背景与目标

HyPoS 是 MPI + OpenMP 混合并行的分布式泊松方程求解器（简历/作品集项目）。深度代码审查发现 6 类实质缺陷：构建必然失败、Halo 交换存在越界读写与面-邻居语义倒置、3D 通信与重叠/可观测性缺失、测试违反 MPI 规则、平台与文档存在占位/夸大。

本 AR 目标：在不改变项目定位（Linux HPC、C++17、MPI+OpenMP）的前提下，使 HyPoS **可构建、通信正确、内存安全、测试合规、文档与实现一致**。

## 2. 需求范围

**In Scope（本 AR 要做的）：**
- F1 构建修复：`CMakeLists.txt` 源清单与仓库一致（`src/io/vtk_io.cpp` 落位）
- F2 Halo 交换内存安全：缓冲区尺寸与 pack/unpack 元素数一致，仅打包内点面范围
- F3 面-邻居语义修正：6 方向发送面与目标邻居物理对应
- F4 功能补全：3D back/front 通信；`--overlap-comm` 真实生效；Profiler/Reporter 通信计时；`lastResidual()` 接口与 main 回填
- F5 测试体系合规：单次 MPI 初始化、ctest+MPIEXEC 多进程测试、非均匀场验证、CI 更新
- F6 平台适配+文档对齐：CMake 编译器分支与悬空引用清理、空目录处理、README/DESIGN/PERFORMANCE/SCALING_REPORT 与实现一致

**Out of Scope（本 AR 不做的）：**
- MSVC/Windows 支持（项目定位保持 Linux/POSIX）
- RMA、HDF5、PAPI 的实际实现（仅在文档标注为路线图）
- 新求解器算法（CG/SOR/Red-Black GS）
- CollectiveExchanger 改为真集合通信（保持 P2P 回退，但文档明确并清理死代码）

## 3. 功能需求

### 3.1 F1 构建可编译

**描述：** 修复 CMake 源清单引用不存在文件导致配置阶段报错的问题。

**触发条件：** 在干净目录执行 `cmake -B build && cmake --build build`。

**期望行为：** VTK 后端代码从 `src/io/binary_io.cpp` 拆分到 `src/io/vtk_io.cpp`（与 README 目录结构一致），CMake 源清单可解析；Debug 与 Release 两种模式均配置、编译成功。

**异常处理：** 编译失败即视为未完成；警告严格处理（`-Wall -Wextra -Wpedantic` 下相对修复前构建日志零新增警告，基线日志在开发阶段留存）。

**验收标准：**
- Given 干净源码树，When 执行 Release 配置与编译，Then 退出码 0 且产出 `hypos` 可执行文件
- Given 干净源码树，When 执行 Debug（ASan/UBSan）配置与编译，Then 退出码 0

### 3.2 F2 Halo 交换内存安全

**描述：** 当前缓冲区按 `local×halo` 分配而 pack/unpack 按 `total×total×halo` 遍历，写入量约为缓冲区 3 倍（未定义行为）。

**触发条件：** 任意多进程 halo 交换调用。

**期望行为：** 每个方向的发送/接收缓冲区元素数与打包循环计数严格一致；打包与解包仅覆盖**内点范围内的面**（左右面：j、k 取内点范围；上下面：i、k 取内点范围；前后面对：i、j 取内点范围）；不越界读写。

**异常处理：** 通过 ASan/UBSan 运行全部测试，任何越界报告即失败。

**验收标准：**
- Given Debug+ASan/UBSan 构建，When 运行 halo 交换测试（2D 4 ranks、3D 8 ranks），Then 无内存错误报告
- Given 任一方向面尺寸 S，When 初始化交换器，Then `sendBuf.size() == recvBuf.size() == pack 循环计数 == S`（单测断言）

### 3.3 F3 面-邻居语义正确

**描述：** 现实现将"本域最右内点列"发给左邻居（应为最左列），均匀平场下无法暴露，非均匀场会污染边界解。

**触发条件：** halo 交换后检查幽灵层数据。

**期望行为（6 方向约定）：**
- 发给左邻居 = 本域最左内点列（对方解入其右 halo）；发给右邻居 = 最右内点列
- 发给下邻居 = 最下内点行；发给上邻居 = 最上行
- 发给后邻居 = 最内层 k-面；发给前邻居 = 最外层 k-面
- MPI tag 配对与 unpack 目标区域保持自洽

**异常处理：** 边界方向使用 `MPI_PROC_NULL` 时跳过对应收发。

**验收标准：**
- Given 2D 4 ranks、每个单元值 = 全局坐标函数（如 `u = 1000*I + J`），When 完成 halo 交换，Then 每个幽灵单元值精确等于对应邻居内点单元值（逐单元断言，非均匀场）
- Given 3D 8 ranks（强制 2×2×2 拓扑），When 完成 halo 交换，Then 6 方向幽灵层均精确匹配

### 3.4 F4.1 3D back/front 通信

**描述：** 3D 运行时 back/front 方向从未通信，z 向分解时边界数据错误。

**触发条件：** `nz > 1` 且 z 方向进程数 > 1。

**期望行为：** P2P 交换器支持 6 方向共 12 个非阻塞请求（发送/接收各 6），缓冲对齐全；z 方向无分解（`dims[2]==1`）时边界仍为物理边界（Dirichlet）。

**验收标准：**
- Given 3D 8 ranks、`UniformPartition(2,2,2)` 拓扑，When halo 交换后，Then back/front 幽灵层与邻居 k 面数据精确一致（3.3 的非均匀场断言覆盖）
- Given 3D 4 ranks、z 不分解，When 求解运行完成，Then 无非物理 z 向通信（`MPI_PROC_NULL` 路径）且结果收敛

### 3.5 F4.2 通信-计算重叠（--overlap-comm）

**描述：** 现有实现仅解析并上报 `--overlap-comm`，solver 仍先同步交换再全量计算。

**触发条件：** 命令行传入 `--overlap-comm`。

**期望行为：** Jacobi 迭代使用 `beginExchange()` → 内点（与 halo 无关区域）计算 → `endExchange()` → 边界点计算的流水线；overlap 开关仅影响执行顺序，不影响数值结果；未开启时不改变现有行为。

**验收标准：**
- Given 同一问题与进程配置，When 分别以 overlap on/off 运行，Then 两者迭代次数相同（残差计算与归约顺序不因重叠分段改变）且解的最大差 ≤ 1e-12
- Given overlap 开启，When 检查 Profiler 报告，Then `halo_exchange`、`stencil_interior`、`stencil_boundary` 区域均被记录（区域命名与 §3.6 一致）

### 3.6 F4.3 性能可观测性补全

**描述：** Profiler 仅记录 `solver_total`/`jacobi_iteration`；Reporter 的 `comm_time_ms`、`overlap_ratio` 永远为 0。

**触发条件：** `--enable-profiling` 或默认运行。

**期望行为：**
- Profiler 增加命名区域：`halo_exchange`、`stencil_interior`、`stencil_boundary`、`residual_allreduce`
- Reporter 输出真实 `comm_time_ms`（交换总耗时）与 `overlap_ratio`（重叠开启时：`1 - wait_time/exchange_total`，无重叠时为 0；语义写入文档）
- 报告时间数据以 rank0（报告写出进程）的 Profiler 统计计算后写入（与 Profiler 打印同源，便于对账）

**验收标准：**
- Given 4 ranks 运行，When 读取 JSON 报告，Then `comm_time_ms > 0` 且与 Profiler `halo_exchange` 总耗时一致（±5%）
- Given overlap on/off 两次运行，Then `overlap_ratio` 分别在 [0,1] 内且关闭时为 0

### 3.7 F4.4 收敛残差回传

**描述：** `main.cpp` 的 `finalResidual` 从未赋值，报告恒为 0。

**触发条件：** solver 完成迭代（收敛或达到 maxIter）。

**期望行为：** `PoissonSolver` 新增 `virtual Real lastResidual() const`（无参、默认返回 0 的基类实现）；`JacobiSolver` 每次迭代后保存全局残差；main 用其填充 `PerformanceMetrics.finalResidual` 并写入报告。

**异常处理：** 未执行过迭代时返回 0。

**验收标准：**
- Given 收敛运行的 JSON 报告，Then `final_residual <= tolerance` 或等于最后一次全局残差且 > 0
- Given `lastResidual()` 未实现子类（接口默认），Then 返回 0 不崩溃

### 3.8 F5 测试体系合规

**描述：** 单个测试二进制内多个 TEST 反复 `MPI_Init/MPI_Finalize`（MPI 规则违规）；ctest 不经 mpirun 运行多进程测试；4 进程通信测试默认被 skip；偶发崩溃或假绿。

**触发条件：** 构建测试并执行 ctest。

**期望行为：**
- 每个测试二进制使用统一入口：`main()` 中一次 `MPI_Init_thread`，测试体不再初始化/终止 MPI，退出前一次 `MPI_Finalize`
- 多进程测试经 `MPIEXEC_EXECUTABLE` + np 参数注册到 ctest（2D halo：np=4；3D halo：np=8；允许 `--oversubscribe`）
- halo 测试改用非均匀场数据（承接 F3 验收）
- solver 多进程一致性测试：串行解 vs 4 ranks 解，最大差 ≤ 1e-10
- CI 增加 GTest 依赖与 `ctest` 步骤；本地 WSL `ctest` 全绿

**验收标准：**
- Given WSL 构建目录，When `ctest --output-on-failure`，Then 所有测试 PASS（含 4/8 rank 测试，无 skip）

### 3.9 F6 平台适配与文档对齐

**描述：** CMake 无条件使用 GCC 专属标志与失败路径；`cmake/`、`benchmarks/reference_results/` 为空；README 宣称 RMA/HDF5/PAPI 等未实现功能；SCALING_REPORT 为空表格模板。

**触发条件：** 构建与文档审阅。

**期望行为：**
- CMake 按编译器（GNU/Clang/Intel）设置警告/优化/`-march`/sanitizer 标志；不支持的编译器给出显式提示
- 清理悬空引用：`vtk_io`（随 F1 解决）、`FindPAPI`（提供最小 `cmake/FindPAPI.cmake` 或移除引用二选一，落实情况写入文档）
- `cmake/` 与 `benchmarks/reference_results/` 不再为空：前者承载编译器选项/查找模块，后者存放修后基准性能 JSON + README 说明
- README/DESIGN/PERFORMANCE：与实际实现一致。重点清理对象：README「Halo Exchange（集合通信）」「通信-计算重叠」表述（L14-15、L237）；DESIGN §6「扩展性预留」表（RMA/HDF5/PAPI 标注"路线图/未实现"）；PERFORMANCE §5.2「通信-计算重叠效果」示意表（标注为参考示意或替换为实测）；CollectiveExchanger 统一标注 P2P 回退
- SCALING_REPORT：填入 WSL 实测 1/2/4 进程强/弱扩展数据（环境表如实填写），未覆盖的组合标注测量范围

**验收标准：**
- Given 全部修改完成，When 通读 README 与 docs，Then 不存在与实际实现不符的功能宣称
- Given WSL 环境，When 运行 scaling 脚本，Then SCALING_REPORT 表格包含实测数值与运行命令

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 数值一致性 | overlap on/off、串行/并行 | ≤ 1e-12（同迭代数）；串并行 ≤ 1e-10 |
| 构建 | 双模式 | Debug(ASan) + Release 均成功，零新增警告 |
| 性能 | 参考基准 | 修复后基准存入 `benchmarks/reference_results/`，供回归对比 |
| 兼容性 | 接口 | 除新增 `lastResidual()` 外不破坏现有公开接口 |

## 5. 约束与假设

**约束：**
- C++17；MPI-3.1 + OpenMP 4.5；Linux/POSIX（GCC/Clang/Intel），不做 MSVC
- 保持现有模块分层与接口风格（策略模式、RAII、AGENT_SPEC.md 约定）
- 验收工具链：WSL Ubuntu（cmake、openmpi、libgtest-dev），由用户授权安装

**假设：**
- WSL 可安装所需软件包（用户提供 sudo 凭据）
- 测试机可运行 4/8 rank 测试（必要时 mpirun `--oversubscribe`）
- 性能数据仅作本机参考（不承诺集群级扩展性）

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| Halo（幽灵层） | 子域边界外侧的缓冲单元，存放邻居内点数据 |
| 面（Face） | 某个方向 halo 层对应的邻居内点切片 |
| Overlap | 通信-计算重叠：内点计算与 halo 通信并发执行 |
| P2P 回退 | CollectiveExchanger 内部委托 PointToPointExchanger 实现交换 |
| FR | 本 AR 功能需求编号前缀（F1–F6） |
