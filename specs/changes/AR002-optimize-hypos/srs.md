# [AR002] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR002 |
| AR 主题 | optimize-hypos（全项目深度优化） |
| 关联 SR | 无（源自 AR001 归档后的遗留问题与优化机会清单） |
| 日期 | 2026-10-07 |
| 状态 | 已确认（用户授权自动化执行） |

> 注：`specs/component-detail-design/` 组件详设仍缺失；以 AR001 归档文档（`specs/archive/AR001-fix-hypos-defects/`）、`AGENT_SPEC.md` 与当前代码为领域基线。

## 1. 背景与目标

AR001 修复了 6 类实质缺陷并归档。其后审查发现一批遗留问题与优化机会：

- **正确性/一致性**：VTK/Binary 输出仅为 rank0 本子域（非全局解），binary 头 offsets 恒 0；`--save-interval` 有文档无实现；`MPIEnv` 未按 AGENT_SPEC 接入 main；`MemoryPool` 从未使用（死代码）；`applyNeumannBC` 未接线；报告 `comm_overhead_ratio` 与 `scaling.*` 恒 0。
- **性能**：残差使用独立全数组扫描（多一遍内存读）；打包逐元素复制；无持久化通信；无 SIMD 对齐提示；FLOPs 估算偏低（2 flops/格点）。
- **算法**：仅有 Jacobi（收敛慢，64²→1e-7 实测需 10996 次）；无 RBGS/CG。
- **工程化**：无 `.gitignore`；ctest 未覆盖 collective/奇数进程/新求解器；scaling 脚本硬编码；基准未含优化前后对比。

本 AR 目标：把 HyPoS 从"修复后可运行"提升为**输出正确、性能有据、算法可选、工程收口**的完整 HPC 作品。

## 2. 需求范围

**In Scope：**
- R1 全局解输出（PVTU 分片 + 索引；Binary 真实 offsets）
- R2 `--save-interval` 中间解保存
- R3 RAII/清理收口（MPIEnv 接入、collective 告警限流、MemoryPool 处置、Logger 注释）
- R4 Neumann 边界接线（`--bc`）
- R5 报告字段与 FLOPs 估算收口
- R6 通信与内核性能优化（残差融合、持久化通信、memcpy 打包、SIMD 提示、退化子域归一化）
- R7 求解器扩展（Red-Black GS 与 CG）
- R8 工程化收口（.gitignore、ctest 扩充、脚本优化、文档同步、基准刷新）

**Out of Scope：**
- 周期边界（periodic）、残差按格点归一化的语义变更、HDF5/RMA、CG 预条件子、Docker/DevContainer、MPI-IO 单文件聚合（PVTU 分片已满足并行输出）、MSVC/Windows 支持、`--tol` 默认值调整

## 3. 功能需求

### 3.1 R1 全局解输出

**描述：** 并行运行时输出可还原的全局解。

**触发条件：** `--output-format vtk|binary` 任意 rank 数运行。

**期望行为：**
- VTK：每 rank 写自身子域 `solution_<step>_r<rank>.vti`（VTK XML ImageData：`WholeExtent` 为全局点范围、`Origin=offset*dx`、正确 `Spacing`），rank0 额外写 `solution_<step>.pvti` 索引（引用全部分片）；最终解以实际迭代数为 step；单 rank 时行为与全局等价。（格式族说明：结构化点数据用 XML `ImageData`/`PImageData`，非 `.vtu` 非结构化格式）
- Binary：每 rank 写 `solution_<step>_r<rank>.bin`，头部 `offsetX/Y/Z` 填真实全局偏移（0 基内点坐标）。

**异常处理：** 无法写文件时告警并跳过（不崩溃）。

**验收标准：**
- Given np=4 运行完成，When 检查输出目录，Then 存在 4 个 `.vti` 分片 + rank0 的 `.pvti`；分片 `Origin`/`Extent` 与 `offsetX+Y+Z` 一致且 4 个 Extent 恰好拼成全局面
- Given np=1 运行，Then 单分片覆盖全局域（Origin=0，Extent=全局）

> 已知边界说明：VTK 与非结构化 `.vtu`/`.pvtu` 族相比，结构化点场使用 `.vti`/`.pvti`（设计预研结论，格式对被 ParaView 正确识别所必需）。

### 3.2 R2 `--save-interval` 中间解保存

**描述：** 每 N 次迭代保存一次中间结果。

**触发条件：** `--save-interval N`（N>0）且配置了 IO 后端（vtk/binary）。

**期望行为：** solver 迭代循环按固定步数回调主流程保存接口；输出沿用 §3.1 命名 `solution_<step>_r<rank>`（step=全局迭代步号）；N=0（默认）不产生任何中间文件且零额外开销。

**异常处理：** 单次保存失败仅告警，不中断求解。

**验收标准：**
- Given `--save-interval 50 --output-format vtk --max-iter 120`，When 运行完成，Then 存在 `solution_50_r*` 与 `solution_100_r*` 分片（step=实际迭代数的最终解另行写出）；未到 50 步不产生中间文件
- Given 不传 `--save-interval`，Then 输出目录除最终解外无任何 `solution_<k>_r*` 中间分片

### 3.3 R3 RAII 与清理收口

**描述：** 消除死代码与规范偏离。

**期望行为：**
- `main.cpp` 改用 `MPIEnv`（RAII）管理 MPI 生命周期（AGENT_SPEC 对齐；线程级别 FUNNELED 不变）
- `CollectiveExchanger` 回退告警仅 rank0 打印一次
- `MemoryPool`：**接入真实用途或整体移除**（design 阶段二选一，不留死代码；若移除需同步 CMake/README/DESIGN）
- `Logger` 头注释更正（不实称 thread-safe；或补互斥）——design 决策
- 不新增通信路径异常

**验收标准：**
- Given 代码检视 + 运行，Then main 无裸 `MPI_Init/MPI_Finalize` 调用，行为与 AR001 一致（e2e 通过）
- Given np=4 `--comm-mode collective`，Then 回退告警恰 1 条
- Given 全仓库 grep，Then MemoryPool 无"定义但零使用"状态；文档与实现一致

### 3.4 R4 Neumann 边界接线

**描述：** 把既有 `applyNeumannBC` 接入 CLI。

**触发条件：** `--bc neumann`（默认 `dirichlet`）。

**期望行为：** neumann 模式下物理边界 halo 按零通量镜像填充（内部边界 halo 仍由交换覆盖）；dirichlet 模式行为与 AR001 完全一致。

**异常处理：** `--bc` 非法值 → 退出码 1 + 错误信息。

**验收标准：**
- Given `--bc neumann` np=1 与 np=4 运行 32²固定迭代，Then 正常完成不崩溃、边界 halo 为镜像值（单元测试断言镜像语义）；纯 Neumann 且净通量非零的泊松问题本身无解（数学性质），不做收敛性要求，须在文档说明
- Given 未传 `--bc`，Then 行为与 AR001 基线逐位一致（同配置 ctest 通过）

### 3.5 R5 报告字段与估算收口

**描述：** 报告不再输出无意义空值，FLOPs 估算真实。

**期望行为：**
- `comm_overhead_ratio = comm_time_ms / iter_time_ms`（含零保护；on/off 均有真实值）
- `compute_time_ms` 回填 = Profiler `stencil_interior + stencil_boundary` 总秒 ÷ 迭代数 × 1000（口径与 comm_time_ms 一致，rank0）
- `memory_bandwidth_gbps` 从 JSON 移除（当前无测量手段，不留假数据）
- `scaling` 块：仅在提供外部输入时输出（否则从 JSON 移除，避免恒 0 假数据）
- FLOPs 估算按实际离散格式（2D 五点/3D 七点更新+残差≈6/10 flops 每格点）修正，并在 README 性能示例同步

**验收标准：**
- Given np=4 运行报告，Then `comm_overhead_ratio>0` 且与 comm/iter 比值一致（±5%）；`compute_time_ms>0` 且与 Profiler `stencil_interior+stencil_boundary` 总耗时一致（±5%）
- Given 默认运行，Then JSON 无恒 0 的 scaling 数值字段、无 memory_bandwidth_gbps 字段
- Given 256² 单 rank，Then `flops_per_sec` 与解析估算（flops_per_iter×iterations/total_time）偏差 ≤5%

### 3.6 R6 通信与内核性能优化

**描述：** 以有量化证据的性能优化提升迭代吞吐。

**期望行为：**
- **残差融合**：取消独立残差扫描；在更新 pass 中按**固定区域顺序**（内点盒→6 slab）累积部分和，overlap on/off 结果仍 bit 级一致
- **退化子域归一化**：`nxLocal ≤ 2·hw` 时 slab 端点归一化，消除同单元双写
- **持久化通信**：exchanger initialize 时 `MPI_Send_init/Recv_init`，迭代内 `MPI_Startall/Waitall`（PROC_NULL 方向不注册）
- **打包 memcpy 化**：缓冲布局调整为最内维连续，行级 `memcpy` 替代逐元素复制
- **SIMD 向量化**：`__restrict` 指针 + `#pragma omp simd`（不用 `aligned` 断言：内点带偏移使 64B 对齐不可保证，错误断言属 UB）；`-fopt-info-vec` 验证核心循环向量化

**异常处理：** 持久化请求在 exchange 未初始化即调用时返回错误日志（防御性）。

**验收标准：**
- Given 256² 固定 10000 迭代，When 优化后运行 np=1 与 np=4（OMP=1，3 次取中位数），Then `iter_time` 相比 AR001 基准（np1 0.0451ms、np4 0.0204ms）改善 **≥10%**
- Given 全局 4×4 np=4（1×1 退化子域），Then 解有限、与串行参考一致（≤1e-12），且 Debug 下区域划分自检断言通过（各 slab 覆盖单元数之和 = 内点单元数，无交叠双写）
- Given overlap on/off 一致性测试，Then 仍 ≤1e-12 且迭代数相同
- Given Debug+ASan 全量 ctest，Then 零报告；`-fopt-info-vec` 日志显示 stencil 内层循环已向量化

### 3.7 R7 求解器扩展（RBGS + CG）

**描述：** 新增 Red-Black Gauss-Seidel 与共轭梯度求解器，沿用 `PoissonSolver` 接口与 `HaloExchanger`。

**触发条件：** `--solver red_black_gs | cg`。

**期望行为：**
- RBGS：红/黑两次半扫，每半扫一次 halo 交换；带残差（`lastResidual()` 语义同 Jacobi：全局 L2 步长范数）
- CG：matvec 复用 stencil + 每迭代一次 halo 交换；2 次全局内积归约；`lastResidual()` 上报残差范数
- 二者 `solve()` 返回实际迭代次数；与 IO/overlap 选项兼容（overlap 仅对 Jacobi 生效，其他 solver 忽略并告警一次）

**异常处理：** 未知 solver 别名仍退出码 1。

**验收标准：**
- Given 制造解 64²（沿用 M3a 问题与解析误差基准），When RBGS/CG 求解 tol=1e-7，Then 收敛、全局 L2 误差 <1e-3，且迭代数显著低于 Jacobi（RBGS ≤50%、CG ≤10%；Jacobi 基线 ~11000）
- Given 与 Jacobi 收敛解（同问题、tol=1e-12，Jacobi 充分迭代）比较，Then 逐点差 ≤1e-8
- Given 3D 冒烟（64×64×8 np=4），Then 正常完成、解有限

### 3.8 R8 工程化收口

**描述：** 补齐工程化缺口并刷新基准。

**期望行为：**
- `.gitignore`：忽略 `build*/`、`scaling_results/`、`scaling_plots/`、`output/`、`*.o` 等生成物
- ctest 新增：`collective_mpi`（np=4 collective 回退运行）、`ranks3_smoke`（np=3 任意分解）、`solver_alt_mpi`（RBGS/CG np=4 收敛用例）
- `run_scaling_tests.py`：新增 `--max-iter/--tol` 参数；`strong_scaling.json/weak_scaling.json` 写入效率数据供报告引用
- `plot_scaling.py`：移除未用 import
- 文档同步：README（输出格式、求解器、BC、save-interval 行为）、DESIGN（新组件行为）、PERFORMANCE（优化手段与实测对比）
- 基准刷新：优化后重跑 256² np=1/2/4 与弱扩展，更新 `benchmarks/reference_results/` 与 `docs/SCALING_REPORT.md`（含优化前后对比表）

**验收标准：**
- Given Release/Debug ctest，Then 全绿（含 ≥3 个新用例）
- Given 文档走查，Then 无失实宣称；基准目录含优化前后数据
- Given `git status`（若初始化）或文件清单核查，Then 生成物被忽略规则覆盖

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 性能 | iter_time 改善 | 256² np=1/4 相比 AR001 基准 ≥10%（中位数） |
| 数值一致性 | overlap on/off；求解器间 | ≤1e-12；≤1e-8（同问题收敛解） |
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 构建 | 双模式 | 0 新增警告 |
| 兼容性 | CLI/接口 | 默认行为与 AR001 一致；`PoissonSolver::solve` 如需扩展须保持现有调用兼容 |

## 5. 约束与假设

**约束：**
- C++17 + MPI + OpenMP；POSIX/Linux（MSVC 仍显式拒绝）
- 保持 AGENT_SPEC.md 约定与既有分层/风格；不引入外部依赖
- 测试环境：WSL；测试 OMP=4、性能基准 OMP=1（AR001 环境结论延续）
- 性能对比必须基于同一环境、同一配置、≥3 次中位数

**假设：**
- AR001 基准数据（`benchmarks/reference_results/`）为对比基线有效
- CG/RBGS 收敛行为在制造解与 -1 右端项两类问题上稳定
- PVTU/VTU 文本格式足以被 ParaView 解析（无外部验证工具，按规范格式自检 + 结构断言）

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| PVTU/VTU | ParaView 并行数据集：`.pvtu` 为索引，`.vtu` 为单分片 |
| 残差融合 | 把残差累积合并进 stencil 更新循环，消除独立扫描 |
| 持久化通信 | MPI_Send_init/Recv_init 预注册请求，迭代内 Startall/Waitall |
| RBGS | Red-Black Gauss-Seidel，双色半扫并行迭代法 |
| CG | 共轭梯度法（无预条件） |
