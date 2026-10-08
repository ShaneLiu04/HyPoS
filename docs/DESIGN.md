# HyPoS 架构设计文档

> 本文档描述 HyPoS 的系统架构、模块划分和关键设计决策。

---

## 1. 系统概述

HyPoS 采用**分层架构**，将数学求解、并行通信、网格管理和 I/O 解耦为独立模块。这种设计允许：
- 求解器算法（Jacobi → Gauss-Seidel → CG）的独立替换
- 通信模式（P2P → Collective → RMA）的独立替换
- 网格划分策略（Uniform → Hilbert → METIS）的独立替换

---

## 2. 模块架构

### 2.1 核心层（core/）

- **types.hpp**: 类型别名（`Real = double`, `Index = std::size_t`），统一类型体系
- **aligned_buffer.hpp**: 64 字节对齐的 RAII 内存分配器，防止 false sharing
- **exception.hpp**: 异常层次结构（HyPoSException → MPIException → ConfigException）

### 2.2 网格层（grid/）

```
Grid (全局描述) --> GridPartition (策略) --> Subgrid (本地子域 + Halo)
```

- **Grid**: 全局网格元数据（nx, ny, nz, dx, dy, dz, haloWidth）
- **GridPartition**: 划分策略接口
  - `UniformPartition`: 使用 `MPI_Cart_create` + `MPI_Dims_create` 的自动拓扑划分
  - 支持余数均匀分配（remainder 进程优先获得额外行/列）
- **Subgrid**: 本地子域管理
  - 内存：三个 `AlignedBuffer<Real>`（u, u_next, rhs）
  - 索引：全局索引到本地索引的偏移映射（含全局 offsets，用于并行输出与 RBGS 全局着色）
  - 邻居：通过 `MPI_Cart_shift` 预计算的邻居 rank
  - 物理边界：Dirichlet/Neumann（`applyPhysicalBoundary`，逐迭代仅刷新无邻居的 halo 面）

### 2.3 求解器层（solver/）

- **PoissonSolver**: 抽象基类，定义 `solve()` 和 `iterate()` 接口
- **JacobiSolver**: 经典 Jacobi 实现
  - 2D 五点 stencil：`u_new = 0.25 * (u_left + u_right + u_down + u_up - dx^2 * f)`
  - 3D 七点 stencil：除以 6.0
  - OpenMP 并行：`#pragma omp parallel for` + `#pragma omp simd`
  - 残差计算：内点 L2 范数，通过 `MPI_Allreduce` 全局聚合
  - 通信-计算重叠（可选，`JacobiSolver(true)`）：`beginExchange` → 内点盒计算 → `endExchange` → 边界 6 slab；残差使用**融合累积**（固定区域顺序），保证 overlap on/off 的迭代数与解 bit 级一致（`lastResidual()` 暴露最近残差）
- **RedBlackGSSolver**: 双色半扫（全局奇偶着色，依赖子域 offsets），每半扫一次 halo 交换；迭代数约为 Jacobi 的一半（模型问题理论比 ≈1/(1+ρ)）
- **CGSolver**: 无预条件共轭梯度（matvec 复用 stencil 并交换辅助向量；每迭代两次全局内积归约）；一般右端项下迭代数比 Jacobi 低 1-2 个数量级

### 2.4 通信层（comm/）

- **HaloExchanger**: 抽象基类，定义 `beginExchange()` / `endExchange()` 接口
- **PointToPointExchanger**: 持久化非阻塞 MPI 实现
  - `MPI_Send_init` / `MPI_Recv_init` + `MPI_Startall` / `MPI_Waitall`（请求在 `initialize` 注册、析构释放）
  - 6 个方向（左/右/下/上/后/前）共 12 个请求；tag 按发送面编号 0-5；支持任意 padded 缓冲（`exchange(sg, data)`）
  - 每个方向的独立 send/recv 缓冲区；仅交换内点面范围（2D 数据位于 k=0 平面）；包/解包为行级 `memcpy`（缓冲布局最内维连续）
- **CollectiveExchanger**: P2P 回退实现（委托 `PointToPointExchanger` 并记录告警）；真 `MPI_Neighbor_allgatherv` 为路线图

### 2.5 性能层（perf/）

- **Timer**: `std::chrono::steady_clock` 高精度计时
- **Profiler**: 分层区域计时器，支持嵌套区域
  - `ProfileScope`: RAII 自动 begin/end
- **Reporter**: 性能报告生成器
  - JSON 格式：结构化、可扩展
  - CSV 格式：便于批量数据分析

### 2.6 I/O 层（io/）

- **IOBackend**: 抽象接口（含可选的并行索引接口 `writeParallelIndex`）
- **BinaryIOBackend**: 原始二进制输出，最高性能（每 rank 分片，头部含真实全局 offsets）
- **VTKIOBackend**: VTK XML ImageData 分片（每 rank `.vti`，带全局 Origin/Extent 偏移）+ rank0 `PImageData` 索引（`.pvti`），直接支持 ParaView
- **MPIIOBinaryBackend**: MPI-IO 单文件二进制输出（`--output-format mpibin`）。`write()` 为 collective：全 rank 在 `subgrid.comm()` 上以相同文件名调用，经 `MPI_File_set_view`（文件侧 subarray=全局网格内本 rank 内点盒）+ `MPI_File_write_all`（内存侧 subarray=padded 缓冲内点盒）零拷贝直写全局位置。文件布局：72 字节自描述头（magic/version/全局尺寸/间距/边界类型/数据区 offset，小端、固定宽度字段）+ 全局行主序数据区（x 最快）；`MPI_File_set_size` 先行截断旧文件。I/O 失败不抛出：WARN（每 rank 限一次）并跳过本次输出，依赖文件默认 error handler `MPI_ERRORS_RETURN` 保证 open 失败全 rank 一致返回，无 collective 挂死。与分片方案（Binary/VTK）的取舍：单文件自描述、产物简洁，适合归档与第三方读取；分片方案无聚合依赖、天然并行可扩展——两者并存，`writeParallelIndex` 在单文件后端为 no-op。

### 2.7 工具层（utils/）

- **MPIEnv**: MPI 初始化/终止的 RAII 封装，线程安全级别管理
- **Logger**: 分级日志（TRACE → FATAL），带时间戳和 rank 前缀
- **CommandLineParser**: 现代 C++ 参数解析，支持 `--key value` 和 `--key=value`

---

## 3. 内存布局

### 3.1 子域数组布局（SoA）

```
Subgrid 内存结构：
  u      [AlignedBuffer<Real>] — 当前解
  u_next [AlignedBuffer<Real>] — 下一迭代解
  rhs    [AlignedBuffer<Real>] — 右端项 f

每个数组大小 = (nx_local + 2*halo) * (ny_local + 2*halo) * (nz_local + 2*halo)
对齐：64 字节（缓存行边界）
```

### 3.2 索引映射

```cpp
// 全局网格坐标 (I, J, K) -> 本地子域坐标
i_local = I - offsetX + haloWidth
j_local = J - offsetY + haloWidth
k_local = K - offsetZ + haloWidth

// 平坦索引
index(i, j, k) = (k * ny_total + j) * nx_total + i
// 其中 nx_total = nx_local + 2*haloWidth
```

### 3.3 Halo 区域定义

```
Left   halo: i = [0, halo-1]
Right  halo: i = [nx_local+halo, nx_total-1]
Bottom halo: j = [0, halo-1]
Top    halo: j = [ny_local+halo, ny_total-1]
Back   halo: k = [0, halo-1]     (3D only)
Front  halo: k = [nz_local+halo, nz_total-1] (3D only)

Interior: i = [halo, nx_local+halo-1]
          j = [halo, ny_local+halo-1]
          k = [halo, nz_local+halo-1]
```

---

## 4. 通信流程

### 4.1 单次迭代通信流程（P2P）

```
1. Pack left face   -> sendBufLeft
2. Pack right face  -> sendBufRight
3. Pack bottom face -> sendBufDown
4. Pack top face    -> sendBufUp
5. MPI_Startall                        — 启动 12 个持久化请求（Send_init/Recv_init 预注册）
6. [可选] 内点计算（与通信重叠）
7. MPI_Waitall                         — 等待所有通信完成
9. Unpack recvBuf -> halo cells
```

### 4.2 通信-计算重叠

```
迭代开始
├── beginExchange()          // 打包 6 面 + 12 个非阻塞请求（halo_exchange）
├── 内点盒计算                // 不依赖 halo 的区域（stencil_interior）
├── endExchange()            // MPI_Waitall + 解包（halo_wait，仅 overlap 模式单独记录）
├── 边界带计算                // 6 个 slab 补齐（stencil_boundary）
├── 融合残差 + Allreduce       // 更新循环内固定区域序累积，与 overlap 无关
└── swapU()，完成一次迭代
```

---

## 5. 关键设计决策

| 决策 | 选择 | 理由 |
|---|---|---|
| 内存分配 | 对齐分配（AlignedBuffer） | 防止 false sharing；关键路径零动态分配 |
| 求解器接口 | 策略模式 | 支持 Jacobi → CG/SOR 的平滑替换 |
| 通信接口 | 策略模式 | 支持 P2P → Collective → RMA 的平滑替换 |
| 数据布局 | SoA（结构体数组） | Stencil 访问局部性好，利于 SIMD |
| 索引计算 | 行主序平坦索引 | 与 C/C++ 数组语义一致，编译器优化友好 |
| 异常安全 | RAII + 关键路径 noexcept | 避免 MPI 通信路径中的异常 |
| 日志系统 | 单例 + 分级过滤 | 避免日志开销影响性能 |

---

## 6. 扩展性预留

> 本表所列能力均为**路线图 / 未实现**；当前版本的实际实现范围以上文 §2 模块架构为准。

| 预留扩展 | 接口位置 | 实现复杂度 |
|---|---|---|
| CG/SOR 求解器（CG 已实现） | `solver/solver.hpp` | 中（SOR 需松弛因子；CG 预条件子） |
| Red-Black GS（已实现） | `solver/solver.hpp` | 低 |
| RMA 通信 | `comm/halo_exchanger.hpp` | 中（MPI_Win_create + 同步） |
| Hilbert 曲线划分 | `grid/partition.hpp` | 高（空间填充曲线算法） |
| HDF5 输出 | `io/io_backend.hpp` | 低（链接 HDF5 库） |
| PAPI 计数器 | `perf/` | 低（可选编译） |
| 自适应网格 | `grid/` | 高（八叉树/AMR） |

---

## 7. 构建系统

CMake 配置支持：
- 多编译器（GCC/Clang/Intel；MSVC 显式报错，项目为 POSIX-only）
- 多模式（Release / Debug）
- 可选依赖（OpenMP、GoogleTest）
- Release 使用 `-march=native`（GNU/Clang）/ `-xHost`（Intel）；Debug 自动启用 Address/UB sanitizers（GNU/Clang）

---

*Last updated: 2025*
