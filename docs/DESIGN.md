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
- **memory_pool.hpp**: Bump Allocator，用于临时网格分配，避免重复 malloc/free
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
  - 索引：全局索引到本地索引的偏移映射
  - 邻居：通过 `MPI_Cart_shift` 预计算的邻居 rank

### 2.3 求解器层（solver/）

- **PoissonSolver**: 抽象基类，定义 `solve()` 和 `iterate()` 接口
- **JacobiSolver**: 经典 Jacobi 实现
  - 2D 五点 stencil：`u_new = 0.25 * (u_left + u_right + u_down + u_up - dx^2 * f)`
  - 3D 七点 stencil：除以 6.0
  - OpenMP 并行：`#pragma omp parallel for` + `#pragma omp simd`
  - 残差计算：内点 L2 范数，通过 `MPI_Allreduce` 全局聚合
  - 通信-计算重叠（可选，`JacobiSolver(true)`）：`beginExchange` → 内点盒计算 → `endExchange` → 边界 6 slab；残差使用**独立全内点扫描**，保证 overlap on/off 的迭代数与解 bit 级一致（`lastResidual()` 暴露最近残差）

### 2.4 通信层（comm/）

- **HaloExchanger**: 抽象基类，定义 `beginExchange()` / `endExchange()` 接口
- **PointToPointExchanger**: 非阻塞 MPI 实现
  - `MPI_Isend` / `MPI_Irecv` + `MPI_Waitall`
  - 6 个方向（左/右/下/上/后/前）共 12 个非阻塞请求；tag 按发送面编号 0-5
  - 每个方向的独立 send/recv 缓冲区；仅交换内点面范围（2D 数据位于 k=0 平面）
- **CollectiveExchanger**: P2P 回退实现（委托 `PointToPointExchanger` 并记录告警）；真 `MPI_Neighbor_allgatherv` 为路线图

### 2.5 性能层（perf/）

- **Timer**: `std::chrono::steady_clock` 高精度计时
- **Profiler**: 分层区域计时器，支持嵌套区域
  - `ProfileScope`: RAII 自动 begin/end
- **Reporter**: 性能报告生成器
  - JSON 格式：结构化、可扩展
  - CSV 格式：便于批量数据分析

### 2.6 I/O 层（io/）

- **IOBackend**: 抽象接口
- **BinaryIOBackend**: 原始二进制输出，最高性能
- **VTKIOBackend**: VTK 结构化点格式，支持 ParaView 可视化

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
5. MPI_Isend (left, right, down, up)   — 非阻塞发送
6. MPI_Irecv (left, right, down, up)   — 非阻塞接收
7. [可选] 内点计算（与通信重叠）
8. MPI_Waitall                         — 等待所有通信完成
9. Unpack recvBuf -> halo cells
```

### 4.2 通信-计算重叠

```
迭代开始
├── beginExchange()          // 打包 6 面 + 12 个非阻塞请求（halo_exchange）
├── 内点盒计算                // 不依赖 halo 的区域（stencil_interior）
├── endExchange()            // MPI_Waitall + 解包（halo_wait，仅 overlap 模式单独记录）
├── 边界带计算                // 6 个 slab 补齐（stencil_boundary）
├── 独立残差扫描 + Allreduce  // 与 overlap 无关，保证 on/off 结果一致
└── swapU()，完成一次迭代
```

---

## 5. 关键设计决策

| 决策 | 选择 | 理由 |
|---|---|---|
| 内存分配 | 对齐分配 + 内存池 | 防止 false sharing，NUMA 友好 |
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
| CG/SOR 求解器 | `solver/solver.hpp` | 中（需要预处理子、Krylov 子空间） |
| Red-Black GS | `solver/solver.hpp` | 低（修改迭代顺序即可） |
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
