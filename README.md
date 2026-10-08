# HyPoS — Hybrid Poisson Solver

> **面向大厂 AI Infra 的高性能分布式求解器**
>
> HyPoS 是一个基于 Jacobi 迭代的混合并行（MPI + OpenMP）泊松方程求解器，展示分布式并行计算、域分解、通信优化、内存管理和性能可观测性等核心 AI Infra 技能。

---

## 项目亮点

| 技术模块 | 对应 AI Infra 技能 | 实现细节 |
|---|---|---|
| **MPI 域分解** | 数据并行 (Data Parallelism) | 2D/3D Cartesian 拓扑，自动负载均衡 |
| **Halo Exchange** | 集合通信 (Collective Communication) | 持久化非阻塞 MPI（6 方向、12 请求、memcpy 打包）；`collective` 模式为 P2P 回退（见路线图） |
| **通信-计算重叠** | 异步流水线 | 内点计算与边界通信重叠 |
| **OpenMP 加速** | 多线程 Kernel | SIMD 向量化、First-touch 策略 |
| **内存管理** | 高效内存管理 | 64 字节对齐 RAII（AlignedBuffer）、关键路径零动态分配 |
| **性能可观测性** | Profiling & Tracing | 分层计时、JSON/CSV 报告、自动化扩展性测试 |
| **求解器可扩展** | 算法与框架解耦 | Jacobi / Red-Black GS / CG（`--solver` 切换，策略模式接口） |

---

## 快速开始

### 依赖

- **CMake** ≥ 3.16
- **MPI** 实现：OpenMPI、MPICH 或 Intel MPI（注意：Ubuntu 24.04 的 `mpich 4.2.0` 包存在 PMI/PMIx 不匹配缺陷，应用会静默退化为单进程运行；该环境请使用 OpenMPI，可用 `mpirun -n 4 <app>` 验证 `MPI_Comm_size` 是否为 4）
- **OpenMP** 支持（GCC/Clang/Intel 编译器）
- **C++17** 编译器：GCC 10+、Clang 14+、Intel oneAPI 2023+
- 可选：GoogleTest（测试）、PAPI（硬件计数器）

### 编译

```bash
# 标准编译（Release，自动检测 OpenMP）
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# 调试模式（GCC/Clang 下启用 Address/UB sanitizers）
cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug -j$(nproc)
```

### 运行

```bash
# 串行运行（单进程）
mpirun -np 1 ./build/hypos --nx 256 --ny 256 --max-iter 10000

# 并行运行（4 进程，每进程 2 线程）
mpirun -np 4 ./build/hypos --nx 512 --ny 512 --omp-threads 2 --max-iter 10000

# 大规模运行（16 进程，自动线程）
mpirun -np 16 ./build/hypos --nx 2048 --ny 2048 --max-iter 10000 --enable-profiling
```

### 运行参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--nx, --ny, --nz` | 1024, 1024, 1 | 全局网格尺寸 |
| `--halo-width` | 1 | 幽灵层宽度 |
| `--solver` | jacobi | 求解器类型（jacobi / red_black_gs / cg） |
| `--max-iter` | 10000 | 最大迭代次数 |
| `--tol` | 1e-6 | 收敛容差 |
| `--bc` | dirichlet | 物理边界条件（dirichlet / neumann） |
| `--omp-threads` | 系统核心数 | OpenMP 线程数 |
| `--comm-mode` | p2p | 通信模式（p2p / collective；collective 当前为 P2P 回退） |
| `--enable-profiling` | false | 启用详细性能分析 |
| `--overlap-comm` | false | 启用通信-计算重叠 |
| `--output-format` | json | 输出格式（json / csv / vtk / binary / mpibin；vtk=每 rank `.vti` 分片+rank0 `.pvti` 索引，binary=每 rank `.bin` 含真实 offsets，mpibin=MPI-IO 单文件 `solution_<step>.bin`：72 字节自描述头+全局行主序数据区，全 rank 集体写，小端） |
| `--output-dir` | ./output | 输出目录（启动时自动创建） |
| `--save-interval` | 0 | 每 N 次迭代输出中间解（文件名含步号；0=不保存） |

---

## 项目结构

```
HyPoS/
├── CMakeLists.txt              # 主构建配置
├── src/
│   ├── main.cpp                # 程序入口
│   ├── core/
│   │   ├── types.hpp           # 基础类型（Real, Index）
│   │   ├── aligned_buffer.hpp  # 64 字节对齐内存分配器
│   │   └── exception.hpp       # 异常体系
│   ├── grid/
│   │   ├── grid.hpp            # 全局网格定义
│   │   ├── subgrid.hpp         # 子域网格（含 halo 管理）
│   │   ├── partition.hpp       # 划分策略接口
│   │   ├── subgrid.cpp         # 子域实现
│   │   └── partition.cpp       # 均匀划分实现
│   ├── solver/
│   │   ├── solver.hpp          # 求解器策略接口（含进度回调/lastResidual）
│   │   ├── jacobi_solver.cpp   # Jacobi（融合残差、重叠、SIMD）
│   │   ├── red_black_gs_solver.cpp # Red-Black GS
│   │   └── cg_solver.cpp       # 共轭梯度（CG）
│   ├── comm/
│   │   ├── halo_exchanger.hpp  # 通信抽象接口
│   │   ├── p2p_exchanger.cpp   # 点对点非阻塞通信
│   │   └── collective_exchanger.cpp # 集体通信
│   ├── perf/
│   │   ├── timer.hpp           # 高精度计时器
│   │   ├── profiler.hpp        # 分层性能分析器
│   │   ├── reporter.hpp        # 报告生成器
│   │   ├── timer.cpp
│   │   ├── profiler.cpp
│   │   └── reporter.cpp
│   ├── io/
│   │   ├── io_backend.hpp      # I/O 抽象接口
│   │   ├── binary_io.cpp       # 二进制输出（每 rank 分片）
│   │   ├── vtk_io.cpp          # VTK 格式输出（.vti/.pvti）
│   │   └── mpiio_binary.cpp    # MPI-IO 单文件二进制输出（mpibin）
│   └── utils/
│       ├── mpi_env.hpp         # MPI 环境 RAII 封装
│       ├── logger.hpp          # 分级日志
│       ├── cmdline_parser.hpp  # 命令行解析
│       ├── mpi_env.cpp
│       ├── logger.cpp
│       └── cmdline_parser.cpp
├── tests/                       # 单元测试与集成测试
├── scripts/
│   ├── run_scaling_tests.py    # 自动化扩展性测试
│   ├── plot_scaling.py         # 生成 Scaling 曲线图
│   └── benchmark.sh            # 一键 benchmark 脚本
├── docs/
│   ├── DESIGN.md               # 架构设计文档
│   ├── PERFORMANCE.md          # 性能调优指南
│   └── SCALING_REPORT.md       # 扩展性报告模板
└── .github/workflows/ci.yml     # CI/CD 配置
```

---

## 测试

```bash
# 构建并运行全部测试（含 4/8 进程 MPI 用例，经 ctest 自动用 mpirun 注册）
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure

# 小规模集成冒烟
mpirun -np 4 ./build/hypos --nx 128 --ny 128 --max-iter 100
```

---

## 性能分析

### 自动化扩展性测试

```bash
# 强扩展性测试（固定网格，增加进程数）
python3 scripts/run_scaling_tests.py --output-dir ./scaling_results --max-procs 16

# 弱扩展性测试（固定每进程网格，增加进程数）
python3 scripts/run_scaling_tests.py --output-dir ./scaling_results --max-procs 16 --weak-only

# 生成曲线图
python3 scripts/plot_scaling.py --input-dir ./scaling_results --output-dir ./scaling_plots
```

### 一键 Benchmark

```bash
chmod +x scripts/benchmark.sh
./scripts/benchmark.sh
```

### 性能输出示例

```json
{
  "run_id": "hypo_20250115_143052",
  "config": {
    "nx": 1024,
    "ny": 1024,
    "nz": 1,
    "mpi_procs": 4,
    "omp_threads": 4,
    "solver": "jacobi",
    "max_iter": 10000,
    "tolerance": 1e-06,
    "overlap_comm": false
  },
  "performance": {
    "total_time_sec": 1.2345,
    "iter_time_ms": 0.123,
    "iterations": 2847,
    "flops_per_sec": 1.85e+10
  }
}
```

---

## 架构设计

详见 [docs/DESIGN.md](docs/DESIGN.md)。

核心架构图：

```
+-------------+     +----------------+     +------------------+
|   main.cpp  | --> |   PoissonSolver  | --> |   Subgrid        |
+-------------+     |   (JacobiSolver) |     |   (local + halo) |
     |              +----------------+     +------------------+
     |                     |                      |
     v                     v                      v
+-------------+     +----------------+     +------------------+
|  MPIEnv     |     |  HaloExchanger |     |  GridPartition   |
|  (RAII)     |     |  (P2P/Collect) |     |  (Uniform)       |
+-------------+     +----------------+     +------------------+
     |                     |                      |
     v                     v                      v
+-------------+     +----------------+     +------------------+
|  Logger     |     |  Profiler      |     |  IOBackend       |
|  Reporter   |     |  Timer         |     |  (VTK/Binary/    |
|             |     |                |     |   MPI-IO 单文件) |
+-------------+     +----------------+     +------------------+
```

---

## 性能调优

详见 [docs/PERFORMANCE.md](docs/PERFORMANCE.md)。

关键优化点：
- **Cache 优化**：Stencil 访问按 x 方向优先（连续内存），提高预取效率
- **SIMD 向量化**：`#pragma omp simd` 确保编译器生成 AVX2/AVX-512 指令
- **NUMA 感知**：First-touch 策略确保内存分配在计算核心本地
- **False Sharing 避免**：64 字节对齐 + 缓存行对齐
- **通信优化**：非阻塞 MPI 实现通信-计算重叠

---

## 扩展性报告

详见 [docs/SCALING_REPORT.md](docs/SCALING_REPORT.md)。

---

## 技术栈

| 层级 | 技术 | 版本 |
|---|---|---|
| 语言 | C++ | 17 |
| 并行 | MPI + OpenMP | 3.1 / 4.5 |
| 构建 | CMake | ≥ 3.16 |
| 测试 | GoogleTest | 可选 |
| 可视化 | matplotlib / ParaView | Python 3 |

---

## 路线图（当前版本未实现）

以下能力在早期文档中被提及，但**当前版本未实现**，在此明示避免误导：

- RMA（单边通信）Exchanger
- 真集合通信（`MPI_Neighbor_allgatherv`）——`--comm-mode collective` 当前委托 P2P 实现
- HDF5 输出、PAPI 硬件计数器
- SOR 求解器、CG 预条件子、Hilbert 曲线分区

---

## 许可证

MIT License — 详见 LICENSE 文件。

---

> **HyPoS** = **Hy**brid **Po**isson **S**olver
> 
> 构建用于展示大厂 AI Infra 工程能力的分布式高性能计算项目。
