# HyPoS 性能调优指南

> 本文档描述 HyPoS 的性能优化策略和调优参数。

---

## 1. 性能瓶颈分析

### 1.1 Roofline 模型

泊松求解器的主要计算模式是 stencil 更新：
- **计算强度**（Arithmetic Intensity）≈ 0.5 FLOP/byte（内存密集型）
- **理论瓶颈**：内存带宽，而非 CPU 计算
- **目标**：达到理论内存带宽的 60-80%

### 1.2 热点分析

通过 Profiler 输出，预期热点分布（定性指引，非本仓库实测）：

| 区域 | 占比 | 优化方向 |
|---|---|---|
| Stencil 内核 | 60-80% | SIMD、cache tiling |
| Halo 通信 | 15-30% | 重叠、消息聚合、减少进程数 |
| 残差计算 | 5-10% | 并行 reduction |
| I/O | <5% | 异步输出、减少频率 |

---

## 2. CPU Cache 优化

### 2.1 Stencil 访问模式

```cpp
// 优化前：列优先访问（不连续）
for (i = 0; i < nx; ++i)
  for (j = 0; j < ny; ++j)
    u_new[j][i] = 0.25 * (u[j][i-1] + u[j][i+1] + ...); // BAD

// 优化后：行优先访问（连续）
for (j = 0; j < ny; ++j)
  for (i = 0; i < nx; ++i)
    u_new[j*nx + i] = 0.25 * (u[j*nx + i-1] + u[j*nx + i+1] + ...); // GOOD
```

### 2.2 Loop Tiling（分块）

当网格大于 L2 cache 容量时，启用分块：
```cpp
// 分块大小：匹配 L1/L2 cache
constexpr Index BLOCK = 64; // 64x64 doubles ≈ 32KB, 可放入 L1

for (Index jj = jBegin; jj < jEnd; jj += BLOCK)
  for (Index ii = iBegin; ii < iEnd; ii += BLOCK)
    for (Index j = jj; j < std::min(jj+BLOCK, jEnd); ++j)
      for (Index i = ii; i < std::min(ii+BLOCK, iEnd); ++i)
        // stencil update
```

### 2.3 Halo 数据局部性

将 halo 数据与内点数据连续存储，减少 TLB miss：
- 当前实现已使用单一 `AlignedBuffer` 存储整个子域（含 halo）
- 访问模式：`u[(j+halo)*nx_total + (i+halo)]`

---

## 3. SIMD 向量化

### 3.1 编译器指令

```cpp
#pragma omp parallel for schedule(static)
for (Index j = jBegin; j < jEnd; ++j) {
    #pragma omp simd safelen(8) aligned(u, u_new: 64)
    for (Index i = iBegin; i < iEnd; ++i) {
        // stencil update
    }
}
```

### 3.2 验证向量化

使用编译器标志确认向量化成功：
```bash
g++ -O3 -fopenmp -fopt-info-vec-all ... 2>&1 | grep "loop vectorized"
```

### 3.3 手动 AVX-512 Intrinsic（可选）

对于极致性能，可手动编写 AVX-512 版本：
```cpp
#include <immintrin.h>
// 使用 _mm512_loadu_pd, _mm512_add_pd, _mm512_mul_pd
// 每次处理 8 个 double（512 bits / 64 bits）
```

---

## 4. NUMA 与内存亲和性

### 4.1 First-Touch 策略

确保内存分配和初始化由同一个线程执行：
```cpp
#pragma omp parallel for schedule(static)
for (Index j = jBegin; j < jEnd; ++j) {
    for (Index i = iBegin; i < iEnd; ++i) {
        u[idx] = 0.0; // 该线程 touch 该内存页
    }
}
```

### 4.2 numactl 绑核

```bash
# 绑定到 NUMA 节点 0（内存 + CPU）
numactl --cpunodebind=0 --membind=0 ./hypos ...

# 交错模式（适合内存带宽受限）
numactl --interleave=all ./hypos ...

# 查看 NUMA 拓扑
numactl --hardware
```

### 4.3 OpenMP 线程绑定

```bash
# Intel 编译器
export KMP_AFFINITY=granularity=fine,compact,1,0

# GCC
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

---

## 5. 通信优化

### 5.1 消息聚合

减少小消息数量，将多个 halo 面合并为单个大消息：
- 当前实现：每个方向独立通信
- 优化：将 left + right 合并为一条消息（使用 `MPI_Pack`）

### 5.2 通信-计算重叠效果

重叠效果与网格规模、网络延迟强相关：内点计算占比越高，可隐藏的 halo 等待越多。本仓库在单机 WSL 小规模下的通信占比很小（见 SCALING_REPORT），overlap 的主要价值在大规模多节点场景；此处不给出未实测的示意数值。

复测方法：`--overlap-comm --enable-profiling`，Profiler 的 `halo_wait` 与 `halo_exchange` 区域可量化等待时间，报告中的 `overlap_ratio` 给出隐藏比例。

### 5.3 减少通信频率

对于多步迭代算法（如 Jacobi），每步都需要 halo 交换：
- 无法减少频率（每步依赖边界值）
- 对于多网格或异步方法，可以异步推进 halo

---

## 6. 多线程竞争消除

### 6.1 False Sharing 避免

确保线程私有数据间隔至少 64 字节：
```cpp
// 错误：多个线程写入同一缓存行的不同变量
struct Bad { double x; double y; }; // 可能同一缓存行

// 正确：显式对齐
struct Good {
    alignas(64) double x;
    alignas(64) double y;
};
```

当前实现已使用 `AlignedBuffer`（64 字节对齐）避免此问题。

### 6.2 Lock-free 设计

- OpenMP 区域使用 `reduction` 而非 `critical`
- 全局残差使用 `MPI_Allreduce`（硬件优化）
- 无互斥锁、无原子操作（除日志外）

---

## 7. 编译器优化标志

### 7.1 GCC

```bash
-O3 -march=native -fopenmp -ffast-math
# -ffast-math: 允许浮点重排序，可能略微改变精度但提升性能
```

### 7.2 Intel oneAPI

```bash
-O3 -xHost -qopenmp -fma -ftz
# -xHost: 针对当前 CPU 生成指令
# -fma: 启用 FMA 指令（融合乘加）
# -ftz: 刷新非规格化数到零（避免性能陷阱）
```

### 7.3 Clang

```bash
-O3 -march=native -fopenmp -ffp-contract=fast
```

---

## 8. 性能调优检查清单

- [ ] 编译 Release 模式（`-O3 -DNDEBUG`）
- [ ] 启用 OpenMP（`-fopenmp`）
- [ ] 使用 `-march=native` 针对当前 CPU
- [ ] 验证 SIMD 向量化（`-fopt-info-vec`）
- [ ] 使用 `numactl` 绑核
- [ ] 设置 `OMP_PROC_BIND=close`
- [ ] 检查无 false sharing（`perf c2c`）
- [ ] 运行 strong scaling 验证效率 > 60%（8 进程）
- [ ] 运行 weak scaling 验证效率 > 80%（16 进程）

---

## 9. AR002 已实施优化与实测（256²×10000 迭代，OMP=1，WSL 单机）

| 优化 | 手段 | 说明 |
|------|------|------|
| 残差融合 | 取消独立残差扫描，按固定区域序在更新循环内累积部分和 | 省一遍全数组读；on/off 残留位级一致 |
| 持久化通信 | `MPI_Send_init/Recv_init` + `MPI_Startall/Waitall` | 降低每迭代通信开销 |
| 打包 memcpy 化 | 缓冲布局调整为最内维连续，行级 `memcpy` | 降低打包/解包成本 |
| SIMD 向量化 | `__restrict` + `#pragma omp simd`；`-fopt-info-vec` 确认 stencil 循环向量化（AVX2 32B） | 内核吞吐提升 |

**实测（`iter_time_ms`，3 次中位数，对比 AR001 基线）**

| 进程 | AR001 | AR002 | 变化 |
|------|-------|-------|------|
| 1 | 0.0451 | 0.02818 | **-37.5%** |
| 4 | 0.0204 | 0.01245 | **-39.0%** |

**求解器对比（64²，np=4，tol≈1e-7/1e-8）**

| 求解器 | 迭代数 | 说明 |
|--------|--------|------|
| Jacobi | 基线（1e-8 停止时 12967） | — |
| Red-Black GS | 52.3% × Jacobi | 理论比 ≈1/(1+ρ)≈0.5 |
| CG | 一般右端项 215 步；制造解为特征向量时 1 步精确收敛 | 一般问题降低 1-2 个数量级 |

---

*Last updated: 2026*
