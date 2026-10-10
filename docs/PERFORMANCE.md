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
    #pragma omp simd   // 注意：内点带偏移，不保证 64B 对齐，勿用 aligned 断言
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

## 10. AR003 MPI-IO 单文件输出基准（io_write 剖面）

取数口径：`--enable-profiling` 报告的 `io_write` 区段 `total_sec`（rank 0）；256² 网格、`--max-iter 100 --save-interval 10`（累计 10 次写入，每次 512 KB）、OMP=1、WSL 单机、Release 构建、每组合 3 次取中位数。

| 输出格式 | np=1 | np=4 |
|----------|------|------|
| `mpibin`（单文件，collective） | **0.0070 s** | 0.0125 s |
| `binary`（分片 + 索引） | 0.0126 s | **0.0049 s** |

- 单次写入量仅 512 KB，耗时以固定开销（open/view/close、collective 同步）为主，绝对差异在毫秒级；
- np=1 时 mpibin 略快（单文件直写，无分片索引）；np=4 时分片写各 rank 独立文件更快，单文件 collective 需额外的聚合同步，此为单文件自描述便利性的已知代价；
- 复测方法：`scripts/bench_t007.sh`（WSL 内执行），原始数据见 `specs/archive/AR003-mpiio-single-file-output/evidence/T007-bench.log`。

---

## 11. AR004 诚实性修复 + CG 并行化基准（before/after）

取数口径：256²、`--max-iter 500 --tol 0.0`（固定迭代数）、OMP=1、WSL 单机、Release 构建、每组合 3 次取中位数；before = 改动前基线（HEAD aa5c282），after = T001-T007 落地后（HEAD 28685e8）。详见 `specs/changes/AR004-honesty-performance-pack/evidence/`。

| case | before ms/iter | after ms/iter | Δ | 说明 |
|------|---------------|---------------|---|------|
| jacobi np=1 | 0.0701 | 0.0729 | +4.0% | 真实残差判据≈零开销（融合换算 r=D·diff，无显式扫描） |
| jacobi np=1 k=10 | — | 0.0705 | — | `--residual-check-interval 10`：省 9/10 Allreduce |
| rbgs np=1 k=1 | 0.1060 | 0.1422 | +34% | 诚实成本：每 iter 真实残差扫描（exchange+核+Allreduce） |
| rbgs np=1 k=5 | — | 0.1093 | +3.1%（vs before k=1） | interval 摊薄扫描成本 |
| rbgs np=1 k=10 | — | 0.1089 | +2.7%（vs before k=1） | k=10 基本收回诚实化成本 |
| cg np=1 | 0.1882 | 0.1938 | +3.0% | 噪声量级 |
| cg np=4 | 0.0565 | 0.0524 | −7.3% | B2a：updatePInterior/初始化循环并行化收益 |

- 收敛口径变化（64² tol=1e-6/1e-7）：jacobi 16141→17327、rbgs 8368→8812（迭代数——before 为递推残差口径，未真正达到 tol；after 为真实残差口径）；cg 198→198（不变，佐证数值一致性）。
- np=4 数据在本机跨会话漂移 ±15%（mpirun oversubscribe + 共享 VM），仅作趋势参考；门槛判定以 np=1 为准（见 evidence/ar004-bench.md）。
- 复测方法：`scripts/bench_ar004.sh`（含门槛判定），before 侧 `scripts/bench_ar004_baseline.sh`。

---

## 12. AR005 I/O 写出路径基准（ASCII→appended binary / 逐元素→块写）

取数口径：`--max-iter 10 --save-interval 10`（io_write 剖面区段，循环内写出）、np=1、OMP=1、WSL 单机、Release 构建、每组合 3 次取中位数；baseline = 改动前（ASCII vtk + 逐元素 binary，HEAD 69305cb 前），after = T002/T003 落地后。门槛（design D6）：median(after) ≤ median(baseline)。复测：`scripts/bench_ar005.sh baseline|after`；详见 `specs/changes/AR005-vtk-binary-ci-trust/evidence/`。

| nx | format | baseline sec / bytes | after sec / bytes | 提速 | 体积 |
|----|--------|---------------------|-------------------|------|------|
| 256 | vtk | 0.019120 / 1,245,572 B | 0.000528 / 524,750 B | **36.2×** | **-57.9%** |
| 256 | binary | 0.001253 / 524,344 B | 0.000374 / 524,344 B | 3.35× | 0（布局不变） |
| 512 | vtk | 0.065241 / 4,981,124 B | 0.001339 / 2,097,614 B | **48.7×** | **-57.9%** |
| 512 | binary | 0.005386 / 2,097,208 B | 0.001422 / 2,097,208 B | 3.79× | 0（布局不变） |

- 四组合全部满足 D6 门槛；`.vti` 体积 = 数据镜像 + ~406B XML 头（与 binary 后端一致），`.bin` 字节数逐位不变（U2 布局守护测试佐证纯性能重构）。
- 提速来源：vtk 每元素 ~24 字符文本格式化 → 每 (k,j) 行一次块写；binary 每 1 元素 1 次 write syscall（256² = 65,536 次）→ 每行 1 次（256 次）。
- ParaView 兼容口径：本地无 GUI 环境，以 U3 最小解析器测试（位级还原）+ VTK XML 规范走查替代（srs §4），未做真实加载演示。

## 13. AR006 halo 通信范式基准（pack p2p vs datatype 直传 vs 真集合 alltoallw）

取数口径：Jacobi 256²、`--tol 0.0 --residual-check-interval 10000`（残差开销旁路）、`halo_exchange` 剖面区段（默认非 overlap 路径：begin+end 均在区段内）、np=1/4（np4 = 2×2 切分，OMP=1 per rank）、WSL 单机、Release 构建。复测：`scripts/bench_ar006.sh`（200 迭代口径）。

**测量稳定性声明（诚实优先）**：本负载单次交换仅 ~5-15µs，与 WSL 单机运行间漂移同量级——三轮测量（下表）排序不稳定，**不下「谁更快」的结论**；2000 迭代大负载 ×3 中位为最稳定口径。

第 1 轮（T006 首测，`bench_ar006.sh` 200 迭代 ×3 中位）：

| np | comm_mode | halo_exchange sec（200 迭代累计） | 相对 p2p |
|----|-----------|-----------------------------------|----------|
| 1 | p2p | 0.000030 | 1.00× |
| 1 | datatype | 0.000025 | 0.83× |
| 1 | collective | 0.000148 | 4.93× |
| 4 | p2p | 0.000936 | 1.00× |
| 4 | datatype | 0.001340 | 1.43× |
| 4 | collective | 0.002195 | 2.34× |

第 2 轮（ST 复测，同脚本同口径）——np4 排序与第 1 轮不同：p2p 0.001494 / datatype 0.001895（1.27×）/ collective 0.001303（**0.87×**）；np1 三模式与第 1 轮一致（collective 空图固定开销 149µs 两轮稳定）。

第 3 轮（ST 稳定性口径，2000 迭代 ×3 中位，np4）：

| comm_mode | halo_exchange sec（2000 迭代累计） | 相对 p2p |
|-----------|-----------------------------------|----------|
| p2p | 0.013856（raw 9.9-15.4ms） | 1.00× |
| datatype | 0.014757（raw 14.1-21.4ms） | 1.06× |
| collective | 0.015803（raw 11.2-16.9ms） | 1.14× |

- **结论（记录式，D6）**：本机该负载（小面 128 列 × hw=1）下三范式 halo 交换耗时**同量级且无稳定优劣**——范式间差异（≤~15%）与单机运行间漂移（rep 间 ±30%）不可区分，第 1 轮 2.34× 触发的 D6「劣化超 2×」评审经复测判定为噪声，无回退动作（默认 `p2p` 本就未动，D1）；datatype/collective 定位为**范式可选与正交性验证**（`--comm-mode` 三值），不作收益宣称。GUIDE C1「小面 datatype 可能更慢」在本机未获稳定印证，也未反转出稳定优势——更大面/更大 halo 宽度未测，不外推。
- np=1 为无面对照（全 PROC_NULL）：p2p 与 datatype 接近零开销；collective 的 ~150µs/200 迭代为每迭代对空图 alltoallw 的固定调用开销（两轮一致，稳定可测）。
- 正确性对照：U4（3D hw=2，P2P vs Datatype 双 exchanger halo 带逐元素一致）、B1-B4（np4 非均匀/np8 3D 两新范式）、B5/E1/E2（collective 与 p2p 等价、三 exchanger np=1 安全、begin/end 拆分）全绿——三范式行为等价，性能差异纯属实现路径。

---

## 14. AR007 两层多重网格（mg2）vs 基线迭代法（256² 制造解）

取数口径：256² 制造解负载（`-Laplace u = 2π²·sin(πx)sin(πy)`，Dirichlet 0）、`--tol 1e-2`（绝对真残差 L2，main 语义同既有求解器）、np=1/4（np4 = 2×2 切分，OMP=1 per rank）、WSL 单机、Release 构建、3 次取中位。复测：`scripts/bench_ar007.sh`。

| np | solver | 迭代数 | 平滑步数（口径见注） | 耗时 sec |
|----|--------|--------|---------------------|----------|
| 1 | red_black_gs | 85540 | 85540（×2 sweep/iter） | 11.63 |
| 1 | jacobi | 166440 | 166440（×1 sweep/iter） | 9.97 |
| 1 | **mg2** | **14 cycles** | **56**（=cycles×(ν1+ν2)=×4） | **0.102** |
| 4 | red_black_gs | 83021 | 83021（×2） | 4.85 |
| 4 | jacobi | 161403 | 161403（×1） | 4.11 |
| 4 | **mg2** | **14 cycles** | **56** | **0.103** |

- **结论（记录式，GUIDE B1a 验收①）**：同 tol 下 mg2 总平滑步数 56 « red_black_gs 85540（≈1500×），耗时 0.10s vs 11.6s（np1 ≈114×）。GUIDE 验收①「迭代数（平滑步计）与 2 层理论定性一致」成立：mg2 收敛所需 cycle 数与网格规模基本无关（同 tol 1e-6 实测 64²/128²/256² = 20/22/23 cycles），而基线迭代数随 n² 增长（RBGS 256² 在 tol 1e-3 下超过 10 万迭代未收敛）。
- **收敛率**：渐近每 cycle 残差 ~0.3-0.4×（64² np1 实测前 12 cycles 逐 cycle ~0.30；tol 扫描 1e-3→1e-8 每数量级约 +2 cycles），与两层 MG 理论定性一致（ν1=ν2=2、粗层 CG 精解）。
- **适用 tol 范围注记（诚实边界）**：粗层容差取外层 tol（绝对口径），当粗 rhs 量级 ‖4·R(r)‖ 逼近 tol 时粗层 CG 提前停止（"饥饿"），格式退化为纯平滑——实测 256² tol 1e-9、64² tol 1e-11 时 5000 cycles 不收敛。有效区间为 tol ≳ ~1e-8（256²；×4 尺度补偿实际把该区间扩大了 4 倍）。后续 AR 可考虑粗层相对容差。
- **实现要点（正确性前提，非收益宣称）**：粗层方程为 rediscretized 五点模板 + 尺度补偿因子 4（=(H/h)²，源码注释有推导）——缺该因子时粗解 e_H≈e/4，cycle 退化为纯平滑（实测 cycles 随 n² 增长、率 ~0.75/cycle，探针验证 0.258 系数后修复）。
- 口径注：red_black_gs 每「iteration」= 红+黑 2 sweeps；jacobi 每 iteration = 1 sweep；mg2 每 cycle = ν1+ν2 = 4 sweeps + 1 次粗层 CG（256² np1 实测逐 cycle 粗层 CG 迭代 205→8 渐降、中位 ~180、14 cycles 合计 ~2300）。三者的「迭代数」不可直接互换，比较以「平滑步数」与「耗时」为准。
- np4 与 np1 的 mg2 cycle 数一致（14 vs 14），符合复制式粗层设计（每 rank 冗余串行 CG，粗层无跨 rank 通信差异）。

## 15. AR008 多层多重网格（mgv W-cycle / mgcg MG-CG 预条件）vs mg2/裸 CG（256²/512² 制造解）

取数口径：制造解负载（同 §14）、`--tol 1e-6`（绝对真残差 L2）、np=1/4（np4 = 2×2 切分，OMP=1 per rank）、WSL 单机、Release 构建、3 次取中位。复测：`scripts/bench_ar008.sh`。

| np | 规模 | solver | 迭代数 | 平滑步数 | 粗根 CG 迭代 | 耗时 sec |
|----|------|--------|--------|----------|--------------|----------|
| 1 | 256² | cg | 700 | — | — | 0.122 |
| 1 | 256² | mg2 | 23 cycles | 92 | — | 0.232 |
| 1 | 256² | **mgv** | **15 cycles** | 300 | 14355 | **0.132** |
| 1 | 256² | **mgcg** | **8 iter** | — | — | **0.076** |
| 4 | 256² | cg | 834 | — | — | 0.037 |
| 4 | 256² | mg2 | 23 cycles | 92 | — | 0.306 |
| 4 | 256² | **mgv** | **15 cycles** | 300 | 14357 | 0.114 |
| 4 | 256² | **mgcg** | **8 iter** | — | — | **0.061** |
| 1 | 512² | cg | 1378 | — | — | 1.303 |
| 1 | 512² | mg2 | 25 cycles | 100 | — | 2.260 |
| 1 | 512² | **mgv** | **14 cycles** | 336 | 26815 | **0.389** |
| 1 | 512² | **mgcg** | **8 iter** | — | — | **0.254** |
| 4 | 512² | cg | 1697 | — | — | 0.445 |
| 4 | 512² | mg2 | 24 cycles | 96 | — | 2.874 |
| 4 | 512² | **mgv** | **14 cycles** | 336 | 26758 | 0.338 |
| 4 | 512² | **mgcg** | **8 iter** | — | — | **0.212** |

- **结论（记录式，GUIDE B1b 验收②）**：迭代数口径——mgcg 8 迭代在 256²/512² **恒定**（规模无关），vs 裸 cg 700/1378（≈87×/172× 迭代数之差）；mgv cycle 数 15/14（4× 细化反降 1，与 U4 判据 |cyc512−cyc256|≤4 一致）。耗时口径——np1 512²：mgcg 0.254s vs cg 1.303s（≈5.1×）、mgv 0.389s vs mg2 2.260s（≈5.8×）；256² np1 mgcg 0.076s 为全场最优。
- **W-cycle 结构（D11）**：mgv 每层 2 次粗修正（每 cycle 平滑步=4×层数：256² 24 步/cycle×15=300；粗根 CG 调用 2^(层数−1)=32 次/cycle，合计 14355/26815 迭代——粗根 8² 极廉价）。**结构注记（诚实记录）**：单 V-cycle（γ=1）下逐层 ×4 补偿对插值杂散过补偿，深层链发散（256² 齐次 [1,1] 稳态率 ×2.01/cycle、512² ×2.69）；W 的第二次粗修正基于刷新残差清洗杂散污染，恢复层数无关收敛（齐次稳态率 ~1e-4）。代价：粗层工作量翻倍（总耗时仍优于 mg2，见上表）。
- **mg2 vs mgv 对比叙事**：mg2 cycle 数随规模缓增（23→25），mgv 恒定（15→14）；mgv 平滑步总量高于 mg2（300 vs 92，W 结构性翻倍）但每步更「值」（cycle 率 ~0.1 vs 0.3-0.4），净耗时 mgv 在 256²/512² 均优（0.132/0.232、0.389/2.260）。
- **mgcg 预条件子质量**：8 迭代中每迭代含 1 个 W-cycle 预条件步——预条件后条件数 O(1)（迭代数规模无关即证）；np4 与 np1 迭代数一致（8 vs 8），粗层复制式冗余对迭代轨迹零影响（U9 轨迹一致性）。
- **np4 观察**：mg2 np4 比 np1 慢（复制式粗层冗余，AR007 §14 已记录）；mgv/mgcg np4 与 np1 耗时相当（细层并行收益抵消粗层冗余）——多层链的粗层占总工作量比例更低。
- 口径注：cg 每「iteration」=1 次 matvec+2 次内积；mg2 每 cycle=4 sweeps+1 次粗层 CG（两层）；mgv 每 cycle=4×层数 sweeps+2^(层数−1) 次粗根 CG；mgcg 每 iteration=1 次 matvec+3 次内积+1 个完整 W-cycle 预条件步。四者「迭代数」不可互换，比较以耗时为准。

## 16. AR009 管线化共轭梯度（pcg）vs 裸 CG（256²/512² 制造解）

取数口径：制造解负载（同 §14/§15）、`--tol 1e-6`（绝对真残差 L2）、np=1/4（np4 = 2×2 切分，OMP=1 per rank）、WSL 单机、Release 构建、3 次取中位。复测：`scripts/bench_ar009.sh`。

| np | 规模 | cg 迭代数 | pcg 迭代数 | cg 耗时 sec | pcg 耗时 sec | 耗时差 |
|----|------|-----------|------------|-------------|--------------|--------|
| 1 | 256² | 700 | **700** | 0.185 | 0.221 | +19% |
| 4 | 256² | 834 | **834** | 0.097 | 0.101 | +4% |
| 1 | 512² | 1378 | **1378** | 1.587 | 1.689 | +6% |
| 4 | 512² | 1697 | **1697** | 0.784 | 0.790 | +0.7% |

- **结论（记录式，B1/R1①）**：**迭代数 4 组全部精确相同**（700/700、834/834、1378/1378、1697/1697）——CG-1（Chronopoulos–Gear）与标准 CG 同 Krylov 空间、同收敛迭代数在实测中零差异（±5% 判据含 512² 半项以 0 差异通过）；np4 vs np1 迭代数差（834 vs 700）系归约分块序不同的轨迹分岔（同 §E3 evidence），非算法差异。
- **耗时（诚实记录，单机下界）**：单机共享内存下 pcg 无 wall 收益（+0.7%~+19%）——同步次数减半（2 阻塞→1 非阻塞融合，µs 级）的收益被 **5 向量工作集（vs cg 3）**与 ν-递推开销抵消；np4 大规模下差距收窄（512²np4 +0.7%）。**收益叙事按 D1 诚实声明**：每迭代 MPI 集合调用 2→1（消息数同减半）+ 非阻塞语义（MPI 进度引擎不被独占）；真正的归约-计算重叠需 CG-2 两步前瞻（数值稳定性劣化，列后续工作）。跨节点收益未测（单机环境限制，不外推）。
- **归约计数（P3 实测佐证）**：pcg 每迭代 `pcg_iallreduce` region 计数=iters+1（init 1+每迭代 1）、`cg_blocking_allreduce`==0；cg 阳性对照==2·iters+1（含 init rho）。GUIDE B2b 原「3→1」口径修正为「**2 阻塞→1 非阻塞融合**」（原文 3 系把一次性 init rho 计入）。
- 口径注：pcg 每「iteration」=1 次 matvec（A·r）+1 次打包双标量 Iallreduce + 4 组 axpy/方向更新（q=A p 递推免二次 matvec）；cg 每 iteration=1 次 matvec+2 次阻塞 dot 归约。二者迭代数语义相同（Krylov 步数）。

---

*Last updated: 2026*
