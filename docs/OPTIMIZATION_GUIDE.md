# HyPoS 深度优化指南（OPTIMIZATION_GUIDE）

| 字段 | 内容 |
|------|------|
| 文档定位 | 后续优化 AR 的选题依据、设计输入与验收基线（活文档） |
| 基线版本 | commit `4cd7024`（AR003 归档后） |
| 创建日期 | 2026-10-08 |
| 维护规则 | 每完成一个优化 AR：勾选对应条目、回填实测数据、标注完成 AR 编号；行号证据随代码漂移，使用前须复核 |

---

## 1. 本文档的用途与使用方式

1. **选题**：新优化 AR 从 §5 路线表选取，或按 §6 决策规则从 §4 中新提；
2. **设计输入**：每个条目的「现状与证据」「优化方案」「验收判据」直接作为 srs.md/design.md 的素材起点（仍走完整 SDD 流程：srs → design → develop → review → ST → 归档）；
3. **验收对照**：AR 完成时对照条目验收判据 + §6 通用底线；
4. **纪律**：本文档自身不得含失实宣称——所有「已完成」必须有对应 AR 归档与测试/基准证据。

---

## 2. 项目定位与优化评估标尺

HyPoS 的本质定位是**「HPC 工程教学/参考样板」**：用一个小而完整的泊松求解器，示范 MPI+OpenMP 混合并行的全套工程实践（域分解、halo 通信、通信-计算重叠、多求解器策略、性能剖析、并行 I/O）。它不是生产级 PDE 库。

因此优化不是单一的 FLOP/s 最大化，而是四个维度的加权：

| 维度 | 判断问题 | 权重说明 |
|------|---------|---------|
| **示范价值** | 是否清晰示范一个 HPC 领域的重要技术？ | 定位主线 |
| **诚实性** | 宣称与实现是否零缺口？ | 教学样板最忌"假示范"，缺口优先级最高 |
| **扩展故事** | 是否强化"能 scale"的叙事与数据？ | 单机 WSL → 多节点就绪 |
| **复杂度预算** | 代码是否保持可读？ | 小项目中复杂度是净负债，每项须评估 |

---

## 3. 现状基线诊断（基于 commit 4cd7024）

### 3.1 已有优势（任何优化不得破坏）

- 持久化非阻塞 halo（6 方向 12 请求，`MPI_Send_init/Recv_init` + `Startall/Waitall`），memcpy 打包布局按最内维连续优化（src/comm/p2p_exchanger.cpp）
- Jacobi 通信-计算重叠流水线，overlap on/off **位级一致**（每 7 区域固定累积序，tests/test_alt_solvers.cpp 有专项断言）
- SIMD 化 stencil（`__restrict__` + `omp simd`，AR002 实测 iter_time -37.5%(np1)/-39%(np4)，见 PERFORMANCE §9）
- SoA + 64B 对齐 AlignedBuffer；swapU 指针交换零拷贝
- MPI-IO 单文件输出（AR003）：双 subarray filetype 零拷贝直写 + 自描述头
- 完整测试面：16 项 ctest（Debug+Release 双模式、ASan/UBSan 零报告）+ SDD 门控流程
- 分层架构与 AGENT_SPEC 约定（库代码禁 MPI_COMM_WORLD、MPI 路径不抛异常、row-major x 最快）

### 3.2 宣称 vs 实现缺口清单（G 表）

> 教学样板的诚实性问题。修复一项勾选一项。

| ID | 缺口 | 证据 | 状态 |
|----|------|------|------|
| G1 | **First-touch 宣称无实现**：README:16 与 AGENT_SPEC.md:54 宣称 first-touch 初始化，但 `zeroInitialize()` 是单线程 `std::fill`（subgrid.cpp:31-35），首次触摸全部由主线程完成——NUMA 语境下是教科书级反例 | subgrid.cpp:31-35 | ✅ 已完成（AR004/T005：外层维度 omp parallel for 三场单遍历全缓冲；applyDirichletBC 决策保持串行——setup 期单次、仅 halo 带） |
| G2 | **`--comm-mode collective` 为纯 P2P 委托** + WARN（collective_exchanger.cpp:6-27），真集合通信不存在 | collective_exchanger.cpp:12-13 | 未开始（对应 C2） |
| G3 | **tol 语义失实**：`--tol` 判据用**递推残差**（更新量平方和，jacobi_solver.cpp:47-48,64-65），非 ‖Au−f‖₂。SCALING_REPORT.md:92 记录恶果：tol 未达时 final_residual=204.6，"没收敛但停了"无从解释 | jacobi_solver.cpp:47-48 | ✅ 已完成（AR004/T001/T002：三求解器判据与 lastResidual 均为真实残差口径；Jacobi 融合换算≈零开销 +4%，RBGS 每 k 步扫描） |
| G4 | **性能叙事数据缺口**：强扩展 256² np4 加速比 2.21（效率 55%）、np4 通信占比 13%——"通信-计算重叠"的价值在当前规模下数据支撑不足 | SCALING_REPORT.md:86-105 | 持续项（B/C 系列改善后复测） |
| G5 | **RBGS 剖面区名误用** `"jacobi_iteration"`**（复制粘贴残留），观测数据失真 | red_black_gs_solver.cpp:101 | ✅ 已完成（AR004/T002 rbgs_iteration、T004 cg_iteration；U11 断言三求解器区名互斥） |
| G6 | **Profiler 非线程安全且与注释不符**：hpp:13 注释称 "atomic counters"，实现是 mutex 保护栈（profiler.hpp:54-56），无 per-thread 聚合，不可在 OpenMP 并行区使用 | perf/profiler.cpp:8-56 | ✅ 已完成（AR004/T007：per-thread ThreadData 注册表、热路径零锁、聚合时加锁；hpp 注释对齐真实机制；U8 多线程精确聚合） |

### 3.3 性能级实现缺口（P 表）

| ID | 问题 | 证据 | 状态 |
|----|------|------|------|
| P1 | **CG 的 p/r/u 更新循环是单线程标量循环**（无 OpenMP、无 SIMD）——CG 是"快 1-2 个量级"的招牌（tests/test_alt_solvers.cpp:449），但每迭代最规整的 axpy 循环未并行 | cg_solver.cpp:194-216, 251-273 | ✅ 已完成（AR004/T004：updatePInterior 并行化 + 初始化循环并行；np=4 实测 −7.3%，U6 golden 位级守护） |
| P2 | **CG 每迭代 3 次全局 Allreduce**（dot(r,r)/dot(p,Ap)/dot(r,r)），无 pipelined 化、无 `MPI_Iallreduce` 重叠 | cg_solver.cpp:159,169,178 | 未开始（对应 B2b） |
| P3 | **RBGS 内层 stride-2 循环无 `omp simd`**；每 sweep 2 次 halo 交换且无重叠——通信翻倍抵消部分收益 | red_black_gs_solver.cpp:30,50,71-81 | 未开始（对应 B4） |
| P4 | **残差 Allreduce 每迭代一次，无检查频率参数**——归约开销管理是 HPC 经典课题，项目把它写死 | jacobi_solver.cpp:169-175 | ✅ 已完成（AR004/T003：`--residual-check-interval N`（Jacobi/RBGS），RunConfig/JSON 透出；Jacobi k=10 np4 实测 −27%/iter） |
| P5 | **VTK 输出为 ASCII Float64**（512² 一次输出数百 MB 文本）；binary 后端逐元素 `ofstream.write` 非块写 | vtk_io.cpp:50-60；binary_io.cpp:44-52 | ✅ 已完成（AR005/T002-T004：.vti 改 appended raw binary（vtk 36-49×/体积 -58%）、.bin 块写（3.4-3.8×，字节布局不变）；U2/U3 位级守护 + pvti 零变更） |
| P6 | **拓扑重复创建**：`UniformPartition::partition()` 内建 cart comm 用后即 free（partition.cpp:36-38），main.cpp:168 再建一次；两次均 `reorder=1`，理论上可给出**不同的 rank→coord 映射**（当前靠实现巧合保持一致） | partition.cpp:38 + main.cpp:168 | ✅ 已完成（AR004/T006：SubgridInfo.cartComm 所有权移交 + cart rank 查 coords（修 reorder 错位）+ main 删除二次创建；B3 拓扑一致性测试） |
| P7 | **halo 打包用 memcpy 中间缓冲**而非 MPI 派生数据类型直传（pack/unpack 与面缓冲均可省） | p2p_exchanger.cpp:139-242 | 未开始（对应 C1） |

### 3.4 架构层"深度天花板"

- **算法复杂度天花板**：Jacobi/GS 类迭代法对低频误差收敛慢（O(N²) 型迭代数），项目没有任何 O(N) 级算法（多重网格）。当前全部优化压"每步多快"，MG 改变"需要多少步"。
- **重叠深度 = 1**：流水线只藏一次通信；无时间分块（temporal blocking）。
- **单机叙事**：实测全部在 WSL 单机；跨节点网络行为（大消息聚合、RMA、真集合的收益曲线）无数据。
- **`--overlap-comm` 仅 Jacobi 支持**（main.cpp:198-199 打 WARN 忽略）——RBGS/CG 无重叠路径。

---

## 4. 优化空间分层详述

> 成本记号：S（<半天/单文件级）、M（1-2 天/跨模块）、L（>2 天/新子系统）。
> 每条含「验收判据」——新 AR 的 srs 验收标准建议直接引用。

### A. 诚实性修复层（最高优先）

#### A1 First-touch 真实现（修复 G1）

- **方案**：`zeroInitialize()` 三个场改 `#pragma omp parallel for` 按与计算循环相同的分区（j/k 外层）初始化；`applyDirichletBC` 同理评估。对齐 AGENT_SPEC.md:54 既有约定。
- **收益**：消除最大宣称缺口；多节点 NUMA 下避免首迭代页迁移（WSL 单机不可见，多节点叙事必备）。
- **成本/风险**：S / 极低（初始化语义不变）。
- **验收判据**：①初始化循环与 stencil 循环分区一致（代码走查）；②全量 ctest 绿（含 overlap 位级一致）；③README/AGENT_SPEC 宣称与实现对齐，无新的失实。

#### A2 真实残差语义 + 检查频率（修复 G3、P4）

- **方案**：收敛判据改真实残差 ‖Au−f‖₂（可在 stencil 循环内顺手累积 `Au-f`，单遍无额外扫），新增 `--residual-check-interval k`（默认 1）：每 k 步做一次归约判定，k 步间用递推残差外推或直接跳过判定。帮助文本明示语义。
- **收益**：`--tol` 语义可信（制造解报告收敛时 final_residual ≤ tol）；归约次数 1/iter → 1/k，通信占比可测下降；"归约开销管理"示范点落地。
- **成本/风险**：M / 中（k>1 时可能多迭代 ≤k 步才停；需文档注明）。
- **验收判据**：①制造解（sin）在 `--tol 1e-6` 下报告的 final_residual ≤ 1e-6（真实残差口径）；②`--residual-check-interval 10` 时 profiler 的 residual_allreduce 调用数/迭代 = 0.1；③PERFORMANCE 补充 k=1/5/10 的耗时对比；④既有测试全部绿。

#### A3 拓扑一致性（修复 P6）

- **方案**：`UniformPartition::partition()` 直接返回/复用 main 的 cart comm（或统一 `reorder=0`），消灭双重 `Cart_create` 的映射隐患。
- **成本/风险**：S / 低（接口小改，注意 AGENT_SPEC「接口变更须审批」——属于内部装配路径，建议 design 内记录）。
- **验收判据**：进程拓扑仅创建一次；新增一个 np=4 下 offsets/邻居/coords 一致性测试。

#### A4 观测自愈（修复 G5、G6）

- **方案**：RBGS 剖面区名独立（`rbgs_red_sweep`/`rbgs_black_sweep` 或统一 `solver_iteration`）；Profiler 改 per-thread 计数器（thread_local 或 per-thread 数组）+ report 时聚合，移除 mutex 栈或仅保护活跃区进出；hpp:13 注释与实现对齐。
- **成本/风险**：S-M / 低。
- **验收判据**：①RBGS 报告区名正确；②Profiler 在 OpenMP 并行区调用无数据竞争（TSan 走查或测试）；③既有 profile 报告字段兼容（不破坏 performance_report.json 消费脚本）。

### B. 算法层（示范价值最高的深度空间）

#### B1 几何多重网格（旗舰，拆两步）

- **现状**：项目已有 RBGS（现成 smoother）、均匀分解（现成 restriction/prolongation 权重结构）、CG（现成粗网格解法器）——离 MG 只差粗网格生成与限制/延拓算子。
- **B1a 两层校正格式（验证正确性）**：细网格 RBGS 预平滑 → 残差限制到粗网格 → 粗网格 CG 精解 → 延拓校正 → 后平滑。收敛理论（2 层即消除全部可分辨高频）可直接对比实测迭代数。
- **B1b 完整 V-cycle + MG-CG 预条件**：多层粗化（到 ≤8³ 粗根），粗层间通信用现有 exchanger 泛化；把 V-cycle 作为 CG 预条件子（对标 DESIGN §6 路线图"CG 预条件子"）。
- **收益**：算法叙事质变（O(N) vs O(N²) 迭代数，泊松 HPC 科研第一课）；256² 场景迭代数预计从数千降至数十 cycle。
- **成本/风险**：B1a M-L / 中（新通信模式：粗层 gather/broadcast 与残差归约）；B1b L / 中。
- **验收判据**：①B1a：同 tol 下总迭代数（平滑步计）与 2 层理论定性一致，制造解收敛；②B1b：256²/512² 下 MG-CG 迭代数 vs 裸 CG 对比入 PERFORMANCE；③RBGS smoother 复用不改其位级一致测试；④新增粗层通信的 np=1/4 正确性测试。

#### B2 CG 修补（最紧急性能修复 + 深化）

- **B2a 向量循环并行化（修复 P1）**：p/r/u/ap 更新全部 `omp parallel for + omp simd`（结构同 `axpyInterior`，cg_solver.cpp:85-105 已有正确范式可复制）。
- **B2b Pipelined CG（P2）**：Gropp 管线化（每迭代 1 次全局同步），`MPI_Iallreduce` 与向量运算重叠。与"每步 Allreduce 的 Jacobi"形成同步成本对照实验。
- **收益**：B2a 立竿见影（CG 每 iter 的串行部分占比高，预期 256² np=1 迭代时间下降，待实测入 PERFORMANCE）；B2b 在 np>1 时同步开销 3→1。
- **成本/风险**：B2a S / 极低；B2b M / 中（管线 CG 数值稳定性略降，需文档注明）。
- **验收判据**：①B2a：CG 256² np=1/np=4 iter_time 中位数对比入 PERFORMANCE（无回退门槛：np=1 不劣于现状）；②B2b：3 归约/iter → 1 归约/iter（profiler 计数断言），制造解收敛迭代数与裸 CG 一致（±5%）。

#### B3 Chebyshev 半迭代

- **方案**：Jacobi + Chebyshev 加速，系数由谱半径解析给定（λ_max/λ_min 由网格尺寸估计），**零全局同步**。
- **收益**：与"每步 Allreduce"的收敛率相近算法对比同步成本——SCALING_REPORT 的通信占比数据即刻有解释力；教学上覆盖"迭代加速不必同步"这一常被忽略的设计维度。
- **成本/风险**：M / 低-中。
- **验收判据**：制造解收敛；与 Jacobi 同 tol 迭代数/耗时/归约次数三方对比入 PERFORMANCE。

#### B4 RBGS 通信减半与重叠（修复 P3）

- **方案**：红黑两色边界合并为一次交换（pack 红+黑一起发，sweep 通信 2→1）；内层 stride-2 循环评估 `omp simd`（或改写为连续索引 + 位掩码选色以利向量化）；接入 `--overlap-comm` 流水线。
- **成本/风险**：M / 中（位级一致测试需同步更新）。
- **验收判据**：halo_exchange 计数/sweep 2→1（profiler 断言）；overlap on/off 位级一致测试延续；RBGS vs Jacobi 收敛比维持 ≈0.5。

### C. 通信层（扩展故事的三范式叙事）

#### C1 MPI 派生数据类型直传 halo（修复 P7）

- **方案**：用 `MPI_Type_create_subarray` 描述面数据（x/y/z 面各一型），`MPI_Send_init/Recv_init` 直接绑定计算缓冲，消掉 pack/unpack 与 `std::vector` 面缓冲。AR003 已积累 subarray 经验（含 MPI_ORDER_C 维度序陷阱，见 AR003 归档 evidence）。
- **收益**：消 memcpy 与中间缓冲；"pack vs 派生数据类型谁快"的经典问题获得本项目实测答案（大面大概率 datatype 赢，小面相反——数据本身就是 PERFORMANCE 好章节）。
- **成本/风险**：M / 中（halo 正确性测试面已厚，风险可控）。
- **验收判据**：①既有全部 halo 正确性测试（含 np4 非均匀、np8 3D）不变绿转；②pack/unpack 路径删除或经 flag 切换保留 A/B 能力；③256² np=1/4 halo_exchange 耗时对比入 PERFORMANCE。

#### C2 真集合 halo（修复 G2）

- **方案**：`MPI_Dist_graph_create_adjacent` 建邻居图 + `MPI_Neighbor_alltoallw` 收发面。与 C1 组合后形成"P2P vs 真集合"完整对比；WARN 与委托回退路径保留为 `--comm-mode p2p`。
- **成本/风险**：M / 中。
- **验收判据**：①collective 模式 np=1/4/8 halo 正确性测试全绿（不再 WARN 委托）；②同规模 P2P vs collective 的 halo_exchange 耗时/等待时间对比入 PERFORMANCE；③README:264 路线图条目移除。

#### C3 RMA halo（路线图条目）

- **方案**：`MPI_Win_create`（计算缓冲）+ epoch 模型（post/start/complete/wait 主动目标同步最贴合 stencil），每邻居 `MPI_Put` 面数据。
- **收益**：与 C1/C2 构成**同一 halo、三种 MPI 范式**的三方对比——本项目定位下含金量最高的通信叙事（近年 HPC 会议常客：P2P vs RMA vs NBC 在 stencil 上的行为差异）。
- **成本/风险**：M-L / 中（epoch 语义易错，正确性测试面已厚）。
- **验收判据**：三范式 halo 正确性等价（同一测试矩阵跑三遍）；np=1/4（及多节点若可）耗时对比入 PERFORMANCE + SCALING_REPORT。

#### C4 重叠深度 2（时间分块/wavefront）

- **方案**：halo=2 双缓冲，边界两层按依赖序流水计算，通信藏两层深。依赖 C1（面数据类型让多层边界更好写）。
- **成本/风险**：L / 高（位级一致约束可能需要重新设计累积序）。
- **定位**：远期；先完成 C1-C3 与 B1。

### D. I/O 层

#### D1 VTK 二进制格式 + binary 块写（修复 P5）✅ 已完成（AR005/T002-T004）

- **方案**：`.vti` 从 ASCII 改 appended binary（raw offset 模式，XML 头写 `format="appended"` + offset）；`BinaryIOBackend` 改块写（攒一行 memcpy 进缓冲或 `ofs.write(ptr, bytes)` 整行）。
- **收益**：512² VTK 输出体积/耗时预计降 1-2 个数量级（数据入 PERFORMANCE）；ParaView 加载时间实测。
- **成本/风险**：S-M / 低（错误处理路径沿用 mpibin 的 WARN-限流模式）。
- **验收判据**：①ParaView 正确加载 np=4 拼装结果（手工验证记录 evidence）；②256²/512² 输出耗时与文件体积 ASCII vs binary 对比入 PERFORMANCE；③既有 VTK 测试更新后全绿。

#### D2 VTI 单文件 MPI-IO（路线图条目）

- **方案**：rank0 写固定长度 XML 头（预留空间）或两遍写（数据后补头），数据区复用 AR003 mpibin 的双 subarray 直写模式。I/O 层随后形成「分片 / 单文件自定义 / 单文件标准格式」完整三件套。
- **成本/风险**：M / 中。
- **验收判据**：np=4 单 `.vti` 全局拼装 + ParaView 可读；与分片 `.pvti` 输出耗时对比入 PERFORMANCE。

#### D3 异步输出双缓冲

- **定位**：**仅设计注记，不建议近期实现**——输出线程引入后 `MPI_THREAD_FUNNELED` 级别与线程安全叙事需重审，复杂度/收益比不佳。

### E. 可观测性层

#### E1 PAPI 真接入（路线图条目）

- **方案**：可选依赖（CMake `find_package(PAPI)`，缺席优雅降级）：FLOPS + 内存带宽计数器入 performance_report；回答 PERFORMANCE §1 roofline "目标 60-80% 带宽"的无实测背书问题。
- **验收判据**：①PAPI 缺席时行为不变（CI 无 PAPI 路径全绿）；②有 PAPI 环境下 stencil 区段实测带宽 vs 理论值入 PERFORMANCE；③README:265 路线图条目移除。

#### E2 残差历史输出

- **方案**：`--residual-history <file>` 每步（或每 k 步）残差写 CSV；配 plot 脚本画收敛曲线。
- **收益**：成本极小的"放大器"——B 系列任何算法改动的效果即刻可视化；也是 E1/基准叙事的基础设施。
- **验收判据**：Jacobi/RBGS/CG 三求解器均可输出；plot 脚本出 PNG；AR002 §9 的对比图可复现。

#### E3 性能回归门禁

- **方案**：`run_scaling_tests.py` 增加与基线快照（JSON，入库）比对、超阈值（如 >10%）CI 报警；只看 ≥3 次中位数。
- **收益**：守护 AR002（-37.5%）与后续 AR 成果不被滚动开发悄悄回退。
- **验收判据**：CI job 演示一次"人为回退被抓住"的证据（evidence 记录）。

### F. 工程配套层

#### F1 CI 矩阵 + CI 可信度修复（**含一项紧急排查**）——探针+runner 固定已完成（AR005/T001/T005），矩阵扩展未开始

- **紧急项**：CI 的 Debug(MPICH) job 跑在 ubuntu-latest（24.04）上——**该平台 mpich 4.2.0 存在 PMI/PMIx 不匹配缺陷，应用静默退化为单进程**（AR003 开发期实测确认，README:28 已记录）。这意味着 **CI 里 MPICH 路径的 np4/np8 多 rank 测试（halo_2d_mpi、solver_mpi 等）很可能一直在"4 个单进程假通过"**。修复：ctest 增加进程数探针测试（`MPI_Comm_size != 预期 np` 即 FAIL），并 pin MPICH 到已修复版本或换发行版。
- **矩阵**：GCC/Clang × OpenMPI/MPICH、ccache 缓存、覆盖率 job（见 F2）。
- **验收判据**：①size 探针入 ctest 且在 np=2/4/8 测试中激活；②CI 矩阵全绿；③缓存使 CI 时长不显著回退。

#### F2 覆盖率基线

- **方案**：gcovr/lcov 接入，CI 出基线数字；不设硬门槛（先有数据再谈阈值，避免指标驱动开发）。
- **验收判据**：CI 产出覆盖率报告；数字记入本文档 §3.1。

#### F3 CLI 类型校验统一

- **方案**：CommandLineParser 为数值型选项显式校验并输出可读错误（当前 `--save-interval abc` 依赖 `std::stoi` 抛异常，ctest parse_error 测的正是这个副作用）。
- **验收判据**：全部数值选项非法值 → 统一错误消息 + 退出码；parse_error 测试更新为显式语义。

---

## 5. 优先级路线图（AR 切分建议）

排序原则：**诚实性 → 快赢 → 深度叙事**；依赖先行；每个 AR 保持项目惯例（小而完整、带基准与文档）。

| 优先 | 建议 AR | 内容（对应条目） | 规模 | 核心价值 |
|------|---------|----------------|------|---------|
| ★1 | AR004 | 诚实性修复包：A1 first-touch + A2 真实残差/检查频率 + A3 拓扑 + A4 观测自愈 + B2a CG 并行化 | M | 消灭全部宣称缺口；CG 性能修复立竿见影 |
| ★2 | AR005 | I/O 快赢 + CI 可信度：D1 VTK 二进制/块写 + F1（**含 CI mpich 假通过紧急排查**） | S-M | 输出提速 1-2 量级；CI 数据可信 |
| ★3 | AR006 | 通信范式三部曲 I：C1 派生数据类型直传 + C2 真集合 halo | M | 收尾半成品示范；pack vs datatype 实测 |
| ★4 | AR007 | 算法深水区第一步：B1a 两层 MG 校正 | M-L | O(N) 算法叙事开局 |
| 后续 | AR008+ | B1b V-cycle/MG-CG、B2b pipelined CG、B3 Chebyshev、B4 RBGS 通信、C3 RMA、D2 VTI 单文件、E1 PAPI、E2 残差历史、E3 性能门禁、F2/F3；工程卫生：test_alt_solvers.cpp 拆分（AR005 review 遗留，IO 用例拆至 test_io_layout.cpp） | — | 按依赖与资源排入 |

依赖关系（简化）：

```
A1/A2/A3/A4 ──(无依赖，可并行)──────────────┐
B2a ──→ B2b（pipelined 建立在并行化之上）      │
C1 ──→ C2 ──→ C3（三范式共用面描述基建）       ├──→ E3（门禁守护一切成果）
B1a ──→ B1b（B1b 复用 B1a 的限制/延拓）       │
E2（放大器，随任意算法 AR 顺带落地）           │
D1 ──→ D2（单文件 VTI 复用 AR003 模式）      ┘
```

---

## 6. 决策规则（新优化项的进出准则）

1. **进**：新条目须能用 §2 四维标尺说清价值，且给出可测量的验收判据；无验收判据的"感觉优化"不进本指南。
2. **出**：完成（附 AR 编号与证据）、或经论证定位不符（记录理由）。
3. **通用验收底线**（所有优化 AR 继承）：
   - 全量 ctest Debug+Release 绿，ASan/UBSan 零报告；
   - overlap on/off 位级一致约束不被破坏（或显式走 SDD 变更流程修订该约束）；
   - AGENT_SPEC 约定不破坏；接口变更须在 design.md 记录授权；
   - 性能敏感改动必须附 PERFORMANCE 数据（≥3 次中位数口径），无失实宣称；
   - 路线图条目实现后从 README「未实现」清单移除（防 G 表再生）。
4. **复杂度预算**：单 AR 净新增代码建议 ≤ ~800 行（AR003 约为基准）；超出须拆分或论证。

---

## 附录 A：测量与环境注意事项

- **环境**：WSL2 单机（Ultra 9 185H，OpenMPI，g++）；构建目录建议放 WSL ext4（`/root/build-*`）保证 I/O 剖面真实；
- **口径**：性能数据一律 ≥3 次取中位数（延续 AR002/AR003 惯例）；iter_time 取 `--enable-profiling` 剖面区段；
- **已知环境陷阱**：Ubuntu 24.04 mpich 4.2.0 PMI/PMIx 不匹配 → 多 rank 静默退化单进程（AR003 实测；验证法：`mpirun -n 4 <app>` 打印 `MPI_Comm_size` 是否为 4）。**任何多 rank 测试"通过"前先确认不是单进程假通过**；
- **行号漂移声明**：本文档 file:line 证据基于基线 commit `4cd7024`，代码演进后以符号名/语义检索为准。

---

*Last updated: 2026-10-08（基线 4cd7024）*
