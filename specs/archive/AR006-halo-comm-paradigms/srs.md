# [AR006] 需求设计说明书

| 字段 | 内容 |
|------|------|
| AR 编号 | AR006 |
| AR 主题 | halo-comm-paradigms（通信范式三部曲 I：C1 派生数据类型直传 + C2 真集合 halo） |
| 关联 SR | 无（源自 docs/OPTIMIZATION_GUIDE.md §5 路线 ★3，条目 C1 + C2；修复 §3 P7 + G2） |
| 日期 | 2026-10-09 |
| 状态 | 已确认（用户授权自动化执行，依据 OPTIMIZATION_GUIDE.md §5 ★3） |

## 1. 背景与目标

OPTIMIZATION_GUIDE.md §3 诊断 P7：halo 打包用 memcpy 中间缓冲（`std::vector` 面缓冲 + pack/unpack 三重循环），而非 MPI 派生数据类型直传；G2：`--comm-mode collective` 为纯 P2P 委托 + WARN，真集合通信不存在。本 AR 落地路线 ★3 通信范式三部曲的前两部：
- C1：`MPI_Type_create_subarray` 描述面数据，`MPI_Send_init/Recv_init`（或等价非阻塞收发）直接绑定计算缓冲，消掉 pack/unpack 与面缓冲——"pack vs 派生数据类型谁快"获得本项目实测答案
- C2：`MPI_Dist_graph_create_adjacent` 邻居图 + `MPI_Neighbor_alltoallw` 真集合收发面（复用 C1 的面 datatype），collective 模式不再委托

## 2. 需求范围

**In Scope（对应指南条目）：**
- C1：新增派生数据类型直传 halo 实现（x/y/z 三对面各一 subarray 型；消 pack/unpack 中间缓冲）；pack 路径与 datatype 路径可切换（A/B 对比能力保留，GUIDE C1 验收②；具体切换机制与默认值由 design 决定并记录）
- C2：`CollectiveExchanger` 真实现（Dist graph + Neighbor_alltoallw；WARN 委托路径退位为回退/对照，`--comm-mode p2p` 保留既有行为）
- A/B 与性能数据：pack vs datatype、p2p vs collective 的 halo_exchange 剖面对比（256²，np=1/4），入 PERFORMANCE
- 文档同步：README:264 路线图「真集合通信」条目移除、GUIDE P7/G2 勾选、AGENT_SPEC 类层级更新、SCALING_REPORT（如涉及宣称）

**Out of Scope：**
- C3 RMA halo（三部曲第三部，路线后续 AR）
- B 系列算法改动（通信范式不变更求解器语义）
- 多节点实测（WSL 单机环境限制，叙事口径沿既有 PERFORMANCE/SCALING_REPORT 惯例）

## 3. 功能需求

### 3.1 R1 派生数据类型直传 halo（C1，修复 P7）

**描述：** 新增 HaloExchanger 实现：以 `MPI_Type_create_subarray` 描述六方向面在 padded 缓冲中的位置与形状，非阻塞收发直接以计算缓冲为源/目的，消除 pack/unpack 与 `std::vector` 面缓冲。

**触发条件：** 选用 datatype 路径的 halo 交换（切换机制由 design 决定）。

**期望行为：**
- x 面（Left/Right）：subarray 描述 [k][j][hw 列]（x 维连续段）；y 面（Down/Up）：[k][hw 行][i]；z 面（Back/Front）：[hw 层][j][i]——三对面在 padded 缓冲（nxT×nyT×nzT，row-major）内均为单一 subarray 可描述
- 发送型 starts 指向 interior 边带，接收型 starts 指向 halo 带；两侧对称构造，wire 布局由 datatype 自述（y 面从现 pack 的 [h][k][i] 变为 [k][h][i]，自洽无兼容性负担——收发两侧同型）
- 2D（nzLocal==1）与 3D k 域（[kBegin,kEnd)）语义沿 pack 版 `faceLayerRange` 口径
- MPI_PROC_NULL 邻居方向跳过（不建请求/不通信）；沿既有错误处理约定（MPI 路径不抛异常，HYPOS_ERROR/WARN 日志）
- datatype 与持久请求在 initialize 建立、析构释放（含 MPI_Finalized 防御，沿 PointToPointExchanger 析构惯例）

**异常处理：** initialize 前/空 data 指针调用 → HYPOS_ERROR + 忽略（沿 p2p 版契约）。

**验收标准：**
- Given 既有 halo 正确性测试矩阵（2D 方向语义/多 halo/PROC_NULL/np4 非均匀、3D 方向语义/np8），When 以 datatype 路径运行，Then 全部通过（GUIDE C1 验收①「不变绿转」）
- Given 同一 Subgrid 与数据场，When pack 路径与 datatype 路径各自完成 exchange，Then halo 带内容一致（等价性测试）
- Given np=1（全 PROC_NULL 边界），When exchange，Then 无 MPI 通信调用且行为安全

### 3.2 R2 真集合 halo（C2，修复 G2）

**描述：** `CollectiveExchanger` 从 P2P 委托改为真集合实现：`MPI_Dist_graph_create_adjacent` 建邻居有向图，`MPI_Neighbor_alltoallw`（begin/end 语义用非阻塞变体）收发面数据。

**触发条件：** `--comm-mode collective`。

**期望行为：**
- 邻居图由六方向邻居关系构造（MPI_PROC_NULL 方向不入图）；图通信子在 initialize 建立、析构释放
- 邻居序（sources/destinations 顺序）与 alltoallw 的 sendbuf/displs/sendtypes 数组严格对应（design 明确约定，如按方向枚举序）
- sendtypes/recvtypes 复用 R1 的面 subarray 型（逐邻居数组）；数据直传计算缓冲（无 pack）
- np=1 无邻居：空图 + 零邻居 alltoallw 调用，行为安全
- 不再输出「falls back to PointToPointExchanger」WARN

**异常处理：** MPI 路径不抛异常；初始化失败 HYPOS_ERROR。

**验收标准：**
- Given collective 模式 np=1/4/8（含非均匀切分、3D），When 跑 halo 正确性测试，Then 全绿且无委托 WARN（GUIDE C2 验收①）
- Given 既有 `CollectiveExchangerDelegatesToP2P` 测试，When 真集合落地，Then 该测试更新为真集合语义断言（委托断言退位，不得静默删除覆盖）
- Given `--comm-mode p2p`，When 运行，Then 既有 P2P 行为零变化（回退路径保留）

### 3.3 R3 A/B 切换与性能对比

**描述：** pack 与 datatype 路径可切换运行（同规模同负载 A/B）；性能数据如实入 PERFORMANCE。

**期望行为：**
- 切换机制不引入求解器语义变化；两路径可被同一测试矩阵覆盖
- PERFORMANCE 新小节：256² halo_exchange 剖面耗时，np=1/4，pack vs datatype、p2p vs collective（≥3 次中位数，口径说明）
- 不预设倍数门槛（诚实优先，沿 AR004/AR005 D6 惯例）：设「datatype 不劣于 pack 超可解释幅度」「collective 不劣于 p2p 超可解释幅度」的记录式判定，实测结论如实记录（含"小面 datatype 可能更慢"的预期）

**验收标准：**
- Given 切换机制，When 同一 ctest 矩阵分别以两路径运行，Then 均全绿
- Given PERFORMANCE，Then 含 pack vs datatype 与 p2p vs collective 两组对比数据及口径说明

### 3.4 R4 文档同步

**描述：** 文档与实现零缺口（延续 AR002-AR005 惯例）。

**期望行为：** README:264 路线图「真集合通信」条目移除（连同 `--comm-mode collective` 委托注记）；GUIDE 勾选 P7/G2（附 AR 编号与证据）；AGENT_SPEC 类层级/通信模块描述更新；halo_exchanger.hpp 中 CollectiveExchanger 的过时 doc 注释（「delegates…roadmap item」）同步更新；README `--comm-mode` 帮助文本与新行为一致。

**验收标准：**
- Given 文档走查，Then 无失实宣称；委托相关的旧表述全部更新
- Given `--help` 输出，Then comm-mode 说明与实际行为一致

## 4. 非功能需求

| 类型 | 指标 | 要求 |
|------|------|------|
| 正确性 | halo 等价性 | 同一场两路径 halo 带内容一致；collective 与 p2p 结果一致 |
| 兼容性 | 既有行为 | `--comm-mode p2p` 默认路径零变化；HaloExchanger 接口不变（新增实现不改基类签名） |
| 内存安全 | ASan/UBSan | 全部测试零报告 |
| 资源 | MPI 对象泄漏 | datatype/comm/request 全部释放（析构防御 MPI_Finalized） |
| 构建 | 双模式 | 0 新增警告 |
| 性能 | 数据如实 | 对比数据不预设结论，含反直觉结果也如实入档 |
| 回归 | 全量测试 | Debug+Release ctest 全绿；既有 halo/solver/overlap 位级守护测试不变绿转 |
| 复杂度 | 预算 | 单 AR 净新增代码 ≤ ~800 行（GUIDE §6） |

## 5. 约束与假设

**约束：**
- 延续 AGENT_SPEC 全部约定（row-major、MPI 路径不抛异常、snake_case）
- 不引入外部依赖；C++17 + MPI（MPI-3 起的 Dist graph/Neighbor_alltoallw，与既有 MPI-IO 假设一致）
- exchange 接口的 data 指针由调用方传入、可能逐次不同——持久请求不得绑定首次调用指针（design 决定请求策略并记录）
- 测试环境：WSL 单机 OpenMPI；np>1 条目受 AR005 探针守卫（HYPOS_EXPECT_NP）

**假设：**
- OpenMPI 支持 MPI_Dist_graph_create_adjacent（MPI-3）与 MPI_Ineighbor_alltoallw（MPI-4.0 新增，OpenMPI ≥4.0 系；design 第 3 轮 P2 勘误）
- AR003 的 subarray 经验（MPI_ORDER_C 维度序 [z][y][x]）适用——面描述以「慢维到快维」顺序给出
- halo 测试矩阵（既有 8 用例 + np4 非均匀 + np8 3D）对 wire 布局变化不敏感（它们断言的是 halo 带内容语义）

## 6. 术语说明

| 术语 | 定义 |
|------|------|
| 派生数据类型直传 | 以 MPI datatype 描述非连续数据在内存中的形状，通信直接以原缓冲为源/目的，无显式打包 |
| subarray 型 | `MPI_Type_create_subarray` 描述的 N 维规则子块（sizes/subsizes/starts，MPI_ORDER_C 慢维在前） |
| 邻居集合通信 | `MPI_Neighbor_alltoall*` 族：仅与_dist_graph 声明的邻居交换，无全局同步语义 |
| pack/unpack 路径 | 现状 P2P 实现：memcpy 面数据进中间 `std::vector` 缓冲再以 MPI_DOUBLE 连续块收发 |
| wire 布局 | 面数据在传输字节流中的元素顺序；pack 版 y 面为 [h][k][i]，datatype 版为 [k][h][i]（自洽即可） |
