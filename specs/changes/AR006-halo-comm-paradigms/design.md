# [AR006] 详细设计文档

| 字段 | 内容 |
|------|------|
| AR 编号 | AR006 |
| AR 主题 | halo-comm-paradigms（C1 派生数据类型直传 + C2 真集合 halo） |
| 关联 srs.md | ./srs.md |
| 日期 | 2026-10-09 |
| 状态 | 已确认（req 门控 PASS 第 1 轮） |

## 1. 设计目标（从 srs 追溯）

- R1：halo 交换消 pack/unpack 中间缓冲（subarray 直传）
- R2：collective 模式真集合实现（Dist graph + Neighbor_alltoallw）
- R3：三实现 A/B 可切换 + 性能数据如实入档
- R4：文档零缺口（README 路线图/GUIDE 勾选/AGENT_SPEC/hpp 注释）

## 2. 功能点分解

| # | 功能点 | 源需求 | 说明 |
|---|--------|--------|------|
| FP1 | 面 subarray datatype 家族（六方向 × 发/收 12 型） | R1 | 单一 `MPI_Type_create_subarray` 描述每面（D3 参数表） |
| FP2 | `DatatypeExchanger`（p2p 范式 datatype 打包） | R1/R3 | Isend/Irecv 直传计算缓冲；非持久请求（D2） |
| FP3 | `CollectiveExchanger` 真实现 | R2/R3 | Dist graph + `MPI_(I)Neighbor_alltoallw`；复用 FP1 类型 |
| FP4 | CLI `--comm-mode` 扩展 `datatype` + 接线 | R3/R4 | 三值 p2p\|datatype\|collective，默认 p2p 不变（D1） |
| FP5 | 测试矩阵（U/B/E 系列 + ctest 条目） | R1/R2/R3 | 语义等价 + 非均匀 np4/np8 + begin/end + np=1 安全 |
| FP6 | bench + 文档同步 | R3/R4 | halo_exchange/halo_wait 剖面对比；README :267 移除等 |

## 3. 影响范围（文件清单）

| 文件 | 变更 |
|------|------|
| src/comm/halo_exchanger.hpp | 新增 `DatatypeExchanger` 类声明；更新 `CollectiveExchanger` 声明（去 delegate_）与 doc 注释（R4 项） |
| src/comm/datatype_exchanger.cpp | 新增（FP1/FP2） |
| src/comm/collective_exchanger.cpp | 重写（FP3） |
| src/main.cpp | comm-mode 接线扩展 `datatype`（:225-230）；`--help` 文本 |
| CMakeLists.txt | 新源文件入 target + 新 ctest 条目（FP5） |
| tests/test_halo_exchange.cpp | 新增 U/B/E 用例；`CollectiveExchangerDelegatesToP2P` 改写 |
| scripts/bench_ar006.sh | 新增（FP6） |
| docs/PERFORMANCE.md、docs/OPTIMIZATION_GUIDE.md、README.md、AGENT_SPEC.md | 文档同步（FP6） |

不改：`PointToPointExchanger`（pack 版原样保留 = A/B 基线）、solvers、`HaloExchanger` 基类签名。

## 4. 实现设计

### 4.1 FP1 面 subarray datatype（D3 参数表——本设计核心）

padded 缓冲 sizes（`MPI_Type_create_subarray` 3 维，`MPI_ORDER_C` 慢维在前，AR003 经验）：`(nzT, nyT, nxT)`，oldtype=`MPI_DOUBLE`，count=1。

k 域：`k0 = (nzLocal==1) ? 0 : kBegin`；`kCount = (nzLocal==1) ? 1 : (kEnd-kBegin)`（沿 pack 版 `faceLayerRange` 口径，p2p_exchanger.cpp:66-74）。

每方向两个型（send 指 interior 边带，recv 指 halo 带）：

| 方向 | subsizes | send starts | recv starts |
|------|----------|-------------|-------------|
| Left | (kCount, nyL, hw) | (k0, jBegin, **iBegin**) | (k0, jBegin, **0**) |
| Right | (kCount, nyL, hw) | (k0, jBegin, **iEnd-hw**) | (k0, jBegin, **iEnd**) |
| Down | (kCount, hw, niL) | (k0, **jBegin**, iBegin) | (k0, **0**, iBegin) |
| Up | (kCount, hw, niL) | (k0, **jEnd-hw**, iBegin) | (k0, **jEnd**, iBegin) |
| Back | (hw, nyL, niL) | (**kBegin**, jBegin, iBegin) | (**0**, jBegin, iBegin) |
| Front | (hw, nyL, niL) | (**kEnd-hw**, jBegin, iBegin) | (**kEnd**, jBegin, iBegin) |

（`niL = iEnd-iBegin`；z 面仅 3D 且邻居非 PROC_NULL 时构造；全局 nz==1 时 Back/Front 邻居恒 PROC_NULL，自然跳过；degenerate 3D 切分（nzLocal==1 但存在真实 z 邻居）下以 neighbor!=PROC_NULL 为门构造，与 pack 版逐面语义等价不受影响（design 门控第 1 轮 G1 观察项注记）。）

wire 布局自洽性：三对面均由单一 subarray 自述（y 面 wire 序 [k][h][i]，与 pack 版 [h][k][i] 不同——收发两侧对称同型，无兼容性负担；srs §6 术语表已声明）。tag 沿 pack 版方向序（tag = direction）。

正确性关键点：subsizes/starts 与 packFace/unpackFace（p2p_exchanger.cpp:129-243）逐面语义对照——x 面读 [k][j][srcCol..srcCol+hw) 列段、写 [k][j][dstCol..) 列段等，设计评审时逐行核对。

### 4.2 FP2 `DatatypeExchanger`（伪码）

```
initialize(subgrid):
  comm_ = subgrid.comm()
  for direction in 0..5:
    if neighbor(direction) == MPI_PROC_NULL: continue        # 不建型不通信
    activeDirs_.push(direction)
    sendTypes_[d] = makeFaceType(subgrid, d, SEND);  MPI_Type_commit
    recvTypes_[d] = makeFaceType(subgrid, d, RECV);  MPI_Type_commit
  initialized_ = true

beginExchange(subgrid, data):     # 非持久请求（D2）：data 指针逐次可变
  for d in activeDirs_: MPI_Irecv(data, 1, recvTypes_[d], neighbor(d),
                                  opposite(d), comm_, &reqs[2i])
  for d in activeDirs_: MPI_Isend(data, 1, sendTypes_[d], neighbor(d),
                                  d, comm_, &reqs[2i+1])

endExchange(subgrid, data):
  MPI_Waitall(2n, reqs, MPI_STATUSES_IGNORE)     # 无 unpack：数据已直落 halo 带

~DatatypeExchanger():  # 防御 MPI_Finalized（沿 p2p 析构惯例 :82-95）
  for type in committed types: MPI_Type_free
  for req in pendingReqs_ (非 MPI_REQUEST_NULL): MPI_Request_free   # begin 未 end 防御（沿 :88-92 惯例）
```

`exchange` = begin + end。错误处理：未 initialize/空 data → `HYPOS_ERROR` + return（沿 p2p 契约 p2p_exchanger.cpp:251-255）。

### 4.3 FP3 `CollectiveExchanger` 真实现

```
initialize(subgrid):
  comm_ = subgrid.comm()
  activeDirs_ = [d for d in 0..5 if neighbor(d) != MPI_PROC_NULL]
  # sources 对称构造（与 destinations 同序同值）：入邻居多重集恒等（推导见下）；不可整体取 opposite 序（第 2 轮复审 Important）
  destinations = [neighbor(d) for d in activeDirs_]
  sources      = [neighbor(d) for d in activeDirs_]
  MPI_Dist_graph_create_adjacent(comm_, n, sources, MPI_UNWEIGHTED,
                                 n, destinations, MPI_UNWEIGHTED,
                                 MPI_INFO_NULL, 0 /*reorder*/, &graphComm_)
  for i, d_i in enumerate(activeDirs_):
    sendTypes_[i] = send 型 of d_i + commit
    f_i = (neighbor(d_i) == self_rank) ? opposite(d_i) : d_i   # 分情形落位（D4）
    recvTypes_[i] = recv 型 of f_i + commit;  recvFace_[i] = f_i
  # alltoallw 参数数组按 activeDirs_ 枚举序
  sendcounts_[i] = recvcounts_[i] = 1
  sdispls_[i] = faceSendFirst(d_i) 相对 data 的字节偏移
  rdispls_[i] = faceRecvFirst(recvFace_[i]) 同理    # 每次 exchange 时按当次 data 重算

exchange(subgrid, data):
  MPI_Neighbor_alltoallw(data, sendcounts_, sdispls_, sendTypes_,
                         data, recvcounts_, rdispls_, recvTypes_, graphComm_)

beginExchange: MPI_Ineighbor_alltoallw(同参数, &req_)     # MPI-3 非阻塞变体
endExchange:    MPI_Wait(&req_)

~CollectiveExchanger(): 防御 finalized；在途 req_ 非 NULL 则 MPI_Request_free；
                       MPI_Type_free ×2n；MPI_Comm_free(&graphComm_)
```

np=1：`activeDirs_` 空 → `MPI_Dist_graph_create_adjacent` 建空图（0 进 0 出，合法）→ alltoallw 无操作。不再输出委托 WARN（collective_exchanger.cpp:12-13 删除）。

**邻居配对与落位定律（design 门控第 2 轮 G3 修复，D4 核心依据）**：
- 数据流恒等式（pack 版语义）：`interior_A(d) → halo_{n_A(d)}(opposite(d))`（p2p 由 tag=opposite(d) 保证，p2p_exchanger.cpp:119）；等价地 `halo_A(e) ← interior_{n_A(e)}(opposite(e))`
- **sources 对称构造合法**：X 是我的入邻居 ⟺ 存在 e 使 n_X(e)=我 ⟺ opposite(e)∈activeDirs_ 且 n(opposite(e))=X——入邻居多重集 = {n(d) : d∈activeDirs_} 严格成立。**非周期笛卡尔边界 rank（activeDirs 不对合自反封闭，如 2×2×2 角点 active={Right,Up,Front}）同样成立**：封闭性来自「共享面被两端各自声明」，而非本 rank 方向集自反；对称序的 sources=[n(Right),n(Up),n(Front)] 恰为全部真实入邻居
- **recv 块 i 分情形落位**（本地判定，不依赖对端枚举序——唯一边按邻居身份配对，天然免疫两侧 activeDirs 差异）：
  - 唯一边（n(d_i)≠self）：对端 X=n(d_i) 面向我的一面是 X 的 opposite(d_i)，X 发来 interior_X(opposite(d_i))，落 **halo(d_i)**
  - 平行边（n(d_i)==self；非周期笛卡尔下平行边仅自环，即该维 np=1）：MPI 按出现序配对——我 destinations 中 self 第 a 次出现（send 块=interior(d_i)）↔ sources 中 self 第 a 次出现（recv 块 i）；由恒等式落 **halo(opposite(d_i))**
  - 统一式：`f_i = (n(d_i)==self) ? opposite(d_i) : d_i`，recvTypes_[i]/rdispls_[i] 取 f_i 面。自环时 n(d_i)=self ⟹ n(opposite(d_i))=self ⟹ opposite(d_i)∈activeDirs_，建型无缺失
- 三情形复算：①唯一边（含边界角点）：rank(0,0,0) recv 块 0 收 n(Right) 的 interior(Left) 落 halo(Right) ✓；②自环（全平行边）：send 块 i=interior(d_i) 按位回到 recv 块 i 落 halo(opposite(d_i))——镜像语义 ✓；③混合（2×1 切分）：rank(0,0) active={Right,Down,Up}（n(Down)=n(Up)=self）：recv 块 0 落 halo(Right)，recv 块 1（self 第 1 次出现↔interior(Down)）落 halo(Up)，recv 块 2 落 halo(Down) ✓
- 前置假设：非周期笛卡尔（本仓库 decomposition 现状）；周期 np=2 维会产生对同一真邻居的平行边（按出现序同法推导，本 AR 不涉及）
- 教训记录：第 1 轮「recv 侧整体 opposite 序」修法只修复自环、破坏笛卡尔边界——sources 退化为 [PROC_NULL×n]，与对端 destinations 声明不一致且 halo 永不填充（第 2 轮复审 Important，B1-B4/collective_mpi 必挂）

### 4.4 FP4 CLI 接线

main.cpp :225-230 扩展：`commMode == "datatype"` → `std::make_unique<DatatypeExchanger>()`。`--help`（:49）comm-mode 说明改 `p2p, datatype, collective (default: p2p)`，datatype 注记「p2p 范式 + MPI 派生数据类型直传（无打包缓冲）」。RunConfig 透出不变（commMode 字段已有）。

### 4.5 FP5 测试设计

单测（tests/test_halo_exchange.cpp，gtest，MPI_COMM_SELF 可单进程跑）：

| ID | 用例 | 断言要点 |
|----|------|---------|
| U1 | Datatype2DDirectionSemanticsSelfLoop | 沿 DirectionSemanticsSelfLoop2D（:27-58）模式：2D 自环（左右=自身），已知编码场 exchange 后 halo 带 == 对侧 interior 带（逐元素 `sg.at` 断言） |
| U2 | DatatypeMultiHaloSelfLoop | hw=2（沿 FaceBuffersSafeSelfLoop2DMultiHalo :61 模式）：多 halo 层全部正确 |
| U3 | Datatype3DDirectionSemanticsSelfLoop | 沿 SelfLoop3DDirectionSemantics（:255）模式：六方向（3D 自环） |
| U4 | DatatypeEquivalentToPack（**3D hw=2 配置**——同时覆盖 3D 多 halo，design 门控第 1 轮 G5-Minor 补强） | 同场两 exchanger 各自 exchange 后 halo 带逐元素一致（wire 布局不同但落位相同——等价性铁证） |
| B1 | MpiNonUniform4RanksHalosDatatype | 沿 :151 模式 np4 非均匀切分，DatatypeExchanger |
| B2 | MpiNonUniform8RanksHalos3DDatatype | 沿 :371 模式 np8 3D |
| B3 | MpiNonUniform4RanksHalosCollective | np4 非均匀，CollectiveExchanger 真集合；沿 :151/:371 模式（含 `if (size!=np) GTEST_SKIP` 守卫——unit 全量跑时不误执行，第 2 轮复审 Minor 补注） |
| B4 | MpiNonUniform8RanksHalos3DCollective | np8 3D 真集合；同 :371 模式（含 np 守卫） |
| B5 | CollectiveExchangerIsTrueCollective | 改写 :129 委托断言：name()=="collective"、真集合路径下 halo 带与 p2p 等价（自环场景）——委托 WARN 退位不以日志断言（脆弱），以行为等价断言 |
| E1 | NoNeighborNp1Safe（参数化三 exchanger） | np=1 全 PROC_NULL：exchange/begin+end 安全，halo 不变 |
| E2 | BeginEndSplitEquivalent | beginExchange 后 endExchange 前 halo 未就绪不读（仅 end 后断言）；三 exchanger begin/end 与 exchange 结果一致 |

ctest 条目（CMakeLists，np>1 带 `HYPOS_EXPECT_NP`，gtest filter 追加 `:MpiEnvTest.*`；**专用点名 filter**——沿既有 halo_2d_mpi（:153 `HaloExchangeTest.MpiNonUniform4RanksHalos:MpiEnvTest.*`）与 halo_3d_mpi 惯例，避免 np 硬编码用例在错误规模下执行）：
- `datatype_np4`（np=4，filter `HaloExchangeTest.MpiNonUniform4RanksHalosDatatype:MpiEnvTest.*`）
- `datatype_np8`（np=8，filter `HaloExchangeTest.MpiNonUniform8RanksHalos3DDatatype:MpiEnvTest.*`）
- `collective_np4`（np=4，filter `HaloExchangeTest.MpiNonUniform4RanksHalosCollective:MpiEnvTest.*`）
- `collective_np8`（np=8，filter `HaloExchangeTest.MpiNonUniform8RanksHalos3DCollective:MpiEnvTest.*`）
- `datatype_e2e`（np=4 直跑 hypos `--comm-mode datatype`，沿 collective_mpi :182-185 模式）
- U1-U4/B5/E1/E2 为单进程（MPI_COMM_SELF/世界大小无关）用例，由 `unit` 条目（:87，无 filter 全量）自然覆盖
- 既有 `halo_2d_mpi`/`halo_3d_mpi`（p2p 路径）filter 专用不扩展；p2p 对照侧证据由 U4 等价性用例 in-test 承担（同场 pack vs datatype 对比）；`collective_mpi` e2e 保留

### 4.6 FP6 bench 与文档

`scripts/bench_ar006.sh`（沿 bench_ar005.sh 模式，mawk 兼容）：
- 负载：Jacobi 256² `--max-iter 200 --tol 0.0 --residual-check-interval 10000`（残差开销旁路）、np=1/4、OMP=1（np4 加 OMP=1 per-rank，与 CI 语义一致）、`--comm-mode p2p|datatype|collective` ×3 中位
- 取数：profiler `halo_exchange` + `halo_wait` 区段（jacobi_solver.cpp:118/131）
- 判定（记录式，D6）：median(datatype) vs median(p2p)、median(collective) vs median(p2p) 如实记录；np4 本机漂移 ±15% 惯例注记；「劣化超 2× 且可解释」触发默认值回退评审（D1 回退条款）
- 文档：PERFORMANCE §13；README :267 路线图真集合条目移除 + :49 帮助文本；GUIDE P7/G2 ✅；AGENT_SPEC 通信模块类层级（DatatypeExchanger）；halo_exchanger.hpp:86-88 注释重写（req 门控 G8 项）

## 5. 接口描述

无库内公共接口变更（`HaloExchanger` 基类签名不动）；新增子类 `DatatypeExchanger`（实现全部 4 纯虚方法 + name()=="datatype"）；`CollectiveExchanger` 公共签名不变、行为变更（委托→真集合）。CLI `--comm-mode` 新增值 `datatype`（默认 `p2p` 不变）。

## 6. 测试设计（§4.5 全表）与交付物映射

| srs 验收 | 设计落点 |
|----------|---------|
| R1① 不变绿转 | U1-U3/B1/B2 + 既有 halo_2d_mpi/halo_3d_mpi 复跑 |
| R1② 等价性 | U4 |
| R1③ np=1 安全 | E1 |
| R2① np=1/4/8 全绿无 WARN | B3/B4/E1 + collective_np4/np8 + 委托 WARN 代码删除 |
| R2② 委托测试更新 | B5 |
| R2③ p2p 零变化 | 不改 PointToPointExchanger + 全量回归 |
| R3 双路径全绿矩阵 | datatype_np4/np8/collective_np4/np8 + datatype_e2e |
| R3 PERFORMANCE 数据 | FP6 bench |
| R4 文档 | FP6 文档清单 |

覆盖率口径：沿 AR005 声明（无 gcovr 环境声明，走查 + 全量绿佐证）。

## 7. 决策记录

| # | 决策 | 备选与否决理由 |
|---|------|---------------|
| D1 | 切换机制：`--comm-mode` 扩展三值 `p2p\|datatype\|collective`，默认 `p2p` 不变 | 备选 a) 新 flag `--halo-pack`：两 flag 交叉语义（collective 恒 datatype）徒增解释成本，否决；备选 b) 默认改 datatype：违背「p2p 默认路径零变化」NFR 且 bench 未证，否决。回退条款：若 bench 实测 256² datatype 劣化超 2× 且可解释（小面 type setup 开销），记录数据但默认保持 p2p（本 AR 默认即 p2p，无需回退动作，未来 AR 依数据再议） |
| D2 | DatatypeExchanger 用非持久 Isend/Irecv（type 持久、请求每次新建） | exchange 的 data 指针逐次可变（CG 传独立 p 缓冲，cg_solver.cpp:44），Send_init 绑定首指针是悬垂级缺陷；持久请求收益（省请求创建）相对 type 直传的通信收益是次要项，正确性优先 |
| D3 | 三对面各用单一 subarray（参数表 §4.1） | 备选 hvector/hindexed 组合型描述 y 面：构建复杂、实现可能内部仍 pack，否决；y 面 wire 序 [k][h][i] 自洽（收发同型），无跨版本兼容负担（读方不存在——wire 只在本对收发间存在） |
| D4 | **sources 对称构造**（sources[i]=destinations[i]=n(d_i)，activeDirs_ 方向枚举序）+ **recv 块分情形落位** `f_i=(n(d_i)==self) ? opposite(d_i) : d_i`（唯一边落 halo(d_i)、自环平行边按出现序落 halo(opposite(d_i))）；in/out 数组独立传入（标准要求） | 入邻居多重集恒等（含 activeDirs 不对合自反的边界 rank——封闭性来自共享面两端各自声明）；落位全本地判定不依赖对端枚举序（第 1 轮对称+固定 halo(d_i) 在自环错配；第 1 轮修法「整体 opposite 序」又破坏笛卡尔边界——第 2 轮复审 Important，教训见 §4.3 推导节；三情形复算通过） |
| D5 | 类型/图/请求生命周期：initialize 建立、析构释放、全部防御 `MPI_Finalized` | 沿 PointToPointExchanger 析构惯例（p2p_exchanger.cpp:82-95）；alltoallw 非阻塞 request 为成员、end 后归 MPI_REQUEST_NULL |
| D6 | 性能判定记录式（无倍数门槛）：数据如实 + 劣化超 2× 触发默认值回退评审 | GUIDE C1 明言「小面 datatype 可能更慢」是领域已知；诚实优先（§6 底线） |
| D7 | `CollectiveExchangerDelegatesToP2P` 改写为 `CollectiveExchangerIsTrueCollective`（行为等价断言），不以日志文本断言委托 WARN 消失 | 日志字符串断言脆弱（重构即碎）；行为等价 + WARN 代码物理删除（diff 可见）是更强证据 |

## 8. 风险与缓解

| 风险 | 缓解 |
|------|------|
| subarray 维度序错误（AR003 陷阱） | §4.1 参数表 + U1-U3 自环测试逐方向暴露错位；与 packFace 逐面语义对照评审 |
| alltoallw displs 字节偏移算错 | U/B 系列 halo 带内容断言（非均匀切分下错位必现）；sdispls 相对 data 首址（不依赖 MPI_BOTTOM） |
| Ineighbor_alltoallw 平台支持 | MPI-3 标准（OpenMPI 1.7+/MPICH 3.0+），WSL OpenMPI 实测即证；不支持则编译期无该符号会立即暴露 |
| np4 本机漂移 | bench 沿 ±15% 惯例注记；np=1 作锚 |
| 行数预算 | 净增估算 ~500 行（实现 ~300 + 测试 ~200 + 脚本/文档），≤800 预算 |
