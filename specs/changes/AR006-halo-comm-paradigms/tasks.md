# [AR006] 任务跟踪

| 字段 | 内容 |
|------|------|
| AR 编号 | AR006 |
| 关联 srs.md | ./srs.md |
| 关联 design.md | ./design.md（待生成）|
| 创建日期 | 2026-10-09 |

## 任务列表

| ID | 任务描述 | 依赖 | 状态 | 备注 |
|----|---------|------|------|------|
| T001 | C1 Red：datatype halo 测试先行——2D/3D 方向语义自环 + 多 halo + PROC_NULL + 与 pack 路径的 halo 带等价性断言（含 begin/end 拆分等价用例 E2，design §4.5；测试矩阵以「datatype 路径」跑；实现缺席时编译挂或 FAIL 即 Red） | - | pending | |
| T002 | C1 Green：DatatypeP2P halo 实现——六方向 subarray 型（发送型指向 interior 边带/接收型指向 halo 带）、非阻塞收发直传计算缓冲、datatype/请求生命周期管理（data 指针逐次不同的请求策略按 design 决策）；既有 halo 测试矩阵以 datatype 路径复跑全绿 | T001 | pending | |
| T003 | C2 Red：真集合测试先行——collective 语义断言（np4 非均匀/np8 3D 用例以 collective 跑 + 无委托 WARN 断言 + `CollectiveExchangerDelegatesToP2P` 改写为真集合断言） | T002 | pending | |
| T004 | C2 Green：CollectiveExchanger 真实现——`MPI_Dist_graph_create_adjacent` 邻居图 + `MPI_(I)Neighbor_alltoallw`（复用 T002 面 datatype），begin/end 语义走非阻塞变体；图通信子/alltoallw 资源生命周期 | T003 | pending | |
| T005 | 接线与测试矩阵：切换机制落地（design 决策）+ CLI 帮助文本 + ctest 条目（datatype/collective 路径的 halo 正确性与既有矩阵覆盖；np>1 条目带 HYPOS_EXPECT_NP） | T002, T004 | pending | |
| T006 | 收尾：bench（256² halo_exchange 剖面，np=1/4，pack vs datatype、p2p vs collective，≥3 次中位 → PERFORMANCE 新小节）+ 文档同步（README:264 路线图移除 + comm-mode 帮助、GUIDE P7/G2 勾选、AGENT_SPEC 类层级）+ 全量回归（Debug/Release + ASan/UBSan 零报告） | T001-T005 | pending | |

## 状态说明

- `pending`：待开始
- `in_progress`：进行中（当前会话）
- `passing`：开发完成，测试通过
- `failed`：测试失败，需修复

## 进度记录

> 每个开发会话结束后追加，记录完成情况。

## 阶段门控记录

> 由 sdd-phase-gate skill 在阶段门控审查后追加，记录每轮审查结果（PASS/FAIL + 轮次）。格式见 sdd-phase-gate SKILL.md Step 6。

### 2026-10-09 req 门控记录

- 门控结果：PASS（第 1 轮）
- 审查项数：10 项（0 NO；1 项 Minor WARN——已顺手修复：srs R4 补 halo_exchanger.hpp CollectiveExchanger 过时 doc 注释的同步项）
- 交叉验证亮点：G6 实核 solver 调用方（CG exchange 传独立 p 缓冲 cg_solver.cpp:44），「持久请求不得绑定调用方指针」约束确认为真；G7 抽查全部 8 halo 用例均断言 halo 带内容而非 wire 布局，布局变化不敏感假设成立
- 审查代理：sdd-gate-reviewer

### 2026-10-09 design 门控记录

- 门控结果：FAIL（第 1 轮）
- 失败项：G3/G5（Important）——①FP3 对称 sources 构造在平行边（B5 自环）下按边序错配（halo(d) 错收 interior(d)，应收对侧）；②np8 filter（`*Datatype3D*`）与 B2/B4 用例名不匹配，两用例无处执行、collective_np8 空洞绿
- 另 5 项 Minor（析构在途请求防御/3D hw=2 覆盖/四处行号勘误/E2 未点名/2D 措辞）——全部顺手修复
- 修复：①recv 侧改 opposite 序（sources[i]=n(opposite(d_i))、recv 块落 halo(opposite(d_i))），笛卡尔与自环两情形推导验证入 §4.3 + D4 更新；②ctest filter 改专用点名式（沿 halo_2d_mpi 惯例）；③U4 定为 3D hw=2 配置兼覆盖 3D 多 halo；④DatatypeExchanger 析构补在途请求防御；⑤行号勘误（collective_mpi :182-185/unit :87/README :267）；⑥G1 六面参数表逐面核对全 OK（审查者逐行对照 packFace/unpackFace）
- 审查代理：sdd-gate-reviewer
