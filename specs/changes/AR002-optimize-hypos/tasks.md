# [AR002] 任务跟踪

| 字段 | 内容 |
|------|------|
| AR 编号 | AR002 |
| 关联 srs.md | ./srs.md |
| 关联 design.md | ./design.md |
| 创建日期 | 2026-10-07 |

## 任务列表

| ID | 任务描述 | 依赖 | 状态 | 备注 |
|----|---------|------|------|------|
| T001 | R1 全局输出：VTK 改为每 rank `.vti` 分片 + rank0 `.pvti` 索引（Origin/Extent 偏移）；Binary 每 rank 文件 + 真实 offsets；Subgrid 携带 offsets | - | pending | 验收：np=4 四片+索引；np=1 全局等价 |
| T002 | R2 `--save-interval` 回调机制与中间输出（文件名含步号） | T001 | pending | 验收：50/100 步中间文件；默认零开销 |
| T003 | R3 收口：MPIEnv 接入 main；collective 告警限 rank0；Logger 注释更正；MemoryPool 接入或移除 | - | pending | 验收：无裸 MPI init；告警 1 条；无死代码 |
| T004 | R4 Neumann 接线：`--bc dirichlet\|neumann` + setupProblem 分支 + 非法值报错 | - | pending | 验收：镜像 halo 断言；默认行为不变 |
| T005 | R5 报告收口：comm_overhead_ratio 与 compute_time_ms 回填；scaling/bandwidth 字段移除；FLOPs 估算修正 | - | pending | 验收：比值一致；无恒 0 字段；估算 ±5% |
| T006 | R6a 残差融合（固定区域顺序）+ 退化子域 slab 归一化 | - | pending | 验收：on/off 一致 ≤1e-12；ASan 零报告 |
| T007 | R6b 持久化通信（Send_init/Recv_init + Startall/Waitall） | T006 | pending | 验收：全量测试绿；行为与一次性请求一致 |
| T008 | R6c 打包 memcpy 化（布局最内维连续）+ SIMD 对齐提示 + `-fopt-info-vec` 验证 | T007 | pending | 验收：向量化日志；行为一致 |
| T009 | R6 验收：256² np=1/4 优化前后基准对比（≥10%，3 次中位数）与证据存档 | T006,T007,T008 | pending | 验收：改善 ≥10%；基准更新 |
| T010 | R7a Red-Black GS 求解器（双色半扫 + 交换 + lastResidual） | T006 | pending | 验收：收敛、迭代数 ≤50% Jacobi、解一致 ≤1e-8 |
| T011 | R7b CG 求解器（matvec 复用 stencil、双内积归约） | T010 | pending | 验收：迭代数 ≤10% Jacobi、解一致 ≤1e-8 |
| T012 | R8a 工程化：`.gitignore`；ctest 新增 collective_mpi/ranks3_smoke/solver_alt_mpi；scaling 脚本参数化；plot 清理 | T010,T011 | pending | 验收：ctest 全绿含新用例 |
| T013 | R8b 文档同步：README/DESIGN/PERFORMANCE（新特性与优化说明） | T001-T012 | pending | 验收：无失实；使用方式齐全 |
| T014 | R8c/总验收：端到端 + 基准刷新 + SCALING_REPORT 优化对比 + 弱扩展重测 | T009,T012,T013 | pending | 验收：I 检查全过；报告含前后对比 |

## 状态说明

- `pending`：待开始
- `in_progress`：进行中（当前会话）
- `passing`：开发完成，测试通过
- `failed`：测试失败，需修复

## 进度记录

> 每个开发会话结束后追加，记录完成情况。

## 阶段门控记录

> 由 sdd-phase-gate skill 在阶段门控审查后追加，记录每轮审查结果（PASS/FAIL + 轮次）。

### 2026-10-07 req 门控记录

- 门控结果：PASS（第 1 轮）
- 审查项数：9 项（8 YES、1 WARN；无 NO）
- 修复的问题：4 项 Minor 已同步修订（R2 命名/glob 统一、R6 退化子域归一化判据补充、R5 字段收口范围扩展至 compute_time/bandwidth、§1 背景数据更正为 64² 实测）；另设计预研修正 3 处规格（VTK XML 格式族改为 `.vti/.pvti`、Neumann 纯边界可解性说明、SIMD aligned 断言可行性、求解器对比收敛口径）
- 审查代理：sdd-gate-reviewer
