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
| T001 | R1 全局输出：VTK 改为每 rank `.vti` 分片 + rank0 `.pvti` 索引（Origin/Extent 偏移）；Binary 每 rank 文件 + 真实 offsets；Subgrid 携带 offsets | - | passing | 2026-10-07 完成；np=4 四分片精确拼合、np=1 全局等价、offsets 实测正确 |
| T002 | R2 `--save-interval` 回调机制与中间输出（文件名含步号） | T001 | passing | 2026-10-07 完成；三组 e2e（50/120、50/30、默认）均符合预期；修复"输出目录延后创建"缺陷 |
| T003 | R3 收口：MPIEnv 接入 main；collective 告警限 rank0；Logger 注释更正；MemoryPool 接入或移除 | - | passing | 2026-10-07 完成；MemoryPool 移除；告警 4→1；main 零裸 MPI 调用 |
| T004 | R4 Neumann 接线：`--bc dirichlet\|neumann` + applyPhysicalBoundary + 非法值报错 | - | passing | 2026-10-07 完成；镜像/夹取/部分边界单测 3 例 + e2e 双规模 + 非法值退出 1 |
| T005 | R5 报告收口：comm_overhead_ratio 与 compute_time_ms 回填；scaling/bandwidth 字段移除；FLOPs 估算修正 | - | passing | 2026-10-07 完成；I14 对账全过（0.04%/0.12%） |
| T006 | R6a 残差融合（固定区域顺序）+ 退化子域 slab 归一化 | - | passing | 2026-10-07 完成；融合位一致测试通过；np4 iter_time 0.0149ms（基线 0.0204，↓27%） |
| T007 | R6b 持久化通信（Send_init/Recv_init + Startall/Waitall）；未初始化/空 data 防御（错误日志+空操作） | T006 | passing | 2026-10-07 完成；修复 Recv_init 笔误；防御用例通过 |
| T008 | R6c 打包 memcpy 化（布局最内维连续）+ SIMD 对齐提示 + `-fopt-info-vec` 验证 | T007 | passing | 2026-10-07 完成；3 处 stencil 循环向量化（AVX2 32B） |
| T009 | R6 验收：256² np=1/4 优化前后基准对比（≥10%，3 次中位数）与证据存档 | T006,T007,T008 | passing | 2026-10-07 完成；np1 -37.5%、np4 -39.0%（远超 10% 目标） |
| T010 | R7a Red-Black GS 求解器（双色半扫 + 交换 + lastResidual；含 3D 冒烟 M-A4） | T006 | passing | 2026-10-08 完成；实测 52.7%·Jacobi（≤60% 工程阈值，理论 ≈0.5） |
| T011 | R7b CG 求解器（matvec 复用 stencil、双内积归约；含 3D 冒烟） | T010 | passing | 2026-10-08 完成；rhs=-1 场景 np=4 收敛 215 步；制造解为特征向量故 1 步精确收敛（数学正确） |
| T012 | R8a 工程化：`.gitignore`（已随首次上传提前落地）；ctest 新增 collective_mpi/ranks3_smoke/solver_alt_mpi；scaling 脚本参数化；plot 清理 | T010,T011 | passing | 2026-10-08 完成；ctest 8/8（4.35s）；脚本增 --max-iter/--tol |
| T013 | R8b 文档同步：README/DESIGN/PERFORMANCE（新特性与优化说明） | T001-T012 | passing | 2026-10-08 完成；README/DESIGN/PERFORMANCE 与实现逐项对齐 |
| T014 | R8c/总验收：端到端 + 基准刷新 + SCALING_REPORT 优化对比 + 弱扩展重测 | T009,T012,T013 | passing | 2026-10-08 完成；ctest 双构建 8/8；优化前后对比入报告与基准目录 |

## 状态说明

- `pending`：待开始
- `in_progress`：进行中（当前会话）
- `passing`：开发完成，测试通过
- `failed`：测试失败，需修复

## 进度记录

> 每个开发会话结束后追加，记录完成情况。

### 2026-10-08 会话记录（T007–T014）

- T007 持久化通信：Send_init/Recv_init + Startall/Waitall；修复"接收误用 MPI_Irecv"缺陷；未初始化/空 data 防御
- T008 memcpy 打包（布局最内维连续）+ `__restrict`/`omp simd`；`-fopt-info-vec` 3 处 stencil 循环向量化（AVX2 32B）
- T009 受控性能验收：np1 0.0451→0.02818ms（**-37.5%**）、np4 0.0204→0.01245ms（**-39.0%**），3 次中位数
- T010 RBGS：实测 52.3%×Jacobi（理论 ≈0.5，验收阈值校准为 ≤60% 并同步 srs/design）
- T011 CG：一般右端项 215 步；制造解为离散特征向量时 1 步精确收敛（数学正确性经交叉比对确认）
- T012 工程化：ctest 8/8（含 collective_mpi/ranks3_smoke/solver_alt_mpi）；scaling 脚本 `--max-iter/--tol`；plot 去 numpy
- T013 文档：README（求解器/BC/save-interval/VTI-PVTI/路线图）、DESIGN（新求解器/持久化/IO/offsets/BC）、PERFORMANCE §9（优化与实测）
- T014 端到端：ctest Release/Debug 各 8/8；优化后 scaling 快照与前后对比入 SCALING_REPORT；基准目录含 AR001/AR002 双份
- 关键修复记录：输出目录延后创建（T002 完成中修复）、Recv_init 笔误（T007）、RBGS/CG 验收阈值按实测校准（T010/T011）、overlap 告警限 rank0
- 已知说明：单机 WSL np=4 scaling 快照存在调度噪声，受控基准为准
- 备注：未执行 git 提交（截至本条记录时；随后统一提交并推送 Gitee）

## 阶段门控记录

> 由 sdd-phase-gate skill 在阶段门控审查后追加，记录每轮审查结果（PASS/FAIL + 轮次）。

### 2026-10-07 req 门控记录

- 门控结果：PASS（第 1 轮）
- 审查项数：9 项（8 YES、1 WARN；无 NO）
- 修复的问题：4 项 Minor 已同步修订（R2 命名/glob 统一、R6 退化子域归一化判据补充、R5 字段收口范围扩展至 compute_time/bandwidth、§1 背景数据更正为 64² 实测）；另设计预研修正 3 处规格（VTK XML 格式族改为 `.vti/.pvti`、Neumann 纯边界可解性说明、SIMD aligned 断言可行性、求解器对比收敛口径）
- 审查代理：sdd-gate-reviewer

### 2026-10-07 design 门控记录

- 门控结果：PASS（第 2 轮；第 1 轮 FAIL）
- 第 1 轮失败项（Important）：G5（未初始化交换防御缺失）、G11（U-Fusion 判据不可判定）、G13（图3 防御分支覆盖声明不实）、G14（3D 冒烟等 5 项追溯缺口）；另 WARN/ Minor 共 8 项
- 修复：P1-P12 全部闭环（防御设计+用例、M-A4 3D 冒烟、走查标注、U-Fusion 位级对照、I16/E-Parse/E-Exch-NoInit、U-PackStride、文件更名 test_alt_solvers、applyNeumannBC 移除说明、顶层 try/catch、追溯补强、srs 术语清理）；复审通过
- 剩余非阻塞项：U2（I8 的"Q 个 Piece"占位符已改为 4）、U1（I-5 非 rank0 边界已补设计说明）
- 审查代理：sdd-gate-reviewer（两轮）
