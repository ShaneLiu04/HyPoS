# [AR002] ST 验收测试用例

| 字段 | 内容 |
|------|------|
| AR 编号 | AR002 |
| 关联 srs.md | ./srs.md |
| 生成日期 | 2026-10-08 |
| 执行环境 | WSL2 Ubuntu / Intel Ultra 9 185H / OpenMPI 5.0.10 + MPICH(ASan) / gcc 15.2 / OMP=4（基准 OMP=1） |

## 测试用例列表（实际结果均已执行，证据见 logs/）

| ID | 关联需求 | 场景 | 期望 | 实际结果 | 状态 |
|----|---------|------|------|---------|------|
| ST-01 | srs 3.1 | np=4/np=1 `--output-format vtk` | 4×`.vti`+`.pvti`（Extent 精确拼合、Origin=offset*dx）；np=1 单分片全局 | 4 片拼合 32²、Origin 正确；np=1 单片 | PASS |
| ST-02 | srs 3.1 | np=4 `--output-format binary` | 4×`.bin`，头部真实 offsets | 4 文件；r1 头 (16,16,1,1,0,16,0)（logs/I15.log） | PASS |
| ST-03 | srs 3.2 | `--save-interval 50 --max-iter 120/30/默认` | 50/100 中间解；未到步数无中间件；默认仅最终解 | 三组均符合（logs/I9.log；T002 会话） | PASS |
| ST-04 | srs 3.3 | main 无裸 MPI 调用 / collective 告警 / MemoryPool 死代码 | 0/1 条/0 引用 | grep 0/1/0（T003 会话） | PASS |
| ST-05 | srs 3.4 | `--bc neumann` np=1/4；非法值 | 完成且镜像语义（U-BC-N）；非法退出 1 | 冒烟 exit 0（logs/I16-*.log）；非法 exit 1（logs/E-BC.log） | PASS |
| ST-06 | srs 3.5 | np=4 报告字段 | comm_overhead≈comm/iter、compute_time>0、无 scaling/bandwidth、FLOPs 公式一致 | 偏差 0.04%/0.12%，字段收口（T005 会话 + I14） | PASS |
| ST-07 | srs 3.6 | 受控性能对比（256²×10000，3 次中位数） | ≤AR001 基线 90% | np1 -37.5%、np4 -39.0%（logs/T009-perf-summary.log） | PASS |
| ST-08 | srs 3.6 | on/off 位级一致 + 退化归一化 + 向量化 | memcmp=0；Debug 自检；-fopt-info-vec 有向量化 | 单测通过；3 处 vectorized（logs/T008-vec-build.log） | PASS |
| ST-09 | srs 3.7 | RBGS/CG vs Jacobi（64²，vp=1/4） | 收敛、迭代 ≤60%/≤10%、交叉 ≤1e-8、3D 冒烟 | RBGS 52.3%、CG 215 步（特征向量制造解 1 步）；全部断言通过 | PASS |
| ST-10 | srs 3.8 | ctest 全量（Release + Debug+ASan） | 全绿（含 collective/3rank/alt solver/parse_error） | 9/9 + 9/9（logs/T014-ctest-*.log） | PASS |
| ST-11 | srs 3.8 | 文档与 .gitignore 走查 | 无失实 | D1 核查 + 审查修订后通过（审查报告两轮） | PASS |
| ST-12 | srs 3.8 | 端到端异常 | E-Parse/E-BC/E-OutDir 行为 | exit 1/1/0+warn（logs/E-*.log） | PASS |

## 执行摘要

| 总计 | 通过 | 失败 | 阻塞 |
|------|------|------|------|
| 12 | 12 | 0 | 0 |

## ST 执行报告

**需求覆盖：** R1-R8 = 8/8（100%）；**回归：** 全量 ctest 9/9 双构建，无回归。
**遗留问题：** 1 项信息级（退化子域 slab 语义已归一化并在 Debug 自检，无残留）；单机 WSL np=4 scaling 快照存在调度噪声（受控基准为准，已文档化）。
**优化成效（受控口径）：** iter_time np1 -37.5%、np4 -39.0%；RBGS 迭代 ≈ 0.52×Jacobi；CG 一般右端项降低 1-2 个数量级。

### 结论

> **Go** — 覆盖率 100%，无 Critical/Major 缺陷，双构建测试全绿，性能目标（≥10%）大幅超额达成，文档与实现一致。
