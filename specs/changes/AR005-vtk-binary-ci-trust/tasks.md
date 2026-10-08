# [AR005] 任务跟踪

| 字段 | 内容 |
|------|------|
| AR 编号 | AR005 |
| 关联 srs.md | ./srs.md |
| 关联 design.md | ./design.md（待生成）|
| 创建日期 | 2026-10-08 |

## 任务列表

| ID | 任务描述 | 依赖 | 状态 | 备注 |
|----|---------|------|------|------|
| T001 | F1-探针：新增 `MpiEnvTest.SizeProbeMatchesExpected`（`HYPOS_EXPECT_NP` 未设 SKIP/不匹配或非法值 FAIL）+ ctest：np>1 全部 12 条目 ENVIRONMENT 注入 + 8 条 gtest 条目 filter 追加 `:MpiEnvTest.*`（消费方）+ 新增 `mpi_size_probe` 条目（np=4）。测试先行：Red 期以单进程直跑探针二进制（HYPOS_EXPECT_NP=4）验证「假通过会被咬」+ 非法值 fail-safe 态 | - | passing | 2026-10-08 完成；计数勘误：np>1 条目实为 11（gtest 7 + hypos 4，design 第 2 轮门控 M1）；四态验证落 evidence/T001-probe-states.log（Red-1 咬缺陷 FAIL 含诊断/Red-2 abc fail-safe/未设 SKIP/Green PASS）；Release 26/26 + Debug(ASan+UBSan) 26/26 |
| T002 | D1-c BinaryIOBackend 块写：数据区 i 行连续段一次 `ofs.write(ptr, rowBytes)`；头部不变。测试先行：Red 期固化「同一 u 场输出逐位一致」fixture（块写重构前后文件 hash 相等）；块写耗时不劣化为走查口径 | T001 | passing | 2026-10-08 完成；U2 以「从零构造期望字节」固化布局（比对比旧实现更硬）——对现状 GREEN 基线 + 块写实现守护；halo 毒化 sentinel 防 halo 字节混入；Release 26/26 + Debug 26/26 零警告；after 基准并入 T006 统一跑 |
| T003 | D1-a VTK `.vti` appended binary：XML 头 `header_type="UInt64"` + `format="appended" offset` + `<AppendedData encoding="raw">` UInt64 长度头 + 行主序镜像；遍历序与旧 ASCII 一致。测试先行：更新 VtkPieceAndParallelIndexFiles + 新增最小 appended 解析器测试（还原数值序列） | T002 | passing | 2026-10-08 完成；U3 VtkAppendedBinaryParsesBack Red 子代理先行（现状 ASCII FAIL 确认）→ Green 重写 write()：binary 流 + appended payload (k,j) 行块写；offset=0 单数组免两遍写（D2）；既有 VtkPiece 测试只断言 Origin/extent 故未红——格式断言补强在 T004；Release 26/26 |
| T004 | D1-b `.pvti` 适配：PDataArray 与 piece DataArray 类型/名称一致；pvti 保持 ASCII。测试先行：pvti 断言更新（格式声明一致性） | T003 | passing | 2026-10-08 完成；VtkPieceAndParallelIndexFiles 补 piece 结构断言（header_type/format=appended/无 format="ascii"）+ pvti 零变更守护（GhostLevel="0"、PDataArray 原样）；payload 逐位校验已在 U3，此处仅结构——职责分离 |
| T005 | F1-CI 修复：ci.yml Debug(MPICH) job 换发行版或 pin 修复版（design 记录依据；无可行方案则回退 Debug(OpenMPI) 并记录）。交付：配置文件 + 本地等价命令验证记录 | T001 | passing | 2026-10-08 完成；runs-on ubuntu-latest→ubuntu-22.04（D4，整 job 随迁含 OpenMPI 段——design 门控第 1 轮已披露副作用）；gitee 主仓库不跑 Actions，本地 /usr/bin/mpicxx.mpich 等价构建即 T001-T005 期间 Debug 26/26 的持续验证；无回退必要 |
| T006 | 收尾：基准（ASCII vs binary 256²/512² 体积/耗时 ≥3 次中位数 + binary 块写前后对比 → PERFORMANCE 新小节）+ ParaView 加载验证（本地无环境则按 srs §4 口径以解析器测试 + 规范合规走查替代，evidence 记录）+ 文档同步（README mpich 记录/宣称、GUIDE 勾选 P5/D1/F1 紧急项、AGENT_SPEC 如有接口变更）+ 全量回归（Debug/Release + ASan/UBSan 零报告） | T001-T005 | passing | 2026-10-08 完成；after 基准四组合全过 D6 门槛（vtk 36.2×/48.7× + 体积 -57.9%，binary 3.35×/3.79× 字节数不变）；ParaView 口径按 srs §4 以 U3 位级解析器测试 + 规范走查替代（局限性声明入 evidence）；PERFORMANCE §12/GUIDE P5+D1+F1/README 探针+vtk 格式宣称同步；AGENT_SPEC 核验无接口变更；终跑 Release 26/26 + Debug(ASan+UBSan) 26/26 |

## 状态说明

- `pending`：待开始
- `in_progress`：进行中（当前会话）
- `passing`：开发完成，测试通过
- `failed`：测试失败，需修复

## 进度记录

> 每个开发会话结束后追加，记录完成情况。

## 阶段门控记录

> 由 sdd-phase-gate skill 在阶段门控审查后追加，记录每轮审查结果（PASS/FAIL + 轮次）。格式见 sdd-phase-gate SKILL.md Step 6。

### 2026-10-08 req 门控记录

- 门控结果：PASS（第 1 轮）
- 审查项数：9 项（0 NO；3 项 Minor WARN）
- 修复的问题：G5 T002-T006 状态列补齐 `pending`（已修复）；G6 T006 聚合判定为可接受（收尾惯例，延续既有 AR 模式，develop 阶段按会话拆分执行）；G9 性能项不设倍数阈值为有意设计（诚实优先，与 GUIDE §6 一致）——design 阶段补「不劣化」底线（binary/vtk 写耗时 ≤ 同场 ASCII 耗时）使其可判定
- 交叉验证：srs 与 OPTIMIZATION_GUIDE §5 ★2/§3 P5/D1/F1 对应准确；README:28 mpich 记录引用属实；gitee 主仓库不跑 Actions 的环境假设与 AGENTS.md 吻合
- 审查代理：sdd-gate-reviewer

### 2026-10-08 design 门控记录

- 门控结果：FAIL（第 1 轮）
- 失败项：G2/G15（Important，同根因）——探针 ENV 注入对 7 个带 `--gtest_filter` 的 gtest 条目为死配置（filter 不含 MpiEnvTest，无消费方）；4 个 hypos 直跑条目无测试体无法断言——与 srs §3.3 期望行为②不符且无偏离决策记录
- 严重性：Important；另 4 项 Minor（§6.5 测试总数算术 37→26；ci.yml 单 job 表述与整 job 迁移副作用未披露；覆盖率口径未声明；HYPOS_EXPECT_NP 非法值行为未定义）
- 审查代理：sdd-gate-reviewer
- 修复：①8 条 gtest 条目 filter 追加 `:MpiEnvTest.*`（消费方）+ 4 条 hypos 直跑条目环境级守卫论证（新增决策 D7），srs §3.3 措辞同步修订；②§6.5 总数改 26；③ci.yml 单 job 整体迁移表述 + OpenMPI 随迁注记；④§6 补覆盖率口径声明；⑤探针非法值按不匹配 FAIL（fail-safe），U1 补第四态断言

### 2026-10-08 design 门控记录（第 2 轮）

- 门控结果：PASS（第 2 轮）
- 审查项数：15 项（G10 SKIP；0 NO；2 项 Minor WARN）
- 第 1 轮修复核验：G2/G15 Important（探针死配置）确认修复——filter 消费方 + D7 环境级守卫 + srs 同步；4 项 Minor 全部到位
- 本轮新发现（2 Minor，已顺手修复）：M1 np>1 条目计数 8/12 → 7/11（枚举名单本就准确，仅计数口径）；M2 4 条 hypos 条目 ENV 注入无消费方 → 注明「声明性自文档，实际守卫由环境级探针承担（D7）」
- 事故记录：M1/M2 修正时误用 PowerShell 文本 cmdlet 处理 UTF-8 中文文件导致编码损坏，design.md 从上下文完整重建（教训：中文文件一律用 edit/write 工具）
- 审查代理：sdd-gate-reviewer
- 亮点：审查者逐行核验 filter 名单与 CMakeLists:139-174 吻合；ci.yml 单 job 结构、subgrid/aligned_buffer 行号引用全部属实
