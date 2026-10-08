# [AR003] 任务跟踪

| 字段 | 内容 |
|------|------|
| AR 编号 | AR003 |
| 关联 srs.md | ./srs.md |
| 关联 design.md | ./design.md（待生成）|
| 创建日期 | 2026-10-08 |

## 任务列表

| ID | 任务描述 | 依赖 | 状态 | 备注 |
|----|---------|------|------|------|
| T001 | 测试先行（红）：MPI-IO 单文件测试骨架——np=1 写出后 POSIX 读回，断言头部字段（magic/version/全局尺寸/间距/边界类型/offset）与数据区 vs 串行参考场逐位一致 | - | passing | 7 用例全绿（unit/mpiio_np1/mpiio_np4） |
| T002 | 实现 `MPIIOBinaryBackend`（绿，np=1）：自描述头写入 + 数据区 `MPI_File_write_at_all`；`writeParallelIndex` no-op；RAII 句柄管理 | T001 | passing | 修复头部写顺序（view 前）与 ORDER_C 维度序两个缺陷后全绿 |
| T003 | np=4 subarray 拼接正确性：`MPI_Type_create_subarray` filetype 各 rank 直写全局位置；测试读回全局场与参考一致（含子域边界无错位/覆盖）；3D 冒烟（64×64×8）数据区尺寸断言 | T002 | passing | FourRanks 两用例在真实 4 rank（OpenMPI）下执行并通过 |
| T004 | 异常路径：输出目录不可写时告警跳过、不崩溃、无 collective 挂死；告警限流；无 MPI_File 泄漏路径 | T002 | passing | X-1/X-2 全绿；泄漏走查记录 evidence/T004-leak-walkthrough.md |
| T005 | main.cpp 集成：`--output-format mpibin` 实例化新后端；`writeSolution` 对单文件后端用无 `_r<rank>` 后缀共享名并跳过 piece gather/索引；`--output-format` 白名单校验（非法退出码 1）；`--save-interval` 中间单文件；`HYPOS_PROFILE("io_write")` 剖面 | T002 | passing | mpiio_e2e（含 --save-interval 10 --enable-profiling）通过，报告含 io_write 区段 |
| T006 | ctest 扩充（`mpiio_np1`/`mpiio_np4`/`mpiio_e2e`）+ 全量回归：Debug/Release + ASan 全绿零报告 | T003, T004, T005 | passing | Debug 12/12 + Release 12/12（evidence/T-ctest-debug.log、T-ctest-release.log） |
| T007 | 基准参考：256² np=1/np=4（OMP=1，≥3 次中位数）mpibin vs binary 分片 io_write 剖面耗时对比 | T006 | passing | 中位数：mpibin 7.0/12.5ms vs binary 12.6/4.9ms（np=1/4）；evidence/T007-bench.log + PERFORMANCE §10 |
| T008 | 文档全套同步：README（格式表 + 路线图节更新）、DESIGN（IO 层新组件）、PERFORMANCE（基准数据）；文档走查无失实宣称 | T007 | passing | README/DESIGN/PERFORMANCE 全部更新；README 补 Ubuntu 24.04 mpich 陷阱提示 |

## 状态说明

- `pending`：待开始
- `in_progress`：进行中（当前会话）
- `passing`：开发完成，测试通过
- `failed`：测试失败，需修复

## 进度记录

> 每个开发会话结束后追加，记录完成情况。

### 2026-10-08 会话记录（代码完成，环境重启待验证）

- 完成任务：T001-T005 代码全部完成（Red 测试 7 用例 + Green 实现 + main 集成）；T006 CMake 注册完成；T008 README/DESIGN 草稿完成
- 修改文件：src/io/io_backend.hpp（新增 MPIIOBinaryBackend 声明）、src/io/mpiio_binary.cpp（新增）、src/main.cpp（白名单/实例化/分流/io_write 剖面）、tests/test_mpiio_output.cpp（新增）、CMakeLists.txt、README.md、docs/DESIGN.md
- 验证状态：**代码未经编译验证**——本机原无 WSL；已安装 WSL 3.0.1.0（内核 6.18.40.1-1），虚拟机平台组件等待系统重启生效
- 环境恢复步骤（重启后）：
  1. `wsl --import hypos D:\wsl\hypos D:\wsl\ubuntu-base.tar.gz`（rootfs 已下载 28.59MB 并校验通过）
  2. `wsl -d hypos -u root -- apt update && apt install -y build-essential cmake mpich libgtest-dev`（apt 慢可换国内镜像源）
  3. WSL 内构建：`cd /mnt/d/hy-po-s && cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug && cmake --build build-debug -j`（Release 同理）
  4. `ctest --test-dir build-debug --output-on-failure`（含新增 mpiio_np1/mpiio_np4/mpiio_e2e）
  5. 通过后回填 tasks.md 状态、补 T007 基准与 PERFORMANCE 数据、留存 evidence、提交代码
- 备注：req/design 门控均已 PASS；下载走公司代理较慢（约 120KB/s）

### 2026-10-08 会话记录 2（环境就绪，全绿 + 基准完成）

- 环境恢复完成：WSL 导入 hypos（ubuntu-base 24.04.3）→ apt 安装构建链 → Debug/Release 构建
- **环境级问题（重要）**：Ubuntu 24.04 的 `mpich 4.2.0` 包存在 PMI/PMIx 不匹配缺陷——`libmpich.so.12` 链接 `libpmix.so.2`，而同包 hydra 只下发 PMI-1 环境（PMI_FD/PMI_SIZE），应用静默退化为**单进程**（4 个独立 singleton），所有 np=4 测试假通过。已切换 OpenMPI（`-DMPI_CXX_COMPILER=/usr/bin/mpicxx.openmpi`），以最小 size 探针验证 `MPI_Comm_size==4` 后重建；README 已加陷阱提示
- 修复缺陷 3 个：
  1. 头部写在 `set_view` 之后——`MPI_File_write_at` 的 offset 是 view 相对，落点错到数据区首 72 字节 → 头部移至 set_view 之前（默认视图）
  2. `MPI_Type_create_subarray` 的 `MPI_ORDER_C` 维度序按 C 声明序（最慢维在前）传入，应为 `{z,y,x}` 而非 `{nx,ny,nz}` → 内存 filetype 选错元素（表现为数据区错位/全零）
  3. 测试 `fill2D` 用 `index(i,j)`（k 默认 0 平面）而 2D 内点在 k=1 平面 → 改用 `index(i,j,kBegin())`
- 验证：Debug 12/12 + Release 12/12（OpenMPI 真实多 rank，FourRanks 用例实跑非跳过）；证据 evidence/T-ctest-debug.log、T-ctest-release.log、T-fourranks-np4.log
- T007 基准完成（scripts/bench_t007.sh，evidence/T007-bench.log，PERFORMANCE.md §10）
- T004 泄漏走查记录：evidence/T004-leak-walkthrough.md
- 文档：README/DESIGN/PERFORMANCE 全部同步
- WSL 构建目录约定：`/root/build-debug`、`/root/build-release`（ext4，保证 IO 测试真实）
- 证据说明：`T-debug-ctest.log` 为调试中间态（红阶段）记录，保留供过程追溯；最终绿证据以 `T-ctest-debug.log`/`T-ctest-release.log`（16/16）为准

### 2026-10-08 审查记录（第 1 轮：FAIL → 修复）

- 审查结果：FAIL（无功能缺陷；1 Major 证据完整性 + 3 Minor 测试覆盖/表述）
- 修复的问题：
  1. **Major（S1）**：绿证据日志因命令行变量转义问题落盘到 WSL 根目录——已补交 evidence/（T-ctest-debug/release/fourranks-np4.log）
  2. **Minor（S1）**：E-3 格式白名单失败无自动化测试 → 新增 ctest `mpibin_bad_format`（WILL_FAIL）
  3. **Minor（S1/S4）**：S-5 仅最终解 / S-4 中间文件存在性无断言 → 新增 `mpiio_e2e_final`（无 --save-interval）与 `mpiio_e2e_verify`/`mpiio_e2e_final_verify`（cmake/VerifyMpiioOutput.cmake 精确文件集断言，fixture 串联）
  4. **Minor（S5）**：T004 走查记录「X-1 断言 WARN 次数」表述失准 → 更正为代码走查口径
- 回归：Debug 16/16 + Release 16/16（新增 4 个 ctest 后）；evidence/T-ctest-debug.log、T-ctest-release.log 已刷新
- 审查代理：sdd-gate-reviewer（独立复跑全部验证后给出结论）

## 阶段门控记录

> 由 sdd-phase-gate skill 在阶段门控审查后追加，记录每轮审查结果（PASS/FAIL + 轮次）。格式见 sdd-phase-gate SKILL.md Step 6。

### 2026-10-08 req 门控记录

- 门控结果：PASS（第 1 轮）
- 审查项数：9 项
- 修复的问题：G2/G3 §3.1 判据表述二义（已改为「np=1 逐位一致；np=4 ≤1e-15」）；G3 §3.3 泄漏检视留痕（T004 备注已注明记入 evidence）
- 审查代理：sdd-gate-reviewer

### 2026-10-08 design 门控记录

- 门控结果：PASS（第 1 轮）
- 审查项数：15 项（4 项 WARN，均 Minor）
- 修复的问题：G8 §6 补覆盖率目标声明；G12 补 S-7（io_write 剖面 e2e 用例）；G13 §6.7 补 rank0 头部写失败分支（走查口径）并显式标注 write_all 注入限制；G15 补 E-1b/E-1c（已带后缀/空文件名边界）并澄清 §4.3.3 空 filename 行为
- 审查代理：sdd-gate-reviewer
