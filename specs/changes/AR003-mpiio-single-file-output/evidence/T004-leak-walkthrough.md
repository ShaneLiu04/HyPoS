# T004 / X-4 泄漏与异常路径走查记录

- 走查对象：`src/io/mpiio_binary.cpp`（最终版，随 12/12 绿提交）
- 走查方式：代码走查（ASan 已启用但 ctest 环境 `ASAN_OPTIONS=detect_leaks=0`——WSL root 下 LeakSanitizer 默认关闭，故按 design §6.7 以走查口径覆盖）
- 日期：2026-10-08

## 结论：无泄漏路径

| 资源 | 创建点 | 释放点 | 提前返回路径 |
|------|--------|--------|--------------|
| `MPI_File fh` | `MPI_File_open`（line 60） | `MPI_File_close`（line 149，无条件） | 唯一提前返回在 open 失败（line 67），此时 `fh == MPI_FILE_NULL`，无资源存在 |
| `MPI_Datatype fileType` | line 94-95 | `MPI_Type_free`（line 147，无条件） | commit 后无提前返回 |
| `MPI_Datatype memType` | line 105-106 | `MPI_Type_free`（line 148，无条件） | 同上 |
| 全局 errhandler | 未修改 | 无需恢复 | MPI-IO 文件默认 errhandler 即 `MPI_ERRORS_RETURN`（MPI 标准），未调用 `MPI_File_set_errhandler`/`MPI_Comm_set_errhandler` |

## collective 一致性

- open 失败：`MPI_ERRORS_RETURN` 下失败在**所有 rank** 一致返回 → 提前 return 不产生集合操作失配（测试 X-1 在 np=1/4 均验证无挂死）
- `set_size`/`set_view`/`write_all`/header `write_at` 失败：仅告警（`warned_` 限流，每 rank 首次），**继续走到 close** → 集合路径各 rank 对齐

## 未检查返回值（走查标注，可接受）

- `MPI_Type_create_subarray`/`MPI_Type_commit`/`MPI_File_set_view`：参数由合法网格推导，失败会随后在 `write_all` 报错并进入告警路径，不产生泄漏或挂死

## 测试证据

- X-1（不可写路径）/X-2（截断）：`unit`、`mpiio_np1`、`mpiio_np4` 全绿（见 T-ctest-debug.log / T-ctest-release.log）
- 告警限流：X-1 断言 WARN 仅出现一次/进程（warn-once 闩锁）
