# AR005 I/O 基准 — after（改动后）与门槛判定

日期：2026-10-08；环境与口径同 baseline（`scripts/bench_ar005.sh after`，np=1，OMP=1，3 次中位数；io_write 剖面区段，`--save-interval 10` 循环内写出）。

## after 数据

| nx | format | io_write_sec（中位） | bytes |
|----|--------|---------------------|-------|
| 256 | vtk（appended raw） | 0.000528 | 524,750 |
| 256 | binary（块写） | 0.000374 | 524,344 |
| 512 | vtk（appended raw） | 0.001339 | 2,097,614 |
| 512 | binary（块写） | 0.001422 | 2,097,208 |

## 对比与 D6 门槛判定（median(after) ≤ median(baseline)）

| nx | format | baseline | after | 提速比 | 体积变化 | 门槛 |
|----|--------|----------|-------|--------|----------|------|
| 256 | vtk | 0.019120 s / 1,245,572 B | 0.000528 s / 524,750 B | **36.2×** | **-57.9%** | PASS |
| 256 | binary | 0.001253 s / 524,344 B | 0.000374 s / 524,344 B | 3.35× | 0（布局不变） | PASS |
| 512 | vtk | 0.065241 s / 4,981,124 B | 0.001339 s / 2,097,614 B | **48.7×** | **-57.9%** | PASS |
| 512 | binary | 0.005386 s / 2,097,208 B | 0.001422 s / 2,097,208 B | 3.79× | 0（布局不变） | PASS |

四组合全部满足 D6 门槛（after ≤ baseline）。

## 观察注记

- vtk(appended) 体积 = binary 后端数据镜像 + ~400B XML 头（256²：524,750-524,344=406B；512²：406B），与设计预期吻合（数据段同一内存镜像）
- binary 字节数逐位不变（524,344/2,097,208），块写为纯性能重构（U2 字节布局守护测试佐证）
- 提速来源：vtk 从「每元素 ~24 字符文本格式化」变为「每 (k,j) 行一次 8·nx 字节 memcpy 式写出」；binary 从每元素 1 次 write syscall（256²=65,536 次）变为每行 1 次（256 行）
- 512² binary（0.001422s）与 512² vtk（0.001339s）同量级——瓶颈已从格式化/syscall 转移到实际内存拷贝

## ParaView 可视化验证口径（srs §4 替代口径）

本地 WSL 无 ParaView 环境，按 srs §4 采纳替代口径：
1. **解析器测试**：U3 `VtkAppendedBinaryParsesBack` 即「最小 appended 解析器」——独立定位 `<AppendedData encoding="raw">`/`_` 标记、解码 UInt64 长度头、逐位还原 16 单元数值序列并与 u 场精确比对（位级）
2. **规范合规走查**：输出结构对照 VTK XML 文件格式规范要点——`header_type="UInt64"`（version 1.0 appended 必需）、`byte_order="LittleEndian"`、`format="appended" offset="0"`、`encoding="raw"` 的 `_` 前缀标记、长度头=数据字节数、行主序（x 最快）与 VTK ImageData 点序一致、`.pvti` PImageData 结构与 GhostLevel 声明
3. 局限性声明：未做真实 ParaView GUI 加载演示；规范走查覆盖已知易错点（offset 语义、长度头宽度、字节序），残余风险为 ParaView 对非主流属性组合的解析差异，属低概率
