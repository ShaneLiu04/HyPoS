# AR005 I/O 基准 — baseline（改动前）

日期：2026-10-08；环境：WSL Release 构建 /root/build-release，np=1，OMP=1，3 次中位数（`scripts/bench_ar005.sh baseline`，提交于 T002 前——此时代码为 ASCII vtk + 逐元素 binary 现状）。

| nx | format | io_write_sec（中位） | bytes |
|----|--------|---------------------|-------|
| 256 | vtk（ASCII） | 0.019120 | 1,245,572 |
| 256 | binary（逐元素） | 0.001253 | 524,344 |
| 512 | vtk（ASCII） | 0.065241 | 4,981,124 |
| 512 | binary（逐元素） | 0.005386 | 2,097,208 |

口径说明：
- 计时取 profiler `io_write` 区段（`--save-interval 10` 触发循环内写出；final 写出发生在 profiler 报告打印之后不进剖面——main.cpp 输出顺序，两模式口径一致可比）
- 体积取 `solution_10_r0.vti/.bin`（interior 数据 + 头部；vtk 含 XML 文本开销，binary 含 56B 头）
- 预期对比（改动后见 ar005-after 数据）：`.vti` appended binary 体积应接近 binary 后端（数据镜像 + 小 XML 头），耗时显著低于 ASCII；binary 块写耗时相对逐元素应有可见改善（syscall 次数 4.2M→4K 级）；门槛判据 median(after) ≤ median(baseline)（design D6）
