# C++ vs Rust 性能对比报告

日期：2026-04-17  
工作区：`/Users/henrykang/Documents/CavalierContours`

## 1. 范围与版本

本报告对比以下两套实现：

- C++ 仓库：当前工作区代码，commit `384b3cf`
- Rust 参考实现：`Downloads/cavalier_contours/cavalier_contours`

说明：

- Rust 目录位于当前工作区内，不是独立 git 仓库，因此没有单独的 commit 可记录。
- Rust crate 版本来自 `Cargo.toml`，为 `0.7.0`。

## 2. 测试与基准方法

### 2.1 C++ 侧

直接使用现有可执行文件：

- 测试：
  - `TEST_sample`
  - `TEST_cavc_pline`
  - `TEST_cavc_pline_function`
  - `TEST_cavc_parallel_offset`
  - `TEST_cavc_combine_plines`
  - `TEST_staticspatialindex`
  - `TEST_cavc_api_regression`
  - `TEST_cavc_offset_islands`
  - `TEST_cavc_internal_slice_view`
- 基准：
  - `offsetbenchmarks`
  - `combinebenchmarks`
  - `extentsbenchmarks`
  - `windingnumberbenchmarks`

### 2.2 Rust 侧

Rust 仓库内没有现成 `benches/` 或 `criterion` harness，因此本次使用了一个**临时本地 benchmark harness**，复现 C++ `tests/benchmarks/benchmarkprofiles.h` 中的 profile 数据与 benchmark 循环体：

- `offset`：与 C++ 相同的 profile、相同的 offset 迭代次数和正负偏移调用
- `combine`：与 C++ 相同的 shifted profile 组装逻辑和 coincident 调用模式
- `extents` / `winding number`：与 C++ benchmark 的输入构造方式保持一致

计时方式：

- Rust 侧通过 `std::time::Instant` 运行到最小持续时间后，取每次迭代的平均耗时
- 因计时框架不同，本报告更适合做**方向性对比**，不应视为严格 publication-grade benchmark

## 3. 测试结果

### 3.1 C++

以下测试全部通过：

- `TEST_sample`
- `TEST_cavc_pline`
- `TEST_cavc_pline_function`
- `TEST_cavc_parallel_offset`
- `TEST_cavc_combine_plines`
- `TEST_staticspatialindex`
- `TEST_cavc_api_regression`
- `TEST_cavc_offset_islands`
- `TEST_cavc_internal_slice_view`

### 3.2 Rust

执行命令：

```bash
cargo test -p cavalier_contours --release
```

结果：

- 存在 1 个失败测试：
  - `tests/test_pline_view.rs` 中的 `attempting_to_wrap_slice_on_open_pline`
- 该失败看起来是 Rust 参考仓库当前自身状态中的现有问题，不是本次 C++ 改动引入

为确认其余测试状态，额外执行：

```bash
cargo test -p cavalier_contours --release -- --skip attempting_to_wrap_slice_on_open_pline
```

结果：

- 其余 unit tests、integration tests、doc tests 全部通过

## 4. 性能对比

表中“倍率”按“快的一方 / 慢的一方”给出。

### 4.1 Offset

| Case | C++ | Rust | 结论 |
|---|---:|---:|---|
| `Profile1` | 0.303 ms | 0.141 ms | Rust 快约 2.14x |
| `Profile2` | 0.628 ms | 0.307 ms | Rust 快约 2.05x |
| `Pathological1/50` | 15.4 ms | 9.85 ms | Rust 快约 1.56x |
| `Pathological1/100` | 58.2 ms | 38.1 ms | Rust 快约 1.53x |
| `Profile1NoArcs` | 4.65 ms | 1.93 ms | Rust 快约 2.41x |
| `Profile2NoArcs` | 10.0 ms | 4.06 ms | Rust 快约 2.46x |

结论：

- Rust 在 `offset` 上仍然明显领先
- 当前 C++ 虽然已经比之前更接近 Rust，但在 `offset` 上仍大约落后 1.5x 到 2.5x

### 4.2 Boolean / Combine

#### 4.2.1 Shifted workload

| Case | C++ | Rust | 结论 |
|---|---:|---:|---|
| `combine16Shifted Profile1` | 159 us | 136 us | Rust 快约 1.17x |
| `combine16Shifted Profile2` | 264 us | 221 us | Rust 快约 1.20x |
| `combine16Shifted Pathological1/50` | 501 us | 433 us | Rust 快约 1.16x |
| `combine16Shifted Pathological1/100` | 614 us | 526 us | Rust 快约 1.17x |
| `combine16Shifted Profile1NoArcs` | 351 us | 306 us | Rust 快约 1.15x |
| `combine16Shifted Profile2NoArcs` | 728 us | 668 us | Rust 快约 1.09x |

结论：

- Rust 在一般 `combine` workload 下仍领先
- 但差距已经缩小到约 1.1x 到 1.2x，属于“接近但仍落后”

#### 4.2.2 Coincident workload

| Case | C++ | Rust | 结论 |
|---|---:|---:|---|
| `combineCoincident Profile1` | 0.106 us | 4.22 us | C++ 快约 39.8x |
| `combineCoincident Profile2` | 0.109 us | 8.09 us | C++ 快约 74.2x |

结论：

- 在“两个输入完全相同”的 coincident case 上，当前 C++ 明显反超
- 主要原因是当前 C++ 代码中对**完全相同输入**加入了直接返回结果的快路径

### 4.3 Extents

| Case | C++ | Rust | 结论 |
|---|---:|---:|---|
| `Profile1` | 119 ns | 108.8 ns | Rust 快约 1.09x |
| `Profile2` | 181 ns | 168.2 ns | Rust 快约 1.08x |
| `Pathological1/100` | 2839 ns | 2535.8 ns | Rust 快约 1.12x |
| `Profile1NoArcs` | 221 ns | 129.2 ns | Rust 快约 1.71x |
| `Profile2NoArcs` | 461 ns | 263.3 ns | Rust 快约 1.75x |

结论：

- Rust 在 `extents` 上仍普遍领先
- 一般 case 差距不大，但在 `NoArcs` case 上 Rust 领先较明显

### 4.4 Winding Number

| Case | C++ | Rust | 结论 |
|---|---:|---:|---|
| `Profile1` | 1.58 us | 1.11 us | Rust 快约 1.43x |
| `Profile2` | 2.64 us | 1.87 us | Rust 快约 1.41x |
| `Pathological1/100` | 44.3 us | 15.2 us | Rust 快约 2.91x |
| `Profile1NoArcs` | 8.09 us | 7.87 us | Rust 略快约 1.03x |
| `Profile2NoArcs` | 14.9 us | 17.0 us | C++ 快约 1.14x |

结论：

- Rust 仍普遍更快，尤其在复杂曲线 case 上优势明显
- 但在个别 `NoArcs` case 上，当前 C++ 已经可以反超

## 5. 总结

这轮对比可以得出几个比较清楚的结论：

1. `offset` 仍然是当前 C++ 相对 Rust 落后最明显的板块
   - 大部分 case Rust 仍领先 1.5x 到 2.5x
   - 如果下一轮继续做性能优化，`offset` 仍应保持最高优先级

2. `boolean/combine` 已经明显逼近 Rust
   - 一般 shifted workload 下差距收敛到约 1.1x 到 1.2x
   - 这说明最近这轮共享 slice-view 的重构是有实际收益的

3. `combine` 的完全 coincident case 已经反超
   - 当前 C++ 在 identical-input 情况下利用快路径获得了数量级优势
   - 这是一个明确的、可归因的正向结果

4. `extents` 和 `winding number` 仍然整体偏慢
   - `extents` 差距中等
   - `winding number` 在复杂 case 上差距仍较大，值得继续优化

## 6. 主要命令

### 6.1 C++ 测试

```bash
./build/mac_clang_release/TEST_sample
./build/mac_clang_release/TEST_cavc_pline
./build/mac_clang_release/TEST_cavc_pline_function
./build/mac_clang_release/TEST_cavc_parallel_offset
./build/mac_clang_release/TEST_cavc_combine_plines
./build/mac_clang_release/TEST_staticspatialindex
./build/mac_clang_release/TEST_cavc_api_regression
./build/mac_clang_release/TEST_cavc_offset_islands
./build/mac_clang_release/TEST_cavc_internal_slice_view
```

### 6.2 C++ benchmark

```bash
./build/mac_clang_release/offsetbenchmarks --benchmark_filter='BM_offset(Profile1|Profile2)$' --benchmark_min_time=1.0s
./build/mac_clang_release/offsetbenchmarks --benchmark_filter='BM_offset(Pathological1/50|Pathological1/100)' --benchmark_min_time=1.0s
./build/mac_clang_release/offsetbenchmarks --benchmark_filter='BM_offset(Profile1NoArcs|Profile2NoArcs)' --benchmark_min_time=1.0s

./build/mac_clang_release/combinebenchmarks --benchmark_filter='BM_combine16Shifted(Profile1|Profile2|Pathological1/50|Pathological1/100)|BM_combineCoincident(Profile1|Profile2)' --benchmark_min_time=0.5s

./build/mac_clang_release/combinebenchmarks --benchmark_filter='BM_combineCoincident(Profile1|Profile2)$' --benchmark_min_time=2.0s

./build/mac_clang_release/extentsbenchmarks --benchmark_filter='BM_extents(Profile1|Profile2|Pathological1/100|Profile1NoArcs|Profile2NoArcs)' --benchmark_min_time=1.0s

./build/mac_clang_release/windingnumberbenchmarks --benchmark_filter='BM_windingNumber100PtGrid(Profile1|Profile2|Pathological1/100|Profile1NoArcs|Profile2NoArcs)' --benchmark_min_time=1.0s
```

### 6.3 Rust 测试

```bash
cargo test -p cavalier_contours --release
cargo test -p cavalier_contours --release -- --skip attempting_to_wrap_slice_on_open_pline
```

### 6.4 Rust benchmark

Rust benchmark 使用临时本地 harness：

```bash
cargo run --release
```

工作目录：

```text
/tmp/cavalier_rust_compare
```

## 7. 注意事项

- 本报告中的 Rust benchmark 并非来自仓库内官方 `benches/`，而是临时复现 C++ benchmark profile 的本地 harness
- 因为 C++ 使用 Google Benchmark，Rust 使用 `Instant` 平均计时，所以绝对值请谨慎解读
- 但在同机、同 profile、同负载结构下，这些结果足以作为本轮优化决策的方向性依据
