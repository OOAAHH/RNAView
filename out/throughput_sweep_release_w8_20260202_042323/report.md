# RNAVIEW 高通量吞吐对比报告（Legacy vs Rust *release*）

本报告基于目录 `out/throughput_sweep_release_w8_20260202_042323/` 的实测数据，目标是用“跑库/高通量”的视角解释：
1) Rust（release）相比 legacy 的吞吐提升有多大；2) “启动/初始化/落盘”开销在不同任务量下占比是多少；3) 哪些输入最拖慢整体吞吐。

---

## 1. 测试怎么跑的（一句话版）

用 `tools/rnaview_throughput.py compare` 在 **8 个 worker 并发**下，分别让：
- **legacy**：跑 `bin/rnaview`
- **rust**：跑 `rust/target/release/rnaview-hotcore`（`from-structure --oracle compute`）

对同一批输入做吞吐对比，并用 `--repeat` 把输入集复制扩增（模拟高通量跑库）。

复现命令（示例）：
```bash
python3 tools/rnaview_throughput.py compare --workers 8 --repeat 20 \
  --legacy-bin bin/rnaview \
  --rust-bin rust/target/release/rnaview-hotcore
```

---

## 2. 指标怎么理解（尽量口语）

每个 repeat 都会生成一个 `throughput_compare.json`，里面最关键的几项：

- `jobs = input_count * repeat`：本轮总共处理多少个结构文件（任务数）。
- `throughput_jobs_per_sec`：吞吐（越大越好），= `executed_jobs / (elapsed_ms/1000)`。
- `elapsed_ms`：`rnaview_batch.py` **真正开始跑 jobs 之后**到结束前的时间（不含进程冷启动等）。
- `wall_ms`：外层 wrapper 看到的 `rnaview_batch.py` 子进程**从启动到退出**的总时间。
- `overhead_ms = wall_ms - elapsed_ms`：可以粗略理解成“启动/初始化 + 收尾落盘”的固定成本。
  - 它不仅包含 Python import/参数解析，也包含 **写 `summary.json`** 等收尾动作，所以 jobs 越多，`summary.json` 越大，`overhead_ms` 也会略增。
- `overhead/job_ms`：把 overhead 摊到每个任务上（jobs 越多，这个数会显著下降）。

---

## 3. 总览结论（最重要的一页）

结果汇总见：`out/throughput_sweep_release_w8_20260202_042323/summary.md`。

结论用一句话说：**在这个回归集（14 个输入）上，Rust release 在 8 并发下吞吐大约是 legacy 的 7.8×～15.9×；当任务量变大时 speedup 稳定在 ~8×～11×。**

表格（来自 `summary.md`）：

| repeat | jobs | speedup | legacy thr(j/s) | rust thr(j/s) | legacy overhead_ms | rust overhead_ms | legacy overhead/job_ms | rust overhead/job_ms |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 14 | 15.942 | 0.734 | 11.696 | 215 | 182 | 15.357 | 13.000 |
| 2 | 28 | 12.259 | 1.434 | 17.577 | 201 | 186 | 7.179 | 6.643 |
| 5 | 70 | 8.940 | 3.136 | 28.034 | 282 | 229 | 4.029 | 3.271 |
| 10 | 140 | 10.767 | 3.182 | 34.263 | 300 | 278 | 2.143 | 1.986 |
| 20 | 280 | 9.300 | 3.857 | 35.865 | 403 | 384 | 1.439 | 1.371 |
| 50 | 700 | 7.801 | 4.288 | 33.450 | 671 | 737 | 0.959 | 1.053 |

怎么读这张表：
- **任务越多（jobs 越大），两边的吞吐都会上升**：并发 worker 更容易被“喂饱”，固定开销被摊薄。
- **overhead/job_ms 会随 jobs 变大迅速下降**：这就是你想量化的“启动/初始化/落盘开销在不同任务量下的影响”。
  - 例如 legacy：从 ~15ms/job（14 个任务）降到 <1ms/job（700 个任务）。
  - rust 也类似：从 ~13ms/job 降到 ~1ms/job。

这些列都是分别对 legacy 引擎和rust 引擎算出来的同一组指标（单位见括号）：

legacy thr(j/s) / rust thr(j/s)：吞吐量（jobs per second）。
公式：thr = executed_jobs / (elapsed_ms / 1000)；越大越快。

legacy overhead_ms / rust overhead_ms：单次运行的“固定开销”（毫秒）。
公式：overhead_ms = wall_ms - elapsed_ms

wall_ms：外层包装脚本看到的子进程从启动到退出的总墙钟时间
elapsed_ms：rnaview_batch.py 内部统计的“开始跑 jobs 之后”的时间
所以这个差值大致代表 Python 启动/import、参数解析、检查二进制、创建目录、写 log/summary 等非核心计算开销。
legacy overhead/job_ms / rust overhead/job_ms：把固定开销摊到每个 job 上（毫秒/任务）。
公式：overhead/job_ms = overhead_ms / executed_jobs（你这次数据里基本等于 / jobs）。

---

## 4. 为什么 Rust 这么快：关键是“最慢的那几个 case”被大幅加速

在高通量场景里，整体耗时经常被“最慢的那几个输入”支配（类似跑大规模任务时的长尾）。

以 repeat=50（每个 input 重复 50 次）为例，按 **每个输入的 median(中位数) 单任务耗时**统计：

| input | legacy median (ms) | rust median (ms) | speedup | legacy p95 (ms) | rust p95 (ms) |
|---|---:|---:|---:|---:|---:|
| `test/mmcif/insertion_code/1VVJ/1VVJ.cif` | 20907 | 1536 | 13.61 | 21463 | 2370 |
| `test/pdb/urx053/urx053.pdb` | 341 | 121 | 2.82 | 405 | 187 |
| `test/pdb/tr0001/tr0001.pdb` | 255 | 96 | 2.66 | 341 | 242 |
| `test/mmcif/other/6pom/6pom.cif` | 300 | 115 | 2.61 | 379 | 183 |
| `test/mmcif/nmr_structure/8if5/8if5.cif` | 285 | 113 | 2.52 | 376 | 180 |
| `test/pdb/test1/test1.pdb` | 220 | 90 | 2.44 | 284 | 175 |
| `test/mmcif/insertion_code/1EFW/1EFW.cif` | 362 | 149 | 2.43 | 481 | 320 |
| `test/pdb/url064/url064.pdb` | 228 | 94 | 2.43 | 272 | 193 |
| `test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif` | 213 | 90 | 2.37 | 268 | 147 |
| `test/mmcif/x-ray/434D/assembly-2/434d-assembly2.cif` | 218 | 93 | 2.34 | 260 | 137 |
| `test/pdb/pdb1nvy/pdb1nvy.pdb` | 210 | 92 | 2.28 | 250 | 168 |
| `test/mmcif/x-ray/4NMG/assembly-1/4nmg-assembly1.cif` | 227 | 101 | 2.25 | 281 | 154 |
| `test/mmcif/insertion_code/4ARC/4ARC.cif` | 294 | 132 | 2.23 | 356 | 220 |
| `test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif` | 196 | 92 | 2.13 | 262 | 142 |

最关键的一行是 `1VVJ.cif`：
- legacy 单个任务中位数约 **20.9 秒**
- rust 单个任务中位数约 **1.5 秒**

这类“超级慢 case”在批处理中会直接决定整体吞吐（尤其是 repeat 很大时），所以你看到整体 speedup 能到 8×～15×，主要原因就是长尾被削掉了。

原始数据位置：
- legacy：`out/throughput_sweep_release_w8_20260202_042323/r50/legacy/summary.json`
- rust：`out/throughput_sweep_release_w8_20260202_042323/r50/rust/summary.json`

---

## 5. “启动/初始化/落盘开销”到底大不大？

看 `overhead_ms` 本身：两边基本都在 **0.2s～0.7s** 这个量级（会随 jobs 增加而略增，主要是写更大的 `summary.json`）。

更有意义的是看 `overhead/job_ms`：
- 小任务量（jobs=14）：每个任务要额外承担 ~13–15ms 的固定开销
- 大任务量（jobs=700）：每个任务额外承担 ~1ms 左右

所以结论很直观：
- **跑库/高通量（几百～几千个输入）时，这个开销几乎可以忽略**；
- **Notebook 里只跑 1–2 个输入时**，你会更容易感知到这部分固定成本 —— 这也是为什么你规划的 “Python 可 import 外壳 +（未来）PyO3 直连库调用”会对交互场景更友好。

---

## 6. 产物清单（你后续引用/贴图用）

- Sweep 汇总表：`out/throughput_sweep_release_w8_20260202_042323/summary.md`
- 每个 repeat 的单独对比小结：`out/throughput_sweep_release_w8_20260202_042323/r20/throughput_compare.md`（以及 r1/r2/r5/r10/r50 同名文件）
- 原始 JSON（可二次分析）：`out/throughput_sweep_release_w8_20260202_042323/r20/throughput_compare.json`（每个 repeat 都有）

