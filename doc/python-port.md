# RNAView 现代化升级（Rust + Python）评估与架构草案

这份文档面向“将当前 C 版 RNAView 逐步现代化升级（高性能部分用 Rust，其余用 Python）”的目标：先搞清楚现状与工作量，再给出可落地的分阶段迁移架构。

规格/契约层建模（Legacy vs 新架构）见：`doc/spec.md`。

## 1. 现有项目快速审核（你现在这个仓库）

### 1.1 代码规模与构成

- 代码量：约 `27k` 行 C（`src/*.c` + `include/*` + `rnaml2ps/*`）。
- 入口程序：`src/rnaview.c`（`main()` + `rna()` + `work_horse()`），构建后生成 `bin/rnaview`。
- 数据文件：`BASEPARS/*`（标准碱基几何、颜色/渲染参数等），运行时通过环境变量 `RNAVIEW` 或当前目录定位。

### 1.2 运行工作流（从入口看真实行为）

`src/rnaview.c` 的高层流水线基本是：

1. 解析 CLI 参数（链选择、VRML、批量、解析 mmCIF 的 label/auth 等）。
2. 判断输入是 RNAML(XML) 还是 PDB/mmCIF（通过扫描文件内容，而不是严格依赖 `-x`）。
3. 对 PDB 做“清洗/过滤”（链过滤、去 ANISOU、NMR 模型、分辨率过滤、alt conformer 处理等），对 mmCIF 则走自带解析器。
4. `rna()`：把输入解析成“扁平数组”数据结构（`AtomName/ResName/ChainID/ResSeq/Miscs/xyz` 等），构建 `seidx`（residue -> atom index range）。
5. `work_horse()`：核心计算与输出编排：
   - `base_info()`：生成每个残基的参考系/原点/关键原子坐标等（用于后续几何判断）。
   - `all_pairs()`：枚举候选残基对，进行几何筛选 + 氢键识别 + Leontis-Westhof 分类，输出 `.out` 主结果并累积统计。
   - `re_ordering()/write_best_pairs()`：做“最佳配对/螺旋段”整理。
   - 若启用 2D：`process_2d_fig()`（布局）→ `write_xml()` 生成 RNAML → `xml2ps()` 读 RNAML 画 PS。
   - 若启用分析：`bp_analyze()` 输出 `analyze.out` 等。
   - 若启用 3D：`process_3d_fig()` 输出 VRML（`.wrl`）。
   - 若启用 motif：`motif()` 基于 `.out` 再输出 pattern 文件。

备注：本仓库已对 legacy 做了一个小 patch：默认不生成 2D/XML/PS（`PS=0`），通过 `-p` 显式开启。这是为了 Phase0/Phase1 先聚焦 “core 一致性回归”，可视化先后置；后续若要兼容上游 legacy 行为，可以在新 CLI 层提供 `--2d/--no-2d` 之类的显式选项，并保留 legacy 兼容模式。

### 1.3 重要输出（建议作为 Python 版的兼容目标）

- `FILEOUT.out`：核心“完全注释的 base-pair 列表”（很多后续模块都以它为输入）。
- `FILEOUT.xml`：RNAML（2D/PS 的中间格式）。
- `FILEOUT.ps`：2D 结构（PostScript）。
- `analyze.out`：结构参数分析（如果启用）。
- `base_pair_statistics.out`：全局统计（`main()` 里总是写）。
- 其他临时/辅助：`*_patt_tmp.out`、`*_tmp.pdb`、`best_pair.out`、`pattern_tmp.out` 等。

历史上的 `test.sh` 以 `.ps` 做字节级 `diff`，容易被 `%%CreationDate` 等非确定字段干扰；在本仓库的迁移/回归阶段：

- Phase 0/1：推荐用 `.out` 的 core 语义回归（`tools/rnaview_out_core.py` / `tools/rnaview_batch.py --regress --regress-mode core`）。
- Phase 2：验收口径提升为 `FILEOUT.out` 逐字节一致（`tools/rnaview_batch.py --regress --regress-mode out`），用它作为 Rust 热点替换的硬门槛。

## 2. C 模块职责 → Rust/Python 组件映射

建议用“职责”而不是“逐文件翻译”来规划迁移，下面是最关键的映射：

- `src/rnaview.c`
  - Python 对应：`cli` + `pipeline` 编排层（参数解析、输入判断、输出路径策略、调用核心分析/渲染）。
- `include/cifparse.c`
  - 对应：优先在 Python 侧调用成熟库；或在 Rust 侧实现/复用解析（后置）。短期要保持一致，可先用 subprocess 调旧二进制做 oracle。
- `src/fpair*.c` + `src/pair_type.c`
  - Rust 对应：`core/pairing/*`（几何筛选、氢键识别、LW 分类、Saenger 对应等）。
  - 这是“最难/最值钱”的部分：决定结果一致性与性能，且属于热点，优先 Rust。
- `src/ps-xy*.c` + `src/xml2ps.c` + `src/rnaxml-new.c`
  - Python 对应：`render/`（2D 布局 + RNAML 写入 + 图形输出）。
  - 现代化建议：以 `SVG/PDF/PNG` 为主，PS 作为兼容输出；RNAML 作为可选中间格式或 legacy 兼容。
- `src/analyze.c`
  - Python 对应：`analysis/helical_params.py`（参数分析输出）。
- `src/pattern.c`
  - Python 对应：`analysis/motifs.py`（基于 base-pair 列表的 motif 搜索）。
- `src/statistics.c`
  - Python 对应：`analysis/stats.py`（统计汇总，最好输出成 JSON/CSV，同时保留 legacy 文本）。
- `BASEPARS/*`
  - Python 对应：`package_data/`（作为包内资源，避免依赖环境变量；同时提供 `--basepars-dir` 覆盖）。

## 3. Python 版建议架构（可落地）

### 3.1 分层

1. **I/O 层**：PDB/mmCIF/RNAML 解析 → 标准化结构对象
2. **Core 层**：碱基识别、参考系、几何筛选、氢键与分类 → `BasePair` 等结果对象
3. **Analysis 层**：螺旋段整理、统计、motif、参数分析
4. **Render 层**：2D/3D/格式导出（SVG/JSON/PS 兼容）
5. **CLI 层**：向后兼容的命令行 + 更现代的子命令/参数

### 3.2 数据模型（建议）

- `Atom`: `name, element, x,y,z, altloc, occupancy, bfactor, charge`
- `ResidueId`: `chain_id, resseq, icode, model`
- `Residue`: `id, resname, atoms, is_nucleic, base_letter, modified_from?`
- `Structure`: `models -> chains -> residues`（或简化为单模型）
- `BasePair`:
  - `i, j`（内部索引，1-based/0-based 可选，但要有映射）
  - `residue_i, residue_j`（原始 ResidueId）
  - `lw: (edge_i, edge_j, cis_trans)`、`saenger`、`is_wc`
  - `hbonds: list[Hbond(atom_i, atom_j, dist)]`
  - `metrics: distances/angles`（与当前 `rtn_val[]` 对齐方便回归）
- `AnalysisResult`:
  - `pairs, helices, loops, isolated_bases, sugar_syn, stats`
  - `artifacts`（输出文件路径、日志信息）

### 3.3 API（示例形态）

- 库接口（给下游脚本/Notebook 用）：
  - `result = rnaview.analyze(path, format="pdb|cif|auto", chains="ABC", cif_ids="auth|label", ...)`
  - `result.pairs` 返回结构化对象，可直接 `to_json()` / `to_dataframe()`
- CLI（兼容 + 现代）：
  - `rnaview pairs <file> [--chains ABC] [--cif-ids auth|label] [--json out.json]`
  - `rnaview draw <file> --svg out.svg`（或 `--ps`）
  - `rnaview convert <file.xml> --ps out.ps`
  - `rnaview validate <file> --golden-dir test/...`（用于回归）

### 3.4 性能策略（Python 版必须提前考虑）

现有算法里大量是 `O(N^2)` 的残基对枚举（`all_pairs()` 里双重循环）。直接 Python 实现会明显慢；策略是：

- 热点（候选对筛选/几何/Hbond/LW 分类）用 Rust 实现并暴露给 Python 调用。
- 先把正确性对齐，再做空间筛选（网格/KD-tree）与并行（线程/向量化）优化。

## 4. 推荐迁移路线（最稳的分阶段）

### Phase 0：定目标与回归基线（1–3 天）

- 明确“权威输出”和“一致性口径”。已确认：以 base-pair 结果为权威，PS/渲染一致性后置。
- 选定 5–10 个代表性结构作为“核心回归集”（PDB + mmCIF + insertion code + NMR）。
- 定义 **核心一致性契约**（用于所有后续验收）：
  - `.out`：只要求 `BEGIN_base-pair … END_base-pair`、`BEGIN_multiplets … END_multiplets`、以及 `The total base pairs = ...` 后的统计表在“科学意义”上相同（允许路径/提示/日志变化）。
  - `pairs.json`：结构化权威产物，要求确定性序列化（可做字节级 diff）。
- 增加一个回归工具：把 `.out` 解析成“规范化的 core 结果”，用于对比（见 `tools/rnaview_out_core.py`）。
  - 抽取：`python3 tools/rnaview_out_core.py extract golden.out > golden.core.json`
  - 对比：`python3 tools/rnaview_out_core.py compare golden.out candidate.out`
  - 批量冻结：`python3 tools/rnaview_out_core.py freeze test`（生成 `test/golden_core/`）
  - `.out` writer（core-only）：`python3 tools/rnaview_pairs_json.py write-out pairs.json > candidate.out`
  - 验证 writer：`python3 tools/rnaview_pairs_json.py validate-golden`

### Phase 1：Python 批处理包装层（1–2 周，最快见效）

目标：面向“跑库/批处理”，不动/少动算法，先把可用性、产物格式和回归体系搭起来。

- 用 Python 做批处理 CLI/库接口，底层调用现有 `bin/rnaview`（subprocess）作为 oracle。
- 解析 `.out` 的 core 段，生成 `pairs.json`（schema v1），并保留 `.out`（用于兼容）。
- 做好并发（按文件并行）、失败隔离、日志、输出目录布局、可重跑（幂等）。

优点：风险最低、很快可交付；可以并行推进 Phase 2/3。

### Phase 2：Rust 核心引擎逐块替换（2–8 周）

- Rust 侧按可回归拆分实现并替换（每一步都用 Phase 1 的 oracle 对照）：
  1. residue/base 识别与编号（对齐 `.out` 的 base index 语义）
  2. `check_pairs()` 几何筛选指标对齐
  3. `Hbond_pair()` 氢键枚举与阈值对齐
  4. `LW_pair_type()` 边类型与 cis/trans 对齐
  5. `all_pairs()` 枚举/去重/tertiary 标注对齐
- Python 仍负责 I/O 编排与输出落盘；Rust 只负责热点计算。
- Phase 2 的验收标准：对同一输入与同一组选项，`FILEOUT.out` 必须与 legacy 逐字节一致（可直接 `diff`）。
  - 基准/剖析建议用 release 构建：`bash tools/build_rnaview_rustcore_release.sh`，并用 `python3 tools/rnaview_bench.py compare --rustcore-bin bin/rnaview_rustcore_release ...` 对比。

### Phase 3：I/O 现代化与可维护性（2–6 周）

- 视需求把 PDB/mmCIF 解析迁移到 Python 或 Rust（以“不破坏 core 一致性”为前提）。
- 引入更清晰的数据模型与错误处理（便于跑库与定位问题）。

### Phase 4：渲染与格式现代化（2–6 周）

- 2D：优先输出 `SVG`（易集成网页/论文），其次 `PDF/PNG`；PS 作为兼容。
- 3D：VRML 可保留，但更现代的是 `glTF` 或直接输出给 PyMOL/ChimeraX 脚本。

## 5. 工作量粗估（以“core/.out 有回归要求”为前提）

如果目标是“现代化 + 保持核心结果一致（base-pair/multiplets/stats 一致，PS/可视化后置）”：

- Phase 1：1–2 周
- Phase 2：2–8 周
- Phase 3：2–6 周
- Phase 4：2–6 周

合计：大约 2–5 个月（取决于一致性要求、性能目标、以及是否引入 numpy/scipy/gemmi 等依赖）。

风险点：

- **一致性**：阈值、边界条件、altloc/NMR 模型选择、modified residues 的处理会导致输出差异。
- **性能**：大 RNA（上千残基）若不做邻域筛选，Python 会明显慢于 C。
- **渲染**：PS 逐字节一致很难；建议把“科学内容一致”与“像素/文本完全一致”分开定义。

## 6. 已确认约束 & 下一步

- Phase 0/1 权威输出：`pairs.json`（确定性序列化）+ `.out`（仅 core 段要求一致）；PS/可视化放到下一阶段。
- Phase 2 验收口径：`FILEOUT.out` 逐字节一致（作为 Rust 核心替换的回归门槛）。
- 技术分工：高性能热点用 Rust，其余编排/批处理/落盘用 Python。
- 阶段目标：第一阶段以“批处理跑库”为主（大量结构、可并发、可重跑、可汇总）。

建议的下一步落地顺序：

1. 先把 `test/**.out` 作为 golden，写一个“抽取 `.out` core 段并做语义对比”的回归脚本（`tools/rnaview_out_core.py`）。
2. 交付 Python 批处理 CLI：并行跑 `bin/rnaview`，生成 `pairs.json`，并用上面的脚本对照 golden。
3. 再开始 Rust 核心引擎替换：每替换一块，就把回归集跑通并对齐 core 结果。
