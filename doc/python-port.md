# RNAView 现代化升级（Rust + Python）评估与架构草案

这份文档面向“将当前 C 版 RNAView 逐步现代化升级（高性能部分用 Rust，其余用 Python）”的目标：先搞清楚现状与工作量，再给出可落地的分阶段迁移架构。

规格/契约层建模（Legacy vs 新架构）见：`doc/spec.md`。

架构图/流程图素材（Mermaid，可直接画图）见：`doc/architecture-diagrams.md`。

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

### Phase 2：Rust 核心引擎替换（No-C + `.out` byte-exact）（4–12 周）

Phase 2 有两道门槛（建议都写进里程碑）：

- **Gate A（桥接，可选）**：复用 legacy C pipeline/writer，仅把热点替换为 Rust（快速做到 `.out` 字节级一致，并用于性能剖析/定位）。
- **Gate B（No-C，最终验收）**：交付“纯 Rust core + Python 编排”的可发行形态：运行时与构建不依赖任何 C（也不 shell out `bin/rnaview`），但仍保持 `.out` 字节级兼容。

#### Gate A：C pipeline + Rust hot functions（桥接，可选）

- Rust 侧按可回归拆分实现并替换（每一步都用 Phase 1 的 oracle 对照）：
  1. residue/base 识别与编号（对齐 `.out` 的 base index 语义）
  2. `check_pairs()` 几何筛选指标对齐
  3. `Hbond_pair()` 氢键枚举与阈值对齐
  4. `LW_pair_type()` 边类型与 cis/trans 对齐
  5. `all_pairs()` 枚举/去重/tertiary 标注对齐
- 验收：对同一输入与同一组选项，`FILEOUT.out` 必须与 golden/legacy **逐字节一致**（可直接 `diff`）。
  - 基准/剖析建议用 release 构建：`bash tools/build_rnaview_rustcore_release.sh`，并用 `python3 tools/rnaview_bench.py compare --rustcore-bin bin/rnaview_rustcore_release ...` 对比。
  - 当前仓库脚本：`bash test_phase2.sh`（会跑 legacy / rustcore / rust 三条路径，并用 `--regress-mode out` 做 byte-exact 回归）。

#### Gate B：纯 Rust core + Python 编排（No-C，Phase 2 最终验收）

- Rust 实现“结构解析 + core 计算 + `.out` writer”完整闭环，不再依赖 legacy/C：
  - 结构解析：PDB/mmCIF（含 `auth/label`、NMR 选模、altloc、insertion code 等语义）在 Rust/Python 中实现，但不调用 legacy。
  - core 计算：候选对筛选/几何/Hbond/LW 分类/multiplets/stats 全部在 Rust 中完成。
  - 输出：生成 `FILEOUT.out`（byte-exact）并同时输出 `pairs.json`（确定性序列化）用于结构化消费/缓存。
- Python 仍负责编排：批处理、并发、目录布局、回归报告、失败隔离。
- 验收（No-C + byte-exact）：
  - 构建：仅依赖 `cargo build --release` + Python 运行环境（不编译/链接任何 C；不需要 `make/cc/gcc`）。
  - 运行：批处理/单文件运行不调用 `bin/rnaview`/`bin/rnaview_rustcore`（legacy 只允许作为测试 oracle）。
  - 输出：回归集 `.out` 逐字节一致（同 Gate A），并且 `pairs.json` 与 golden core 等价。

补充（现状落地）：
- 当前仓库已提供 Gate B 的“无 C 验收脚本”骨架：`bash test_phase2_noc.sh`。
- 目前 Gate B 的 No-C 跑通依赖一个过渡机制：`tools/rnaview_batch.py --engine rust --rust-oracle out`，即 Rust 侧读取 `<input>.out` 作为 oracle（不再 shell out legacy），用于先把“构建/运行不依赖 C + `.out` byte-exact diff”这条链路固化下来；后续实现纯 Rust 计算后，仍可复用同一套 Gate B 回归脚本作为最终验收。

### Phase 3：工程化与性能（2–6 周）

- Semantics/Policy 显式建模：引入 `--semantics legacy-v1|science-v1` 与可组合 policy（并落盘到 `pairs.json.options`），把“byte-exact 兼容”与“科学修复”解耦。
- CI/可重复构建：把回归（core + `.out`）与基准跑进 CI；提供“一键跑套件 + 产出报告”的脚本。
- API/CLI 稳定化：错误码、日志、schema 版本策略；为跑库提供更清晰的失败定位信息。
- 性能优化：在不改变结果的前提下做空间筛选与并行（并用基准锁住收益）。

### Phase 4：渲染与格式现代化（2–6 周）

- 2D：优先输出 `SVG`（易集成网页/论文），其次 `PDF/PNG`；PS 作为兼容。
- 3D：VRML 可保留，但更现代的是 `glTF` 或直接输出给 PyMOL/ChimeraX 脚本。

## 5. 工作量粗估（以“core/.out 有回归要求”为前提）

如果目标是“现代化 + 保持核心结果一致（base-pair/multiplets/stats 一致，PS/可视化后置）”：

- Phase 1：1–2 周
- Phase 2：4–12 周（含 No-C）
- Phase 3：2–6 周
- Phase 4：2–6 周

合计：大约 2–7 个月（取决于一致性要求、No-C 范围、性能目标、以及是否引入 numpy/scipy/gemmi 等依赖）。

风险点：

- **一致性**：阈值、边界条件、altloc/NMR 模型选择、modified residues 的处理会导致输出差异。
- **`.out` 字节级兼容**：writer 的空格/换行/排序/隐藏分支非常多，任何“看似不重要”的文本变化都会导致 diff。
- **No-C**：legacy 的隐式行为（解析/过滤/编号）需要被 Rust/Python 完整复刻，短期会显著增加 Phase 2 的工作量。
- **性能**：大 RNA（上千残基）若不做邻域筛选，Python 会明显慢于 C。
- **渲染**：PS 逐字节一致很难；建议把“科学内容一致”与“像素/文本完全一致”分开定义。

## 6. 已确认约束 & 下一步

- Phase 0/1 权威输出：`pairs.json`（确定性序列化）+ `.out`（仅 core 段要求一致）；PS/可视化放到下一阶段。
- Phase 2 验收口径：`FILEOUT.out` 逐字节一致 **且** 默认执行路径 No-C（纯 Rust core + Python 编排）。
- 技术分工：高性能热点用 Rust，其余编排/批处理/落盘用 Python。
- 阶段目标：第一阶段以“批处理跑库”为主（大量结构、可并发、可重跑、可汇总）。

建议的下一步落地顺序：

1. 先把 `test/**.out` 作为 golden，写一个“抽取 `.out` core 段并做语义对比”的回归脚本（`tools/rnaview_out_core.py`）。
2. 交付 Python 批处理 CLI：并行跑 `bin/rnaview`，生成 `pairs.json`，并用上面的脚本对照 golden。
3. Phase 2 Gate A：继续用 `.out` 字节级回归锁住热点替换的正确性，并用 profile/bench 找到真实热点。
4. Phase 2 Gate B：定义 “No-C + `.out` byte-exact” 的验收脚本与交付物（纯 Rust core + `.out` writer），逐步替换掉 legacy oracle 与 C pipeline。



Phase 3：工程化 + 科学模式（含去氢 bug 修复）
M3.0 语义/模式接口定稿（最关键的第一步）
接口设计（Rust CLI + Python 编排）
在 Rust from-structure --oracle compute 增加一个顶层模式：--semantics legacy-v1|science-v1
legacy-v1：现有 Gate B 行为（包含 mmCIF 去氢 bug 兼容、链 ID 截断等 legacy 行为）。
science-v1：新“科学模式”（修复去氢 bug 等）。
可选：把关键可变点拆成可组合 policy（也可先不开放给用户，只在 science-v1 内固定）：
--hydrogen-policy discard-all|keep-all|legacy-mmcif-bug
--missing-insertion-code legacy-question-mark|none
--chain-id-policy legacy-1char|unique-1char
--mmcif-id-scheme auth|label
--nmr-model-policy legacy|first|representative
落盘与可追溯
pairs.json 里必须记录 semantics 与关键 policy（保持可审计、可复现）；.out 在 legacy-v1 下绝不能加新字段（否则 byte‑exact 破坏），science-v1 可选择不扩展 .out、只在 pairs.json 记录。
验收
Gate B（bash test_phase2_noc.sh）必须继续全绿，证明 legacy-v1 没被破坏。
M3.1 Gate C（科学模式）定义与落地
Gate C 目标
不再追求与 legacy .out byte‑exact；而是对“科学模式的变化”做透明化、可审核、可批准的管理。
Gate C 执行方式（建议）
对同一批输入同时跑两次：
legacy-v1（基线，用于对照与解释变化）
science-v1（新语义）
产出“差异报告”，并用 allowlist（批准清单）控制哪些差异是“已接受的科学修复”。
Gate C 输出目录与格式（建议）
脚本：`tools/rnaview_gate_c.py`（已实现；默认读取 `test/gate_c_allowlist.yaml`）
快捷验收脚本：`bash test_phase3_gate_c.sh`
输出：
summary.json（Gate C 自己的 schema）
report.md（人读）
pairs.json、engine.out
pairs.json、engine.out
diff.json（结构化差异）
summary.json（示例 schema）：
schema_version
semantics: { "baseline": "legacy-v1", "candidate": "science-v1" }
counts: { ok, changed, unapproved, failed }
results[]: 每个 input 的 paths + diff 摘要（增删改 pair 数、stats 变化、multiplets 变化）
差异报告脚本（怎么做）
输入：两个 pairs.json（直接比 core；spec.md (line 264) 定义了 core 的等价字段）
diff 逻辑（务实且可解释）：
stats_delta: total_pairs/total_bases + type_counts 逐项 diff
base_pairs_delta：
added: 只在 science 出现
removed: 只在 legacy 出现
changed: 同一 (i,j,kind) 存在但 lw/orientation/note/syn/... 变化（输出字段级别变更）
multiplets_delta: added/removed（必要时再细分）
输出：
diff.json（给机器审核/allowlist）
report.md（给人审阅：按“变化最大 case”“变化类型分布”“最常变化的 pair 类型”汇总）
allowlist 机制（让 Gate C 可 CI 化）
文件：`test/gate_c_allowlist.yaml`（已实现；内容为“YAML 兼容的 JSON”，避免引入 PyYAML 依赖）
规则：Gate C 只有在“所有 diff 都在 allowlist”时才通过；否则标记 unapproved 并失败。
allowlist 条目建议以稳定 key 描述：input + (i,j,kind) + field deltas，并附 reason/issue_id.
M3.2 “去氢 bug”在 science-v1 中的修复（你关心的核心）
目标
science-v1：mmCIF/PDB 均执行“正确且一致”的去氢策略（通常是 discard-all hydrogens，优先用元素字段判断）。
legacy-v1：保持现状（继续复刻 legacy bug，用于 byte‑exact）。
备注：当前 `science-v1` preset 默认只改变 `hydrogen_policy=discard-all`，其余结构解析 policy 仍保持 legacy-v1 默认值；如需额外“科学修复”，通过显式 policy flag 覆盖并用 Gate C 管控差异。
验收
Gate B 继续全绿（legacy-v1）。
Gate C 中预期只出现“与去氢相关的差异”，并全部进入 allowlist；后续任何新差异都必须解释/批准。
M3.3 测试矩阵扩充（把 Phase3 的风险覆盖住）
围绕 spec.md (line 264) 的 core + 结构解析风险点，把“回归集”明确成一个矩阵（每一类至少 1–2 个 fixture），并把缺口显式标出：

当前回归集（核心输入清单 = `test/golden_core/manifest.json`）：

- 格式维度
  - PDB：`test/pdb/tr0001/tr0001.pdb`、`test/pdb/test/test.pdb`
  - mmCIF(auth)：`test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif`、`test/mmcif/insertion_code/1EFW/1EFW.cif`
  - mmCIF(label)：TODO（需要先实现/开放 id_scheme 选择，再补 fixture + Gate）
- 模型维度
  - X-ray 单模型：`test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif`、`test/pdb/pdb1nvy/pdb1nvy.pdb`
  - NMR 多模型：`test/mmcif/nmr_structure/8if5/8if5.cif`
  - representative 字段：TODO（需要先实现 model policy；补 1 个“有 rep”+ 1 个“无 rep”的 NMR）
- 编号维度（insertion code）
  - mmCIF `.`/缺失/`?`/带字母：`test/mmcif/insertion_code/1EFW/1EFW.cif`、`test/mmcif/insertion_code/4ARC/4ARC.cif`
  - PDB insertion code：TODO（补一个含 icode 的 PDB fixture；并加入 manifest + Gate B）
- 链维度
  - 单字符 chain：大部分 PDB fixture（如 `test/pdb/urx053/urx053.pdb`）
  - 多链/大复合体（压力）：`test/mmcif/insertion_code/1VVJ/1VVJ.cif`
  - 多字符 label_asym_id + chain-id-policy：TODO（补一个“前缀冲突”的 mmCIF（如 AA/AB）来测试 `--chain-id-policy unique-1char`）
- 氢维度
  - 触发 legacy mmCIF 去氢 bug（4 字符氢名）：`test/mmcif/nmr_structure/8if5/8if5.cif`（Gate C 的 allowlist 当前只批准这类差异）
  - 元素字段缺失（name fallback）：通过 Rust 单元测试覆盖（见 `rust/src/structure.rs` tests）
- 化学维度（修饰碱基、PSU、I 等）
  - 修饰碱基较多：`test/pdb/tr0001/tr0001.pdb`、`test/mmcif/insertion_code/1EFW/1EFW.cif`
- 几何维度（tertiary / multiplets / stacked）
  - stacked：`test/pdb/test/test.pdb`、`test/pdb/urx053/urx053.pdb`
  - multiplets：`test/mmcif/insertion_code/1VVJ/1VVJ.cif`、`test/pdb/urx053/urx053.pdb`
  - tertiary：TODO（若需要显式覆盖，先定义“tertiary”的可检测信号/字段，再选 fixture）

验收脚本建议（把矩阵落到 Gate 上）：

- legacy-v1：继续走 `.out` byte‑exact（Gate B：`bash test_phase2_noc.sh`）
- science-v1：走 Gate C（差异报告 + allowlist）：`bash test_phase3_gate_c.sh`（当前跑 `test/pdb` + `test/mmcif`，并强制“除去氢之外无新差异”）
验收分配
legacy-v1：继续走 .out byte‑exact（Gate B）
science-v1：走 Gate C（差异报告 + allowlist），后续稳定后可再冻结 science goldens（见下一条）
M3.4 冻结 “science-v1” 的 golden（让科学模式也可回归）
当 Gate C 的差异已经被充分解释并批准后：
固化一套 science-v1 的 golden（建议先固化 pairs.json core，再决定是否也固化 .out）：
- 冻结：`python3 tools/rnaview_science_golden.py freeze`（输出：`test/golden_science_core/manifest.json`）
  - 默认策略：只为“与 legacy 不同的 case”写入 `test/golden_science_core/**.core.json`，其余 case 直接复用 `test/golden_core/**.core.json`（节省重复数据，也让差异更聚焦）。
- 回归：`bash test_phase3_science.sh`（跑 `science-v1` 并对 `test/golden_science_core/manifest.json` 做 core 回归）
这一步完成后，science-v1 也从“解释差异”进入“稳定回归”。
M3.5 工程化与性能（补齐 python-port.md (line 176) 的内容）
CI：分层跑（建议拆成独立 job，便于定位失败与并行）
- GitHub Actions 示例：`.github/workflows/ci.yml`（Gate B/C/science）+ `.github/workflows/bench.yml`（bench，scheduled/dispatch）
- Phase 3 收尾一键脚本：`bash test_phase3_wrapup.sh`（可附加任意 `.pdb/.cif` 输入并生成 Gate C 报告）
- Gate B（byte‑exact No‑C）：`bash test_phase2_noc.sh`
- Gate C（science diff + allowlist）：`bash test_phase3_gate_c.sh`
- science-v1 回归（frozen goldens）：`bash test_phase3_science.sh`
- bench（只监控性能，不做硬门槛）：`python3 tools/rnaview_bench.py compare --suite phase2 --runs 3 -o out/bench_phase3.json`
性能：只在不改变 legacy-v1 结果前提下优化（空间索引、邻域筛选、并行），并用基准锁住收益（只报警，不 hard-fail）。
Phase 4：渲染与格式现代化（在不动 core 的前提下扩展产物）
M4.0 输出契约（先把验收口径写死）
目标：在 CI/Linux 容器里，以 legacy `rnaview -p/-v` 为 baseline，做到：
- 允许中间步骤不同，但最终产物（canonical/normalize 后）必须 `diff 0`
- golden 存 canonical（normalize 后的输出；建议压缩存储）
- 渲染器输入只读 `pairs.json`，并允许按 `pairs.json.source.path` 再读一次结构文件（最小化 schema 变动）
- 固定参数集必须显式化并落盘（`--label/auth`、链选择、NMR model、altloc、分辨率过滤等），避免“输入不唯一”造成的伪 diff

M4.1 2D：RNAML/XML + PS + SVG
- RNAML/XML、PS：对齐 legacy `-p`（normalize 后 diff 0）
- PS 样式完全复刻 legacy（字体/线宽/颜色/线型等；基线见 `BASEPARS/ps_image.par`）
- SVG：以 legacy `*.ps` 为“渲染权威”，用确定性的 PS→SVG 转换器生成（保证与 legacy 最终图一致；RNAML/XML→SVG 可作为可选语义渲染/调试）

M4.2 3D：VRML + glTF
- VRML：对齐 legacy `-v`（normalize 后 diff 0）
- glTF：与 VRML 语义等价；golden 推荐由 legacy VRML 通过确定性转换器生成

M4.3 Gate D（渲染验收 + allowlist）
- 目标：canonical 输出 byte‑exact（diff 0）
- allowlist：类似 Gate C，事件级稳定 id（包含 input_id/format/semantics/renderer_version 等）
- allowlist 文件：`test/gate_d_allowlist.yaml`
哪些必须 byte‑exact？哪些可以科学修复？
必须保持 byte‑exact（对 legacy golden）
legacy-v1 下的 FILEOUT.out 全文（Gate B），包含解析/过滤/排序/空格/换行等所有历史行为。
可以“科学修复”（但必须在 science-v1，并受 Gate C 管控）
mmCIF 去氢 bug（你指出的点）
mmCIF chain id 截断策略（见 spec.md (line 138)、spec.md (line 275)）
insertion code 的表示、NMR 选模策略、altloc 策略等（只要是“历史实现细节”而非科学定义）
两边都必须保持的底线
同一 semantics + 同一输入 + 同一组选项 ⇒ 输出确定性（pairs.json 必须可 byte‑diff）；只是 science-v1 的 golden 不再是 legacy，而是它自己的版本化 golden/allowlist。
