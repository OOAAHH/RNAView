# RNAView 现代化升级：架构图/流程图素材（Legacy C → Rust + Python）

> 这份文档的目标是“给你画图用”：把 **重构前/后** 的系统边界、模块划分、数据流、回归/验收与未来封装计划，用**可直接复制的 Mermaid** 和“模块关系说明”写清楚。  
> 规格/契约（输出一致性口径、gate）请以 `doc/spec.md` 为准；迁移阶段与工作量评估见 `doc/python-port.md`。

---

## 1. 术语与产物（图里建议统一用这些名字）

**输入**
- `Structure`: `*.pdb/*.ent/*.cif`
- `RNAML`: `*.xml`（legacy 2D 渲染中间格式）

**关键输出（Phase 0–3 的“科学 core”）**
- `FILEOUT.out`：legacy 传统文本主结果（包含 `BEGIN_base-pair` / `BEGIN_multiplets` / 统计表）
- `pairs.json`：结构化权威输出（schema v1；用于 diff/缓存/未来 API）

**渲染输出（Phase 4）**
- 2D：`*.xml`（RNAML）、`*.ps`（PostScript）、`*.svg`（新增）
- 3D：`*.wrl`（VRML）、`*.gltf/.glb`（新增）

---

## 2. 重构前：Legacy C（单体程序）系统上下文图

```mermaid
flowchart LR
  user[User / Pipeline] -->|CLI: bin/rnaview| rnaview[Legacy RNAVIEW (C)]
  rnaview -->|Structure input| in1[PDB / mmCIF]
  rnaview --> out1[FILEOUT.out]
  rnaview -->|optional -p| xml[FILEOUT.xml (RNAML)]
  rnaview -->|optional -p| ps[FILEOUT.ps]
  rnaview -->|optional -v| wrl[FILEOUT.wrl / *_new.wrl]
  rnaview --> aux[other outputs: analyze.out / pattern*.out / statistics]
  rnaview --> basepars[BASEPARS/* (runtime data)\nvia $RNAVIEW]
```

**代码锚点（画图时可做注释）**
- 入口与编排：`src/rnaview.c`
- mmCIF 解析：`include/cifparse.c`
- pairing 核心：`src/fpair.c`、`src/fpair_sub.c`、`src/pair_type.c`
- 2D：`src/rnaxml-new.c`（写 RNAML）→ `src/xml2ps.c`（画 PS）
- 3D：`src/vrml.c`

---

## 3. 重构前：Legacy C 内部流水线（计算 + 渲染）

> 这张图建议你作为“流程图”放在重构前一页：它解释 legacy 为什么难改（计算/输出/渲染高度耦合）。

```mermaid
flowchart TD
  cli[argv/flags] --> detect[detect input type\n(PDB/mmCIF vs RNAML)\nsrc/rnaview.c]

  detect -->|RNAML| xml2ps_only[xml2ps (RNAML→PS)\nsrc/xml2ps.c] --> ps_only[*.ps]

  detect -->|Structure| clean[clean/filter input\nsrc/rnaview.c:clean_inpfile]
  clean --> parse[parse atoms/residues\nPDB: src/rnaview.c\nmmCIF: include/cifparse.c]
  parse --> idx[build residue index / base seq\n(seidx/bseq/...) src/rnaview.c]

  idx --> base_info[base_info / geometry prep\nsrc/fpair_sub.c]
  base_info --> all_pairs[all_pairs: enumerate i<j\nH-bond + LW classify\nsrc/fpair.c]
  all_pairs --> out[FILEOUT.out\n(base-pair + multiplets + stats)]

  out -->|if -p| rnaml[write RNAML\nsrc/rnaxml-new.c] --> xml[*.xml]
  xml --> xml2ps[RNAML→PS\nsrc/xml2ps.c] --> ps[*.ps]

  out -->|if -v| vrml[VRML render\nsrc/vrml.c] --> wrl[*.wrl]

  out -->|optional| analyze[analyze\nsrc/analyze.c] --> analyze_out[analyze.out]
  out -->|optional| motif[pattern/motif\nsrc/pattern.c] --> patt_out[pattern*.out]
```

---

## 4. 重构后：目标架构（Rust + Python 分层）系统上下文图

> 这张图是“重构后”的首页：强调 **Python 编排** 与 **Rust 热点/核心** 的边界，以及 legacy 仅作为 oracle 的角色（最终 No‑C）。

```mermaid
flowchart LR
  user[User / Pipeline / Notebook] --> py[Python entrypoints\n(batch/render/regress; future: importable API)]
  py -->|PDB/mmCIF| rust[Rust core engine\n(rnaview-hotcore)]
  rust --> pairs[pairs.json (deterministic)]
  rust --> out[FILEOUT.out (byte-exact target)]

  py --> render[Render pipeline\n(PS/XML/WRL + SVG/glTF)]
  render --> svg[*.svg]
  render --> gltf[*.gltf/.glb]

  legacy[Legacy bin/rnaview]:::legacy -->|oracle for regress only| py

  classDef legacy fill:#eee,stroke:#999,color:#333;
```

**实现落点（当前仓库）**
- Python 编排/回归：`tools/rnaview_batch.py`、`tools/rnaview_gate_c.py`、`tools/rnaview_gate_d.py`
- Rust core（No‑C 计算 + `.out(full)` writer）：`rust/src/noc_engine.rs`、`rust/src/legacy_pairing.rs`、`rust/src/out_full.rs`
- Rust CLI（面向 Python 调用/未来 API）：`rust/src/main.rs`（`rnaview-hotcore`）

---

## 5. 重构后：运行路径（Phase 0–3）= “引擎”视角流程图

> 这张图适合解释“为什么仓库里有三种 engine”：legacy / rustcore(桥接) / rust(No‑C)。

```mermaid
flowchart TB
  input[PDB/mmCIF] --> choose{Engine / Gate}

  choose -->|legacy| legacy_bin[bin/rnaview\n(C monolith)]
  legacy_bin --> legacy_out[legacy .out]
  legacy_out --> core_extract[extract core\ntools/rnaview_out_core.py]
  core_extract --> pairs[pairs.json\ntools/rnaview_pairs_json.py]

  choose -->|Gate A bridge| rustcore_bin[bin/rnaview_rustcore(_release)\n(C pipeline + Rust staticlib)]
  rustcore_bin --> out_exact[.out (byte-exact gate)]

  choose -->|No-C (target)| hotcore[rnaview-hotcore from-structure\n--oracle compute]
  hotcore --> pairs2[pairs.json]
  hotcore --> out2[.out(full)]

  out_exact --> diff_out[diff .out bytes\nPhase2 Gate A/B]
  out2 --> diff_out
  legacy_out --> diff_out
```

**桥接（Gate A）的实现含义**
- `tools/build_rnaview_rustcore*.sh`：把 legacy C 源码重新编译链接，但通过 `-DRNAVIEW_RUST_*` 宏把热点函数替换为 Rust 实现（`rust/src/legacy_ffi.rs` + 对应 Rust 实现）。
- 好处：先锁定 `.out` byte‑exact，再逐步把 C 依赖清零。

### 5.1 新版本算法图示（v2.0 当前实现：No‑C compute + 2D/3D render）

> 这张图更偏“算法/实现调用链”，用于解释 **新版在做什么**，以及 compute/render 两条路径如何共用结构解析与 `pairs.json` 契约。

```mermaid
flowchart TB
  structure[Structure input\nPDB/mmCIF]:::io

  out_full[FILEOUT.out (full)\n(byte-exact target)]:::artifact
  pairs[pairs.json (schema v1)\n(deterministic)]:::artifact

  subgraph compute["No-C compute (Rust)"]
    cfg[SemanticsConfig + StructurePolicies\n(--semantics + policy overrides)]
    arrays[build_legacy_arrays\n(1-based legacy arrays)]
    base_info[compute_base_info\n(base frames / geometry prep)]
    pair_enum[all_pairs\n(enumerate i<j + classify)]
    multiplets[bp_network_multiplets]
    stats[pair_type_statistics]
    ir[OutFull IR]
    writer[write_out_full]
    core_extract[extract_core_from_out_*]
  end

  structure --> cfg --> arrays --> base_info --> pair_enum
  pair_enum --> multiplets --> ir
  pair_enum --> stats --> ir
  ir --> writer --> out_full
  out_full --> core_extract --> pairs

  subgraph render2d["2D render (Rust)"]
    arrays2[build_legacy_arrays\n(re-parse structure)]
    syn[syn_or_anti]
    layout[compute_layout_2d\n(best-pair → helix → XY)]
    rnaml[write_rnaml_xml\n(RNAML XML)]
    ps[ps_from_rnaml_xml\n(PostScript PS)]
  end

  structure --> arrays2 --> syn
  arrays2 --> layout
  pairs --> layout
  syn --> rnaml
  layout --> rnaml --> ps

  subgraph render3d["3D render (Rust)"]
    wrl_fn[render_vrml_from_pairs_json]
    wrl[VRML .wrl]
  end

  structure --> wrl_fn
  pairs --> wrl_fn --> wrl

  ps --> svg[PS→SVG converter\n(tools/rnaview_ps_svg.py)]
  wrl --> gltf[VRML→glTF converter\n(tools/rnaview_vrml_gltf.py)]

  classDef io fill:#eef,stroke:#99f,color:#000;
  classDef artifact fill:#efe,stroke:#9f9,color:#000;
```

**代码锚点（和上图 1:1 对应）**
- 结构解析 + policy：`rust/src/structure.rs`、`rust/src/semantics.rs`
- No‑C compute 主流程：`rust/src/noc_engine.rs`
- base frames / pairing：`rust/src/legacy_alg.rs`、`rust/src/legacy_pairing.rs`
- `.out(full)` IR + writer：`rust/src/out_full.rs`
- `pairs.json` schema：`rust/src/lib.rs`（`PairsJson`/`Core`）
- 2D：`rust/src/render_2d.rs` → `rust/src/legacy_2d_layout.rs` → `rust/src/legacy_rnaml.rs` → `rust/src/legacy_xml2ps.rs`
- 3D：`rust/src/vrml_render.rs`
- 派生格式 converter：`tools/rnaview_ps_svg.py`、`tools/rnaview_vrml_gltf.py`

---

## 6. 重构后：代码模块图（“画架构图”用的组件清单）

> 你可以把这一节直接映射到 draw.io / PPT 的“盒子 + 箭头”。

### 6.1 Python 层（编排/验收/渲染）

- **Batch & regression**：`tools/rnaview_batch.py`
  - 输入收集（dir/glob/list）
  - 并发调度（`--workers`）
  - engine 选择：`legacy` vs `rust`
  - oracle/语义开关：`--rust-oracle legacy|out|compute`，`--semantics legacy-v1|science-v1`，以及结构 policy（H/链/插入码）
  - 产物落盘：`<out_dir>/<job_id>/pairs.json`、`legacy.out`/`engine.out`、`summary.json`
- **Core diff / golden**：`tools/rnaview_out_core.py`、`test/golden_core/*`
- **语义差异 gate（allowlist）**：`tools/rnaview_gate_c.py`、`test/gate_c_allowlist.yaml`
- **渲染 gate（golden + canonicalize + allowlist）**：`tools/rnaview_gate_d.py`、`test/golden_render/*`、`test/gate_d_allowlist.yaml`
- **渲染入口（可替换 renderer）**：`tools/rnaview_render.py`
  - `--backend legacy|rustcore|rustcore-release|pairs-out|pairs-json|pairs-out-noc3d`
  - `pairs-out*` 会先用 `rnaview-hotcore` 生成 `engine.out`/`pairs.json`，再通过 `RNAVIEW_OUT_PATH` 注入到 renderer（当前默认是 `bin/rnaview_rustcore_release`）
- **格式转换器（派生产物）**
  - `tools/rnaview_ps_svg.py`：PS → SVG（deterministic/canonicalize）
  - `tools/rnaview_rnaml_svg.py`：RNAML → SVG（可作为替代路径/对照）
  - `tools/rnaview_vrml_gltf.py`：VRML → glTF

### 6.2 Rust 层（核心引擎 + writer + 渲染）

- **CLI / 可调用入口**：`rust/src/main.rs`（`rnaview-hotcore`）
  - `from-structure`：Structure →（oracle/compute）→ `pairs.json` + `FILEOUT.out`
  - `from-out`：`.out` → `pairs.json`
  - `write-out`：`pairs.json` → `.out(core)`
  - `render-wrl` / `render-2d`：`pairs.json` + source → `*.wrl` / `*.xml/*.ps`
- **结构解析与“legacy 语义”对齐**：`rust/src/structure.rs`
  - hydrogen 处理策略：`HydrogenPolicy`
  - mmCIF 插入码策略：`MissingInsertionCodePolicy`
  - chain id 映射策略：`ChainIdPolicy`
- **No‑C 计算主流程**：`rust/src/noc_engine.rs`
  - `parse_structure_bases_with_atoms_with_policies` → `build_legacy_arrays`
  - `legacy_alg::compute_base_info`
  - `legacy_pairing::all_pairs` / `bp_network_multiplets` / `pair_type_statistics`
  - `out_full::OutFull`（writer 输入 IR）
- **`.out(full)` IR 与 byte‑exact writer**：`rust/src/out_full.rs`
- **渲染**：`rust/src/vrml_render.rs`、`rust/src/legacy_2d_layout.rs`、`rust/src/legacy_rnaml.rs`、`rust/src/legacy_xml2ps.rs`
- **语义/政策配置**：`rust/src/semantics.rs`

---

## 7. 重构后：数据契约与回归（“输出为什么可验收”）

```mermaid
flowchart LR
  out[FILEOUT.out] --> core[core extract\n(base-pair/multiplets/stats)]
  core --> pairs[pairs.json (schema v1)]
  pairs --> out_core[.out(core)\n(pairs.json→out writer)]

  legacy_out[legacy .out] --> core
  rust_out[rust .out(full)] --> core

  subgraph Gates[Regression Gates]
    g0[Phase0/1: core equivalence]
    g2[Phase2: .out byte-exact]
    g3[Phase3: science diff + allowlist]
    g4[Phase4: render goldens + canonicalize]
  end

  core --> g0
  out --> g2
  core --> g3
  render[ps/xml/wrl/svg/gltf] --> g4
```

**对应脚本（方便你在图旁边标注）**
- Phase 2 Gate A：`bash test_phase2.sh`
- Phase 2 Gate B（No‑C）：`bash test_phase2_noc.sh`
- Phase 3 Gate C：`bash test_phase3_gate_c.sh`
- Phase 4 Gate D：`bash test_phase4_gate_d.sh`

---

## 8. 未来：Notebook 可 import 的 Python 外壳（Conda 发布）建议架构

> 这部分是“计划图”：你可以把它放到 roadmap 一页，说明最终交付形态和升级路径。

```mermaid
flowchart LR
  nb[Jupyter / Python user] --> api[Python package: rnaview\n(high-level API)]
  api --> opts[Options/Semantics/Policies]
  api --> cache[Cache/Artifacts\n(pairs.json, .out, svg, gltf)]

  api --> sel{Engine backend}
  sel -->|preferred| pyo3[Native extension\n(PyO3/maturin)\nRust core as library]
  sel -->|fallback| cli[subprocess\nrnaview-hotcore]
  sel -->|debug/oracle| legacy[subprocess\nbin/rnaview]

  pyo3 --> rustcore[Rust core compute\n(structure→pairs/out)]
  cli --> rustcore

  conda[Conda package/channel] --> api
  conda --> pyo3
```

**实现策略（建议）**
- API 先以“调用现有 CLI + 读写 pairs.json”为主，稳定后再引入 PyO3（可平滑替换 backend）。
- `pairs.json` 保持“确定性序列化 + schema version”，作为 notebook/批处理/缓存共享的中心契约。

---

## 9. 高通量性能测试（模拟“跑库”场景）的建议

> 目标不是做硬门槛，而是提供可重复的 **吞吐量/资源占用** 报告，便于你后续优化 Rust 热点与并发策略。

### 9.1 建议的测试模型（可画成一张图）

```mermaid
flowchart LR
  suite[Input suite\n(many PDB/mmCIF)] --> queue[Job queue]
  queue --> w1[Worker 1\n(engine)]
  queue --> w2[Worker 2\n(engine)]
  queue --> wN[Worker N\n(engine)]
  w1 --> artifacts[out_dir/job_id/*]
  w2 --> artifacts
  wN --> artifacts
  artifacts --> summary[summary.json\n(elapsed_ms per job)]
  summary --> report[Throughput report\n(files/s, p50/p95, failures)]
```

### 9.2 落地工具（当前仓库已有的“最接近”入口）

- 批处理 + 并发：`python3 tools/rnaview_batch.py run ... --workers N`
  - 可分别跑 `--engine legacy` 与 `--engine rust --rust-oracle compute`，对比吞吐与失败率
- 单机对比 benchmark（偏 micro）：`python3 tools/rnaview_bench.py compare --suite phase2 --runs 3`

- 专门的吞吐测试脚本（支持 `--repeat` 扩充 case 数，输出统一 JSON 报告）：
  - 单引擎：`python3 tools/rnaview_throughput.py run --engine rust --workers 8 --repeat 20`
  - 对比：`python3 tools/rnaview_throughput.py compare --workers 8 --repeat 20`
