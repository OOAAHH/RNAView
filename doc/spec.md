# RNAView 现代化升级：Spec 层建模（Legacy C vs Rust+Python）

本文是“规格/契约（spec）层”的建模：描述系统**做什么**、输入输出**是什么**、一致性**如何判定**、模块边界**如何划分**。不追求贴近某种实现细节，目的是让团队在重构/替换时有共同语言与验收口径。

## 0. 目标与非目标

### 0.1 目标（已确认）

- **Phase 0/1（强约束）**：升级前后 `.out` 的 core 内容一致：
  - `BEGIN_base-pair … END_base-pair`
  - `BEGIN_multiplets … END_multiplets`
  - `The total base pairs = ... (from ... bases)` 及其后的统计表
- **Phase 2（强约束）**：
  - 对同一输入与同一组选项，`FILEOUT.out` **逐字节一致**（可直接 `diff`）。
  - **No-C**：最终交付路径不依赖 C（不编译/链接 C，也不 shell out legacy 二进制）；legacy 仅允许作为测试/回归 oracle。
  - 验收脚本（当前仓库）：Gate A `test_phase2.sh`；Gate B `test_phase2_noc.sh`。
- **批处理优先（第一阶段）**：支持“很多结构跑库”的高吞吐执行、可并发、可重跑、可汇总。
- **技术分工**：对性能有要求的区域用 Rust；编排、批处理、落盘与生态集成用 Python。
- **权威结构化输出**：`pairs.json`（确定性序列化，可字节级 diff）。

### 0.2 非目标（第一阶段不做/弱约束）

- PS/VRML 等渲染输出的逐字节一致；它们作为派生产物后置。
- 交互式可视化（Notebook/Web）后置。

## 1. 术语（Domain Vocabulary）

- **Structure**：输入的三维结构（PDB 或 mmCIF），可包含多个模型（NMR）与多个链。
- **Residue**：残基；RNAView 逻辑上以“残基”为单位形成 base index（1-based）。
- **ResidueId**：用于唯一定位残基的一组字段（至少包含 chain、resseq；mmCIF 还涉及 insertion code / label/auth 选择）。
- **Base**：被 RNAView 视为核酸碱基的残基（含常见修饰碱基），以一个字符表示（可能是大写或小写）。
- **BaseIndex**：RNAView 内部编号（输出中 `n1_n2` 的 `n1,n2`），从 1 开始。
- **BasePair**：两碱基之间的配对或堆叠关系（包含 LW 分类、cis/trans、Saenger 对应、tertiary 标注等）。
- **Multiplet**：三联体及更高阶的多重配对组合（从 base-pair 结果推导）。
- **Stats**：配对类型统计汇总（表格形式）。

## 2. Legacy C 系统（现状）的 Spec 模型

### 2.1 系统边界（System Context）

输入：
- 单个结构文件：PDB / mmCIF
- 或 RNAML(XML)：作为“渲染/转换”输入（走 xml2ps）
- 或“文件列表”：触发批处理路径

输出（关键）：
- `*.out`：文本主结果（包含 base-pair、multiplets、统计）
- 以及派生：`*.xml`（RNAML）、`*.ps`、`*.wrl`、`*_torsion.out`、`*_sort.out`、`*_patt.out` 等

### 2.1.1 Legacy 模块地图（用于定位，不是契约）

| Spec 元素 | 主要实现位置（legacy C） | 备注 |
|---|---|---|
| 输入探测（PDB/mmCIF/RNAML） | `src/rnaview.c` | 通过扫描文件内容而非严格依赖 `-x` |
| PDB 清洗（去 ANISOU/选链/分辨率过滤） | `src/rnaview.c:clean_inpfile` | 会影响“进入核心计算”的结构内容 |
| mmCIF 解析（auth/label、NMR 选模） | `include/cifparse.c` + `src/rnaview.c` | legacy 链 ID 常被截断为 1 字符 |
| residue/base 编号（BaseIndex） | `residue_idx`（多处调用） | 输出的 `i/j` 来源于此顺序 |
| 修饰碱基映射（uncommon -> a/g/u/c/t/P/I） | `src/fpair_sub.c:get_seq/identify_uncommon` | 大小写带语义（常见：小写=修饰） |
| base-pair core 输出 | `src/fpair.c:all_pairs` | `BEGIN_base-pair` 块 |
| multiplets core 输出 | `src/fpair_sub.c:multiplets` | `BEGIN_multiplets` 块 |
| total/stats 输出 | `src/rnaview.c:work_horse` + `pair_type_statistics/...` | `The total base pairs...` 及统计表 |

### 2.2 高层流水线（黑盒）

```
CLI/argv
  └─(探测输入类型：PDB/mmCIF vs RNAML)
      ├─ RNAML -> xml2ps -> PS (结束)
      └─ Structure -> (清洗/过滤/选链/选模/去 H/去 water/altloc 等)
            └─ 解析为扁平数组 Atom/Res/Misc/xyz + residue index (seidx)
                └─ 识别核酸残基 + base 序列(bseq) + 修饰碱基映射
                    └─ all_pairs: 枚举 i<j 残基对 -> 几何筛选 -> Hbond -> LW 分类 -> 输出 BEGIN_base-pair
                        └─ write_multiplets + 统计表输出
```

> spec 视角关心的是：哪些输入决定了 `.out` 的科学内容（core），以及在 Phase 2 里如何把 `.out` 的“格式/文本”也纳入逐字节一致的验收口径。

### 2.3 Legacy `.out` 的“core 契约”

#### 2.3.1 Base-pair core 块

core 块由以下结构构成：

```
BEGIN_base-pair
<0..N 行 base-pair/stacked 行>
END_base-pair
```

每行（典型）是：

```
<i>_<j>, <chain_i>: <resseq_i> <base_i>-<base_j> <resseq_j> <chain_j>: <type> <pa_int> <syn_i><syn_j> <saenger_or_note>
```

其中：
- `<i>,<j>`：BaseIndex（1-based）
- `<chain_*>`：单字符链 ID（legacy mmCIF 路径会截断到 1 字符）
- `<resseq_*>`：整数残基号（PDB `ResSeq`；mmCIF 对应 auth/label seq id）
- `<base_*>`：单字符碱基代码（A/G/U/C/T/I/P 或其小写变体，legacy 用小写表示“uncommon/modified”一类）
- `<type>`：LW 或 WC 风格类型字符串（如 `+/+ cis`、`W/H tran`、`S/? cis` 等；legacy 会出现 `tran`）
- `<pa_int>`：tertiary/附加标注字符（legacy 用 `!` 等标记）
- `<syn_i><syn_j>`：是否为 syn 构象（`syn` 字样）
- `<saenger_or_note>`：Saenger 对应（罗马数字或 `n/a`）与/或额外说明（如 `!(b_s)`）

**Spec 一致性口径（core）**：
- 不要求保持行的原始空格对齐，但要求“记录集合”在字段层面一致（见第 4 节的新系统契约）。
- stacked 行是 core 的一部分（它们也属于科学信息）。
- 允许 legacy 中存在无法解析的“异常行”，但升级后应尽量避免；若仍存在，需作为 `unknown` 记录保留原文。

#### 2.3.2 Multiplets core 块

```
BEGIN_multiplets
<0..M 行>
END_multiplets
```

Spec 一致性口径（core）：
- 以“参与的 base index 组合 + 描述文本”等价为准（组合顺序可忽略，文本需可比对或规范化）。

#### 2.3.3 Stats core 块

核心锚点：
- `The total base pairs =  <X> (from <Y> bases)`

其后的统计表（legacy 输出为 ASCII 表格）在 spec 中表示为：
- `total_pairs: X`
- `total_bases: Y`
- `pair_type_counts: { "WW--cis": n, "WW-tran": n, ... }`（键集合随 legacy 输出而定）

### 2.4 Legacy 的关键语义（必须显式进入新 spec）

这些语义会直接影响 base index、配对集合与分类，因此必须被新实现显式建模：

- **BaseIndex 定义**：base index 是“参与计算的核酸残基序列”的 1-based 编号；它与输入文件中 residue 的原始编号不是同一个概念。
- **修饰碱基映射**：legacy 对“uncommon residue”按几何/原子存在性猜测为 `a/g/u/c/t` 或 `P/I` 等；大小写携带语义（常见：小写=修饰/不标准）。
- **mmCIF 的 label/auth 选择**：legacy 通过 `--label` 切换；并且链 ID 被截断为 1 字符（这是兼容性风险点，必须在新 spec 中决定是否保持）。
- **NMR 选模**：mmCIF NMR 通过 `_pdbx_nmr_representative/_pdbx_nmr_ensemble` 选择代表模型；没有就默认模型 1。
- **过滤策略**：去水、去 H、去 ANISOU、altloc 处理、分辨率过滤（`-r`）等会改变输入结构集合，从而改变结果。

## 3. 新系统（Rust + Python）的 Spec 模型

### 3.1 系统分层（Components）

```
Python (Orchestration)
  - CLI / 批处理调度 / 并发 / 重跑 / 日志 / 目录布局
  - I/O glue：调用 Rust engine（legacy oracle 仅用于测试/回归对照）
  - 产物落盘：pairs.json + (可选) .out writer
  - 回归验证：对照 golden 的 core 等价性

Rust (Hot Core Engine)
  - 结构域对象：Residue/Atom/Geometry
  - 核酸识别与修饰碱基映射
  - 候选对枚举 + 几何筛选
  - Hbond 枚举
  - LW 分类与 Saenger 对应
  - （可选）multiplets 推导与统计
```

迁移阶段推荐形态（按 gate）：
- Phase 0/1：允许 Python 调 `bin/rnaview` 做 oracle，确保回归对齐；
- Phase 2 Gate A：允许复用 legacy C writer（`legacy/rustcore`）做 `.out` byte-exact 的桥接；
- Phase 2 Gate B（No-C）：默认执行路径不再调用 legacy/C。

### 3.2 统一的 Domain Model（语言无关）

#### 3.2.1 基础对象

- `ResidueId`
  - `chain: str`（允许多字符；是否降级为 1 字符要作为兼容开关）
  - `resseq: int`
  - `icode: str`（`" "` / `"A"` / `"?"` 等；legacy mmCIF 会出现 `?`）
  - `id_scheme: "auth" | "label"`（mmCIF）
  - `model: int`（NMR）

- `Base`
  - `letter: str`（单字符，保持 legacy 语义：可大写/小写）
  - `resname: str`（3-letter，如 `PSU`）
  - `is_modified: bool`（从 letter 大小写/映射来源推导）

#### 3.2.2 核心结果对象

- `BasePair`
  - `i: int`, `j: int`（BaseIndex，1-based）
  - `res_i: ResidueId`, `res_j: ResidueId`
  - `base_i: str`, `base_j: str`
  - `kind: "pair" | "stacked" | "unknown"`
  - `lw: str | null`（如 `W/H`、`+/+`、`-/-`、`S/?`）
  - `orientation: "cis" | "tran" | null`
  - `syn: { "i": bool, "j": bool }`
  - `note: str | null`（如 `XIX`、`n/a`、`!(b_s)`、`!1H(b_b)` 等 legacy 尾注）

- `Multiplet`
  - `indices: list[int]`
  - `text: str`（保留 legacy 描述，便于一致性验证）

- `Stats`
  - `total_pairs: int`
  - `total_bases: int`
  - `pair_type_counts: dict[str,int]`

### 3.3 输出契约（Output Contracts）

#### 3.3.1 `.out`

`.out` 在新系统中分两类实现路径：
- A) legacy / rustcore（复用 legacy C writer）生成：仅作为桥接/对照工具，不满足 Phase 2（No-C）最终验收。
- B) No-C writer（Rust/Python）生成：需要复刻 legacy 的 `.out` 全文格式，用于 Phase 2（No-C）最终验收。

**验收标准分阶段**：
- Phase 0/1：仅对比 core（base-pair/multiplets/stats），非 core 内容不纳入一致性。
- Phase 2 Gate A：`FILEOUT.out` 逐字节一致（可直接 `diff`；允许 A 路径桥接）。
- Phase 2 Gate B（No-C）：`FILEOUT.out` 逐字节一致 **且** 默认执行路径 No-C（B 路径）。

对比工具：
- core：`tools/rnaview_out_core.py`（抽取 core → 规范化 JSON → compare）
- `.out` 全文：`tools/rnaview_batch.py --regress-mode out`（逐字节 diff）

#### 3.3.2 `pairs.json`（权威）

目标：确定性序列化，便于大规模跑库的回归与缓存。

建议 schema v1（当前实现约定）：

```
{
  "schema_version": 1,
  "source": { "path": "...", "format": "out|pdb|cif", "id_scheme": "auth|label", "model": 1 },
  "options": { "...": "pipeline options (optional)" },
  "core": {
    "base_pairs": [ BasePair ... ],
    "multiplets": [ Multiplet ... ],
    "stats": Stats
  }
}
```

其中 `core` 字段的结构与 `tools/rnaview_out_core.py extract` 输出一致，便于复用回归。

**确定性规则（必须写进 spec）**：
- `base_pairs` 按 `(i,j,kind,...)` 排序；
- `pair_type_counts` 按 key 排序；
- JSON 序列化固定：`sort_keys=true`、固定 `separators`、UTF-8、无随机字段（时间戳等）。

#### 3.3.3 Semantics & Policies（Phase 3+）

为了把 **Phase 2 的 byte-exact 兼容性** 与 **Phase 3 的科学修复** 解耦，新系统引入显式的：

- `--semantics legacy-v1|science-v1`：一个“语义预设（preset）”，定义默认 policy 组合。
- 一组可组合的 `policy`：用户可以按自己的科学假设覆盖 preset 的默认值。

这些参数必须写入 `pairs.json.options`（**落盘 + 可追溯**），以保证批处理结果可复现、可审计。

当前第一批 policy（结构解析层）：

- `structure.hydrogen_policy`
  - `legacy-mmcif-bug`：复刻 legacy mmCIF 去氢 bug（保留部分 4 字符 H 原子名），用于 `legacy-v1` 的 byte-exact 回归。
  - `discard-all`：彻底去氢（`science-v1` 默认）。
  - `keep-all`：保留所有氢（用于研究“氢坐标对分类的影响”，不保证与 legacy 一致）。
- `structure.missing_insertion_code_policy`
  - `legacy-question-mark`：复刻 legacy：mmCIF 缺失/`.` 的 insertion code 被映射为 `?`。
  - `none`：科学模式（可选覆盖）：mmCIF 缺失/`.`/`?` 的 insertion code 视为缺失（`None`）。
- `structure.chain_id_policy`
  - `legacy-1char`：复刻 legacy：链 ID 截断为 1 字符（可能产生碰撞）。
  - `unique-1char`：科学模式（可选覆盖）：将不同链 ID 映射到不同的 1 字符（保持 `.out` 格式要求，同时避免碰撞）。

preset 默认值（便于回归口径）：

- `legacy-v1`：`legacy-mmcif-bug` + `legacy-question-mark` + `legacy-1char`
- `science-v1`：`discard-all` + `legacy-question-mark` + `legacy-1char`（先聚焦去氢 bug；其余差异通过显式 policy flag 启用并用 Gate C 管控）

接口约束：

- 这些参数仅对 **`--oracle compute`（纯 Rust 计算）** 生效；对 `legacy/out` oracle 不允许传入（避免误导性的“伪科学配置”）。
- `legacy-v1` 的 `.out(full)` 必须保持 byte-exact，因此 **不得**在 `.out` 文本中增加任何新字段；追溯信息写入 `pairs.json.options`。

#### 3.3.4 Gate C（science-v1 差异管理）

Gate C 用于把 “science-v1 相对 legacy-v1 的差异” 做可落盘、可审阅、可批准的管理：

- 工具：`python3 tools/rnaview_gate_c.py run ...`（仓库内验收脚本：`bash test_phase3_gate_c.sh`）。
- allowlist：`test/gate_c_allowlist.yaml`（YAML 兼容的 JSON；存放已批准的 diff event id）。
- 退出码：仅当 **无 failed 且无 unapproved** 时返回 `0`；允许存在 “changed（已批准差异）”。
- 输出：`<out_dir>/summary.json`、`<out_dir>/report.md`、`<out_dir>/cases/<job_id>/{legacy-v1,science-v1}/...` 与 `diff.json`（差异详情）。

#### 3.3.5 science-v1 Golden 回归（Phase 3.4）

当 Gate C 的差异已经被解释并批准后，science-v1 需要从“差异解释”进入“稳定回归”：

- 冻结工具：`python3 tools/rnaview_science_golden.py freeze`
  - 输出：`test/golden_science_core/manifest.json`
  - 约定：manifest 的 `core_json` 可以指向 `test/golden_core/**.core.json`（当 science 与 legacy 完全一致时），也可以指向 `test/golden_science_core/**.core.json`（当存在已批准差异时）。
- 回归验收脚本：`bash test_phase3_science.sh`
  - 对 `--semantics science-v1` 的 core 做回归（不比较 `.out` bytes）。

### 3.5 批处理（第一阶段）建议的外部接口契约

> 这是面向“跑库”的接口 spec（不绑定具体实现）。

- 输入选择：支持单文件、目录递归、glob、清单文件（list）
- 并发执行：可配置 `workers`；同一输入应当幂等（可重跑不破坏结果）
- 输出布局（建议）：
  - `<out_dir>/<job_id>/pairs.json`
  - `<out_dir>/<job_id>/legacy.out`（若 engine=legacy，便于审计）
  - `<out_dir>/<job_id>/engine.out`（若 engine=rust；byte-exact 回归可对它做 diff，路径也会写入 `summary.json` 的 `out_path`）
  - `<out_dir>/summary.json`（成功/失败计数、错误列表、耗时）
- 退出码：
  - `0`：全部成功且通过回归（`core`/`.out`，若启用）
  - `1`：存在失败或回归不一致
  - `2`：参数/输入错误
  - `3`：内部异常（bug）

### 3.4 “科学一致性”的正式定义（用于验收）

给定同一输入结构与同一组选项，legacy 与新实现应满足：

1. `Stats.total_pairs/total_bases` 相同
2. `Stats.pair_type_counts` 的键集合与计数相同（允许未来扩展字段，但回归集内必须对齐）
3. `BasePair` 集合等价（顺序无关）：
   - `(i,j,chain_i,resseq_i,base_i, chain_j,resseq_j,base_j, kind, lw, orientation, syn, note)` 等价
   - stacked 视为 `kind="stacked"` 的记录参与等价
4. `Multiplet` 集合等价（顺序无关）：`indices + text` 等价

> 注：如果未来决定“链 ID 截断”不是科学含义，则需要在 spec 中声明 canonicalization（例如：新系统保留完整 chain，但对比时投影到 legacy 的 1 字符）。

## 4. 迁移时的“规格驱动”验收方法（建议）

1. 以 `test/**.out` 为 golden，先用 `tools/rnaview_out_core.py extract` 固化 core 期望。
   - 批量生成：`python3 tools/rnaview_out_core.py freeze test`（输出到 `test/golden_core/`，含 `manifest.json`）
2. 新实现每个阶段只要能产出 `.out` 或 `pairs.json`，都必须能投影成同一份 core，并通过 compare。
3. 当 Rust engine 能直接产出 `pairs.json` 时：
   - 继续保留 `pairs.json -> .out(core)` 作为低风险的纯函数回归工具（用于核心语义验证与最小化 diff）。
   - 若 Phase 2 要求 `.out` **byte-exact** 且 **No-C**，则必须进一步实现 `.out(full)` writer，并把非 core 字段也纳入可验证的数据模型（必要时扩展 `pairs.json`/另加 AST），避免隐藏依赖与非确定性。
   - 现有 core writer/验证工具：`python3 tools/rnaview_pairs_json.py validate-golden`

---

如果你希望下一步更“工程化”：我可以基于这份 spec 再补一份 `pairs.json` 的严格 JSON Schema（draft-07/2020-12）+ 若干 golden 样例（从 `test/**.out` 自动导出），用于后续 CI/回归。 

## 5. Phase 4：渲染与格式现代化（Gate D）

Phase 4 的目标是：**在不改变 core（Phase 0–3 已锁定）** 的前提下，把 legacy 的 2D/3D 产物做成“可复现、可回归、可对齐”的现代化输出（并新增 SVG/glTF）。

本阶段的“完全一致”定义为：**先做 normalization，再 byte‑exact（diff 0）**；并且 **baseline 是 legacy `rnaview -p/-v`**（在 CI/Linux 容器中）。

### 5.1 输入契约（Render Inputs）

渲染器以 `pairs.json` 为主输入（权威 core），并允许按 `pairs.json.source.path` 再读取一次原始结构文件（PDB/mmCIF），以最小化 `pairs.json` schema 变动。

渲染结果必须是以下输入的纯函数：
- `pairs.json.core`（base_pairs/multiplets/stats）
- `pairs.json.source.path` 指向的结构文件内容（用于坐标/链信息/残基名称等）
- 固定的一组渲染选项（必须落盘到 `pairs.json.options`，至少包括）：
  - mmCIF：`id_scheme=auth|label`
  - chain 选择（如 `chains=ABC` 或 `all`）
  - NMR model 选择策略
  - altloc 策略
  - 分辨率过滤策略（若启用）
- `semantics` + `policies`（legacy-v1/science-v1；渲染要能在两种语义下工作）

### 5.2 输出契约（Render Outputs）

对每个输入 job（与 Phase 2/3 相同的 job_id 规则），渲染器至少产出：

2D：
- `*.xml`（RNAML；legacy 兼容）
- `*.ps`（PostScript；legacy 兼容）
- `*.svg`（新增；与 legacy 2D 视图语义等价）

3D：
- `*.wrl`（VRML；legacy 兼容）
- `*.gltf` 或 `*.glb`（新增；与 VRML 语义等价）

> 注：SVG/glTF 的 baseline 仍然来自 legacy：2D 推荐以 legacy 的 `*.ps` 作为“渲染权威”（最终画出来的结果），用确定性的 PS→SVG 转换器生成 SVG golden；3D 仍以 legacy `*.wrl` 作为语义权威，再生成 glTF golden。

### 5.3 Normalization（仅用于验收的规范化）

Gate D 只对 **canonical（normalize 后）** 的产物做 byte‑exact 比较；golden 也存 canonical。

Normalization 的目标是去掉“与渲染语义无关但会导致 diff”的字段；规则应当最小化且按格式分开定义（便于审计）。

通用规则：
- 统一换行符为 `LF`
- 去掉行尾空白
- 文件末尾保证单个 `LF`

PS（`.ps`）：
- 去掉 `%%CreationDate:` 等时间戳行（其余保持不变）
- 字体/线宽/颜色/版面参数必须复刻 legacy（基线见 `BASEPARS/ps_image.par`）

VRML（`.wrl`）：
- 去掉 `# Creation Date:` 与 `# UserName:` 行（其余保持不变）

RNAML/XML（`.xml`）与 SVG（`.svg`）：
- 保守起见，默认只做通用规则（不做重排/重格式化），避免破坏文本节点中的表格布局
- 若后续发现“仅空白差异”高频，可再引入更强的 canonicalization（需写入 spec 并产出迁移说明）

glTF（`.gltf/.glb`）：
- 若输出为 `.gltf`（JSON），canonical 形式应满足：键排序稳定、浮点格式稳定、无随机 id
- 若输出为 `.glb`，建议回归时先投影为 canonical `.gltf`（或自定义 IR）再比较

Normalization 是否作为“最终落盘输出”：
- 默认只用于验收（保留 raw 输出便于 debug）
- 若 canonical 与 raw 在可视化上无可见差异，可选择直接落盘 canonical 作为最终输出

### 5.4 Gate D（渲染验收：legacy 基线 + allowlist）

Gate D 的目标是把渲染产物纳入与 Phase 2/3 同等级别的回归体系：

- 回归集：与 Gate C 一致（`test/pdb` + `test/mmcif`）
- baseline：legacy `rnaview -p`（2D）与 `rnaview -v`（3D）
- golden：存储 canonical 输出（建议压缩存储，例如 `*.canon.gz`）
- 允许差异：与 Gate C 类似，使用 allowlist 管控（例如 mmCIF 去氢 bug 在 `science-v1` 下导致的可视化差异）

工具与目录约定（当前仓库）：
- 冻结（生成 canonical goldens）：`python3 tools/rnaview_gate_d.py freeze`（默认输出到 `test/golden_render/manifest.json`）
- 对比（生成报告与 diff 事件）：`python3 tools/rnaview_gate_d.py compare --out-dir out_phase4_gate_d`
  - candidate renderer（默认）：`python3 tools/rnaview_render.py render`（可用 `--candidate-cmd ...` 替换成 Rust/新渲染端；或用 `--candidate-engine legacy` 强制直接跑 legacy；注意 `--candidate-cmd` 需要放在命令行最后）
  - new-renderer→golden（推荐入口）：`CANDIDATE_BACKEND=rustcore bash test_phase4_gate_d.sh`（或 `rustcore-release`）
  - 视觉 sanity（人工）：对少量代表 case，把 legacy `*.ps` 与 `out_phase4_gate_d/cases/*/candidate.svg` 并排打开确认“画得对”。必要时先用 `python3 tools/rnaview_render.py render --input <file> --out-dir out_legacy_render --formats ps,xml,wrl` 生成一份干净的 legacy 输出用于对照。

allowlist：
- 文件：`test/gate_d_allowlist.yaml`
- 粒度：事件级（稳定 id）
- id 必须至少包含：`input_id` + `format(ps|xml|svg|wrl|gltf)` + `semantics` + `renderer_version` + `hash(payload)`

退出码建议：
- `0`：无 failed 且无 unapproved diff
- `1`：存在 failed 或 unapproved diff
- `2`：参数/输入错误

### 5.5 样式映射（legacy style map）

为了保证 PS/SVG 在“视觉语义”上一致，样式映射必须进入规格：
- 字体：使用 legacy 的内建字体设置（当前基线见 `BASEPARS/ps_image.par`，例如 `/Times-Bold`）
- 线宽：`W1/W2/W3/W4`（`1 / 1.5 / 2 / 3`）
- 颜色：按 legacy 的 HSB/RGB 配置复刻（`BASEPARS/ps_image.par` 中的 `Al/Cl/Gl/Ul/...` 以及 minor/major saturation）
- 线型：`LINE/DASHLINE` 等（legacy dash 模式 `[2 4]`）

> 这些参数在 Phase 4 的实现中应当被抽象成可版本化的“style preset”（例如 `style=legacy-ps-v1`），并作为 `renderer_version` 的一部分参与 Gate D 的稳定 id 计算。
