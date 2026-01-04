# RNAView 现代化升级：Spec 层建模（Legacy C vs Rust+Python）

本文是“规格/契约（spec）层”的建模：描述系统**做什么**、输入输出**是什么**、一致性**如何判定**、模块边界**如何划分**。不追求贴近某种实现细节，目的是让团队在重构/替换时有共同语言与验收口径。

## 0. 目标与非目标

### 0.1 目标（已确认）

- **科学一致性（强约束）**：升级前后 `.out` 的 core 内容一致：
  - `BEGIN_base-pair … END_base-pair`
  - `BEGIN_multiplets … END_multiplets`
  - `The total base pairs = ... (from ... bases)` 及其后的统计表
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

> spec 视角关心的是：输出 `.out` 的 core 内容由哪些输入决定，以及哪些“环境/实现细节”不应纳入一致性口径。

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
  - I/O glue：调用 legacy oracle 或 Rust engine
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

第一阶段推荐形态：
- Python 默认仍可调用 `bin/rnaview` 作为 oracle，确保回归对齐；
- Rust engine 逐块替换热点逻辑，替换一块就跑回归集。

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

#### 3.3.1 `.out`（只比较 core）

`.out` 在新系统中分两种来源：
- A) legacy oracle 直接生成
- B) 新实现由 `pairs.json`/内存结构重新渲染生成（推荐：保证确定性与一致性）

**验收标准**：仅对比 core（base-pair/multiplets/stats），非 core 内容不纳入一致性。

对比工具：`tools/rnaview_out_core.py`（抽取 core → 规范化 JSON → compare）。

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

### 3.5 批处理（第一阶段）建议的外部接口契约

> 这是面向“跑库”的接口 spec（不绑定具体实现）。

- 输入选择：支持单文件、目录递归、glob、清单文件（list）
- 并发执行：可配置 `workers`；同一输入应当幂等（可重跑不破坏结果）
- 输出布局（建议）：
  - `<out_dir>/<job_id>/pairs.json`
  - `<out_dir>/<job_id>/legacy.out`（若 engine=legacy，便于审计）
  - `<out_dir>/summary.json`（成功/失败计数、错误列表、耗时）
- 退出码：
  - `0`：全部成功且通过 core 回归（若启用）
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
3. 当 Rust engine 能直接产出 `pairs.json` 时，把 `.out` writer 作为纯函数（`pairs.json -> .out core`），从根上减少不确定性。
   - writer/验证工具：`python3 tools/rnaview_pairs_json.py validate-golden`

---

如果你希望下一步更“工程化”：我可以基于这份 spec 再补一份 `pairs.json` 的严格 JSON Schema（draft-07/2020-12）+ 若干 golden 样例（从 `test/**.out` 自动导出），用于后续 CI/回归。 
