# 2003–2025 RNA basepair 认知与方法学进化 → RNAView2.0 科学路线图（草案）

本文的目标不是“写历史综述”，而是把 **2003–2025 这二十多年关于 RNA basepair 的认知/方法学进步**，转译成 RNAView2.0 的 **可落地工程任务**：输出契约怎么扩、算法怎么演进、怎么验收、如何分版本发布。

> 时间轴与里程碑在这里先以“组织材料”为主（不追求精确年份/完整引用）。项目启动后建议用 WP0 补齐：关键论文/工具/数据库的引用、确切时间点、以及它们对应的工程决策。

关联文档：
- 规格/验收：`doc/spec.md`
- 架构/算法图：`doc/architecture-diagrams.md`（尤其 §5.1）
- 交付路线：`doc/delivery-plan.md`

---

## 1. 项目目标（Outcome）

1) **兼容性不倒退**：`legacy-v1` 继续保持 Gate B（`.out(full)` byte‑exact）与 Gate D（渲染 diff‑0）全绿。  
2) **科学改进可控发布**：把“修复/增强”作为新的 semantics（`science-v*`）发布，并在 Gate C 中用 allowlist + 解释把差异纳入回归。  
3) **把进步吸收到“可消费的契约”**：以 `pairs.json` 为中心，逐步扩展 schema（版本化、确定性序列化、可回滚），让 Notebook/批处理/下游工具能稳定消费新信息。

---

## 2. 2003–2025 的材料组织方式（建议）

建议把材料按“**问题域**”拆成 5 条主线，每条主线都对应一个可落地的工程模块：

1) **交互类型的本体（Ontology）**：base–base（canonical + non‑canonical）、stacking、base‑phosphate/base‑ribose、triples/multiplets 等，到底要不要进入“权威输出”。  
2) **分类体系与互操作**：Leontis‑Westhof（边类型 + cis/trans）仍是中心，但需要补齐哪些扩展信息（Saenger 对应、isostericity 组、注释/证据字段）。  
3) **结构数据的现实约束**：mmCIF 语义、链 ID、多模型（NMR）、altloc、不完整原子、修饰核苷（CCD/化学组分）等。  
4) **算法与性能**：候选对枚举（O(n²) → 空间筛选）、氢键判定、几何阈值、确定性与可重复性。  
5) **验证范式**：从“和 legacy 一致”走向“和外部共识/基准一致”，以及如何在 gate 体系里逐步引入新的 truth‑set。

---

## 3. 时间轴（按“时代/变化/对 RNAView2.0 的动作”）

> 下表的“时代划分”是为了组织工作包，不代表严格的历史断点；具体引用与精确年份请在 WP0 补齐。

| 时代（粗分） | 变化与方法学趋势（要吸收的点） | 对 RNAView2.0 的动作（可落地） |
|---|---|---|
| 2003–2007（RNAView/非经典配对进入主流） | 非 canonical basepair 与 RNA motif 的系统化描述加速；“注释工具”开始成为基础设施 | 锁住 `legacy-v1` 兼容；把 `pairs.json` 作为权威输出的雏形（已完成 v1） |
| 2008–2012（规模化注释 + family/等构概念） | 对 basepair family、等构（isostericity）与 motif 级别的组织更系统 | 规划 `pairs.json` schema v2：在不破坏 v1 的前提下引入可选字段（family/isostericity/证据） |
| 2013–2016（工具生态成熟 + 更多交互类型） | 以“可互操作的注释”为目标：除 base–base 外，更多工作会显式输出 stacking、base‑phosphate 等 | 把交互类型升级为可扩展枚举（不仅 `pair/stacked`）；新增交互类型先走 `science-v*` 语义试点 |
| 2017–2019（大体系结构 + cryo‑EM/低分辨率现实） | 更大更复杂的 RNA/RNP；结构不完整、链/编号更复杂，鲁棒性与性能更重要 | 把结构解析与归一化做成“政策化”（policy）并落盘到 `pairs.json.options`；引入空间筛选但保持确定性 |
| 2020–2025（预测结构/高通量/不确定性） | 预测结构与多源数据更多；“输出不确定性/证据”变重要 | 在 `pairs.json` 中引入可选的 evidence/confidence 字段；允许缺原子/低置信结构的“降级推断”但必须标注 |

---

## 4. 融合到 RNAView2.0 的技术抓手（推荐设计）

### 4.1 用 semantics 做“科学演进的版本轨道”

- `legacy-v1`：只做兼容（byte‑exact + render diff‑0）。  
- `science-v1`：当前已落地的第一个科学修复轨道（示例：mmCIF hydrogen 处理修复）。  
- 后续按主题推进：`science-v2`（修饰核苷映射/更严格的 donor/acceptor）、`science-v3`（新增交互类型）、……

原则：
- **每个 semantics 版本必须可解释、可回归**：差异进入 Gate C allowlist（带原因），并写入冻结回归集。
- **不把“新信息”硬塞进 legacy 输出格式**：优先放在 `pairs.json` 的新字段里（schema 版本化），legacy `.out` 只在必要时做兼容附注。

### 4.2 用 `pairs.json` schema 做“可消费的知识承载层”

当前 `schema_version=1` 的 basepair 记录主要承载：索引/链与残基号、`kind`、LW/方向、`out_index`（渲染用）、`note/text` 等。

面向 2003–2025 的吸收，建议把未来扩展拆成 3 级（避免一次性重构）：

1) **Schema v1.x（仅加可选字段，不破坏下游）**
   - residue id 更完整：`icode`、`model`、`id_scheme`（mmCIF auth/label）
   - geometry 摘要：距离/角度/平面偏移等（只要确定性、可复现）
2) **Schema v2（交互类型与注释结构升级）**
   - 交互类型扩展：`basepair/stacking/base-phosphate/base-ribose/...`
   - 证据字段：H-bond 列表、donor/acceptor、阈值来源、置信度（可选）
3) **Schema v3（面向 motif/network 的上层对象）**
   - triples/multiplets 的结构化模型
   - motif/graph 表达（可选）

---

## 5. 工作包（WPs）与交付物

### WP0：文献/工具/数据库地图（“把 2003–2025 搞清楚”）

交付物：
- `references/`：关键论文与工具的清单（带年份、链接、一句话影响）
- “里程碑时间轴”补齐到可引用状态（把本文件 §3 的表格落地为可追溯条目）

### WP1：数据集与基准（Truth‑set / Regression‑set）

交付物：
- 一个可复现的结构集合清单（PDB/mmCIF、NMR、多链、插入码、修饰核苷、大体系）
- 对每个语义版本的冻结回归集（Gate B/C/D 的 manifest 扩展）

### WP2：契约与 schema 演进（pairs.json v2+）

交付物：
- `pairs.json` schema 文档（版本化 + 兼容策略）
- schema migration/validation 工具（CI 中可跑）

### WP3：结构解析与归一化（Reality handling）

交付物：
- policy 的完整矩阵（hydrogen/chain id/icode/model/altloc/…）
- mmCIF/PDB 行为差异的测试与文档（何时保持 legacy，何时进入 science）

### WP4：配对/交互算法增强（在 science 轨道里推进）

交付物：
- 每次只推进一小块“可解释差异”（比如：修饰核苷 donor/acceptor 修复）
- 对应的 Gate C diff 报告 + allowlist 更新 + 说明文档

### WP5：渲染与展示（把新注释变成可看见的价值）

交付物：
- 2D/3D 输出在保持 legacy 兼容的同时，能够把新增注释以“可选层”的方式显示出来（例如：用不同颜色/线型区分新交互类型）

---

## 6. 验收口径（建议写进项目章程）

1) **不会破坏 legacy 交付**：任何科学改动不得影响 `legacy-v1` 的 Gate B/D。  
2) **科学改动必须进入 Gate C**：每个差异要么修回，要么解释并进入 allowlist（并冻结回归）。  
3) **schema 变更必须有版本与迁移**：下游可以选择“只看 v1 字段”，同时逐步引入 v2/v3。  
4) **确定性优先**：无论新字段还是新算法，都必须保证在同一输入/同一 policy 下稳定可重复。

---

## 7. 建议的推进节奏（不绑定具体人力的粗排期）

- Sprint 1（1–2 周）：WP0 + WP1（先把材料与基准集准备好，避免后续“改了但不知道对不对”）
- Sprint 2（2–4 周）：WP2（schema v1.x → v2 的设计与最小落地）+ 一条最小的 `science-v2` 试点
- Sprint 3+（滚动）：WP3/WP4/WP5 按主题小步推进，每次只引入一个“可解释差异”

