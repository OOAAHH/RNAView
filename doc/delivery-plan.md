# RNAView 现代化升级：交付视角 Roadmap（基于当前实测）

本文聚焦“交付/发布”而不是算法细节：在现有 gate 与实测数据基础上，总结当前完成度、缺口与下一步里程碑。

- Spec/契约口径：`doc/spec.md`
- 架构图与模块关系：`doc/architecture-diagrams.md`
  - 新版本算法图示（No‑C compute + 2D/3D render）：见 `doc/architecture-diagrams.md` §5.1

## 1. 当前完成度（截至 2026-03-03）

**计算端（Phase 0–3）**

- ✅ Gate B（No‑C + `.out(full)` byte‑exact，legacy-v1）：`bash test_phase2_noc.sh`
- ✅ Gate C（science-v1 diff + allowlist）：`bash test_phase3_gate_c.sh`（当前回归集：ok=14, changed=1(已批准), failed=0）
- ✅ Phase 3 的性能与工程化证据：`out/throughput_sweep_release_w8_20260202_042323/report.md`

**渲染端（Phase 4 / Gate D）**

- ✅ Gate D（legacy baseline + canonicalize + diff‑0）：`bash test_phase4_gate_d.sh`（CI candidate backend：`pairs-out-noc`；当前回归集 `ok=15, unapproved=0, failed=0`）
- ✅ 2D No‑C（XML/PS）：Rust `render-2d` 已与 legacy golden 收敛并通过 Gate D；最后 3 个差异 case 的根因是 legacy `xml2ps.c` 在 `k1==99/999` 时不会输出序号 label（“缺口逻辑”），Rust 端已保持一致。
- ✅ 3D No‑C（VRML）：Rust `render-wrl` 已通过 Gate D（同 `pairs-out-noc` backend）。
- 🆕 Gate NA（DNA-inclusive，science-v1 reviewed regression）：`bash test_phase4_gate_na.sh`（基于 `test/science_dna_cases.json`；golden 在 `test/golden_na/**`；当前 reviewed 回归集包含 pure DNA + hybrid；不要求 legacy C 具备对 DNA case 的分析能力）

## 2. 交付版（v1）Definition of Done（建议）

为避免“gate 都绿了但交付形态不完整”，建议把 v1 的 DoD 明确成三条：

1) **No‑C 运行时闭环**
- 默认路径不编译/链接任何 C，也不 shell out legacy 二进制（legacy 仅允许作为测试 oracle）。

2) **一致性验收口径**
- legacy-v1：Gate B（`.out(full)` byte‑exact）与 Gate D（渲染产物 canonical diff‑0）都必须在 CI/Linux 容器中全绿。
- science-v1：通过 Gate C 的 allowlist 管控差异，差异解释后进入冻结回归（Phase 3.4）。

3) **Python 生态交付（PyPI / Notebook）**
- `pip install rnaview` 后可 `import rnaview`（Notebook 场景成立）。
- 提供可复现的高层 API（建议以 `pairs.json` 为中心契约），并提供 CLI entrypoint。
- `BASEPARS/*` 作为包内资源（不依赖 `$RNAVIEW` 环境变量）。

## 3. 近期里程碑（按交付优先级）

### M4.x：补齐 2D No‑C renderer（已收敛）

目标：实现真正的 `pairs.json (+source.path) → byte‑exact XML/PS`，并切换 Gate D 到 **new-renderer→golden**。

落地拆分（建议顺序）：

1. **Rust CLI/IR 定稿**：`rnaview-hotcore render-2d <pairs.json> --source <pdb/cif> --semantics legacy-v1|science-v1 --out-xml ... --out-ps ...`
2. **复刻 RNAML writer**：按 legacy `rnaxml-new.c` 逐字节复刻 XML（包含缩进/换行/浮点格式/输出顺序；base-pair 顺序用 `BasePair.out_index`）
3. **复刻 xml2ps**：按 legacy `xml2ps.c` 逐字节复刻 PS（并遵循 Gate D 的 canonicalization：仅剔除不稳定 header；注意 legacy 在 `k1==99/999` 时不输出序号 label，Rust 端必须保留该“缺口逻辑”以保持 byte‑exact）
4. **接入 Python 调度**：在 `tools/rnaview_render.py` 增加新 backend（建议命名 `pairs-out-noc`），让 Gate D 的 candidate 不再调用 `bin/rnaview_rustcore_release`
5. **CI 切换**：`.github/workflows/ci.yml` 的 Gate D `CANDIDATE_BACKEND` 从 `pairs-out-noc3d` 切到 `pairs-out-noc`

### M4.y：3D No‑C 交付收尾（已基本完成）

- 现状：Rust `render-wrl` 已通过 Gate D；保持 determinism 与输出稳定即可。
- glTF：目前走确定性 converter（`tools/rnaview_vrml_gltf.py`）；可选增强是“直接输出 glTF/GLB”（非阻塞项）。

### M5：PyPI/Notebook 交付（import rnaview）

建议先做“可用、可发布、可维护”，再做“最优性能”：

1. **包结构**：引入 `pyproject.toml`，把 `tools/` 中可复用逻辑沉淀为 `rnaview/` 包（CLI 脚本保留但改为 thin wrapper）。
2. **后端策略（推荐）**：
   - v0：wheel 内置 `rnaview-hotcore`（subprocess 调用；实现快、风险低）
   - v1：可选 PyO3/maturin（把 Rust core 暴露成 extension，减少 subprocess 开销，Notebook 体验更好）
3. **资源打包**：把 `BASEPARS/*` 变成 package data，并提供 `rnaview.data_dir()`/`rnaview.get_resource_path()` 之类的稳定 API。
4. **最小可用 API**（示例）：
   - `rnaview.analyze(path, *, semantics="legacy-v1", formats=("out","pairs","ps","svg","wrl","gltf"))`
   - `rnaview.render(pairs_json, *, source_path, formats=...)`
5. **发布与 CI**：增加 wheel 构建与 smoke test（至少验证 `import rnaview` + 跑一个小 case；建议用 `.github/workflows/wheels.yml`）。

## 4. 关键风险与决策点（现在就要想清楚）

- **2D byte‑exact 的工作量**：PS/XML 都是“格式敏感”的 legacy 输出，任何浮点/空白/排序差异都会导致 diff；建议把 2D renderer 当成一个独立子项目推进。
- **SVG 的可信口径**：当前 Gate D 的 SVG 是从 PS 经 converter 生成；若要保证“SVG 视觉等价于 PS”，需要额外的视觉/光栅化 sanity（不建议一开始就做成 hard gate）。
- **跨平台 wheel**：v0 先落地 Linux x86_64 + macOS arm64 的平台专用 wheels；Windows / 其他平台后续补齐（建议走 CI/cibuildwheel/auditwheel 的标准发布链路）。

## 5. Legacy → RNAView2.0：代码升级说明（旧到新的“讲解版”）

> 这节不是 spec（验收口径看 `doc/spec.md`），而是“为什么要这么分层、怎么一步步升级”的叙事摘要，方便对外讲清楚/对内 onboarding。

### 5.1 总体变化（一句话）

从 **C 单体（解析+计算+输出+渲染强耦合）** 升级为 **Python 编排 + Rust 热核（可验证、可替换、可逐步科学修复）**，并用 gate 把“兼容性”和“科学改进”分开管理。

### 5.2 迁移策略（按 Gate/Phase 的“可回归切片”）

- **Phase 0–1（oracle + 契约先行）**：把 legacy `.out` 抽成 `pairs.json`（schema v1、确定性序列化），先锁住“科学 core 等价”的回归体系。
- **Phase 2（byte‑exact + No‑C）**：
  - Gate A（可选桥接）：`bin/rnaview_rustcore(_release)` 继续跑 C pipeline/writer，但把热点函数替换为 Rust，快速锁定 `.out` 字节级一致。
  - Gate B（最终交付）：`rnaview-hotcore --oracle compute` 走纯 Rust 解析+计算+writer，默认运行时不依赖任何 C（legacy 仅做测试 oracle）。
- **Phase 3（科学差异可控）**：引入 `--semantics legacy-v1|science-v1` + policy，允许在 Gate C 中用 allowlist 管控“被解释过的差异”。
- **Phase 4（渲染收敛）**：渲染侧以 `pairs.json (+source)` 为输入，通过 Gate D 把 XML/PS/WRL 的输出稳定性锁住（必要时做 canonicalize）。

### 5.3 旧 → 新模块映射（定位代码用）

| Legacy C（主要职责） | RNAView2.0（当前仓库实现落点） |
|---|---|
| `src/rnaview.c`（CLI 编排 + 解析/清洗 + work_horse） | Python 编排：`tools/rnaview_batch.py` / gates；Rust CLI：`rust/src/main.rs` |
| `include/cifparse.c`（legacy mmCIF 解析） | Rust 结构解析：`rust/src/structure.rs`（并在 CLI 层提供 mmCIF parser 策略） |
| `src/fpair*.c`/`pair_type.c`（配对核心） | Rust port：`rust/src/legacy_pairing.rs`、`rust/src/legacy_alg.rs` |
| `.out` writer（legacy 输出编排） | Rust IR+writer：`rust/src/out_full.rs` |
| 2D：`ps-xy*.c` + `rnaxml-new.c` + `xml2ps.c` | Rust 2D：`rust/src/legacy_2d_layout.rs` + `rust/src/legacy_rnaml.rs` + `rust/src/legacy_xml2ps.rs`（入口 `rust/src/render_2d.rs`） |
| 3D：`vrml.c` | Rust 3D：`rust/src/vrml_render.rs` |
| `$RNAVIEW/BASEPARS/*`（运行时数据） | 交付形态：作为包内资源 + 稳定 API（规划见本文件 §3 M5） |

### 5.4 新架构的“最小心智模型”（写代码时盯这三个点）

1) **`pairs.json` 是中心契约**：批处理、回归、渲染、未来 API 都优先围绕它组织（而不是围绕 legacy `.out` 文本）。
2) **兼容与科学分叉用 semantics 管**：`legacy-v1` 追求 byte‑exact；`science-v*` 追求科学修复，但必须在 Gate C 中可解释、可回归。
3) **render 是独立产品线**：尽量做到 `pairs.json (+source)` → 图形产物；渲染差异不要污染 compute 的核心回归口径。

## 6. 2003–2025 basepair 认知进化吸收进 RNAView2.0：项目计划（草案）

这部分不属于“短期交付 v1”的 hard gate，但建议尽早启动（先做材料/基准集），否则后续每一次科学增强都会变成“无锚点的争论”。

- 计划文档：`doc/basepair-evolution-2003-2025.md`
- 推荐落地抓手（和当前工程体系对齐）：
  - 用 `Semantics`（`science-v*`）承载科学演进，用 Gate C allowlist 把差异解释进回归
  - 用 `pairs.json` schema 版本化承载新注释（v1.x → v2 → v3），保持确定性与可迁移
  - 用 WP0/WP1 把 2003–2025 的关键材料与 truth/regression‑set 固化下来
