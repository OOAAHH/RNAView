# RNAView 现代化升级：交付视角 Roadmap（基于当前实测）

本文聚焦“交付/发布”而不是算法细节：在现有 gate 与实测数据基础上，总结当前完成度、缺口与下一步里程碑。

- Spec/契约口径：`doc/spec.md`
- 架构图与模块关系：`doc/architecture-diagrams.md`

## 1. 当前完成度（截至 2026-02-07）

**计算端（Phase 0–3）**

- ✅ Gate B（No‑C + `.out(full)` byte‑exact，legacy-v1）：`bash test_phase2_noc.sh`
- ✅ Gate C（science-v1 diff + allowlist）：`bash test_phase3_gate_c.sh`（当前回归集：ok=14, changed=1(已批准), failed=0）
- ✅ Phase 3 的性能与工程化证据：`out/throughput_sweep_release_w8_20260202_042323/report.md`

**渲染端（Phase 4 / Gate D）**

- ✅ Gate D（legacy baseline + canonicalize + diff‑0）：`bash test_phase4_gate_d.sh`
- ✅ 3D No‑C（VRML）：`CANDIDATE_BACKEND=pairs-out-noc3d` 下，`.wrl` 由 Rust `render-wrl` 生成并通过 Gate D。
- ⚠️ 2D 仍非 No‑C：`pairs-out-noc3d` 目前的 `.xml/.ps` 仍由 `bin/rnaview_rustcore_release`（C）生成，只是通过 `RNAVIEW_OUT_PATH` 注入新的 `.out`；因此 Gate D 目前还不能证明“2D 渲染彻底 No‑C”。

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

### M4.x：补齐 2D No‑C renderer（关键阻塞项）

目标：实现真正的 `pairs.json (+source.path) → byte‑exact XML/PS`，并切换 Gate D 到 **new-renderer→golden**。

落地拆分（建议顺序）：

1. **Rust CLI/IR 定稿**：`rnaview-hotcore render-2d <pairs.json> --source <pdb/cif> --out-xml ... --out-ps ...`
2. **复刻 RNAML writer**：按 legacy `rnaxml-new.c` 逐字节复刻 XML（包含缩进/换行/浮点格式/输出顺序；base-pair 顺序用 `BasePair.out_index`）
3. **复刻 xml2ps**：按 legacy `xml2ps.c` 逐字节复刻 PS（并遵循 Gate D 的 canonicalization：仅剔除不稳定 header）
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
5. **发布与 CI**：增加 wheel 构建与 smoke test（至少验证 `import rnaview` + 跑一个小 case）。

## 4. 关键风险与决策点（现在就要想清楚）

- **2D byte‑exact 的工作量**：PS/XML 都是“格式敏感”的 legacy 输出，任何浮点/空白/排序差异都会导致 diff；建议把 2D renderer 当成一个独立子项目推进。
- **SVG 的可信口径**：当前 Gate D 的 SVG 是从 PS 经 converter 生成；若要保证“SVG 视觉等价于 PS”，需要额外的视觉/光栅化 sanity（不建议一开始就做成 hard gate）。
- **跨平台 wheel**：如果 v0 先只支持 Linux（CI 容器），文档里要写清楚；后续再补 macOS/Windows wheels。
