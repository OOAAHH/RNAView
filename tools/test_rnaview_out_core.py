import contextlib
import importlib.util
import io
import sys
import unittest
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_rnaview_out_core():
    root = _repo_root()
    mod_path = root / "tools" / "rnaview_out_core.py"
    spec = importlib.util.spec_from_file_location("rnaview_out_core", mod_path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


class TestRnaViewOutCore(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _repo_root()
        cls.mod = _load_rnaview_out_core()

    def test_extract_core_tr0001(self) -> None:
        out_path = self.root / "test" / "pdb" / "tr0001" / "tr0001.pdb.out"
        core = self.mod.extract_core(out_path)

        self.assertEqual(core["stats"]["total_pairs"], 30)
        self.assertEqual(core["stats"]["total_bases"], 76)
        self.assertIn("pair_type_counts", core["stats"])
        self.assertGreater(len(core["stats"]["pair_type_counts"]), 0)

        # A canonical stacked record should parse.
        self.assertIn(
            {"i": 30, "j": 31, "chain_i": "A", "resseq_i": 30, "base_i": "G", "base_j": "A", "resseq_j": 31, "chain_j": "A", "kind": "stacked"},
            core["base_pairs"],
        )

        # A syn + tran record should parse (lw/orientation/syn/note).
        self.assertIn(
            {
                "i": 9,
                "j": 23,
                "chain_i": "A",
                "resseq_i": 9,
                "base_i": "A",
                "base_j": "A",
                "resseq_j": 23,
                "chain_j": "A",
                "kind": "pair",
                "lw": "H/H",
                "orientation": "tran",
                "syn": 1,
                "note": "II",
            },
            core["base_pairs"],
        )

        # Multiplets should preserve indices.
        self.assertEqual(
            [m["indices"] for m in core["multiplets"]],
            [[9, 12, 23], [13, 22, 46]],
        )

    def test_compare_same_file_exit_zero(self) -> None:
        out_path = self.root / "test" / "pdb" / "tr0001" / "tr0001.pdb.out"
        code = self.mod.main(["compare", str(out_path), str(out_path)])
        self.assertEqual(code, 0)

    def test_compare_different_files_exit_one(self) -> None:
        left = self.root / "test" / "pdb" / "tr0001" / "tr0001.pdb.out"
        right = self.root / "test" / "pdb" / "tr0001" / "tr0001.pdb_sort.out"
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            code = self.mod.main(["compare", str(left), str(right), "--max-diffs", "3"])
        self.assertEqual(code, 1)
        self.assertIn("/base_pairs", stderr.getvalue())

    def test_extract_known_outputs_no_unknown_base_pairs(self) -> None:
        candidates: list[Path] = []
        for out_file in (self.root / "test").rglob("*.out"):
            name = out_file.name
            if name.endswith(("_patt.out", "_sort.out", ".out.out")):
                continue
            with out_file.open("r", encoding="utf-8", errors="replace") as handle:
                if not any(line.strip() == "BEGIN_base-pair" for line in handle):
                    continue
            candidates.append(out_file)

        self.assertGreater(len(candidates), 0)

        for out_path in candidates:
            core = self.mod.extract_core(out_path)
            unknown = sum(1 for bp in core["base_pairs"] if bp.get("kind") == "unknown")
            self.assertEqual(unknown, 0, msg=f"unexpected unknown base-pair lines in {out_path}")
            self.assertIn("total_pairs", core["stats"], msg=f"missing total_pairs in {out_path}")
            self.assertIn("total_bases", core["stats"], msg=f"missing total_bases in {out_path}")
            self.assertIn("pair_type_counts", core["stats"], msg=f"missing pair_type_counts in {out_path}")
