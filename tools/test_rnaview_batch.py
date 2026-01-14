import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


class TestBatchRunner(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo = _repo_root()
        cls.batch_mod = _load_module("rnaview_batch", cls.repo / "tools" / "rnaview_batch.py")

    def test_run_single_with_regress(self) -> None:
        if not (self.repo / "bin" / "rnaview").exists():
            self.skipTest("missing legacy binary: bin/rnaview")

        inp = self.repo / "test" / "pdb" / "tr0001" / "tr0001.pdb"
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td)
            code = self.batch_mod.main(
                [
                    "run",
                    str(inp),
                    "--out-dir",
                    str(out_dir),
                    "--workers",
                    "1",
                    "--regress",
                    "--regress-mode",
                    "out",
                ]
            )
            self.assertEqual(code, 0)
            summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["counts"]["failed"], 0)
            self.assertEqual(summary["counts"]["regress_failed"], 0)
            self.assertEqual(len(summary["results"]), 1)
            res = summary["results"][0]
            job_dir = Path(res["job_dir"])
            self.assertTrue((job_dir / "pairs.json").exists())
            self.assertTrue((job_dir / "legacy.out").exists())

    def test_regress_manifest_supports_explicit_input(self) -> None:
        inp = self.repo / "test" / "pdb" / "tr0001" / "tr0001.pdb"
        core_json = self.repo / "test" / "golden_core" / "pdb" / "tr0001" / "tr0001.pdb.core.json"

        manifest = {
            "schema_version": 1,
            "entries": [
                {
                    "input": "test/pdb/tr0001/tr0001.pdb",
                    "core_json": "test/golden_core/pdb/tr0001/tr0001.pdb.core.json",
                }
            ],
        }

        with tempfile.TemporaryDirectory() as td:
            mp = Path(td) / "manifest.json"
            mp.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

            idx = self.batch_mod._build_regress_index(mp)
            golden = self.batch_mod._lookup_golden_entry(inp.resolve(), idx)
            self.assertIsNotNone(golden)
            assert golden is not None
            self.assertEqual(golden.core_path.resolve(), core_json.resolve())
            self.assertIsNone(golden.out_path)
