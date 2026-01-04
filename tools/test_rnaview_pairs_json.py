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


class TestPairsJson(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo = _repo_root()
        cls.pairs_mod = _load_module("rnaview_pairs_json", cls.repo / "tools" / "rnaview_pairs_json.py")
        cls.core_mod = _load_module("rnaview_out_core", cls.repo / "tools" / "rnaview_out_core.py")

    def test_pairs_json_out_writer_roundtrip(self) -> None:
        core_path = self.repo / "test" / "golden_core" / "pdb" / "tr0001" / "tr0001.pdb.core.json"
        golden_core = json.loads(core_path.read_text(encoding="utf-8"))
        pairs = self.pairs_mod.pairs_json_from_core(
            golden_core,
            source_path="test/pdb/tr0001/tr0001.pdb.out",
            source_format="out",
            options={},
        )
        out_text = self.pairs_mod.out_core_from_pairs_json(pairs)
        with tempfile.NamedTemporaryFile("w", suffix=".out", encoding="utf-8", delete=True) as tmp:
            tmp.write(out_text)
            tmp.flush()
            extracted = self.core_mod.extract_core(Path(tmp.name))
        self.assertEqual(extracted, golden_core)

    def test_validate_golden_manifest(self) -> None:
        manifest = self.repo / "test" / "golden_core" / "manifest.json"
        code = self.pairs_mod.main(["validate-golden", "--manifest", str(manifest)])
        self.assertEqual(code, 0)

