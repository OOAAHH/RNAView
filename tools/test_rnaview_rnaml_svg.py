import importlib.util
import re
import sys
import unittest
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_rnaml_svg():
    root = _repo_root()
    mod_path = root / "tools" / "rnaview_rnaml_svg.py"
    spec = importlib.util.spec_from_file_location("rnaview_rnaml_svg", mod_path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


class TestRnamlSvg(unittest.TestCase):
    def test_multi_molecule_parsing_urx053(self) -> None:
        mod = _load_rnaml_svg()
        xml_path = _repo_root() / "test" / "pdb" / "urx053" / "urx053.pdb.xml"
        xml_text = xml_path.read_text(encoding="utf-8", errors="replace")

        bases_by_mol, bps = mod._parse_rnaml(xml_text)

        self.assertEqual(sorted(bases_by_mol.keys()), ["1", "2"])
        self.assertEqual(len(bases_by_mol["1"]), 158)
        self.assertEqual(len(bases_by_mol["2"]), 158)

        inter = [bp for bp in bps if bp.a.mol != bp.b.mol]
        self.assertEqual(len(inter), 4)

    def test_svg_output_contains_both_molecules(self) -> None:
        mod = _load_rnaml_svg()
        xml_path = _repo_root() / "test" / "pdb" / "urx053" / "urx053.pdb.xml"
        xml_text = xml_path.read_text(encoding="utf-8", errors="replace")
        svg = mod.svg_from_rnaml(xml_text)

        self.assertIn('id="backbone-mol-1"', svg)
        self.assertIn('id="backbone-mol-2"', svg)

        base_texts = len(re.findall(r'<text class="base"', svg))
        self.assertEqual(base_texts, 316)
