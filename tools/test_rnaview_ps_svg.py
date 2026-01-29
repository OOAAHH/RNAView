import importlib.util
import re
import sys
import unittest
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_ps_svg():
    root = _repo_root()
    mod_path = root / "tools" / "rnaview_ps_svg.py"
    spec = importlib.util.spec_from_file_location("rnaview_ps_svg", mod_path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


class TestPsSvg(unittest.TestCase):
    def test_viewbox_matches_bounding_box_urx053(self) -> None:
        mod = _load_ps_svg()
        ps_path = _repo_root() / "test" / "pdb" / "urx053" / "urx053.pdb.ps"
        ps_text = ps_path.read_text(encoding="utf-8", errors="replace")

        m = re.search(r"^%%BoundingBox:\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*$", ps_text, re.M)
        self.assertIsNotNone(m)
        llx, lly, urx, ury = (float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4)))

        svg = mod.svg_from_ps(ps_text)
        vm = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg)
        self.assertIsNotNone(vm)
        width = float(vm.group(1))
        height = float(vm.group(2))
        self.assertAlmostEqual(width, urx - llx, places=2)
        self.assertAlmostEqual(height, ury - lly, places=2)

        self.assertIn("<text", svg)
        self.assertIn("<line", svg)
