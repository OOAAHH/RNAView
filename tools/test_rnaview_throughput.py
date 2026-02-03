import importlib.util
import sys
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


class TestThroughputHelpers(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        repo = _repo_root()
        cls.mod = _load_module("rnaview_throughput", repo / "tools" / "rnaview_throughput.py")

    def test_pct(self) -> None:
        pct = self.mod._pct
        self.assertIsNone(pct([], 50))
        self.assertEqual(pct([10], 0), 10)
        self.assertEqual(pct([10], 50), 10)
        self.assertEqual(pct([10], 100), 10)
        self.assertEqual(pct([10, 20], 0), 10)
        self.assertEqual(pct([10, 20], 100), 20)
        self.assertEqual(pct([10, 20], 50), 15)

    def test_latency_summary_has_expected_keys(self) -> None:
        s = self.mod._summarize_latencies([5, 10, 15, 20, 25])
        for k in ("n", "min_ms", "max_ms", "median_ms", "mean_ms", "p90_ms", "p95_ms", "p99_ms"):
            self.assertIn(k, s)

