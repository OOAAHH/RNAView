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


class TestGateNa(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo = _repo_root()
        cls.mod = _load_module("rnaview_gate_na", cls.repo / "tools" / "rnaview_gate_na.py")

    def test_cases_manifest_filters_to_approved_inputs(self) -> None:
        manifest = self.repo / "test" / "science_dna_cases.json"
        approved = self.mod._approved_inputs_from_cases_manifest(self.repo, manifest)
        self.assertEqual(
            approved,
            {
                "test/pdb/pdb1nvy/pdb1nvy.pdb",
                "test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif",
                "test/mmcif/nmr_structure/8if5/8if5.cif",
                "test/hybrid/pdb1nvy_hybrid/pdb1nvy_dna_rna_hybrid.pdb",
            },
        )

    def test_allowlist_loads_ids(self) -> None:
        allowlist = self.repo / "test" / "gate_na_allowlist.yaml"
        self.assertEqual(self.mod._load_allowlist(allowlist), set())

        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "allowlist.json"
            path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "approved": [
                            {"id": "case|core|abc"},
                            {"id": "case|out|def"},
                        ],
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            self.assertEqual(self.mod._load_allowlist(path), {"case|core|abc", "case|out|def"})

    def test_read_golden_canonical_supports_raw_text(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            raw = Path(td) / "candidate.ps"
            raw.write_text("%%CreationDate: today\nline-1\n", encoding="utf-8")
            canon = self.mod._read_golden_canonical(fmt="ps", path=raw)
            self.assertEqual(canon.decode("utf-8"), "line-1\n")

    def test_golden_na_manifest_matches_approved_inputs(self) -> None:
        manifest = json.loads((self.repo / "test" / "golden_na" / "manifest.json").read_text(encoding="utf-8"))
        approved = self.mod._approved_inputs_from_cases_manifest(self.repo, self.repo / "test" / "science_dna_cases.json")
        inputs = {entry["input"] for entry in manifest["entries"]}
        self.assertEqual(inputs, approved)

    def test_dna_core_goldens_carry_residue_annotations(self) -> None:
        golden_paths = [
            self.repo / "test" / "golden_na" / "pdb" / "pdb1nvy" / "pdb1nvy.pdb.core.json",
            self.repo / "test" / "golden_na" / "mmcif" / "x-ray" / "3P4J" / "assembly-1" / "3p4j-assembly1.cif.core.json",
            self.repo / "test" / "golden_na" / "mmcif" / "nmr_structure" / "8if5" / "8if5.cif.core.json",
            self.repo / "test" / "golden_na" / "hybrid" / "pdb1nvy_hybrid" / "pdb1nvy_dna_rna_hybrid.pdb.core.json",
        ]
        for path in golden_paths:
            payload = json.loads(path.read_text(encoding="utf-8"))
            residues = payload.get("residues") or []
            self.assertTrue(residues, msg=str(path))
            self.assertTrue(all("polymer_class" in r for r in residues), msg=str(path))


if __name__ == "__main__":
    unittest.main()
