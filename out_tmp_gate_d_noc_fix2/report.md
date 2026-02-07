# Gate D report
- semantics: `legacy-v1`
- renderer_version: `dev`
- candidate_engine: `command`
- candidate_cmd: `['python3', 'tools/rnaview_render.py', 'render', '--backend', 'pairs-out-noc']`
- baseline manifest: `test/golden_render/manifest.json`
- allowlist: `test/gate_d_allowlist.yaml`
- counts: `{'ok': 1, 'changed': 0, 'unapproved': 14, 'failed': 0}`

## Unapproved

- unapproved: `test/mmcif/insertion_code/1EFW/1EFW.cif` job_id=`1EFW__26f53222` events=3
- unapproved: `test/mmcif/insertion_code/1VVJ/1VVJ.cif` job_id=`1VVJ__c7a40134` events=3
- unapproved: `test/mmcif/insertion_code/4ARC/4ARC.cif` job_id=`4ARC__ab9c7e43` events=3
- unapproved: `test/mmcif/nmr_structure/8if5/8if5.cif` job_id=`8if5__9ce4cf9f` events=3
- unapproved: `test/mmcif/other/6pom/6pom.cif` job_id=`6pom__82d3711c` events=3
- unapproved: `test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif` job_id=`3p4j-assembly1__b2610a38` events=2
- unapproved: `test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif` job_id=`434d-assembly1__b2a692ef` events=3
- unapproved: `test/mmcif/x-ray/4NMG/assembly-1/4nmg-assembly1.cif` job_id=`4nmg-assembly1__2f46a8f9` events=3
- unapproved: `test/pdb/pdb1nvy/pdb1nvy.ent` job_id=`pdb1nvy__e152e1d2` events=3
- unapproved: `test/pdb/pdb1nvy/pdb1nvy.pdb` job_id=`pdb1nvy__93035b1b` events=3
- unapproved: `test/pdb/test1/test1.pdb` job_id=`test1__704911ff` events=3
- unapproved: `test/pdb/tr0001/tr0001.pdb` job_id=`tr0001__561f5292` events=3
- unapproved: `test/pdb/url064/url064.pdb` job_id=`url064__bc7ec225` events=3
- unapproved: `test/pdb/urx053/urx053.pdb` job_id=`urx053__93caf9b1` events=3

