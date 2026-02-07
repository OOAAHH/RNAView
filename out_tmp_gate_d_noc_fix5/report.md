# Gate D report
- semantics: `legacy-v1`
- renderer_version: `dev`
- candidate_engine: `command`
- candidate_cmd: `['python3', 'tools/rnaview_render.py', 'render', '--backend', 'pairs-out-noc']`
- baseline manifest: `test/golden_render/manifest.json`
- allowlist: `test/gate_d_allowlist.yaml`
- counts: `{'ok': 10, 'changed': 0, 'unapproved': 5, 'failed': 0}`

## Unapproved

- unapproved: `test/mmcif/insertion_code/1EFW/1EFW.cif` job_id=`1EFW__26f53222` events=3
- unapproved: `test/mmcif/insertion_code/1VVJ/1VVJ.cif` job_id=`1VVJ__c7a40134` events=3
- unapproved: `test/mmcif/other/6pom/6pom.cif` job_id=`6pom__82d3711c` events=3
- unapproved: `test/pdb/tr0001/tr0001.pdb` job_id=`tr0001__561f5292` events=3
- unapproved: `test/pdb/urx053/urx053.pdb` job_id=`urx053__93caf9b1` events=3

