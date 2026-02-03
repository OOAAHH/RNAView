# RNAVIEW throughput compare (release rust)

- workers: 8
- repeat: 20 (jobs = input_count * repeat = 14 * 20 = 280)
- rust binary: rust/target/release/rnaview-hotcore
- legacy binary: bin/rnaview
- speedup (rust/legacy): 8.391

| engine | throughput (jobs/s) | executed_jobs | elapsed_ms | wall_ms | overhead_ms | median_ms | p95_ms | ok | failed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 3.939 | 280 | 71089 | 71626 | 537 | 276 | 20775 | 280 | 0 |
| rust | 33.050 | 280 | 8472 | 8973 | 501 | 110 | 1456 | 280 | 0 |

## Artifacts

- throughput_compare.json: out/throughput_compare_release_w8_r20_20260202_042055/throughput_compare.json
- legacy summary: /workspaces/RNAView/out/throughput_compare_release_w8_r20_20260202_042055/legacy/summary.json
- rust summary: /workspaces/RNAView/out/throughput_compare_release_w8_r20_20260202_042055/rust/summary.json
- legacy log: /workspaces/RNAView/out/throughput_compare_release_w8_r20_20260202_042055/legacy/batch.log
- rust log: /workspaces/RNAView/out/throughput_compare_release_w8_r20_20260202_042055/rust/batch.log
