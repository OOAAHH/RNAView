# RNAVIEW throughput compare

- workers: 8  
- repeat: 20  
- input_count (unique): 14  
- staged inputs: ['/workspaces/RNAView/out/throughput_compare_w8_r20_20260202_024518/_inputs']  
- speedup (rust/legacy): 1.389

| engine | throughput (jobs/s) | executed_jobs | elapsed_ms | ok | failed | median_ms | p95_ms |
|---|---:|---:|---:|---:|---:|---:|---:|
| legacy | 3.953 | 280 | 70833 | 280 | 0 | 263 | 20100 |
| rust | 5.491 | 280 | 50988 | 280 | 0 | 175 | 14355 |

## Artifacts

- legacy out dir: /workspaces/RNAView/out/throughput_compare_w8_r20_20260202_024518/legacy/summary.json
- rust out dir: /workspaces/RNAView/out/throughput_compare_w8_r20_20260202_024518/rust/summary.json
- legacy log: /workspaces/RNAView/out/throughput_compare_w8_r20_20260202_024518/legacy/batch.log
- rust log: /workspaces/RNAView/out/throughput_compare_w8_r20_20260202_024518/rust/batch.log
