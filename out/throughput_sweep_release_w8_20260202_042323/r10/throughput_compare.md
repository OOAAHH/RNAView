# throughput compare (workers=8 repeat=10)

- jobs: 14 * 10 = 140
- speedup (rust/legacy): 10.767

| engine | throughput (jobs/s) | elapsed_ms | wall_ms | overhead_ms | overhead/job_ms | median_ms | p95_ms | ok | failed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 3.182 | 43995 | 44295 | 300 | 2.143 | 244 | 21519 | 140 | 0 |
| rust | 34.263 | 4086 | 4364 | 278 | 1.986 | 96 | 1251 | 140 | 0 |
