# throughput compare (workers=8 repeat=1)

- jobs: 14 * 1 = 14
- speedup (rust/legacy): 15.942

| engine | throughput (jobs/s) | elapsed_ms | wall_ms | overhead_ms | overhead/job_ms | median_ms | p95_ms | ok | failed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 0.734 | 19083 | 19298 | 215 | 15.357 | 314 | 6992 | 14 | 0 |
| rust | 11.696 | 1197 | 1379 | 182 | 13.000 | 100 | 508 | 14 | 0 |
