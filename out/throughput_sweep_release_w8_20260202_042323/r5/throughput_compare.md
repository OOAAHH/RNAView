# throughput compare (workers=8 repeat=5)

- jobs: 14 * 5 = 70
- speedup (rust/legacy): 8.940

| engine | throughput (jobs/s) | elapsed_ms | wall_ms | overhead_ms | overhead/job_ms | median_ms | p95_ms | ok | failed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 3.136 | 22323 | 22605 | 282 | 4.029 | 279 | 19609 | 70 | 0 |
| rust | 28.034 | 2497 | 2726 | 229 | 3.271 | 104 | 1424 | 70 | 0 |
