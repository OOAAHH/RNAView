# throughput compare (workers=8 repeat=2)

- jobs: 14 * 2 = 28
- speedup (rust/legacy): 12.259

| engine | throughput (jobs/s) | elapsed_ms | wall_ms | overhead_ms | overhead/job_ms | median_ms | p95_ms | ok | failed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 1.434 | 19528 | 19729 | 201 | 7.179 | 248 | 12530 | 28 | 0 |
| rust | 17.577 | 1593 | 1779 | 186 | 6.643 | 93 | 953 | 28 | 0 |
