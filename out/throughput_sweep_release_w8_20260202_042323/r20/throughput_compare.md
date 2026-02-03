# throughput compare (workers=8 repeat=20)

- jobs: 14 * 20 = 280
- speedup (rust/legacy): 9.300

| engine | throughput (jobs/s) | elapsed_ms | wall_ms | overhead_ms | overhead/job_ms | median_ms | p95_ms | ok | failed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 3.857 | 72602 | 73005 | 403 | 1.439 | 264 | 21496 | 280 | 0 |
| rust | 35.865 | 7807 | 8191 | 384 | 1.371 | 100 | 1322 | 280 | 0 |
