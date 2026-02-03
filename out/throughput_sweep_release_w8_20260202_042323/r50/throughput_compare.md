# throughput compare (workers=8 repeat=50)

- jobs: 14 * 50 = 700
- speedup (rust/legacy): 7.801

| engine | throughput (jobs/s) | elapsed_ms | wall_ms | overhead_ms | overhead/job_ms | median_ms | p95_ms | ok | failed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 4.288 | 163253 | 163924 | 671 | 0.959 | 258 | 20740 | 700 | 0 |
| rust | 33.450 | 20927 | 21664 | 737 | 1.053 | 108 | 1337 | 700 | 0 |
