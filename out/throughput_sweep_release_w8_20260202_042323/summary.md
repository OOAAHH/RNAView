# Throughput sweep (release rust)

- workers: 8
- engines: legacy=bin/rnaview, rust=rust/target/release/rnaview-hotcore

| repeat | jobs | speedup | legacy thr(j/s) | rust thr(j/s) | legacy overhead_ms | rust overhead_ms | legacy overhead/job_ms | rust overhead/job_ms |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 14 | 15.942 | 0.734 | 11.696 | 215 | 182 | 15.357 | 13.000 |
| 2 | 28 | 12.259 | 1.434 | 17.577 | 201 | 186 | 7.179 | 6.643 |
| 5 | 70 | 8.940 | 3.136 | 28.034 | 282 | 229 | 4.029 | 3.271 |
| 10 | 140 | 10.767 | 3.182 | 34.263 | 300 | 278 | 2.143 | 1.986 |
| 20 | 280 | 9.300 | 3.857 | 35.865 | 403 | 384 | 1.439 | 1.371 |
| 50 | 700 | 7.801 | 4.288 | 33.450 | 671 | 737 | 0.959 | 1.053 |
