window.BENCHMARK_DATA = {
  "lastUpdate": 1775492820883,
  "repoUrl": "https://github.com/jungdaesuh/vmecpp",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "78460559+jungdaesuh@users.noreply.github.com",
            "name": "jungdaesuh",
            "username": "jungdaesuh"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c88cbf748da963c61da7229cda192a3ac7163090",
          "message": "Merge branch 'proximafusion:main' into main",
          "timestamp": "2026-04-07T01:22:40+09:00",
          "tree_id": "822a97757f797559878c6a25e283b92420d4c609",
          "url": "https://github.com/jungdaesuh/vmecpp/commit/c88cbf748da963c61da7229cda192a3ac7163090"
        },
        "date": 1775492819383,
        "tool": "pytest",
        "benches": [
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_cli_startup",
            "value": 2.979786109594375,
            "unit": "iter/sec",
            "range": "stddev: 0.002037023716779428",
            "extra": "mean: 335.59455719998823 msec\nrounds: 5"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_cli_invalid_input",
            "value": 2.9389003752281067,
            "unit": "iter/sec",
            "range": "stddev: 0.0026803002943520298",
            "extra": "mean: 340.2633204000267 msec\nrounds: 5"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_fixed_boundary_w7x",
            "value": 0.23265516346520992,
            "unit": "iter/sec",
            "range": "stddev: 0.03526371714174215",
            "extra": "mean: 4.298206775666661 sec\nrounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_fixed_boundary_cma",
            "value": 0.570947570177336,
            "unit": "iter/sec",
            "range": "stddev: 0.006086665743540889",
            "extra": "mean: 1.7514743073333345 sec\nrounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_response_table_from_coils",
            "value": 0.27865079474266613,
            "unit": "iter/sec",
            "range": "stddev: 0.009551108490183174",
            "extra": "mean: 3.588721147999953 sec\nrounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_free_boundary",
            "value": 0.10189962591012883,
            "unit": "iter/sec",
            "range": "stddev: 0.032480149235080635",
            "extra": "mean: 9.813578715999975 sec\nrounds: 3"
          }
        ]
      }
    ]
  }
}