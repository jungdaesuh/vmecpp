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
      },
      {
        "commit": {
          "author": {
            "email": "166746189+jurasic-pf@users.noreply.github.com",
            "name": "Philipp Jura\u0161i\u0107",
            "username": "jurasic-pf"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0bd147f270a43f3fe52d0aac8049dce97d76668f",
          "message": "Show grid metadata and progress bars in vmecpp mgrid (#988)\n\n* Show grid metadata and progress bars in vmecpp mgrid\n\n* Use a single bar color in vmecpp mgrid",
          "timestamp": "2026-10-15T09:21:14+02:00",
          "tree_id": "c775056b24de62ef1bbca9548c69d7f677701c55",
          "url": "https://github.com/proximafusion/vmecpp/commit/0bd147f270a43f3fe52d0aac8049dce97d76668f"
        },
        "date": 1792049313407,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_cli_startup",
            "value": 0.667670561400007,
            "range": "stddev: 0.009982232933967793",
            "unit": "seconds",
            "extra": "rounds: 5"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_fixed_boundary_w7x",
            "value": 2.6134478723333245,
            "range": "stddev: 0.06167122028916781",
            "unit": "seconds",
            "extra": "rounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_fixed_boundary_cma",
            "value": 1.1063102596666379,
            "range": "stddev: 0.012198433448357394",
            "unit": "seconds",
            "extra": "rounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_fixed_boundary_cma_6x8",
            "value": 1.6723256190000104,
            "range": "stddev: 0.02199715101635629",
            "unit": "seconds",
            "extra": "rounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_response_table_from_coils",
            "value": 1.4783067523333482,
            "range": "stddev: 0.006700434580871672",
            "unit": "seconds",
            "extra": "rounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_free_boundary",
            "value": 6.772897599666673,
            "range": "stddev: 0.03707761759870263",
            "unit": "seconds",
            "extra": "rounds: 3"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_simsopt_adjoint_gradient",
            "value": 2.5429080330000033,
            "range": "stddev: 0",
            "unit": "seconds",
            "extra": "rounds: 1"
          },
          {
            "name": "benchmarks/test_benchmarks.py::test_bench_simsopt_finite_difference_gradient",
            "value": 0.5064535330000126,
            "range": "stddev: 0",
            "unit": "seconds",
            "extra": "rounds: 1"
          }
        ]
      }
    ],
    "C++ Microbenchmarks": [
      {
        "commit": {
          "author": {
            "email": "166746189+jurasic-pf@users.noreply.github.com",
            "name": "Philipp Jura\u0161i\u0107",
            "username": "jurasic-pf"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0bd147f270a43f3fe52d0aac8049dce97d76668f",
          "message": "Show grid metadata and progress bars in vmecpp mgrid (#988)\n\n* Show grid metadata and progress bars in vmecpp mgrid\n\n* Use a single bar color in vmecpp mgrid",
          "timestamp": "2026-10-15T09:21:14+02:00",
          "tree_id": "c775056b24de62ef1bbca9548c69d7f677701c55",
          "url": "https://github.com/proximafusion/vmecpp/commit/0bd147f270a43f3fe52d0aac8049dce97d76668f"
        },
        "date": 1792049708744,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "DeAliasConstraintForce/4x4",
            "value": 2.4547138709055486e-05,
            "unit": "seconds",
            "extra": "iterations: 48903\ncpu: 2.4545978978794765e-05 seconds\nthreads: 1"
          },
          {
            "name": "DeAliasConstraintForce/7x1",
            "value": 3.195886219593765e-05,
            "unit": "seconds",
            "extra": "iterations: 43777\ncpu: 3.195796642072321e-05 seconds\nthreads: 1"
          },
          {
            "name": "DeAliasConstraintForce/12x12",
            "value": 0.000559026050567627,
            "unit": "seconds",
            "extra": "iterations: 2500\ncpu: 0.0005590082787999998 seconds\nthreads: 1"
          },
          {
            "name": "DeAliasConstraintForce/16x18",
            "value": 0.0014787943445231544,
            "unit": "seconds",
            "extra": "iterations: 959\ncpu: 0.0014786980604796667 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalFourierToReal/4x4",
            "value": 0.00015896787641239758,
            "unit": "seconds",
            "extra": "iterations: 8698\ncpu: 0.00015896257323522647 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalForcesToFourier/4x4",
            "value": 0.00013900196146271387,
            "unit": "seconds",
            "extra": "iterations: 10071\ncpu: 0.00013899350044682755 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalFourierToReal/6x8",
            "value": 0.00032192375510931015,
            "unit": "seconds",
            "extra": "iterations: 4352\ncpu: 0.00032190436167279417 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalForcesToFourier/6x8",
            "value": 0.0002833452697229337,
            "unit": "seconds",
            "extra": "iterations: 4945\ncpu: 0.00028332736056622867 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalFourierToReal/12x12",
            "value": 0.0005220221442337406,
            "unit": "seconds",
            "extra": "iterations: 2683\ncpu: 0.000521996677226985 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalForcesToFourier/12x12",
            "value": 0.0004694178425027277,
            "unit": "seconds",
            "extra": "iterations: 3138\ncpu: 0.00045753651975780736 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalFourierToReal/12x13",
            "value": 0.0017673441426670968,
            "unit": "seconds",
            "extra": "iterations: 804\ncpu: 0.0017667460024875634 seconds\nthreads: 1"
          },
          {
            "name": "ToroidalForcesToFourier/12x13",
            "value": 0.0019351120311689575,
            "unit": "seconds",
            "extra": "iterations: 723\ncpu: 0.0019347344661134168 seconds\nthreads: 1"
          },
          {
            "name": "LaplaceSolve/5x4",
            "value": 3.289537765357036e-05,
            "unit": "seconds",
            "extra": "iterations: 42149\ncpu: 3.293676110939876e-05 seconds\nthreads: 1"
          },
          {
            "name": "LaplaceSolve/8x6",
            "value": 0.00017490860473001624,
            "unit": "seconds",
            "extra": "iterations: 8026\ncpu: 0.00017492462397208618 seconds\nthreads: 1"
          },
          {
            "name": "LaplaceSolve/12x8",
            "value": 0.0008571565735573862,
            "unit": "seconds",
            "extra": "iterations: 1632\ncpu: 0.0008571201948529372 seconds\nthreads: 1"
          },
          {
            "name": "LaplaceDecompose/5x4",
            "value": 2.73204261367592e-05,
            "unit": "seconds",
            "extra": "iterations: 50919\ncpu: 2.7343865099472986e-05 seconds\nthreads: 1"
          },
          {
            "name": "LaplaceDecompose/8x6",
            "value": 0.00015498681177162607,
            "unit": "seconds",
            "extra": "iterations: 9077\ncpu: 0.00015493948220778254 seconds\nthreads: 1"
          },
          {
            "name": "LaplaceDecompose/12x8",
            "value": 0.0007803843855928934,
            "unit": "seconds",
            "extra": "iterations: 1677\ncpu: 0.0007803847358377816 seconds\nthreads: 1"
          },
          {
            "name": "TransformGreensFunctionDerivative/5x4",
            "value": 0.0002994006061472375,
            "unit": "seconds",
            "extra": "iterations: 4676\ncpu: 0.00029939385479041956 seconds\nthreads: 1"
          },
          {
            "name": "TransformGreensFunctionDerivative/8x6",
            "value": 0.0011864792400745137,
            "unit": "seconds",
            "extra": "iterations: 1179\ncpu: 0.0011864380084817638 seconds\nthreads: 1"
          },
          {
            "name": "TransformGreensFunctionDerivative/12x8",
            "value": 0.005262451960628194,
            "unit": "seconds",
            "extra": "iterations: 266\ncpu: 0.005262131368421057 seconds\nthreads: 1"
          },
          {
            "name": "ComputeOutputQuantities/cma",
            "value": 0.005806990201464929,
            "unit": "seconds",
            "extra": "iterations: 287\ncpu: 0.005723704954703833 seconds\nthreads: 1"
          }
        ]
      }
    ]
  }
}
