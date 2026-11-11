#!/usr/bin/env bash
# Full pipeline: build, test, identify, tune, evaluate, analyze. About 2 minutes on 2 cores.
set -euo pipefail
cd "$(dirname "$0")/.."
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
(cd build && ctest --output-on-failure)
mkdir -p results
./build/actlab_sim chirp_plant --out results          # 1. excite the plant
python3 scripts/identify.py --results results         # 2. fit two-mass model -> identified.txt
./build/actlab_sim lqr_info | tee results/lqr_info.txt
./build/actlab_sim tune --out results                 # 3. grid-search every controller family
./build/actlab_sim all --out results --mc 200         # 4. track, loop gain, sweep, Monte Carlo, adapt
./build/rt_loop_bench --seconds 10 --out results/rt_bench.csv | tee results/rt_bench.txt
./build/rt_loop_bench --seconds 10 --stress 2 --out results/rt_bench_stress.csv | tee results/rt_bench_stress.txt
python3 scripts/analyze.py --results results          # 5. figures + results/summary.md
