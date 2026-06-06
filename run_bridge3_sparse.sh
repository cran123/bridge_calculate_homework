#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

cmake -S ./make -B ./build_sparse -DCMAKE_BUILD_TYPE=Release -DSTAPPP_USE_EIGEN_SOLVER=ON
cmake --build ./build_sparse --config Release --target stap++ -j "$(nproc)"

export STAPPP_SPARSE_ASSEMBLY=1
unset STAPPP_SOLVER
export STAPPP_CHECK_RESIDUAL=1
./build_sparse/stap++ ./data/data-3/Bridge-3.rcm.refrot.norot.generated.dat \
  > ./data/data-3/Bridge-3.rcm.refrot.norot.generated.sparse.run.log 2>&1
