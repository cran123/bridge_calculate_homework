#!/usr/bin/env bash
set -euo pipefail

cd "$HOME/bridge3_stappp_sparse"

rm -f \
  data/data-3/Bridge-3.rcm.refrot.norot.generated.out \
  data/data-3/Bridge-3.rcm.refrot.norot.generated_lc1.vtu \
  data/data-3/Bridge-3.rcm.refrot.norot.generated.sparse.run.log \
  data/data-3/Bridge-3.rcm.refrot.norot.generated.time.log \
  server_bridge3_driver.log

nohup bash -lc '
  cd "$HOME/bridge3_stappp_sparse"
  export STAPPP_SPARSE_ASSEMBLY=1
  unset STAPPP_SOLVER
  export STAPPP_CHECK_RESIDUAL=1
  /usr/bin/time -v ./build_sparse/stap++ ./data/data-3/Bridge-3.rcm.refrot.norot.generated.dat \
    > ./data/data-3/Bridge-3.rcm.refrot.norot.generated.sparse.run.log \
    2> ./data/data-3/Bridge-3.rcm.refrot.norot.generated.time.log
' > server_bridge3_driver.log 2>&1 &

echo $!
