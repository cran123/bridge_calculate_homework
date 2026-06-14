# Bridge-3 sparse run notes

This package contains a modified STAP++ build for Bridge-3.

Key points:

- `Bridge-3.generated.dat` needs about 121.7 GiB with the original skyline matrix.
- `Bridge-3.rcm.generated.dat` uses RCM node numbering and direct Eigen sparse assembly.
- Sparse matrix check reported `NEQ = 809838` and `nonzeros = 55974804`.
- Run with `STAPPP_SPARSE_ASSEMBLY=1` to bypass skyline storage.

Server command:

```bash
cd ~/bridge3_stappp_sparse
bash run_bridge3_sparse.sh
```

Expected outputs:

```text
data/data-3/Bridge-3.rcm.generated.out
data/data-3/Bridge-3.rcm.generated_lc1.vtu
data/data-3/Bridge-3.rcm.generated.sparse.run.log
```
