# STAP++ Edit

## Competition Run Command

These commands run the solver with Intel oneMKL PARDISO, compact result output, and VTU output enabled.

PowerShell:

```powershell
cd C:\Users\huawei\Desktop\有限元\STAPpp_edit

$env:ONEAPI_ROOT = "D:\Program Files (x86)\Intel\oneAPI"
$env:MKLROOT = "D:\Program Files (x86)\Intel\oneAPI\mkl\latest"
$env:Path = "$env:MKLROOT\bin;$env:ONEAPI_ROOT\compiler\latest\bin;$env:Path"

$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:MKL_NUM_THREADS = "4"
$env:OMP_NUM_THREADS = "4"

build_pardiso\Release\stap++.exe data\data-3\Bridge-1.final_check.dat
```

For another input file, replace the last path, for example:

```powershell
build_pardiso\Release\stap++.exe data\data-3\Bridge-2.eigen_mpc.dat
```

Expected output:

- `.out`: result tables and timing log only; input echo tables are skipped.
- `.vtu`: generated for visualization.
- Console: short progress messages only.

PARDISO is active when the console prints:

```text
Intel oneMKL PARDISO solver prepared
```

## Build With PARDISO

If `build_pardiso` does not exist or the executable is outdated:

```powershell
cd C:\Users\huawei\Desktop\有限元\STAPpp_edit

$env:ONEAPI_ROOT = "D:\Program Files (x86)\Intel\oneAPI"
$env:MKLROOT = "D:\Program Files (x86)\Intel\oneAPI\mkl\latest"
$env:Path = "$env:MKLROOT\bin;$env:ONEAPI_ROOT\compiler\latest\bin;$env:Path"

cmake -S make -B build_pardiso -DMKL_DIR="$env:MKLROOT\lib\cmake\mkl"
cmake --build build_pardiso --config Release
```

On the i5-1135G7 machine, start with 4 threads. You can also test 8 threads:

```powershell
$env:MKL_NUM_THREADS = "8"
$env:OMP_NUM_THREADS = "8"
```
