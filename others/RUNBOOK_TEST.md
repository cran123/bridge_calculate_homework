# Test Runbook

Open PowerShell in this repository root.

## 1. Prepare Environment

```powershell
$env:ONEAPI_ROOT = "D:\Program Files (x86)\Intel\oneAPI"
$env:MKLROOT = "D:\Program Files (x86)\Intel\oneAPI\mkl\latest"
$env:Path = "$env:MKLROOT\bin;$env:ONEAPI_ROOT\compiler\latest\bin;$env:Path"
$PY = "$env:USERPROFILE\miniconda3\python.exe"

$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
```

For a quick local test, set both thread counts to `4`.

## 2. Build

```powershell
cmake -S make -B build_pardiso -DMKL_DIR="$env:MKLROOT\lib\cmake\mkl"
cmake --build build_pardiso --config Release
```

## 3. Convert INP To DAT

```powershell
& $PY tools\abaqus_to_stappp.py data\data-3\Bridge-1.inp -o data\data-3\Bridge-1.rcm.generated.dat --renumber rcm --mode 1
```

If gravity must be supplied manually:

```powershell
& $PY tools\abaqus_to_stappp.py data\data-3\Bridge-1.inp -o data\data-3\Bridge-1.rcm.generated.dat --renumber rcm --mode 1 --gravity 0 0 -9.81
```

## 4. Solve DAT

```powershell
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "0"

build_pardiso\Release\stap++.exe data\data-3\Bridge-1.rcm.generated.dat
```

Use Task Manager for memory peak. The program prints only the competition solve time:

```text
COMPETITION SOLVE WALL TIME AFTER INPUT BEFORE OUTPUT, SEC = ...
```

Do not export CSV during the official timing run. Keep it as a post-processing step.

If memory is too high, switch to PARDISO out-of-core before another run:

```powershell
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "32000"
$env:MKL_PARDISO_OOC_PATH = "."
```

The console should print `ooc=1`. This usually reduces RAM pressure but costs time.

For a deliberate low-memory Bridge-1/Bridge-2 run, use:

```powershell
$env:MKL_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "300"
$env:MKL_PARDISO_OOC_PATH = "data\conversion_timing\pardiso_ooc"

build_pardiso\Release\stap++.exe data\data-3\Bridge-1.rcm.generated.dat
build_pardiso\Release\stap++.exe data\data-3\Bridge-2.eigen_mpc.dat
```

If 1 thread is too slow, use the middle low-memory setting:

```powershell
$env:MKL_NUM_THREADS = "4"
$env:OMP_NUM_THREADS = "4"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "500"
```

Measured local low-memory notes:

```text
Bridge-1: 12T in-core 60.3 MB, 1T OOC 100MB 56.2 MB
Bridge-2: 12T in-core 831.6 MB, 1T OOC 300MB 559.5 MB
Bridge-3: 12T in-core solve time 109.505 s; OOC caps 4000/6000/7000MB failed, 8000MB worked
```

Bridge-3 low-memory command:

```powershell
$env:MKL_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "8000"
$env:MKL_PARDISO_OOC_PATH = "data\conversion_timing\pardiso_ooc"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-3.fast.rcm.dat
```

Measured local Bridge-3:

```text
4T in-core:       solve 113.194 s, peak 8718.5 MB, passed
1T OOC 8000MB:   solve 231.832 s, peak 8298.7 MB, passed
12T in-core:     solve 109.505 s passed earlier, but a later monitored run failed with PARDISO error -2
```

So Bridge-3 OOC 8000MB is mainly a fallback. It saves only about 420MB versus the successful 4T in-core run, while taking about twice as long.

## 5. Optional VTU Post-Processing

```powershell
& $PY tools\export_vtu_displacements.py data\data-3\Bridge-1.rcm.generated_lc1.vtu --csv data\data-3\Bridge-1.node_displacements.csv
& $PY tools\export_vtu_displacements.py data\data-3\Bridge-1.rcm.generated_lc1.vtu --stappp-out data\data-3\Bridge-1.displacements.out
```

CSV columns:

```text
node,x,y,z,ux,uy,uz,u_mag
```

The `--stappp-out` file contains only the STAP++-style `D I S P L A C E M E N T S` table.
