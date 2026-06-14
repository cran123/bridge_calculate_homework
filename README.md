# STAP++ Edit

This repository builds the course-design STAP++ solver, converts Abaqus `inp` bridge models to STAP++ `dat`, writes VTU files for ParaView, and can export nodal displacement results from VTU after the timing run.

## 1. Environment

Open PowerShell in this repository root.

```powershell
$env:ONEAPI_ROOT = "D:\Program Files (x86)\Intel\oneAPI"
$env:MKLROOT = "D:\Program Files (x86)\Intel\oneAPI\mkl\latest"
$env:Path = "$env:MKLROOT\bin;$env:ONEAPI_ROOT\compiler\latest\bin;$env:Path"
$PY = "$env:USERPROFILE\miniconda3\python.exe"
```

On another machine or server, keep paths configurable:

```powershell
$PY = "python"
$OOC_DIR = "D:\pardiso_ooc"
$BRIDGE4_DAT = "data\data-3\Bridge-4.dat"
```

Use relative paths from the repository root whenever possible. Only `$OOC_DIR` and `$BRIDGE4_DAT` usually need to change on the test server.

Official 12-thread setting:

```powershell
$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
```

Common solver output setting:

```powershell
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
```

## 2. Build

```powershell
cmake -S make -B build_pardiso -DMKL_DIR="$env:MKLROOT\lib\cmake\mkl"
cmake --build build_pardiso --config Release
```

If the server uses a different build directory, keep the same convention and replace `build_pardiso` in the run commands with that directory.

Linux/bash equivalent, if needed:

```bash
export MKL_NUM_THREADS=12
export OMP_NUM_THREADS=12
export MKL_DYNAMIC=FALSE
export STAPPP_FAST_OUTPUT=1
export STAPPP_VTK_OUTPUT=1
export STAPPP_SOLVER=pardiso
export STAPPP_PARDISO_OOC=1
export MKL_PARDISO_OOC_MAX_CORE_SIZE=48000
export MKL_PARDISO_OOC_PATH=/scratch/pardiso_ooc

./build_pardiso/stap++ data/data-3/Bridge-4.dat
```

PARDISO is active when the console prints:

```text
Intel oneMKL PARDISO solver prepared
```

## 3. Convert INP To DAT

Use the original converter:

```powershell
& $PY tools\abaqus_to_stappp.py data\data-3\Bridge-1.inp -o data\data-3\Bridge-1.rcm.generated.dat --renumber rcm --mode 1
```

Use the faster converter for large files. It preserves the original converter and only adds cached tie processing:

```powershell
& $PY tools\abaqus_to_stappp_fast.py data\data-3\Bridge-3.inp -o data\conversion_timing\Bridge-3.fast.rcm.dat --renumber rcm --mode 1
```

If gravity must be supplied manually:

```powershell
& $PY tools\abaqus_to_stappp_fast.py data\data-3\Bridge-1.inp -o data\data-3\Bridge-1.rcm.generated.dat --renumber rcm --mode 1 --gravity 0 0 -9.81
```

The converter prints node count, element count, boundary DOFs, tie constraints, MPC equations, renumbering, and gravity flag.

## 4. Official Run Commands

The `.out` file is intentionally kept as a small log and prints only:

```text
COMPETITION SOLVE WALL TIME AFTER INPUT BEFORE OUTPUT, SEC = ...
```

Memory is read from Task Manager during the run. Do not export CSV during the official timing run.

### Bridge-1

Recommended speed run:

```powershell
$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "0"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-1.fast.rcm.dat
```

Low-memory run. Bridge-1 is already small, so this only saves a few MB:

```powershell
$env:MKL_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "100"
$env:MKL_PARDISO_OOC_PATH = "data\conversion_timing\pardiso_ooc"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-1.fast.rcm.dat
```

Measured locally:

```text
12T in-core:      solve 0.729 s, peak 60.3 MB
1T OOC 100MB:    solve 0.674 s, peak 56.2 MB
```

### Bridge-2

Recommended speed run:

```powershell
$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "0"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-2.fast.rcm.dat
```

Low-memory run:

```powershell
$env:MKL_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "300"
$env:MKL_PARDISO_OOC_PATH = "data\conversion_timing\pardiso_ooc"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-2.fast.rcm.dat
```

Measured locally:

```text
12T in-core:      solve 5.178 s, peak 831.6 MB
1T in-core:       solve 3.781 s, peak 822.1 MB
1T OOC 500MB:    solve 31.562 s, peak 708.0 MB
1T OOC 300MB:    solve 33.923 s, peak 559.5 MB
```

### Bridge-3

Recommended stable speed run:

```powershell
$env:MKL_NUM_THREADS = "4"
$env:OMP_NUM_THREADS = "4"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "0"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-3.fast.rcm.dat
```

12-thread speed attempt. This passed once locally with `109.505 s`, but a later monitored run failed with PARDISO error `-2`, so treat it as less stable on this machine:

```powershell
$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "0"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-3.fast.rcm.dat
```

Out-of-core fallback. Caps `4000/6000/7000MB` failed locally; `8000MB` worked:

```powershell
$env:MKL_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "8000"
$env:MKL_PARDISO_OOC_PATH = "data\conversion_timing\pardiso_ooc"

build_pardiso\Release\stap++.exe data\conversion_timing\Bridge-3.fast.rcm.dat
```

Measured locally:

```text
12T in-core:        solve 109.505 s, passed once
4T in-core:         solve 113.194 s, peak 8718.5 MB, passed
1T OOC 8000MB:     solve 231.832 s, peak 8298.7 MB, passed
```

Accuracy check at the Abaqus maximum-displacement node:

```text
Bridge-3 max Abaqus displacement node percent error = 3.36056%, pass < 5%
```

### Bridge-4

Bridge-4 must stay below 64GB. Use PARDISO out-of-core by default. Confirm the OOC path is on an SSD with enough free space before running. Recommended free space: at least `100GB`, preferably `200GB`.

Recommended official attempt:

```powershell
$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "48000"
$env:MKL_PARDISO_OOC_PATH = $OOC_DIR

build_pardiso\Release\stap++.exe $BRIDGE4_DAT
```

More conservative memory attempt:

```powershell
$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "32000"
$env:MKL_PARDISO_OOC_PATH = $OOC_DIR

build_pardiso\Release\stap++.exe $BRIDGE4_DAT
```

If the 12-thread attempts approach 64GB or fail, reduce threads but keep OOC:

```powershell
$env:MKL_NUM_THREADS = "8"
$env:OMP_NUM_THREADS = "8"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "32000"
$env:MKL_PARDISO_OOC_PATH = $OOC_DIR

build_pardiso\Release\stap++.exe $BRIDGE4_DAT
```

Last-resort low-memory attempt:

```powershell
$env:MKL_NUM_THREADS = "4"
$env:OMP_NUM_THREADS = "4"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "32000"
$env:MKL_PARDISO_OOC_PATH = $OOC_DIR

build_pardiso\Release\stap++.exe $BRIDGE4_DAT
```

Bridge-4 DAT sanity check from the current file:

```text
nodes = 1,910,327
free DOF estimate = 5,837,778
H8 = 1,760,000
Plate4 = 40,000
Beam = 5,420
Truss = 20
MPC = 432
```

Do not use `STAPPP_PARDISO_OOC=0` for Bridge-4 unless you have already confirmed it stays below 64GB.

## 5. Suggested Test Order

For Bridge-1/2/3, use the speed run first. If memory score matters and time still passes, use the low-memory run as another measured attempt.

For Bridge-4:

```text
1. 12T OOC 48000MB
2. 12T OOC 32000MB
3. 8T  OOC 32000MB
4. 4T  OOC 32000MB
```

The console must print `ooc=1` for Bridge-4.

## 6. Post-Process VTU

The VTU contains nodal coordinates and `Displacement = ux uy uz`. Export every node to CSV after the official run:

```powershell
& $PY tools\export_vtu_displacements.py data\data-3\Bridge-1.rcm.generated_lc1.vtu --csv data\data-3\Bridge-1.node_displacements.csv
```

CSV columns:

```text
node,x,y,z,ux,uy,uz,u_mag
```

To generate a STAP++-style displacement table from the VTU:

```powershell
& $PY tools\export_vtu_displacements.py data\data-3\Bridge-1.rcm.generated_lc1.vtu --stappp-out data\data-3\Bridge-1.displacements.out
```

Open the `.vtu` in ParaView for displacement and stress cloud plots. Open the CSV only if you need per-node numeric displacement values.
