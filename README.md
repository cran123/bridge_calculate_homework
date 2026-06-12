# STAP++ Edit

This repository builds the course-design STAP++ solver, converts Abaqus `inp` bridge models to STAP++ `dat`, writes VTU files for ParaView, and can export nodal displacement results from VTU to CSV.

## Environment

Open PowerShell in this repository root, then run:

```powershell
$env:ONEAPI_ROOT = "D:\Program Files (x86)\Intel\oneAPI"
$env:MKLROOT = "D:\Program Files (x86)\Intel\oneAPI\mkl\latest"
$env:Path = "$env:MKLROOT\bin;$env:ONEAPI_ROOT\compiler\latest\bin;$env:Path"
$PY = "$env:USERPROFILE\miniconda3\python.exe"
```

For the official test, use 12 threads:

```powershell
$env:MKL_NUM_THREADS = "12"
$env:OMP_NUM_THREADS = "12"
$env:MKL_DYNAMIC = "FALSE"
```

For a local quick trial, set 4 threads instead:

```powershell
$env:MKL_NUM_THREADS = "4"
$env:OMP_NUM_THREADS = "4"
$env:MKL_DYNAMIC = "FALSE"
```

## Build

```powershell
cmake -S make -B build_pardiso -DMKL_DIR="$env:MKLROOT\lib\cmake\mkl"
cmake --build build_pardiso --config Release
```

PARDISO is active when the configure/build or solver output includes:

```text
Using Intel oneMKL PARDISO solver by default
Intel oneMKL PARDISO solver prepared
```

## Convert INP To DAT

Use RCM renumbering for bridge cases:

```powershell
& $PY tools\abaqus_to_stappp.py data\data-3\Bridge-1.inp -o data\data-3\Bridge-1.rcm.generated.dat --renumber rcm --mode 1
```

If the `inp` file has no explicit gravity load and you need to add one:

```powershell
& $PY tools\abaqus_to_stappp.py data\data-3\Bridge-1.inp -o data\data-3\Bridge-1.rcm.generated.dat --renumber rcm --mode 1 --gravity 0 0 -9.81
```

The converter prints model statistics such as node count, element count, boundary DOFs, tie constraints, renumbering, and gravity flag.

## Solve DAT

Use compact log output plus full VTU output for the official test:

```powershell
$env:STAPPP_FAST_OUTPUT = "1"
$env:STAPPP_VTK_OUTPUT = "1"
$env:STAPPP_SOLVER = "pardiso"
$env:STAPPP_PARDISO_OOC = "0"

build_pardiso\Release\stap++.exe data\data-3\Bridge-1.rcm.generated.dat
```

Expected outputs:

```text
data\data-3\Bridge-1.rcm.generated.out
data\data-3\Bridge-1.rcm.generated_lc1.vtu
```

The `.out` file is intentionally kept as a small log. It prints one timing metric:

```text
COMPETITION SOLVE WALL TIME AFTER INPUT BEFORE OUTPUT, SEC = ...
```

Memory consumption is read from Task Manager during the run.

Do not export CSV during the official timing run unless it is explicitly needed. CSV export is a separate post-processing step and can be done after the solver finishes.

## PARDISO Memory Fallback

Default PARDISO is in-core and fastest:

```powershell
$env:STAPPP_PARDISO_OOC = "0"
```

If Task Manager shows memory pressure, switch to oneMKL PARDISO out-of-core before running the solver:

```powershell
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "32000"
$env:MKL_PARDISO_OOC_PATH = "."
```

`STAPPP_PARDISO_OOC=1` trades speed for lower RAM pressure by allowing PARDISO to spill factorization data to disk. Use an SSD path with enough free space. The console prints `ooc=1` when this mode is active.

For an intentional low-memory run on Bridge-1 or Bridge-2, use fewer threads plus out-of-core:

```powershell
$env:MKL_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "12000"
$env:MKL_PARDISO_OOC_PATH = "."
```

If this is too slow, try the middle setting:

```powershell
$env:MKL_NUM_THREADS = "4"
$env:OMP_NUM_THREADS = "4"
$env:MKL_DYNAMIC = "FALSE"
$env:STAPPP_PARDISO_OOC = "1"
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE = "16000"
$env:MKL_PARDISO_OOC_PATH = "."
```

This is a time-for-space profile. It is useful for comparing memory scores, but it is usually slower than the official 12-thread in-core run.

## Post-Process VTU

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

You can also write both formats in one post-processing call:

```powershell
& $PY tools\export_vtu_displacements.py data\data-3\Bridge-1.rcm.generated_lc1.vtu --csv data\data-3\Bridge-1.node_displacements.csv --stappp-out data\data-3\Bridge-1.displacements.out
```

Open the `.vtu` in ParaView for displacement and stress cloud plots. Open the CSV if you need per-node numeric displacement values.
