# STAP++ Program Structure and Implemented Elements

## Project layout

The repository is a C++ finite-element program derived from STAP++/STAP90 style input files.

| Path | Purpose |
| --- | --- |
| `src/h`, `src/cpp` | Core solver, domain model, element implementations, output modules |
| `make/CMakeLists.txt` | CMake project file |
| `data/` | Original and large bridge input/output files |
| `tests/` | Existing unit, patch, and convergence examples |
| `cases/` | Course-design examples, generated STAP++ outputs, and ParaView files |
| `tools/` | Local diagnostics and course-design batch scripts |
| `abaqus_scripts/` | Abaqus reference-job and CSV comparison scripts |
| `doc/` | Generated documentation and report material |

## Main solution flow

`src/cpp/main.cpp` accepts one input stem or `.dat` file:

```text
stap++ model.dat
```

The program writes:

```text
model.out
model_lc1.vtu
model_lc2.vtu
...
```

The main phases are:

1. `CDomain::ReadData`: read title, control line, nodes, loads, optional MPCs, element groups.
2. `CalculateEquationNumber`: assign global equation numbers to unconstrained DOFs.
3. `AllocateMatrices`: compute skyline column heights and allocate the global matrix/vector.
4. `AssembleStiffnessMatrix`: call every element stiffness routine and assemble to skyline storage.
5. `ApplyDisplacementPenalty`: impose prescribed displacement constraints by penalty terms.
6. `CLDLTSolver`: factor and solve each load case.
7. `COutputter` and `CVTKOutputter`: write text results and ParaView `.vtu` output.

## Implemented element groups

`ElementTypes` is defined in `src/h/ElementGroup.h`.

| Input code | Class | Abaqus analogue | Main DOFs | Status |
| --- | --- | --- | --- | --- |
| `1` | `CBar` | `T3D2` | `ux uy uz` | 3D two-node bar/truss |
| `4` | `CHex8` | `C3D8` / `C3D8R` | `ux uy uz` | 8-node solid |
| `5` | `CBeam3D2` | `B31` | `ux uy uz rx ry rz` | 3D Euler-Bernoulli beam |
| `6` | `CPlate4` | `S4R` | `ux uy uz rx ry rz` | 4-node engineering plate/shell |
| `8` | `CBbarHex8` | B-bar brick | `ux uy uz` | 8-node solid with B-bar volumetric strain treatment |

The global node has six DOFs. Elements override `GenerateLocationMatrix()` when they use only a subset of those DOFs. Solid and bar elements use translational DOFs only; beam and plate elements use translations and rotations.

## Input format

The common input layout is:

```text
title
NUMNP NUMEG NLCASE MODEX
node bcx bcy bcz bcrx bcry bcrz x y z
...
loadcase_id
nloads [gravity_flag gx gy gz [ndisp]]
node dof value
...
[prescribed_node prescribed_dof value]
...
[optional MPC block]
element_type NUME NUMMAT
material lines
element lines
...
```

Boundary code `1` means constrained and `0` means active. DOF order is:

```text
ux uy uz rx ry rz
```

Material formats:

| Element | Material line |
| --- | --- |
| Bar | `set E A rho` |
| Hex8 / BbarHex8 | `set E nu rho` |
| Beam3D2 | `set E nu A Iy Iz J rho` |
| Plate4 | `set E nu thickness rho` |

Element connectivity formats:

```text
Bar:
ele n1 n2 mat

Hex8 / BbarHex8:
ele n1 n2 n3 n4 n5 n6 n7 n8 mat

Beam3D2:
ele n1 n2 mat ref_vx ref_vy ref_vz

Plate4:
ele n1 n2 n3 n4 mat
```

## Output and ParaView fields

Each solved load case writes one `.vtu`. ParaView can directly open the files under `cases/*/paraview`.

Point fields:

| Field | Meaning |
| --- | --- |
| `Displacement` | `ux uy uz` |
| `Rotation` | `rx ry rz`, zero for elements/nodes without active rotations |
| `DisplacementMagnitude` | Euclidean norm of translational displacement |

Cell fields:

| Field | Meaning |
| --- | --- |
| `ElementType` | STAP++ element type code |
| `AxialForce` | Bar/beam axial force |
| `BeamMy`, `BeamMz`, `BeamTorque` | Beam internal force resultants |
| `SXX SYY SZZ SXY SYZ SXZ`, `VonMises` | Present as common fields; current small-case VTU output keeps most solid/plate stress fields as zero placeholders unless stress export is extended |
| `MX MY MXY` | Plate moment placeholders for uniform ParaView schema |

For report figures, use `Displacement` with ParaView `Warp By Vector`, color by `DisplacementMagnitude`, and optionally color beam/bar cells by `AxialForce` or beam moment fields.

## Batch workflow

Run all course-design STAP++ examples:

```powershell
powershell -ExecutionPolicy Bypass -File tools\run_course_design_suite.ps1
python tools\summarize_course_results.py
```

Optional Abaqus validation:

```powershell
abaqus cae noGUI=abaqus_scripts\run_validation_jobs.py -- --repo-root .
abaqus python abaqus_scripts\compare_stappp_abaqus.py --repo-root .
```

For the course-design comparison workflow, `tools/generate_abaqus_inp_from_stappp.py` parses the optional STAP++ MPC block and writes each constraint as an Abaqus `*Equation`. This keeps the STAP++ tie model and the Abaqus reference model using the same explicit linear constraints.
