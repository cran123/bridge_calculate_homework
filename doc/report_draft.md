# STAP++ Course Design Report Draft

## 1. Introduction

This project extends and validates a C++ finite-element solver compatible with STAP90-style input. The implemented functions include static linear elastic analysis, multiple load cases, concentrated loads, optional gravity loads, prescribed displacement penalties, multi-point constraint support, skyline/global sparse assembly, text result output, and ParaView `.vtu` post-processing.

Implemented element capabilities:

| Element | Function |
| --- | --- |
| Bar | 3D two-node axial bar/truss |
| Beam3D2 | 3D two-node Euler-Bernoulli beam with axial, bending, and torsional stiffness |
| Plate4 | Four-node engineering plate/shell element |
| Hex8 | Eight-node solid element |
| BbarHex8 | Eight-node B-bar solid for reduced volumetric locking |

## 2. Algorithm summary

The solver follows the standard displacement finite-element procedure:

1. Discretize the structure into element groups.
2. Compute each element stiffness matrix in local/global DOF order.
3. Assemble the global equilibrium equation `K u = f`.
4. Apply displacement constraints and MPC treatment.
5. Factorize the global matrix and solve each load case.
6. Recover nodal displacement and element response quantities.
7. Export `.out` and `.vtu` files for checking and visualization.

The current implementation uses skyline-oriented storage and an LDLT/sparse direct solve path. For large bridge models, factorization dominates the elapsed time.

## 3. Implementation plan and structure

Core classes:

| Class | Responsibility |
| --- | --- |
| `CDomain` | Global model data, input parsing, equation numbering, assembly |
| `CElementGroup` | Element/material group allocation and dispatch |
| `CElement` | Base interface for element read/write/stiffness/stress/body force |
| `CNode` | Coordinates and boundary/equation codes |
| `CMaterial` | Base material interface and derived material/section data |
| `CLDLTSolver` | Global equation factorization and back substitution |
| `COutputter` | Text output |
| `CVTKOutputter` | ParaView `.vtu` output |

See `doc/program_structure.md` for full input format, element codes, and output fields.

## 4. Usage

Build the project with CMake or Visual Studio, then run:

```powershell
build_eigen_run\Release\stap++.exe cases\mixed_beam_plate_hex\input\mixed_beam_plate_hex.dat
```

For the packaged course-design examples:

```powershell
powershell -ExecutionPolicy Bypass -File tools\run_course_design_suite.ps1
python tools\summarize_course_results.py
```

Outputs are written to:

```text
cases/<case>/stappp/*.out
cases/<case>/paraview/*.vtu
cases/result_summary.csv
```

## 5. Verification and validation

Verification includes:

| Test type | Covered elements | Files |
| --- | --- | --- |
| Patch tests | Bar, Plate4, Hex8 | `cases/bar_patch`, `cases/plate4_membrane`, `cases/hex8_patch` |
| Beam benchmark | Beam3D2 | `cases/beam_cantilever` |
| Locking comparison | BbarHex8 | `cases/bbarhex8_locking` |
| Mixed use | Hex8 + Plate4 + Beam3D2 | `cases/mixed_beam_plate_hex` |
| Existing convergence | Plate4 and generated suite | `tests/convergence`, `others/tools/run_validation_suite.ps1` |
| Abaqus reference | T3D2, B31, S4R, C3D8 | `abaqus_scripts/` |

The mixed case verifies that the program can simultaneously read, assemble, solve, and export beam, plate, and hexahedral solid elements in one model.

ParaView verification uses the generated `.vtu` files. Use `DisplacementMagnitude` for contours and `Displacement` for deformation. Use `ElementType` to distinguish mixed element groups.

## 6. Scale and efficiency

Small validation cases solve in less than 0.02 s on the current machine. Bridge cases demonstrate the larger scale:

| Case | Equations | Matrix storage items | Factorization time |
| --- | ---: | ---: | ---: |
| Bridge-1 RCM | 15336 | 9161622 | 7.533 s |
| Bridge-2 MPC/Eigen | 121578 | 74046702 | 119.931 s |

The main performance bottleneck is matrix factorization. Assembly, load vector construction, and output are secondary for these examples. For competition scoring, recommended improvement directions are node reordering, reducing half bandwidth, replacing skyline storage with a modern sparse format, and supporting iterative/preconditioned solves for very large models.

## 7. Team work

Fill this section with the final team allocation. A suggested split:

| Member | Work |
| --- | --- |
| Member A | Beam3D2 formulation, tests, Abaqus B31 comparison |
| Member B | Plate4 formulation, patch/convergence tests, report figures |
| Member C | Hex8/BbarHex8 formulation, large bridge tests, efficiency discussion |
| Member D | Input conversion, ParaView output, documentation and integration |

## 8. Conclusion

The program now supports bar, beam, plate, Hex8 solid, and BbarHex8 solid elements, including simultaneous mixed-element analysis. The generated examples provide text results and ParaView visualization files. Abaqus scripts and existing bridge comparison data provide independent validation paths. The current solver is adequate for course-scale and bridge-scale examples, while factorization and bandwidth control remain the main opportunities for further efficiency improvements.

## References

1. K. J. Bathe, *Finite Element Procedures*, Prentice Hall, 1996.
2. T. J. R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*, Dover, 2000.
3. O. C. Zienkiewicz, R. L. Taylor, and J. Z. Zhu, *The Finite Element Method: Its Basis and Fundamentals*, 7th ed., Butterworth-Heinemann, 2013.
4. Dassault Systemes, *Abaqus Analysis User's Guide*, SIMULIA.
5. The Visualization Toolkit documentation, VTK UnstructuredGrid XML file format.
