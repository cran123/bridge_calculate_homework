# Verification and Course-Design Results

## Generated course-design cases

The command below was used to generate `.out` and `.vtu` files:

```powershell
powershell -ExecutionPolicy Bypass -File tools\run_course_design_suite.ps1
python tools\summarize_course_results.py
```

Current generated results are summarized in `cases/result_summary.csv`.

| Case | Elements tested | Nodes | Elements | VTU files | Status |
| --- | --- | ---: | ---: | ---: | --- |
| `bar_patch` | Bar | 2 | 1 | 1 | OK |
| `beam_cantilever` | Beam3D2 | 2 | 1 | 1 | OK |
| `plate4_membrane` | Plate4 | 4 | 1 | 1 | OK |
| `hex8_patch` | Hex8 | 8 | 1 | 1 | OK |
| `bbarhex8_locking` | BbarHex8 | 8 | 1 | 1 | OK |
| `mixed_beam_plate_hex` | Hex8 + Plate4 + Beam3D2 + MPC tie | 14 | 3 | 1 | OK |

The mixed case confirms simultaneous use of several element groups and MPC/tie constraints in a single `.dat` file. It assembles one H8 solid, one Plate4 shell/plate, one Beam3D2 member, and explicit interface MPC equations that are translated to Abaqus `*Equation`.

## Patch and convergence checks

Patch experiments:

| Element | Patch case | Expected behavior |
| --- | --- | --- |
| Bar | `cases/bar_patch` | Constant axial strain/stress; tip displacement `FL/EA` |
| Beam3D2 | `cases/beam_cantilever` | Cantilever tip deflection and end-force recovery under transverse load |
| Plate4 | `cases/plate4_membrane` | Membrane tension patch with uniform in-plane response |
| Hex8 | `cases/hex8_patch` | Single cube under face load with compatible translational DOFs |
| BbarHex8 | `cases/bbarhex8_locking` | Near-incompressible cantilever-style case for locking discussion |

Existing convergence examples are also available under `tests/convergence`, especially the Plate4 cantilever gravity sequence `1x1`, `2x2`, and `4x4`. The existing validation suite in `others/tools/run_validation_suite.ps1` generates additional truss, H8, and BbarH8 convergence data.

## Abaqus comparison

Abaqus support is provided by `abaqus_scripts/`:

| Script | Role |
| --- | --- |
| `run_validation_jobs.py` | Builds and runs compact Abaqus reference jobs |
| `extract_odb_results.py` | Extracts nodal displacement and stress CSV files from ODB |
| `compare_stappp_abaqus.py` | Compares Abaqus nodal CSV files with STAP++ `.out` displacement tables |

The existing bridge reference comparison in `data/data_ref/bridge1_abaqus_vs_stappp_summary.txt` reports:

| Quantity | Value |
| --- | --- |
| Abaqus rows | 4163 |
| STAP nodes | 4091 |
| Exact/near coordinate matches | 3051 / 3051 |
| Match distance max/mean/rms | 0.000467 / 2.1784915e-05 / 7.9723596e-05 |
| Abaqus max UMAG | 0.55737034 |
| Matched STAP UMAG at that point | 0.55575202 |
| Vector error at Abaqus max UMAG | 0.0017907761 |

The larger component-wise differences in the same summary should be discussed as a bridge-model conversion/reference alignment issue rather than a single-element formulation error. For the report, use the compact Abaqus scripts for direct element-level validation and the bridge summary for large-model comparison.

The compact course-design Abaqus comparisons have also been run for all generated cases. The final nodal displacement comparison is stored in:

```text
cases/abaqus_comparison.csv
cases/abaqus_comparison_details.csv
doc/abaqus_comparison.md
```

| Case | Max absolute displacement error | Main conclusion |
| --- | ---: | --- |
| `bar_patch` | `1.0e-36` | Matches Abaqus to roundoff |
| `plate4_membrane` | `5.0e-33` | Matches Abaqus to roundoff |
| `bbarhex8_locking` | `3.8e-08` | Excellent agreement with Abaqus `C3D8H` |
| `beam_cantilever` | `1.542838e-05` | Same response trend; small beam-section convention difference |
| `hex8_patch` | `1.0986328e-04` | Acceptable absolute error; relative error inflated by near-zero components |
| `mixed_beam_plate_hex` | `1.5914e-02` | Same explicit MPC/`*Equation` tie; residual from element formulation/section conventions |

## ParaView outputs

Each course-design case has a ParaView-ready copy under:

```text
cases/<case>/paraview/*_lc1.vtu
```

Recommended ParaView operations:

1. Open the `.vtu` file.
2. Color by `DisplacementMagnitude`.
3. Apply `Warp By Vector` using `Displacement`.
4. For line elements, color by `AxialForce`, `BeamMy`, `BeamMz`, or `BeamTorque`.
5. For mixed models, color by `ElementType` first to confirm H8, Plate4, and Beam3D2 are all present.

## Scale and efficiency

Small validation cases solve in milliseconds and are intended for correctness checks. Bridge-1 and Bridge-2 provide scale data:

| Case | NEQ | NWK | Max half bandwidth | Input s | Stiffness s | Factorization s | Load solve s | Output s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Bridge-1 RCM | 15336 | 9161622 | 1554 | 0.430 | 0.049 | 7.533 | 0.054 | 0.471 |
| Bridge-2 MPC/Eigen | 121578 | 74046702 | 103609 | 7.782 | 2.088 | 119.931 | 0.301 | 3.261 |

The runtime is dominated by sparse/skyline factorization. Assembly and load-case back substitution are much smaller. Therefore the most effective efficiency improvements are bandwidth reduction, better ordering, more compact sparse storage, and robust sparse direct factorization for large mixed models.
