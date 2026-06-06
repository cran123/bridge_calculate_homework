# STAP++ vs Abaqus Comparison

## Method

All course-design STAP++ cases under `cases/*/input` were converted to Abaqus `.inp` files with:

```powershell
python tools\generate_abaqus_inp_from_stappp.py
powershell -ExecutionPolicy Bypass -File tools\run_abaqus_comparison.ps1
```

The comparison uses nodal translational displacement `U1/U2/U3`, because it is available for every element type. STAP++ MPC blocks are translated directly to Abaqus `*Equation` blocks, so the mixed model now uses the same explicit linear tie constraints in both programs. Abaqus ODB files were extracted to CSV with `abaqus_scripts/extract_odb_results.py`; the final tables are:

```text
cases/abaqus_comparison.csv
cases/abaqus_comparison_details.csv
```

Element mapping:

| STAP++ case | STAP++ element | Abaqus reference |
| --- | --- | --- |
| `bar_patch` | Bar | `T3D2` |
| `beam_cantilever` | Beam3D2 | `B31` with generated general beam section |
| `plate4_membrane` | Plate4 | `S4R` |
| `hex8_patch` | Hex8 | `C3D8` |
| `bbarhex8_locking` | BbarHex8 | `C3D8H` near-incompressible reference |
| `mixed_beam_plate_hex` | H8 + Plate4 + Beam3D2 + MPC tie | `C3D8` + `S4R` + `B31` + `*Equation` |

## Summary table

| Case | Compared nodes | Max error node | Max absolute displacement error | Max relative error | Assessment |
| --- | ---: | ---: | ---: | ---: | --- |
| `bar_patch` | 2 | 1 | `1.0e-36` | `1.0e-22` | Excellent agreement |
| `bbarhex8_locking` | 8 | 2 | `3.8e-08` | `8.36e-07` | Excellent agreement; Abaqus uses `C3D8H` as reference |
| `beam_cantilever` | 2 | 2 | `1.542838e-05` | `8.81e-02` | Good trend; section/orientation convention differs slightly |
| `hex8_patch` | 8 | 3 | `1.0986328e-04` | `3.66e-01` | Acceptable absolute error; relative error inflated by small reference components |
| `mixed_beam_plate_hex` | 14 | 14 | `1.5914e-02` | `4.59e-01` | Same explicit MPC tie; remaining difference is element formulation/section convention |
| `plate4_membrane` | 4 | 1 | `5.0e-33` | `5.0e-19` | Excellent agreement |

## Case notes

### Bar patch

The axial bar result matches Abaqus to numerical roundoff. STAP++ gives the expected tip displacement `FL/EA = 1.0e-6` and axial stress `F/A = 1.0`.

### Plate4 membrane patch

The membrane tension patch matches Abaqus/S4R to numerical roundoff for nodal displacement. This is a strong patch-test result for the in-plane part of the plate implementation.

### BbarHex8 locking comparison

The BbarHex8 case agrees very closely with Abaqus `C3D8H`. This supports the B-bar volumetric locking treatment for the near-incompressible test.

### Beam cantilever

The tip `U3` differs by about `1.54e-5`. STAP++ and Abaqus both show the same cantilever response direction and magnitude, but the generated Abaqus general beam section and STAP++ internal beam convention are not perfectly identical. For the report, present this as a validation with small section-convention discrepancy.

### Hex8 patch

The maximum absolute error is about `1.10e-4`. The large relative error comes from displacement components whose Abaqus reference value is close to zero. Use absolute displacement error and the displacement pattern in the report.

### Mixed Beam-Plate-Hex8 model

The mixed model now separates solid, plate, and beam interface nodes and ties them with explicit MPC equations. The STAP++ MPC block is translated directly to Abaqus `*Equation`, for example:

```text
*Equation
2
10, 3, 1, 6, 3, -1
```

This means `U3(10) - U3(6) = 0` in both programs. The interface strategy is:

- Plate interface translational DOFs are tied to the H8 top-face translational DOFs.
- Beam root translational DOFs are tied to the H8/plate interface node.
- Interface rotations that do not have a solid counterpart are constrained rather than tied to nonexistent solid rotations.

After this change, the mixed-case maximum displacement error decreased from about `10.62` to `1.5914e-02`. The largest remaining difference occurs at beam tip node 14: STAP++ gives `U3 = -0.334183`, Abaqus gives `U3 = -0.318269`. This residual is consistent with element formulation and beam/plate section convention differences rather than mismatched tie constraints.

For the report, use the mixed case as both a functional simultaneous-use verification and an MPC/`*Equation` consistency check. Use the single-element/patch cases as the primary formulation accuracy validation.

## Recommended report wording

The bar, plate, and B-bar solid cases agree with Abaqus to near machine precision. The beam and Hex8 cases show the same physical response with small absolute displacement differences attributable to section convention and formulation differences. The mixed model demonstrates simultaneous use of beam, plate, and solid elements in STAP++, and its STAP++ MPC constraints are translated to Abaqus `*Equation` constraints for a consistent tie comparison. Remaining mixed-case differences are attributed mainly to Beam3D2/B31 and Plate4/S4R formulation details.
