# Abaqus validation scripts

Run from the repository root with an Abaqus command prompt:

```powershell
abaqus cae noGUI=abaqus_scripts/run_validation_jobs.py -- --repo-root .
abaqus python abaqus_scripts/compare_stappp_abaqus.py --repo-root .
```

The first command creates and submits small reference jobs for bar, beam, plate,
hex, and mixed-element checks. The second command compares extracted Abaqus CSV
files against STAP++ output files when both are available.

Large bridge models can be checked with the existing converter and reference
CSV files under `data/data_ref`.
