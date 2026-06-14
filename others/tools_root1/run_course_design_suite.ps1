param(
    [string]$RepoRoot = $(Split-Path -Parent $PSScriptRoot),
    [string]$Exe = "",
    [switch]$IncludeBridge,
    [switch]$RunAbaqus
)

$ErrorActionPreference = "Stop"

if (-not $Exe) {
    $candidates = @(
        "build_eigen_run\Release\stap++.exe",
        "build-merge\Release\stap++.exe",
        "build\Debug\stap++.exe"
    ) | ForEach-Object { Join-Path $RepoRoot $_ }
    $Exe = $candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
}

if (-not $Exe -or -not (Test-Path $Exe)) {
    throw "Cannot find stap++.exe. Build the project or pass -Exe <path>."
}

function Copy-InputCase {
    param([string]$CaseDir)
    $inputDir = Join-Path $CaseDir "input"
    $stapDir = Join-Path $CaseDir "stappp"
    $paraDir = Join-Path $CaseDir "paraview"
    New-Item -ItemType Directory -Force -Path $stapDir | Out-Null
    New-Item -ItemType Directory -Force -Path $paraDir | Out-Null
    Get-ChildItem -Path $inputDir -Filter *.dat | ForEach-Object {
        $target = Join-Path $stapDir $_.Name
        Copy-Item -Force -Path $_.FullName -Destination $target
        $target
    }
}

function Invoke-StapCase {
    param([string]$InputPath)
    $caseWork = Split-Path -Parent $InputPath
    Push-Location $caseWork
    try {
        & $Exe $InputPath | Out-Null
        $exit = $LASTEXITCODE
    }
    finally {
        Pop-Location
    }
    if ($exit -ne 0) {
        return [pscustomobject]@{ Input=$InputPath; Status="FAILED"; Out=""; Vtu="" }
    }

    $base = [System.IO.Path]::GetFileNameWithoutExtension($InputPath)
    $out = Join-Path $caseWork ($base + ".out")
    $vtus = Get-ChildItem -Path $caseWork -Filter ($base + "_lc*.vtu") -ErrorAction SilentlyContinue
    $paraDir = Join-Path (Split-Path -Parent $caseWork) "paraview"
    foreach ($vtu in $vtus) {
        Copy-Item -Force -Path $vtu.FullName -Destination (Join-Path $paraDir $vtu.Name)
    }
    return [pscustomobject]@{
        Input = $InputPath
        Status = "OK"
        Out = $out
        Vtu = ($vtus | ForEach-Object { $_.FullName }) -join ";"
    }
}

$caseRoot = Join-Path $RepoRoot "cases"
$caseDirs = Get-ChildItem -Path $caseRoot -Directory | Sort-Object Name
$results = New-Object System.Collections.Generic.List[object]

foreach ($case in $caseDirs) {
    $inputDir = Join-Path $case.FullName "input"
    if (-not (Test-Path $inputDir)) { continue }
    foreach ($input in Copy-InputCase -CaseDir $case.FullName) {
        Write-Host "Running $input"
        $results.Add((Invoke-StapCase -InputPath $input))
    }
}

if ($IncludeBridge) {
    $bridgeInputs = @(
        "data\data-3\Bridge-1.rcm.generated.dat",
        "data\data-3\Bridge-2.eigen_mpc.dat"
    ) | ForEach-Object { Join-Path $RepoRoot $_ } | Where-Object { Test-Path $_ }
    foreach ($input in $bridgeInputs) {
        Write-Host "Running bridge case $input"
        $results.Add((Invoke-StapCase -InputPath $input))
    }
}

$summary = Join-Path $caseRoot "run_summary.csv"
$results | Export-Csv -NoTypeInformation -Encoding UTF8 -Path $summary
Write-Host "Wrote $summary"

if ($RunAbaqus) {
    $script = Join-Path $RepoRoot "abaqus_scripts\run_validation_jobs.py"
    & abaqus cae noGUI=$script -- --repo-root $RepoRoot
}
