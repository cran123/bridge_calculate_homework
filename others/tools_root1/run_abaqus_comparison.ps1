param(
    [string]$RepoRoot = $(Split-Path -Parent $PSScriptRoot),
    [switch]$SkipSubmit
)

$ErrorActionPreference = "Stop"

python (Join-Path $RepoRoot "tools\generate_abaqus_inp_from_stappp.py")

$extractor = Join-Path $RepoRoot "abaqus_scripts\extract_odb_results.py"

Get-ChildItem -Path (Join-Path $RepoRoot "cases") -Recurse -Filter "*_abaqus.inp" | ForEach-Object {
    $inp = $_
    $work = $inp.Directory.FullName
    $job = [System.IO.Path]::GetFileNameWithoutExtension($inp.Name)
    Push-Location $work
    try {
        if (-not $SkipSubmit) {
            Write-Host "Running Abaqus job $job"
            & abaqus job=$job input=$($inp.Name) interactive
            if ($LASTEXITCODE -ne 0) {
                Write-Warning "Abaqus job failed: $job"
                return
            }
        }
        $odb = Join-Path $work ($job + ".odb")
        if (Test-Path $odb) {
            & abaqus python $extractor --odb $odb
        }
    }
    finally {
        Pop-Location
    }
}

python (Join-Path $RepoRoot "abaqus_scripts\compare_stappp_abaqus.py") --repo-root $RepoRoot
