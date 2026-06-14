param(
    [string]$RepoRoot = $(Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Continue'
$Exe = Join-Path $RepoRoot 'build\stap++.exe'
$WorkDir = Join-Path $RepoRoot 'build\validation_suite'
New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null

# =================== helpers ===================
function Write-TextFile { param([string]$Path, [string[]]$Lines)
    $Lines | Set-Content -Encoding ASCII -Path $Path
}

function Invoke-Stappp { param([string]$InputPath)
    & $Exe $InputPath | Out-Null
    if ($LASTEXITCODE -ne 0) { return $null }
    return ([System.IO.Path]::ChangeExtension($InputPath, '.out'))
}

function Get-Disp { param([string]$OutPath, [int]$Node, [int]$Comp)
    $pat = '^\s*' + $Node + '\s+([+-]?(?:\d+\.\d*|\d*\.\d+)[eE][+-]?\d+)'
    foreach ($l in Get-Content $OutPath) {
        if ($l -match $pat) {
            $v = ($l -replace '^\s+','') -split '\s+'
            return [double]$v[$Comp]
        }
    }
    return $null
}

function Get-Stress { param([string]$OutPath, [int]$Elem, [int]$Comp)
    $inStress = $false; $cnt = 0
    foreach ($l in Get-Content $OutPath) {
        if ($l -match 'S T R E S S|ELEMENT.*SXX|FORCE.*STRESS') { $inStress = $true; continue }
        if ($inStress -and $l -match '^\s*\d') {
            $cnt++
            if ($cnt -eq $Elem) {
                $nums = [regex]::Matches($l, '[+-]?\d+\.\d+[eE][+-]\d+')
                if ($Comp -lt $nums.Count) { return [double]$nums[$Comp].Value }
            }
        }
    }
    return $null
}

function Get-H8NodeId { param([int]$Ix,[int]$Iy,[int]$Iz,[int]$Nx,[int]$Ny,[int]$Nz)
    return 1 + $Ix + ($Nx+1)*($Iy + (($Ny+1)*$Iz))
}

# =================== Part 1: CONVERGENCE ===================
Write-Host "`n============================================================"
Write-Host "  PART 1 : CONVERGENCE ANALYSIS"
Write-Host "============================================================"

# -- Truss convergence --
function Write-TrussConv { param([int]$N,[string]$Path)
    $E=1000; $A=1; $Rho=1; $L=1; $dx=$L/$N
    $l=[System.Collections.Generic.List[string]]::new()
    $l.Add("Truss_N$N")
    $l.Add("$($N+1) 1 1 1")
    for ($i=0; $i -le $N; $i++) {
        $bc = if ($i -eq 0) {1} else {0}
        $l.Add(('{0} {1} 1 1 {2} 0 0' -f ($i+1), $bc, ($i*$dx)))
    }
    $l.Add('1')
    $l.Add('0 1 1 0 0')
    $l.Add("1 $N 1")
    $l.Add("1 $E $A $Rho")
    for ($i=1; $i -le $N; $i++) { $l.Add("$i $i $(($i+1)) 1") }
    Write-TextFile -Path $Path -Lines $l
}

# -- H8 cantilever convergence --
function Write-H8Conv { param([int]$N,[string]$Path,[int]$EType)
    $E=1000; $Nu=0.3; $Rho=1; $L=10; $Ny=2; $Nz=2; $dx=$L/$N; $ec=$N*$Ny*$Nz
    $l=[System.Collections.Generic.List[string]]::new()
    $l.Add("H8_ET$EType`_N$N")
    $l.Add("$(($N+1)*($Ny+1)*($Nz+1)) 1 1 1")
    for ($iz=0; $iz -le $Nz; $iz++) { for ($iy=0; $iy -le $Ny; $iy++) { for ($ix=0; $ix -le $N; $ix++) {
        $nd = Get-H8NodeId -Ix $ix -Iy $iy -Iz $iz -Nx $N -Ny $Ny -Nz $Nz
        if ($ix -eq 0) { $l.Add(('{0} 1 1 1 {1} {2} {3}' -f $nd,($ix*$dx),($iy/$Ny),($iz/$Nz))) }
        else           { $l.Add(('{0} 0 0 0 {1} {2} {3}' -f $nd,($ix*$dx),($iy/$Ny),($iz/$Nz))) }
    }}}
    $l.Add('1')
    $l.Add('0 1 0 0 -1')
    $l.Add("$EType $ec 1")
    $l.Add("1 $E $Nu $Rho")
    for ($ix=0; $ix -lt $N; $ix++) { for ($iz=0; $iz -lt $Nz; $iz++) { for ($iy=0; $iy -lt $Ny; $iy++) {
        $n1=Get-H8NodeId -Ix $ix -Iy $iy -Iz $iz -Nx $N -Ny $Ny -Nz $Nz
        $n2=Get-H8NodeId -Ix ($ix+1) -Iy $iy -Iz $iz -Nx $N -Ny $Ny -Nz $Nz
        $n3=Get-H8NodeId -Ix ($ix+1) -Iy ($iy+1) -Iz $iz -Nx $N -Ny $Ny -Nz $Nz
        $n4=Get-H8NodeId -Ix $ix -Iy ($iy+1) -Iz $iz -Nx $N -Ny $Ny -Nz $Nz
        $n5=Get-H8NodeId -Ix $ix -Iy $iy -Iz ($iz+1) -Nx $N -Ny $Ny -Nz $Nz
        $n6=Get-H8NodeId -Ix ($ix+1) -Iy $iy -Iz ($iz+1) -Nx $N -Ny $Ny -Nz $Nz
        $n7=Get-H8NodeId -Ix ($ix+1) -Iy ($iy+1) -Iz ($iz+1) -Nx $N -Ny $Ny -Nz $Nz
        $n8=Get-H8NodeId -Ix $ix -Iy ($iy+1) -Iz ($iz+1) -Nx $N -Ny $Ny -Nz $Nz
        $l.Add("$((($ix*$Nz+$iz)*$Ny+$iy)+1) $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 1")
    }}}
    Write-TextFile -Path $Path -Lines $l
}

function Run-Conv {
    param([string]$Label,[int[]]$Ns,[scriptblock]$Writer,[scriptblock]$TipNode,
          [int]$Comp,[double]$Exact)
    Write-Host "`n$Label"
    Write-Host ('{0,3} {1,14} {2,14} {3,12}' -f 'N','value','rel_err','rate')
    $prevE=$null
    foreach ($n in $Ns) {
        $path = (Join-Path $WorkDir "${Label}_N$n.dat") -replace ' ','_'
        & $Writer $n $path
        $out = Invoke-Stappp $path
        if (-not $out) { Write-Host ('{0,3}  solver failed' -f $n); continue }
        $v = Get-Disp $out (& $TipNode $n) $Comp
        if ($null -eq $v) { Write-Host ('{0,3}  parse error' -f $n); continue }
        $err = [Math]::Abs($v - $Exact) / [Math]::Abs($Exact)
        $rate=''
        if ($prevE -ne $null -and $err -gt 1e-16) { $rate = [Math]::Log($prevE/$err)/[Math]::Log(2) }
        Write-Host ('{0,3} {1,14:G6} {2,14:G6} {3,12:F3}' -f $n,$v,$err,$rate)
        $prevE=$err
    }
    Write-Host ('  Exact = {0:G6}' -f $Exact)
}

Run-Conv 'Truss self-weight' @(1,2,4,8,16,32) `
    ${function:Write-TrussConv} {param($n) $n+1} 1 0.0005

Run-Conv 'H8 cantilever self-weight' @(2,4,8,16,32) `
    {param($n,$p) Write-H8Conv $n $p 4} {param($n) 5*$n+5} 3 -15.0

Run-Conv 'BbarH8 cantilever self-weight' @(2,4,8,16,32) `
    {param($n,$p) Write-H8Conv $n $p 8} {param($n) 5*$n+5} 3 -15.0

# =================== Part 2: PATCH TESTS ===================
Write-Host "`n============================================================"
Write-Host "  PART 2 : PATCH TESTS"
Write-Host "============================================================"

# -- Truss patch: single element axial load --
Write-Host "`n-- Truss Patch (E=1e6,A=1,L=1,F=1.0 at tip) --"
$path = Join-Path $WorkDir 'patch_truss.dat'
@(
'Truss_Patch'
'2 1 1 1'
'1 1 1 1 0.0 0.0 0.0'
'2 0 1 1 1.0 0.0 0.0'
'1'
'1 0 0 0 0'
'2 1 1.0'
'1 1 1'
'1 1.0E6 1.0 0.0'
'1 1 2 1'
) | Set-Content -Encoding ASCII -Path $path
$out = Invoke-Stappp $path
if ($out) {
    $s = Get-Stress $out 1 1
    $ux = Get-Disp $out 2 1
    Write-Host "  Expected: stress=F/A=1.0,  ux=FL/EA=1.0e-6"
    Write-Host ('  FEM:      stress={0:G6}, ux={1:G6}' -f $s, $ux)
}

# -- H8 patch: single element axial --
$label='H8'; $etype=4
Write-Host "`n-- $label Patch (single cube, uniaxial, E=1000,nu=0.3,Fx=1.0) --"
$path = Join-Path $WorkDir "patch_$label.dat"
$l=[System.Collections.Generic.List[string]]::new()
$l.Add("${label}_Patch")
$l.Add('8 1 1 1')
$l.Add('1 1 1 1 0.0 0.0 0.0')
$l.Add('2 0 0 0 1.0 0.0 0.0')
$l.Add('3 0 0 0 1.0 1.0 0.0')
$l.Add('4 1 0 0 0.0 1.0 0.0')
$l.Add('5 1 0 0 0.0 0.0 1.0')
$l.Add('6 0 0 0 1.0 0.0 1.0')
$l.Add('7 0 0 0 1.0 1.0 1.0')
$l.Add('8 1 0 0 0.0 1.0 1.0')
$l.Add('1')
$l.Add('4 0 0 0 0')
$l.Add('2 1 0.25')
$l.Add('3 1 0.25')
$l.Add('6 1 0.25')
$l.Add('7 1 0.25')
$l.Add("$etype 1 1")
$l.Add('1 1000 0.3 0')
$l.Add('1 1 2 3 4 5 6 7 8 1')
Write-TextFile -Path $path -Lines $l
$out = Invoke-Stappp $path
if ($out) {
    $s = Get-Stress $out 1 0
    $ux = Get-Disp $out 2 1
    Write-Host "  Expected: sigma_xx = F/A = 1.0"
    Write-Host ('  FEM:      sigma_xx={0:G6}, tip_ux={1:G6}' -f $s, $ux)
}

# -- BbarH8 patch --
$label='BbarH8'; $etype=8
Write-Host "`n-- $label Patch (single cube, uniaxial, E=1000,nu=0.3,Fx=1.0) --"
$path = Join-Path $WorkDir "patch_$label.dat"
$l=[System.Collections.Generic.List[string]]::new()
$l.Add("${label}_Patch")
$l.Add('8 1 1 1')
$l.Add('1 1 1 1 0.0 0.0 0.0')
$l.Add('2 0 0 0 1.0 0.0 0.0')
$l.Add('3 0 0 0 1.0 1.0 0.0')
$l.Add('4 1 0 0 0.0 1.0 0.0')
$l.Add('5 1 0 0 0.0 0.0 1.0')
$l.Add('6 0 0 0 1.0 0.0 1.0')
$l.Add('7 0 0 0 1.0 1.0 1.0')
$l.Add('8 1 0 0 0.0 1.0 1.0')
$l.Add('1')
$l.Add('4 0 0 0 0')
$l.Add('2 1 0.25')
$l.Add('3 1 0.25')
$l.Add('6 1 0.25')
$l.Add('7 1 0.25')
$l.Add("$etype 1 1")
$l.Add('1 1000 0.3 0')
$l.Add('1 1 2 3 4 5 6 7 8 1')
Write-TextFile -Path $path -Lines $l
$out = Invoke-Stappp $path
if ($out) {
    $s = Get-Stress $out 1 0
    $ux = Get-Disp $out 2 1
    Write-Host "  Expected: sigma_xx = F/A = 1.0"
    Write-Host ('  FEM:      sigma_xx={0:G6}, tip_ux={1:G6}' -f $s, $ux)
}

# =================== Part 3: VERIFICATION ===================
Write-Host "`n============================================================"
Write-Host "  PART 3 : VERIFICATION CASES"
Write-Host "============================================================"

Write-Host "`n-- verify_truss (E=2e11,A=0.01,L=1,F=1000) --"
$out = Invoke-Stappp (Join-Path $RepoRoot 'data\verify_truss.dat')
if ($out) {
    $s=Get-Stress $out 1 1; $ux=Get-Disp $out 2 1
    Write-Host "  Analytical: stress=F/A=1e5, ux=FL/EA=5e-7"
    Write-Host ('  FEM:        stress={0:G6}, ux={1:G6}' -f $s,$ux)
}

Write-Host "`n-- verify_h8 (single hex, 4-point x-load, E=1000) --"
$out = Invoke-Stappp (Join-Path $RepoRoot 'data\verify_h8.dat')
if ($out) {
    $ux=Get-Disp $out 2 1
    Write-Host "  Expected ux(node2) ~ 0.004 (by FL/EA = 4*1/(1000*1) = 0.004)"
    Write-Host ('  FEM:     ux(node2) = {0:G6}' -f $ux)
}

Write-Host "`n-- verify_h8_consistent --"
$out = Invoke-Stappp (Join-Path $RepoRoot 'data\verify_h8_consistent.dat')
if ($out) {
    $ux=Get-Disp $out 2 1
    Write-Host ('  FEM: ux(node2) = {0:G6}' -f $ux)
}

Write-Host "`n============================================================"
Write-Host "  ALL VALIDATION COMPLETE"
Write-Host "============================================================"
