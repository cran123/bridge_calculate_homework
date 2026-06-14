# ============================================================
#  桥梁模型 Abaqus 批量运行脚本
#  用法: 在装有 Abaqus 的机器上运行此脚本
#    powershell -ExecutionPolicy Bypass -File tools\run_bridge_abaqus.ps1
# ============================================================

param(
    [int]$Cpus = 8,
    [switch]$Bridge1,
    [switch]$Bridge2,
    [switch]$All
)

$RepoRoot = Split-Path -Parent $PSScriptRoot
$DataDir = Join-Path $RepoRoot "data\data-3"
$ScriptsDir = Join-Path $RepoRoot "abaqus_scripts"
$CasesDir = Join-Path $RepoRoot "cases"

# 检测 Abaqus
$abaqusCmd = Get-Command abaqus -ErrorAction SilentlyContinue
if (-not $abaqusCmd) {
    Write-Host "错误: 未找到 abaqus 命令，请确认 Abaqus 已安装并加入 PATH。" -ForegroundColor Red
    Write-Host "典型路径: C:\SIMULIA\Commands 或 C:\Program Files\Dassault Systemes\..." 
    exit 1
}

Write-Host "Abaqus: $($abaqusCmd.Source)" -ForegroundColor Green

# ===================== Bridge-1 =====================
function Run-Bridge1 {
    Write-Host ""
    Write-Host "========== Bridge-1 ==========" -ForegroundColor Cyan
    
    $inpFile = Join-Path $DataDir "Bridge-1.inp"
    if (-not (Test-Path $inpFile)) {
        Write-Host "  INP 文件不存在: $inpFile" -ForegroundColor Red
        return
    }

    $workDir = Join-Path $CasesDir "bridge1_run\abaqus"
    New-Item -ItemType Directory -Force -Path $workDir | Out-Null
    Copy-Item $inpFile $workDir -Force

    Push-Location $workDir
    try {
        Write-Host "  运行 Abaqus (cpus=$Cpus)..."
        abaqus job=Bridge-1 input=Bridge-1.inp cpus=$Cpus double interactive
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  Bridge-1 求解完成" -ForegroundColor Green
            
            Write-Host "  提取节点位移..."
            $extractScript = Join-Path $ScriptsDir "extract_odb_results.py"
            Push-Location $RepoRoot
            abaqus python $extractScript --odb (Join-Path $workDir "Bridge-1.odb")
            Pop-Location
            
            # 将 nodes CSV 也复制到 data-3 方便对比脚本找到
            Copy-Item "Bridge-1_nodes.csv" $DataDir -Force
            Write-Host "  节点位移已提取" -ForegroundColor Green
        } else {
            Write-Host "  Bridge-1 求解失败 (exit=$LASTEXITCODE)" -ForegroundColor Red
        }
    }
    finally {
        Pop-Location
    }
}

# ===================== Bridge-2 =====================
function Run-Bridge2 {
    Write-Host ""
    Write-Host "========== Bridge-2 ==========" -ForegroundColor Cyan
    
    $inpFile = Join-Path $DataDir "Bridge-2.inp"
    if (-not (Test-Path $inpFile)) {
        Write-Host "  INP 文件不存在: $inpFile" -ForegroundColor Red
        return
    }

    $workDir = Join-Path $CasesDir "bridge2_run\abaqus"
    New-Item -ItemType Directory -Force -Path $workDir | Out-Null
    Copy-Item $inpFile $workDir -Force

    Push-Location $workDir
    try {
        Write-Host "  运行 Abaqus (cpus=$Cpus)..."
        abaqus job=Bridge-2 input=Bridge-2.inp cpus=$Cpus double interactive
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  Bridge-2 求解完成" -ForegroundColor Green
            
            Write-Host "  提取节点位移..."
            $extractScript = Join-Path $ScriptsDir "extract_odb_results.py"
            Push-Location $RepoRoot
            abaqus python $extractScript --odb (Join-Path $workDir "Bridge-2.odb")
            Pop-Location
            
            Copy-Item "Bridge-2_nodes.csv" $DataDir -Force
            Write-Host "  节点位移已提取" -ForegroundColor Green
        } else {
            Write-Host "  Bridge-2 求解失败 (exit=$LASTEXITCODE)" -ForegroundColor Red
        }
    }
    finally {
        Pop-Location
    }
}

# ===================== 主流程 =====================
Write-Host "仓库根目录: $RepoRoot"
Write-Host "数据目录:   $DataDir"
Write-Host ""

if ($All -or $Bridge1) { Run-Bridge1 }
if ($All -or $Bridge2) { Run-Bridge2 }

if (-not $All -and -not $Bridge1 -and -not $Bridge2) {
    Write-Host "请指定要运行的桥梁模型:" -ForegroundColor Yellow
    Write-Host "  -Bridge1    仅运行 Bridge-1"
    Write-Host "  -Bridge2    仅运行 Bridge-2"  
    Write-Host "  -All        运行全部"
    Write-Host ""
    Write-Host "示例: .\run_bridge_abaqus.ps1 -Bridge2 -Cpus 8"
}

Write-Host ""
Write-Host "完成后运行对比: python tools\compare_bridges.py" -ForegroundColor Green
