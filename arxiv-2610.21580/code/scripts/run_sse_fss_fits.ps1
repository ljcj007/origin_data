param(
    [string]$RunDir = "results/sse_fss_v1",
    [string]$Python = "D:/program/Anaconda3_2025.06/python.exe",
    [string]$FitScript = "D:/Project/harnessing-quantum-fork/submission/quantum.harness/scripts/scaling_fit.py",
    [string]$FitScriptSha256 = "9D6C8ADED7A388AE15AFE84C5B5BF2C3324573162E6C62E731E4B2FE4039FB32"
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
if (-not [System.IO.Path]::IsPathRooted($RunDir)) {
    $RunDir = Join-Path $projectRoot $RunDir
}
$table = Join-Path $RunDir "fss_observables.csv"
$fitDir = Join-Path $RunDir "fits"
if (-not (Test-Path -LiteralPath $FitScript)) {
    throw "缺少冻结的拟合程序：$FitScript"
}
$actualFitHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $FitScript).Hash
if ($actualFitHash -ne $FitScriptSha256) {
    throw "拟合程序哈希不匹配：$actualFitHash"
}

# 独立调用拟合脚本时也必须先通过 14 单元完整性门。
& $Python (Join-Path $PSScriptRoot "analyze_sse_fss.py") --run-dir $RunDir
if ($LASTEXITCODE -ne 0) { throw "正式 SSE 汇总失败，禁止拟合" }

New-Item -ItemType Directory -Force $fitDir | Out-Null

# 采用由二次 Goldstone 启发的预注册能量形式。
& $Python $FitScript $table `
    --form polynomial --degree 6 `
    --obs-col energy_per_site --err-col energy_per_site_error `
    --pin c1=0 --pin c2=0 --pin c3=0 --pin c5=0 `
    --bootstrap 10000 --seed 20260903 `
    --out-csv (Join-Path $fitDir "energy_fit.csv") `
    --out-png (Join-Path $fitDir "energy_fit.png")
if ($LASTEXITCODE -ne 0) { throw "能量外推失败" }

# 序参量采用预注册的偶次修正形式。
& $Python $FitScript $table `
    --form polynomial --degree 4 `
    --obs-col staggered_m2 --err-col staggered_m2_error `
    --pin c1=0 --pin c3=0 `
    --bootstrap 10000 --seed 20260903 `
    --out-csv (Join-Path $fitDir "staggered_m2_fit.csv") `
    --out-png (Join-Path $fitDir "staggered_m2_fit.png")
if ($LASTEXITCODE -ne 0) { throw "序参量外推失败" }

# 刚度按同一预注册形式检验。
& $Python $FitScript $table `
    --form polynomial --degree 4 `
    --obs-col spin_stiffness_cell_average --err-col spin_stiffness_cell_average_error `
    --pin c1=0 --pin c3=0 `
    --bootstrap 10000 --seed 20260903 `
    --out-csv (Join-Path $fitDir "stiffness_fit.csv") `
    --out-png (Join-Path $fitDir "stiffness_fit.png")
if ($LASTEXITCODE -ne 0) { throw "刚度外推失败" }

# 子晶格量为补充外推，并与精确净磁矩约束交叉检查。
foreach ($observable in @("sublattice_plus_m2", "sublattice_minus_m2", "sublattice_cross")) {
    & $Python $FitScript $table `
        --form polynomial --degree 4 `
        --obs-col $observable --err-col "${observable}_error" `
        --pin c1=0 --pin c3=0 `
        --bootstrap 10000 --seed 20260903 `
        --out-csv (Join-Path $fitDir "${observable}_secondary_fit.csv") `
        --out-png (Join-Path $fitDir "${observable}_secondary_fit.png")
    if ($LASTEXITCODE -ne 0) { throw "子晶格补充外推失败：$observable" }
}

# 固定公式，仅改变最小尺寸检查截断敏感性。
$sensitivityDir = Join-Path $fitDir "sensitivity"
New-Item -ItemType Directory -Force $sensitivityDir | Out-Null
foreach ($minimumSize in @(6, 8)) {
    $subset = Join-Path $sensitivityDir "fss_Lmin${minimumSize}.csv"
    Import-Csv $table | Where-Object { [int]$_.L -ge $minimumSize } | Export-Csv $subset -NoTypeInformation

    & $Python $FitScript $subset `
        --form polynomial --degree 6 `
        --obs-col energy_per_site --err-col energy_per_site_error `
        --pin c1=0 --pin c2=0 --pin c3=0 --pin c5=0 `
        --bootstrap 10000 --seed 20260903 `
        --out-csv (Join-Path $sensitivityDir "energy_Lmin${minimumSize}.csv") `
        --out-png (Join-Path $sensitivityDir "energy_Lmin${minimumSize}.png")
    if ($LASTEXITCODE -ne 0) { throw "能量截断检查失败：Lmin=$minimumSize" }

    foreach ($observable in @("staggered_m2", "spin_stiffness_cell_average")) {
        & $Python $FitScript $subset `
            --form polynomial --degree 4 `
            --obs-col $observable --err-col "${observable}_error" `
            --pin c1=0 --pin c3=0 `
            --bootstrap 10000 --seed 20260903 `
            --out-csv (Join-Path $sensitivityDir "${observable}_Lmin${minimumSize}.csv") `
            --out-png (Join-Path $sensitivityDir "${observable}_Lmin${minimumSize}.png")
        if ($LASTEXITCODE -ne 0) { throw "截断检查失败：$observable, Lmin=$minimumSize" }
    }
}

& $Python (Join-Path $PSScriptRoot "summarize_sse_fss_fits.py") --run-dir $RunDir
if ($LASTEXITCODE -ne 0) { throw "拟合汇总失败" }
