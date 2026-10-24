param(
    [string]$RunDir = "results/sse_fss_v1",
    [string]$Python = "D:/program/Anaconda3_2025.06/python.exe"
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
if (-not [System.IO.Path]::IsPathRooted($RunDir)) {
    $RunDir = Join-Path $projectRoot $RunDir
}

# 先执行完整性门，再运行冻结拟合和绘图。
& $Python (Join-Path $PSScriptRoot "analyze_sse_fss.py") --run-dir $RunDir
if ($LASTEXITCODE -ne 0) { throw "正式 SSE 汇总失败" }

# 独立审计只报告提示，不自动排除数据。
& $Python (Join-Path $PSScriptRoot "audit_sse_stationarity.py") --run-dir $RunDir
if ($LASTEXITCODE -ne 0) { throw "SSE 独立平稳性审计失败" }
$stationarityPath = Join-Path $projectRoot "results/sse_fss_stationarity_v1/summary.json"
$stationarity = Get-Content -LiteralPath $stationarityPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($stationarity.provisional -ne $false -or $stationarity.completeness.available_cell_count -ne 14) {
    throw "SSE 独立平稳性审计未覆盖 14/14 单元"
}

& (Join-Path $PSScriptRoot "run_sse_fss_fits.ps1") -RunDir $RunDir -Python $Python
if ($LASTEXITCODE -ne 0) { throw "正式 SSE 拟合失败" }

& $Python (Join-Path $PSScriptRoot "figure3_sse.py") --run-dir $RunDir
if ($LASTEXITCODE -ne 0) { throw "Figure 3 生成失败" }

Push-Location $projectRoot
try {
    & $Python -m unittest discover -s tests -v
    if ($LASTEXITCODE -ne 0) { throw "本地回归测试失败" }
}
finally {
    Pop-Location
}

& $Python (Join-Path $PSScriptRoot "checksum_ledger.py") $RunDir
if ($LASTEXITCODE -ne 0) { throw "正式 SSE 哈希账本生成失败" }

& $Python (Join-Path $PSScriptRoot "checksum_ledger.py") $RunDir --verify
if ($LASTEXITCODE -ne 0) { throw "正式 SSE 哈希账本复核失败" }

Write-Output "正式 SSE 分析、独立平稳性审计、拟合、Figure 3 和测试均已完成。"
