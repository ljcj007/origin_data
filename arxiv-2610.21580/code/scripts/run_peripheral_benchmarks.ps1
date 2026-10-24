param(
    [string]$Output = "results/peripheral_benchmarks_v1"
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

$Python = "D:\program\Anaconda3_2025.06\python.exe"
if (-not (Test-Path -LiteralPath $Python)) {
    $Python = "python"
}

& $Python scripts/peripheral_benchmarks.py --output $Output
if ($LASTEXITCODE -ne 0) { throw "外围基准重建失败" }
& $Python scripts/checksum_ledger.py $Output
if ($LASTEXITCODE -ne 0) { throw "校验账本生成失败" }
& $Python scripts/checksum_ledger.py $Output --verify
if ($LASTEXITCODE -ne 0) { throw "校验账本核对失败" }
