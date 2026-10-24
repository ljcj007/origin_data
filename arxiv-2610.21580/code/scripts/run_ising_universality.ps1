param(
    [string]$Python = "D:/program/Anaconda3_2025.06/python.exe",
    [string]$FitScript = "D:/Project/harnessing-quantum-fork/submission/quantum.harness/scripts/scaling_fit.py",
    [string]$FitScriptSha256 = "9D6C8ADED7A388AE15AFE84C5B5BF2C3324573162E6C62E731E4B2FE4039FB32"
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$runDir = Join-Path $projectRoot "results/ising_tc_wolff_scaling_v1"
$table = Join-Path $runDir "critical_observables.csv"
$subset = Join-Path $runDir "critical_observables_Lmin16.csv"

if ((Get-FileHash -Algorithm SHA256 -LiteralPath $FitScript).Hash -ne $FitScriptSha256) {
    throw "Fit script hash mismatch"
}

# 固定精确临界温度，仅作普适类一致性检查。
& $Python (Join-Path $PSScriptRoot "analyze_ising_universality.py")
if ($LASTEXITCODE -ne 0) { throw "Critical observable reconstruction failed" }

Import-Csv $table | Where-Object { [int]$_.L -ge 16 } | Export-Csv $subset -NoTypeInformation

& $Python $FitScript (Join-Path $projectRoot "results/ising_tc_wolff_audit_v1/audit_observables.csv") `
    --form data-collapse --size-col L --obs-col binder --err-col binder_error_conservative `
    --param-col temperature --pin h_c=1.0041912981341914 --pin gamma_over_nu=0 `
    --bootstrap 1000 --seed 20260904 `
    --out-csv (Join-Path $runDir "binder_collapse.csv") `
    --out-png (Join-Path $runDir "binder_collapse.png")
if ($LASTEXITCODE -ne 0) { throw "Binder collapse failed" }

foreach ($fit in @(
    @{Table=$table; Observable="abs_m"; Error="abs_m_error"; Stem="abs_m_power"; Seed=20260904},
    @{Table=$table; Observable="connected_chi"; Error="connected_chi_error"; Stem="connected_chi_power"; Seed=20260905},
    @{Table=$subset; Observable="abs_m"; Error="abs_m_error"; Stem="abs_m_power_Lmin16"; Seed=20260906},
    @{Table=$subset; Observable="connected_chi"; Error="connected_chi_error"; Stem="connected_chi_power_Lmin16"; Seed=20260907}
)) {
    & $Python $FitScript $fit.Table --form power-law --size-col L `
        --obs-col $fit.Observable --err-col $fit.Error --bootstrap 10000 --seed $fit.Seed `
        --out-csv (Join-Path $runDir ($fit.Stem + ".csv")) `
        --out-png (Join-Path $runDir ($fit.Stem + ".png"))
    if ($LASTEXITCODE -ne 0) { throw "Power-law fit failed: $($fit.Stem)" }
}

& $Python (Join-Path $PSScriptRoot "summarize_ising_universality.py") --fit-script $FitScript
if ($LASTEXITCODE -ne 0) { throw "Universality summary failed" }

& $Python (Join-Path $PSScriptRoot "checksum_ledger.py") $runDir
if ($LASTEXITCODE -ne 0) { throw "Checksum ledger failed" }

Write-Output "Ising universality consistency workflow passed."
