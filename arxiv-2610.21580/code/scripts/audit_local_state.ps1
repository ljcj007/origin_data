param(
    [string]$Python = "D:/program/Anaconda3_2025.06/python.exe",
    [ValidateSet("pre-sse", "final")]
    [string]$Stage = "pre-sse",
    [switch]$BuildManuscript
)

$ErrorActionPreference = "Stop"
$projectRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
Set-Location $projectRoot

if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    throw "Python 不存在: $Python"
}

# 冻结输入不得静默变化。
$expectedHashes = [ordered]@{
    "graphs/super_honeycomb_poscar.dat" = "58BE9C1BC747B0B17F45B31216E096E82CEBFB97874B3AAE5ADCC2DFAE1D65FF"
    "ising_graph.f90" = "879E10F321CABBD588A311F1C6FA3A2CF58056D3AC028E574D357B53DD9A56C4"
    "sse_graph_stiffness_v2.f90" = "D918DF52DD5236D125D707C20FFD39B813526C01DC45409EAD89E8ED8F4A58A0"
}
foreach ($entry in $expectedHashes.GetEnumerator()) {
    $actual = (Get-FileHash -LiteralPath $entry.Key -Algorithm SHA256).Hash
    if ($actual -ne $entry.Value) {
        throw "冻结哈希不一致: $($entry.Key)"
    }
}

$testOutput = & $Python -m unittest discover -s tests -p "test_*.py" 2>&1
if ($LASTEXITCODE -ne 0) {
    $testOutput | Write-Output
    throw "Python 回归测试失败"
}
$testText = $testOutput -join "`n"
$testMatch = [regex]::Match($testText, "Ran\s+(\d+)\s+tests?")
if (-not $testMatch.Success -or $testText -notmatch "(?m)^OK$") {
    throw "无法确认测试数量或 OK 状态"
}
$testCount = [int]$testMatch.Groups[1].Value

$expectedLedgerDirs = @(
    "results/poscar2",
    "results/lswt_super_honeycomb_v1",
    "results/ising_tc_wolff_v1",
    "results/ising_tc_wolff_audit_v1",
    "results/ising_tc_wolff_scaling_v1",
    "results/sse_beta_pilot_v2",
    "results/one_magnon_high_field_v1",
    "results/sse_dimer_validation_v1",
    "results/ising_l2_exact_validation_v1",
    "results/flatband_cls_topology_v1",
    "results/lswt_low_temperature_v1",
    "results/peripheral_benchmarks_v1",
    "results/sse_fss_stationarity_v1"
)
if ($Stage -eq "final") {
    $expectedLedgerDirs += @(
        "results/sse_fss_v1",
        "results/ed_super_honeycomb_v1"
    )
}
foreach ($dir in $expectedLedgerDirs) {
    $ledgerPath = Join-Path $dir "checksums.sha256.json"
    if (-not (Test-Path -LiteralPath $ledgerPath -PathType Leaf)) {
        throw "缺少预期结果账本: $dir"
    }
}

$ledgers = Get-ChildItem results -Recurse -Filter checksums.sha256.json -File | Sort-Object FullName
foreach ($ledger in $ledgers) {
    & $Python scripts/checksum_ledger.py $ledger.Directory.FullName --verify | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "结果账本失败: $($ledger.Directory.FullName)"
    }
}

$jsonCount = 0
Get-ChildItem . -Recurse -Filter *.json -File | ForEach-Object {
    Get-Content -LiteralPath $_.FullName -Raw -Encoding UTF8 | ConvertFrom-Json | Out-Null
    $jsonCount++
}

& gfortran -fsyntax-only -std=f2008 -Wall ising_graph.f90
if ($LASTEXITCODE -ne 0) { throw "Ising Fortran 语法失败" }
& gfortran -fsyntax-only -std=f2008 -Wall sse_graph_stiffness_v2.f90
if ($LASTEXITCODE -ne 0) { throw "SSE Fortran 语法失败" }
& gfortran -fsyntax-only -std=legacy ssebasic.f90
if ($LASTEXITCODE -ne 0) { throw "旧 SSE 兼容语法失败" }

$powershellErrors = 0
Get-ChildItem scripts -Filter *.ps1 -File | ForEach-Object {
    $tokens = $null
    $errors = $null
    [System.Management.Automation.Language.Parser]::ParseFile(
        $_.FullName,
        [ref]$tokens,
        [ref]$errors
    ) | Out-Null
    $powershellErrors += $errors.Count
}
if ($powershellErrors -ne 0) {
    throw "PowerShell AST 错误: $powershellErrors"
}

$shouldBuild = [bool]$BuildManuscript -or $Stage -eq "final"
if ($shouldBuild) {
    & scripts/build_manuscript.ps1
    if ($LASTEXITCODE -ne 0) {
        throw "论文构建失败"
    }
}

$auditPath = Join-Path $projectRoot "results/manuscript_audit.json"
& $Python scripts/audit_manuscript.py --output $auditPath | Out-Null
$auditExit = $LASTEXITCODE
$audit = Get-Content -LiteralPath $auditPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($Stage -eq "pre-sse") {
    $codes = @($audit.findings | ForEach-Object { $_.code } | Select-Object -Unique)
    if ($auditExit -ne 1 -or $audit.summary.error -ne 4 -or $audit.summary.warning -ne 0) {
        throw "SSE 前论文占位状态异常"
    }
    if ($codes.Count -ne 1 -or $codes[0] -ne "PRODUCTION_RESULT_TODO") {
        throw "论文出现非预期审计问题"
    }
} else {
    if ($auditExit -ne 0 -or $audit.summary.error -ne 0 -or $audit.summary.warning -ne 0) {
        throw "最终论文审计未通过"
    }

    $requiredFinalFiles = @(
        "results/sse_fss_v1/direct_batch_status.json",
        "results/sse_fss_v1/fss_observables.csv",
        "results/sse_fss_v1/fss_summary.json",
        "results/sse_fss_v1/fit_results_summary.json",
        "results/sse_fss_stationarity_v1/summary.json",
        "results/sse_fss_stationarity_v1/report.md",
        "figures/figure3_sse.pdf",
        "figures/figure3_sse.svg",
        "figures/figure3_sse.png",
        "figures/figure3_sse.tiff",
        "figures/figure3_qa.md",
        "results/ed_super_honeycomb_v1/result.json",
        "results/ed_super_honeycomb_v1/validation_report.json",
        "paper/manuscript_v1.pdf",
        "paper/manuscript_v1_qa.md"
    )
    foreach ($path in $requiredFinalFiles) {
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
            throw "缺少最终产物: $path"
        }
    }

    $cellDirs = @(Get-ChildItem results/sse_fss_v1/cells -Directory -Filter "cell-*")
    $manifests = @(Get-ChildItem results/sse_fss_v1/cells -Recurse -Filter manifest.json -File)
    if ($cellDirs.Count -ne 14 -or $manifests.Count -ne 14) {
        throw "正式 SSE 不是 14/14 完整单元"
    }

    $fss = Get-Content results/sse_fss_v1/fss_summary.json -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($fss.cell_count -ne 14 -or $fss.size_count -ne 7) {
        throw "正式 SSE 汇总完整性失败"
    }
    if ((@($fss.sizes) -join ",") -ne "4,6,8,10,12,14,16") {
        throw "正式 SSE 尺寸集合失败"
    }
    foreach ($bins in @($fss.bins_per_size)) {
        if ([int]$bins -ne 128) {
            throw "正式 SSE 每尺寸 bin 数不是 128"
        }
    }

    $stationarity = Get-Content results/sse_fss_stationarity_v1/summary.json -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($stationarity.provisional -ne $false -or $stationarity.completeness.available_cell_count -ne 14) {
        throw "SSE 独立平稳性审计尚未覆盖 14/14 单元"
    }

    $edReport = Get-Content results/ed_super_honeycomb_v1/validation_report.json -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($edReport.status -ne "passed") {
        throw "目标晶格 ED 验收未通过"
    }

    $checklist = Get-Content docs/completion_checklist.json -Raw -Encoding UTF8 | ConvertFrom-Json
    $unsatisfied = @($checklist.requirements | Where-Object { $_.blocking -eq $true -and $_.satisfied -ne $true })
    if ($unsatisfied.Count -ne 0 -or $checklist.completion_ready -ne $true) {
        throw "完成清单仍有阻断项"
    }
}

[ordered]@{
    status = "pass"
    stage = $Stage
    python_tests = $testCount
    checksum_ledgers = $ledgers.Count
    required_checksum_ledgers = $expectedLedgerDirs.Count
    json_files = $jsonCount
    frozen_hashes = $expectedHashes.Count
    fortran_sources = 3
    powershell_ast_errors = $powershellErrors
    manuscript_errors = $audit.summary.error
    manuscript_warnings = $audit.summary.warning
    manuscript_built = $shouldBuild
} | ConvertTo-Json
