param(
    [string]$SshConfig = "D:/Project/QMC/.remote/ssh_config",
    [string]$HostName = "qmc-direct",
    [string]$RemoteRoot = "/home/wangjieren/code/1x5_Heisenberg_Pentamer/workflows/super_honeycomb_qmc",
    [string]$RunDir = "results/sse_fss_v1",
    [string]$Python = "D:/program/Anaconda3_2025.06/python.exe"
)

$ErrorActionPreference = "Stop"
$Earliest = [datetimeoffset]"2026-09-04T05:30:00+08:00"
$projectRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
if ([datetimeoffset]::Now -lt $Earliest) {
    throw "SSH is disabled before $Earliest"
}

if (-not [System.IO.Path]::IsPathRooted($RunDir)) {
    $RunDir = Join-Path $projectRoot $RunDir
}
$RunDir = [System.IO.Path]::GetFullPath($RunDir)
if (-not $RunDir.StartsWith($projectRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "RunDir must stay inside the project"
}
$stageValidator = Join-Path $projectRoot "scripts/validate_sse_completion_stage.py"
if (-not (Test-Path -LiteralPath $Python -PathType Leaf) -or -not (Test-Path -LiteralPath $stageValidator -PathType Leaf)) {
    throw "Local stage validator environment is unavailable"
}

# 仅检查完成状态，不会重新提交。
$manifestText = & ssh -F $SshConfig $HostName "find '$RemoteRoot/results/sse_fss_v1/cells' -mindepth 2 -maxdepth 2 -name manifest.json -type f | wc -l"
if ($LASTEXITCODE -ne 0) { throw "Remote manifest check failed" }
$statusText = & ssh -F $SshConfig $HostName "find '$RemoteRoot/results/sse_fss_v1/direct_batch_status.json' -maxdepth 0 -type f | wc -l"
if ($LASTEXITCODE -ne 0) { throw "Remote status check failed" }

$manifestCount = [int](($manifestText | Select-Object -Last 1).Trim())
$statusCount = [int](($statusText | Select-Object -Last 1).Trim())
if ($manifestCount -ne 14 -or $statusCount -ne 1) {
    Write-Output "Remote run is not complete: manifests=$manifestCount, status=$statusCount"
    exit 2
}

$cellsDir = Join-Path $RunDir "cells"
foreach ($cell in @("cell-0013", "cell-0014")) {
    if (Test-Path -LiteralPath (Join-Path $cellsDir $cell)) {
        throw "Local target already exists: $cell"
    }
}
foreach ($name in @("direct_batch_status.json", "direct_batch.log")) {
    if (Test-Path -LiteralPath (Join-Path $RunDir $name)) {
        throw "Local target already exists: $name"
    }
}

$stage = Join-Path $projectRoot ("tmp/sse_fss_fetch_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $stage | Out-Null
foreach ($cell in @("cell-0013", "cell-0014")) {
    & scp -F $SshConfig -r "${HostName}:$RemoteRoot/results/sse_fss_v1/cells/$cell" $stage
    if ($LASTEXITCODE -ne 0) { throw "SCP failed: $cell" }
}
foreach ($name in @("direct_batch_status.json", "direct_batch.log")) {
    & scp -F $SshConfig "${HostName}:$RemoteRoot/results/sse_fss_v1/$name" $stage
    if ($LASTEXITCODE -ne 0) { throw "SCP failed: $name" }
}

foreach ($cell in @("cell-0013", "cell-0014")) {
    foreach ($name in @("manifest.json", "sse_graph.in", "sse_graph_bins.csv", "sse_graph_summary.txt", "stdout.txt", "stderr.txt")) {
        if (-not (Test-Path -LiteralPath (Join-Path $stage "$cell/$name"))) {
            throw "Staged artifact is missing: $cell/$name"
        }
    }
}

& $Python $stageValidator --run-spec (Join-Path $RunDir "run_spec.json") --stage $stage --cell-id cell-0013 --cell-id cell-0014
if ($LASTEXITCODE -ne 0) { throw "Independent staged SSE validation failed" }

# 暂存区通过后才进入正式目录。
$spec = Get-Content -LiteralPath (Join-Path $RunDir "run_spec.json") -Raw -Encoding UTF8 | ConvertFrom-Json
$batch = Get-Content -LiteralPath (Join-Path $stage "direct_batch_status.json") -Raw -Encoding UTF8 | ConvertFrom-Json
$expectedIds = @($spec.cells | ForEach-Object { $_.cell_id } | Sort-Object)
$batchIds = @($batch.cells | ForEach-Object { $_.cell_id } | Sort-Object)
if ($batch.run_id -ne $spec.run_id -or $batchIds.Count -ne 14 -or ($batchIds -join ",") -ne ($expectedIds -join ",")) {
    throw "Staged batch status does not describe the registered 14 cells"
}
if (@($batch.cells | Where-Object { [int]$_.returncode -ne 0 }).Count -ne 0) {
    throw "Staged batch status contains failed cells"
}

foreach ($cell in @("cell-0013", "cell-0014")) {
    $cellSpec = $spec.cells | Where-Object { $_.cell_id -eq $cell }
    $cellDir = Join-Path $stage $cell
    $manifest = Get-Content -LiteralPath (Join-Path $cellDir "manifest.json") -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($manifest.cell_id -ne $cell -or $manifest.status.success -ne $true -or [int]$manifest.execution.returncode -ne 0) {
        throw "Staged manifest failed: $cell"
    }
    if ([int]$manifest.params.L -ne [int]$cellSpec.params.L -or [int]$manifest.params.replica -ne [int]$cellSpec.params.replica) {
        throw "Staged manifest parameters mismatch: $cell"
    }
    if (($manifest.settings | ConvertTo-Json -Depth 20 -Compress) -ne ($spec.settings | ConvertTo-Json -Depth 20 -Compress)) {
        throw "Staged settings mismatch: $cell"
    }
    if (($manifest.provenance | ConvertTo-Json -Depth 20 -Compress) -ne ($spec.provenance | ConvertTo-Json -Depth 20 -Compress)) {
        throw "Staged provenance mismatch: $cell"
    }
    if ($manifest.execution.graph_sha256 -ne $spec.provenance.graph_sha256 -or $manifest.execution.executable_sha256 -ne $spec.provenance.executable_sha256) {
        throw "Staged execution hash mismatch: $cell"
    }
    foreach ($item in @(
        @{ key = "input"; name = "sse_graph.in" },
        @{ key = "bins"; name = "sse_graph_bins.csv" },
        @{ key = "summary"; name = "sse_graph_summary.txt" }
    )) {
        $actualHash = (Get-FileHash -LiteralPath (Join-Path $cellDir $item.name) -Algorithm SHA256).Hash
        $expectedHash = [string]$manifest.artifact_sha256.($item.key)
        if ($actualHash -ne $expectedHash.ToUpperInvariant()) {
            throw "Staged artifact hash mismatch: $cell/$($item.name)"
        }
    }
}

foreach ($cell in @("cell-0013", "cell-0014")) {
    Move-Item -LiteralPath (Join-Path $stage $cell) -Destination $cellsDir
}
foreach ($name in @("direct_batch_status.json", "direct_batch.log")) {
    Move-Item -LiteralPath (Join-Path $stage $name) -Destination $RunDir
}

Write-Output "Downloaded the completed 14/14 SSE run without overwriting local data."
