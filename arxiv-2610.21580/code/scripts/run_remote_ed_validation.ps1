param(
    [string]$SshConfig = "D:/Project/QMC/.remote/ssh_config",
    [string]$HostName = "qmc-direct",
    [string]$RemoteRoot = "/home/wangjieren/code/1x5_Heisenberg_Pentamer/workflows/super_honeycomb_qmc",
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$Earliest = [datetimeoffset]"2026-09-04T05:30:00+08:00"
if ([datetimeoffset]::Now -lt $Earliest) {
    throw "SSH is disabled before $Earliest"
}

$projectRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$source = Join-Path $projectRoot "scripts/ed_graph_heisenberg.py"
$validator = Join-Path $projectRoot "scripts/validate_ed_super_honeycomb.py"
$graph = Join-Path $projectRoot "graphs/super_honeycomb_poscar.dat"
$sseManifest = Join-Path $projectRoot "results/sse_beta_pilot_v2/cells/cell-0013/manifest.json"
$localDir = Join-Path $projectRoot "results/ed_super_honeycomb_v1"
$checksumTool = Join-Path $projectRoot "scripts/checksum_ledger.py"

foreach ($path in @($source, $validator, $graph, $sseManifest, $checksumTool)) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Required local file is missing: $path"
    }
}
if (-not (Test-Path -LiteralPath $localDir -PathType Container)) {
    throw "Local result directory is missing: $localDir"
}
if (-not (Get-Command $Python -ErrorAction SilentlyContinue)) {
    throw "Local Python command is unavailable: $Python"
}

foreach ($name in @(
    "result.json",
    "stdout.txt",
    "stderr.txt",
    "environment.txt",
    "validation_report.json",
    "checksums.sha256.json"
)) {
    if (Test-Path -LiteralPath (Join-Path $localDir $name)) {
        throw "Local artifact already exists: $name"
    }
}

$expectedGraphHash = "58BE9C1BC747B0B17F45B31216E096E82CEBFB97874B3AAE5ADCC2DFAE1D65FF"
$sourceHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $source).Hash
$graphHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $graph).Hash
if ($graphHash -ne $expectedGraphHash) {
    throw "Local graph hash mismatch: $graphHash"
}

# 先只读检查容量，繁忙时不启动 CPU 任务。
$capacityText = & ssh -F $SshConfig $HostName "python3 -c `"import json,os; n=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1); a=os.getloadavg(); print(json.dumps({'nproc':n,'load1':a[0],'load5':a[1],'load15':a[2]}))`""
if ($LASTEXITCODE -ne 0) { throw "Remote capacity probe failed" }
$capacity = ($capacityText | Select-Object -Last 1) | ConvertFrom-Json
if ([int]$capacity.nproc -le 0) { throw "Remote CPU count is invalid" }
$loadFraction = [double]$capacity.load5 / [double]$capacity.nproc
$processText = & ssh -F $SshConfig $HostName "ps -eo pid,comm,%cpu,%mem,args --sort=-%cpu | head -n 12"
if ($LASTEXITCODE -ne 0) { throw "Remote process probe failed" }
$gpuText = & ssh -F $SshConfig $HostName "if command -v nvidia-smi >/dev/null 2>&1; then nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader || echo GPU_PROBE_FAILED; else echo GPU_PROBE_UNAVAILABLE; fi"
if ($LASTEXITCODE -ne 0) { throw "Remote GPU probe failed" }
Write-Output ("Remote capacity: nproc={0}, load1={1}, load5={2}, load15={3}, load5/nproc={4:N3}" -f $capacity.nproc, $capacity.load1, $capacity.load5, $capacity.load15, $loadFraction)
Write-Output "Remote top processes:"
$processText | Write-Output
Write-Output "Remote GPU probe:"
$gpuText | Write-Output
if ($loadFraction -ge 0.80) {
    throw "Remote CPU is busy; the registered ED source has no equivalent GPU backend, so no computation was started"
}

$sourceTag = $sourceHash.Substring(0, 12).ToLowerInvariant()
$remoteBin = "$RemoteRoot/scripts/remote-bin"
$remoteSource = "$remoteBin/ed_graph_heisenberg_$sourceTag.py"
$remoteDir = "$RemoteRoot/results/ed_super_honeycomb_L2_v1"
$remoteResult = "$remoteDir/result.json"
$remoteStdout = "$remoteDir/stdout.txt"
$remoteStderr = "$remoteDir/stderr.txt"
$remoteEnvironment = "$remoteDir/environment.txt"

# 仅选择服务器已有环境。
$remotePythonCandidates = @(
    "/home/wangjieren/code/quantum-harness-ch66/.venv-q66-v1/bin/python",
    "/opt/miniconda3/bin/python",
    "python3"
)
$remotePython = $null
$environmentText = @()
foreach ($candidate in $remotePythonCandidates) {
    $probeText = & ssh -F $SshConfig $HostName "cd '$RemoteRoot' && if command -v '$candidate' >/dev/null 2>&1 || [ -x '$candidate' ]; then env PYTHONDONTWRITEBYTECODE=1 '$candidate' -c `"import sys,numpy,scipy; print('python='+sys.executable); print('version='+sys.version.split()[0]); print('numpy='+numpy.__version__); print('scipy='+scipy.__version__)`"; else exit 127; fi"
    if ($LASTEXITCODE -eq 0) {
        $remotePython = $candidate
        $environmentText = @($probeText)
        break
    }
}
if ($null -eq $remotePython) { throw "Remote SciPy environment is unavailable" }
Write-Output "Remote ED Python: $remotePython"
$environmentText | Write-Output

$remoteGraphLine = & ssh -F $SshConfig $HostName "sha256sum '$RemoteRoot/graphs/super_honeycomb_poscar.dat'"
if ($LASTEXITCODE -ne 0) { throw "Remote graph hash check failed" }
$remoteGraphHash = (($remoteGraphLine | Select-Object -Last 1) -split '\s+')[0].ToUpperInvariant()
if ($remoteGraphHash -ne $expectedGraphHash) {
    throw "Remote graph hash mismatch: $remoteGraphHash"
}

$remoteSourceLine = & ssh -F $SshConfig $HostName "if [ -f '$remoteSource' ]; then sha256sum '$remoteSource'; else echo MISSING; fi"
if ($LASTEXITCODE -ne 0) { throw "Remote source check failed" }
if (($remoteSourceLine | Select-Object -Last 1).Trim() -eq "MISSING") {
    & ssh -F $SshConfig $HostName "mkdir -p '$remoteBin'"
    if ($LASTEXITCODE -ne 0) { throw "Remote source directory creation failed" }
    & scp -F $SshConfig $source "${HostName}:$remoteSource"
    if ($LASTEXITCODE -ne 0) { throw "ED source upload failed" }
}

$remoteSourceLine = & ssh -F $SshConfig $HostName "sha256sum '$remoteSource'"
if ($LASTEXITCODE -ne 0) { throw "Remote source hash check failed" }
$remoteSourceHash = (($remoteSourceLine | Select-Object -Last 1) -split '\s+')[0].ToUpperInvariant()
if ($remoteSourceHash -ne $sourceHash) {
    throw "Remote source hash mismatch: $remoteSourceHash"
}

$remoteState = & ssh -F $SshConfig $HostName "if [ -e '$remoteResult' ]; then echo result.json; fi; if [ -e '$remoteStdout' ]; then echo stdout.txt; fi; if [ -e '$remoteStderr' ]; then echo stderr.txt; fi; if [ -e '$remoteEnvironment' ]; then echo environment.txt; fi"
if ($LASTEXITCODE -ne 0) { throw "Remote ED state check failed" }
$remotePresent = @($remoteState | Where-Object { $_.Trim() -ne "" })
if ($remotePresent.Count -eq 0) {
    & ssh -F $SshConfig $HostName "set -C; mkdir -p '$remoteDir' && cd '$RemoteRoot' && env PYTHONDONTWRITEBYTECODE=1 '$remotePython' -c `"import sys,numpy,scipy; print('python='+sys.executable); print('version='+sys.version.split()[0]); print('numpy='+numpy.__version__); print('scipy='+scipy.__version__)`" > '$remoteEnvironment' && env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 '$remotePython' '$remoteSource' --graph graphs/super_honeycomb_poscar.dat --lx 2 --ly 2 --tol 1e-11 --maxiter 1200 --output '$remoteResult' > '$remoteStdout' 2> '$remoteStderr'"
    if ($LASTEXITCODE -ne 0) { throw "Remote ED run failed; logs remain in $remoteDir" }
}
elseif ($remotePresent.Count -ne 4) {
    throw "Remote ED directory contains a partial artifact set: $($remotePresent -join ', ')"
}

$stage = Join-Path $projectRoot ("tmp/ed_fetch_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $stage | Out-Null
foreach ($name in @("result.json", "stdout.txt", "stderr.txt", "environment.txt")) {
    & scp -F $SshConfig "${HostName}:$remoteDir/$name" $stage
    if ($LASTEXITCODE -ne 0) { throw "ED artifact download failed: $name" }
}

$report = Join-Path $stage "validation_report.json"
& $Python $validator --result (Join-Path $stage "result.json") --graph $graph --sse-manifest $sseManifest --report $report
if ($LASTEXITCODE -ne 0) { throw "ED validation failed" }

foreach ($name in @("result.json", "stdout.txt", "stderr.txt", "environment.txt", "validation_report.json")) {
    Move-Item -LiteralPath (Join-Path $stage $name) -Destination $localDir
}
& $Python $checksumTool $localDir
if ($LASTEXITCODE -ne 0) { throw "ED checksum ledger failed" }
& $Python $checksumTool $localDir --verify
if ($LASTEXITCODE -ne 0) { throw "ED checksum ledger verification failed" }

Write-Output "Remote N=20 ED completed, downloaded, and validated."
