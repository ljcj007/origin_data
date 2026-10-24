param(
    [string]$Python = "D:/program/Anaconda3_2025.06/python.exe",
    [string]$Fortran = "gfortran"
)

$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$runDir = Join-Path $root "results/ising_l2_exact_validation_v1"
$buildDir = Join-Path $root "build/ising_l2_exact_validation_v1"
$graph = Join-Path $root "graphs/super_honeycomb_poscar.dat"
$source = Join-Path $root "ising_graph.f90"
$executable = Join-Path $buildDir "ising_graph_l2_validation.exe"
$acoshArgument = 2.0 + [Math]::Sqrt(3.0)
$temperature = 2.0 / [Math]::Log($acoshArgument + [Math]::Sqrt($acoshArgument * $acoshArgument - 1.0))

if (Test-Path -LiteralPath $runDir) {
    throw "结果目录已存在，拒绝覆盖：$runDir"
}
New-Item -ItemType Directory -Force -Path $runDir, $buildDir | Out-Null

& $Fortran -O3 -std=f2008 -Wall -Wextra -J $buildDir -o $executable $source
if ($LASTEXITCODE -ne 0) { throw "Fortran 编译失败" }

& $Python (Join-Path $root "scripts/ising_l2_exact.py") `
    --graph $graph --lx 2 --ly 2 --temperature $temperature --coupling 1 `
    --output-dir $runDir
if ($LASTEXITCODE -ne 0) { throw "精确枚举失败" }

$hashGraph = (Get-FileHash -LiteralPath $graph -Algorithm SHA256).Hash
$hashSource = (Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash
$hashExecutable = (Get-FileHash -LiteralPath $executable -Algorithm SHA256).Hash
$compiler = (& $Fortran --version | Select-Object -First 1)
$spec = [ordered]@{
    run_id = "ising_l2_exact_validation_v1"
    run_dir = "results/ising_l2_exact_validation_v1"
    settings = [ordered]@{
        model = "classical_antiferromagnetic_ising"
        graph = "graphs/super_honeycomb_poscar.dat"
        boundary = "periodic_xy"
        coupling = 1.0
        field = 0.0
        nbins = 256
        equil_sweeps = 100000
        budget = 5120000
    }
    provenance = [ordered]@{
        protocol = "L=2 exact-enumeration cross-check v1"
        graph_sha256 = $hashGraph
        source_sha256 = $hashSource
        compiler = "$compiler; -O3 -std=f2008 -Wall -Wextra"
        executable_sha256 = $hashExecutable
    }
    validation = [ordered]@{
        comparison_temperature = "2/acosh(2+sqrt(3))"
        conservative_block_sizes_bins = @(1, 2, 4, 8, 16)
        z_limit = 4.0
        metropolis_uniform_susceptibility = "diagnostic split; test 1/n_bin extrapolation"
    }
    cells = @(
        [ordered]@{
            cell_id = "wolff"
            params = [ordered]@{ temperature = $temperature; L = 2 }
            settings = [ordered]@{ algorithm = "wolff" }
        },
        [ordered]@{
            cell_id = "metropolis-short"
            params = [ordered]@{ temperature = $temperature; L = 2 }
            settings = [ordered]@{ algorithm = "metropolis" }
        },
        [ordered]@{
            cell_id = "metropolis"
            params = [ordered]@{ temperature = $temperature; L = 2 }
            settings = [ordered]@{ algorithm = "metropolis"; budget = 51200000 }
        }
    )
}
$specJson = $spec | ConvertTo-Json -Depth 8
[IO.File]::WriteAllText(
    (Join-Path $runDir "run_spec.json"),
    $specJson + [Environment]::NewLine,
    [Text.UTF8Encoding]::new($false)
)

foreach ($algorithm in @("wolff", "metropolis-short", "metropolis")) {
    & $Python (Join-Path $root "scripts/run_ising_cell.py") `
        --run-spec (Join-Path $runDir "run_spec.json") --cell-id $algorithm `
        --executable $executable
    if ($LASTEXITCODE -ne 0) { throw "蒙卡单元失败：$algorithm" }
}

& $Python (Join-Path $root "scripts/validate_ising_l2_exact.py") $runDir `
    --source $source --executable $executable --z-limit 4
if ($LASTEXITCODE -ne 0) { throw "精确对照未通过" }

@(
    "python=$(& $Python -c 'import sys; print(sys.version.split()[0])')"
    "numpy=$(& $Python -c 'import numpy; print(numpy.__version__)')"
    "compiler=$compiler"
    "compile_flags=-O3 -std=f2008 -Wall -Wextra"
) | Set-Content -LiteralPath (Join-Path $runDir "environment.txt") -Encoding utf8

& $Python (Join-Path $root "scripts/checksum_ledger.py") $runDir
if ($LASTEXITCODE -ne 0) { throw "哈希账本生成失败" }
& $Python (Join-Path $root "scripts/checksum_ledger.py") $runDir --verify
if ($LASTEXITCODE -ne 0) { throw "哈希账本核验失败" }
