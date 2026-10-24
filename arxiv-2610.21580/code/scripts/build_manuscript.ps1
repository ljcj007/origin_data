param(
    [string]$Latex = "pdflatex",
    [string]$Bibtex = "bibtex"
)

$ErrorActionPreference = "Stop"
$paperDir = Join-Path (Split-Path -Parent $PSScriptRoot) "paper"
$stem = "manuscript_v1"

Push-Location $paperDir
try {
    # 禁止自动下载，确保离线构建可复现。
    & $Latex -disable-installer -interaction=nonstopmode -halt-on-error "$stem.tex"
    if ($LASTEXITCODE -ne 0) { throw "Initial LaTeX pass failed" }

    & $Bibtex $stem
    if ($LASTEXITCODE -ne 0) { throw "BibTeX pass failed" }

    foreach ($pass in 1..2) {
        & $Latex -disable-installer -interaction=nonstopmode -halt-on-error "$stem.tex"
        if ($LASTEXITCODE -ne 0) { throw "LaTeX convergence pass $pass failed" }
    }

    $log = Get-Content "$stem.log" -Raw
    if ($log -match "Citation.*undefined|There were undefined references|Overfull \\hbox") {
        throw "The manuscript log contains unresolved references or overfull boxes"
    }
    if (-not (Test-Path -LiteralPath "$stem.pdf")) {
        throw "The manuscript PDF was not created"
    }
}
finally {
    Pop-Location
}

Write-Output "Offline manuscript build passed: paper/$stem.pdf"
