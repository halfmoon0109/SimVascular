# scripts/build.ps1 (Windows PowerShell에서 실행)
#
# 실제 빌드는 컨테이너 안에서 별도로(몇 번이든) 수동으로 수행한다.
# 로그는 항상 컨테이너의 /work/logs/ 에 쌓이며, 이는 bind mount이므로
# 호스트에서는 이 저장소의 상위 폴더(D:\sv\logs)로 그대로 보인다.
# 이 스크립트는 그 최신 상태를 저장소 안 logs/ 로 복사하고,
# 실제로 바뀐 내용이 있을 때만 commit+push 한다.

$ErrorActionPreference = "Stop"

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$WorkRoot = Resolve-Path (Join-Path $RepoRoot "..")
$SourceLogs = Join-Path $WorkRoot "logs"
$DestLogs = Join-Path $RepoRoot "logs"

if (-not (Test-Path $SourceLogs)) {
    Write-Host "로그 소스 디렉터리를 찾을 수 없습니다: $SourceLogs" -ForegroundColor Yellow
    exit 0
}

New-Item -ItemType Directory -Force -Path $DestLogs | Out-Null
Copy-Item -Path (Join-Path $SourceLogs "*.log") -Destination $DestLogs -Force -ErrorAction SilentlyContinue

Push-Location $RepoRoot
try {
    git add logs
    git diff --cached --quiet
    if ($LASTEXITCODE -eq 0) {
        Write-Host "새로 바뀐 로그가 없습니다 -- 커밋 생략"
        exit 0
    }

    git commit -m "[Windows] chore: build log $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
    git push
}
finally {
    Pop-Location
}
