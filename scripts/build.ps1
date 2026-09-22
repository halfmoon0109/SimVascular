# scripts/build.ps1 (Windows PowerShell에서 실행)
#
# 실제 빌드는 컨테이너 안에서 별도로(몇 번이든) 수동으로 수행한다.
# 로그는 항상 컨테이너의 /work/logs/ 에 쌓이며, 이는 bind mount이므로
# 호스트에서는 이 저장소의 상위 폴더(D:\sv\logs)로 그대로 보인다.
# 이 스크립트는 그 최신 상태를 저장소 안 logs/ 로 복사하고,
# 실제로 바뀐 내용이 있을 때만 commit+push 한다.
# 메셔가 GUI의 작업 폴더(컨테이너의 /work = 호스트의 D:\sv)에 남기는 교차 진단
# vtp(wall_*_crossings*.vtp)도 같이 복사한다: 웹 세션이 VTK 없이 읽어 교차를
# 국소 재현하는 데 쓴다. 지난 실행의 파일이 그대로 남아 있을 수 있으니, 교차가
# 없는 실행 뒤에 남은 것은 로그의 날짜와 맞지 않는 옛 파일이다.

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
Copy-Item -Path (Join-Path $WorkRoot "wall_*_crossings*.vtp") -Destination $DestLogs -Force -ErrorAction SilentlyContinue

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
