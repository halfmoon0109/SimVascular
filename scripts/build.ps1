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
# 벽 인터페이스(법선·두께·ModelFaceID가 붙은 wall_offset_diagnostics.vtp, 오프셋
# 직전에 메셔가 남김)도 복사한다: 웹 세션이 재빌드 없이 오프셋 면부터 혼합 벽까지
# 재현하는 입력이다. 실행마다 같은 내용으로 다시 써지므로 모델이나 두께가 바뀔
# 때만 커밋에 잡힌다(압축 바이너리, 수십 MB; GitHub 한도 100 MB).

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
Copy-Item -Path (Join-Path $WorkRoot "wall_offset_diagnostics.vtp") -Destination $DestLogs -Force -ErrorAction SilentlyContinue

Push-Location $RepoRoot
try {
    $Branch = (git rev-parse --abbrev-ref HEAD).Trim()

    # 다른 세션(웹의 Claude 등)이 먼저 푸시한 커밋 위로 올라선다. 복사해 둔 로그의
    # 변경은 autostash로 보존된다. 이걸 빼먹으면 push가 "fetch first"로 거부되고,
    # 다음 실행은 새 로그가 없다며 커밋도 푸시도 건너뛰어 로그 커밋이 로컬에 갇힌다
    # (2026-09-30 21:44 실행).
    git pull --rebase --autostash origin $Branch
    if ($LASTEXITCODE -ne 0) {
        Write-Host "git pull --rebase 실패: 충돌을 풀고(git status) 다시 실행하세요" -ForegroundColor Red
        exit 1
    }

    git add logs
    git diff --cached --quiet
    if ($LASTEXITCODE -ne 0) {
        git commit -m "[Windows] chore: build log $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
    }
    else {
        Write-Host "새로 바뀐 로그가 없습니다 -- 커밋 생략"
    }

    # 이번 커밋뿐 아니라 지난 실행에서 푸시가 거부돼 남은 커밋도 함께: 원격보다 앞서 있으면 푸시한다
    $Ahead = (git rev-list --count "origin/$Branch..HEAD").Trim()
    if ($Ahead -ne "0") {
        git push -u origin $Branch
        if ($LASTEXITCODE -ne 0) {
            Write-Host "git push 실패: 스크립트를 다시 실행하면 pull --rebase 뒤 재시도합니다" -ForegroundColor Red
            exit 1
        }
    }
    else {
        Write-Host "원격과 같습니다 -- 푸시 생략"
    }
}
finally {
    Pop-Location
}
