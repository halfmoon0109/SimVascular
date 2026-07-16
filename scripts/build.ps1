# scripts/build.ps1 (Windows PowerShell에서 실행)

$timestamp = Get-Date -Format "yyyyMMdd_HHmm"
$buildLog = "logs/build_$timestamp.log"

New-Item -ItemType Directory -Force -Path "logs" | Out-Null

docker build . *> $buildLog

git add $buildLog
git commit -m "[Windows] build log $(Get-Date -Format 'yyyy-MM-dd')"
git push