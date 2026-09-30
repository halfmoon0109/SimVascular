#!/bin/bash
# scripts/build.sh (git-bash/WSL 등에서 실행)
#
# 실제 빌드는 컨테이너 안에서 별도로(몇 번이든) 수동으로 수행한다.
# 로그는 항상 컨테이너의 /work/logs/ 에 쌓이며, 이는 bind mount이므로
# 호스트에서는 이 저장소의 상위 폴더(D:\sv\logs)로 그대로 보인다.
# 이 스크립트는 그 최신 상태를 저장소 안 logs/ 로 복사하고,
# 실제로 바뀐 내용이 있을 때만 commit+push 한다.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WORK_ROOT="$(cd "$REPO_ROOT/.." && pwd)"
SOURCE_LOGS="$WORK_ROOT/logs"
DEST_LOGS="$REPO_ROOT/logs"

if [ ! -d "$SOURCE_LOGS" ]; then
  echo "로그 소스 디렉터리를 찾을 수 없습니다: $SOURCE_LOGS" >&2
  exit 0
fi

mkdir -p "$DEST_LOGS"
cp -f "$SOURCE_LOGS"/*.log "$DEST_LOGS"/ 2>/dev/null || true
# build.ps1과 같이: 교차 진단 vtp와 벽 인터페이스(wall_offset_diagnostics.vtp)도 담는다
cp -f "$WORK_ROOT"/wall_*_crossings*.vtp "$DEST_LOGS"/ 2>/dev/null || true
cp -f "$WORK_ROOT"/wall_offset_diagnostics.vtp "$DEST_LOGS"/ 2>/dev/null || true

cd "$REPO_ROOT"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"

# build.ps1과 같이: 다른 세션이 먼저 푸시한 커밋 위로 올라선 뒤 커밋하고, 지난 실행에서
# 푸시가 거부돼 남은 커밋까지 원격보다 앞서 있으면 푸시한다
git pull --rebase --autostash origin "$BRANCH"

git add logs
if git diff --cached --quiet; then
  echo "새로 바뀐 로그가 없습니다 -- 커밋 생략"
else
  git commit -m "[Windows] chore: build log $(date +%Y-%m-%d\ %H:%M)"
fi

if [ "$(git rev-list --count "origin/$BRANCH..HEAD")" != "0" ]; then
  git push -u origin "$BRANCH"
else
  echo "원격과 같습니다 -- 푸시 생략"
fi
