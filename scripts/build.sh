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

cd "$REPO_ROOT"
git add logs

if git diff --cached --quiet; then
  echo "새로 바뀐 로그가 없습니다 -- 커밋 생략"
  exit 0
fi

git commit -m "[Windows] chore: build log $(date +%Y-%m-%d\ %H:%M)"
git push
