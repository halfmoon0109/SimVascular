#!/bin/bash
BUILD_LOG="logs/build_$(date +%Y%m%d_%H%M).log"
mkdir -p logs
docker build . > "$BUILD_LOG" 2>&1
git add "$BUILD_LOG"
git commit -m "[Linux] build log $(date +%Y-%m-%d)"
git push
