# 프로젝트 작업 컨벤션

## 코드 스타일 (최초 1회 분석 후 아래 채워 넣기)
- 언어/버전: 대부분 C++ (.cxx/.h, ~515 헤더 + ~514 소스, CMake 빌드), 일부 모듈(`sv4gui` Qt 플러그인)은 C++17 명시(`CMAKE_CXX_STANDARD 17`), 그 외는 표준 미명시. Python 스크립트/바인딩(.py, ~198개, `Python/site-packages/sv*`)은 Python 3 기준 표준 라이브러리 + vtk/numpy 사용. `Code/ThirdParty/*`는 외부 코드이므로 스타일 분석 대상에서 제외.
- 들여쓰기: 탭 미사용, 항상 스페이스. 레거시 코드(`Code/Source/sv`, `sv2` — Tcl 시대의 `cv` 접두사 클래스들)는 2-space, 최신 코드(`sv3`, `sv4gui`, `PythonAPI`, `vtkSV`)는 4-space가 우세. 새 코드 작성 시 해당 파일이 속한 모듈의 기존 들여쓰기를 따를 것(파일 내 일관성 우선). 중괄호는 Allman 스타일(개행 후 여는 중괄호)이 다수.
- 네이밍 규칙:
  - 클래스: PascalCase. 레거시 모듈은 `cv` 접두사(`cvRepositoryData`, `cvMath`), 최신 모듈(`sv3` 네임스페이스)은 접두사 없이 `PathElement` 등으로 사용하되 `sv3::` 네임스페이스로 스코프.
  - 멤버 함수: PascalCase (`GetControlPoint`, `SetSpacing`).
  - 멤버 변수: 최신 코드는 `m_PascalCase` 접두사(`m_ControlPoints`, `m_Spacing`), 레거시 코드는 트레일링 언더스코어(`type_`, `lockCnt_`).
  - 지역 변수/매개변수: camelCase (`controlPoint`, `sampleRate`).
  - 매크로/enum 상수: SCREAMING_SNAKE_CASE, enum 타입명은 `...T` 접미사 관행(`RepositoryDataT`, `CalculationMethod`).
  - 헤더 가드: `__모듈명_파일명_H__` 형태의 이중 언더스코어 패턴(`__SV3_PATHELEMENT_H__`).
  - Python: PEP8 준수 — 함수/변수 snake_case (`read_geo`, `write_geo`).
  - PythonAPI 바인딩 파일: `모듈명_PyClass.cxx` / `모듈명_PyModule.cxx` 네이밍 컨벤션.
- 주석/문서화 스타일:
  - 모든 파일 상단에 Stanford/SimVascular 표준 라이선스 헤더 블록(`/* Copyright (c) Stanford University... */`) 필수.
  - 함수 앞에 대시 배너 주석으로 함수명을 표기하는 관행 다수 사용: `//---------------\n// FunctionName\n//---------------\n// 설명...`.
  - 인라인 주석은 `//`, 블록 설명은 `/* ... */`. 문서화 태그(Doxygen 등) 강제 사용 없음.
  - Python은 함수 docstring에 `"""..."""` 형태로 Args/Returns를 기술(Google 스타일에 가까움).
- 기타 프로젝트 고유 컨벤션:
  - 소스 트리가 시대별로 계층화되어 있음: `sv`(레거시 C/C++, cv 접두사) → `sv2`/`sv3`(과도기, PascalCase+m_ 접두사) → `sv4gui`(현재 Qt/MITK 기반 GUI, C++17) → `PythonAPI`(Python 바인딩 레이어). 새 코드는 작업 대상 레이어의 기존 관행을 우선 따를 것.
  - 익스포트 매크로(`SV_EXPORT_모듈명`, 예: `SV_EXPORT_REPOSITORY`, `SV_EXPORT_PATH`)로 DLL/공유 라이브러리 심볼 가시성 관리.
  - 각 모듈은 자체 `CMakeLists.txt`와 `*Exports.h` 헤더를 가짐.

> 위 항목은 코드베이스를 최초 1회 스캔해서 채워 넣는다.
> 이후 작업에서는 이 섹션만 참고하고, **매번 전체 코드 구조를 다시 스캔하지 않는다.**
> 예외: 이 요약이 최신 코드와 명백히 다르다고 판단될 때만, 관련 파일 1~2개만 국소적으로 확인 후 이 섹션을 업데이트한다.

> 위 항목은 코드베이스를 최초 1회 스캔해서 채워 넣는다.
> 이후 작업에서는 이 섹션만 참고하고, **매번 전체 코드 구조를 다시 스캔하지 않는다.**
> 예외: 이 요약이 최신 코드와 명백히 다르다고 판단될 때만, 관련 파일 1~2개만 국소적으로 확인 후 이 섹션을 업데이트한다.

## 브랜치 전략
- 새 작업은 별도 브랜치에서 진행 (main 직접 커밋 금지)
- 하나의 브랜치를 Claude Code → Codex → Claude Code 순서로 이어받아 작업 (매번 새 브랜치 생성하지 않음)
- 병합은 PR을 통해서만, 최종 merge는 사람이 확인 후 진행

## 커밋 규칙
- 형식: `[작업자] 타입: 요약` 예) `[Claude] feat: 초기 구현`, `[Codex] review: 예외 처리 보강`
- 하나의 커밋에는 하나의 논리적 변경만 포함

## 빌드/실행 환경 분리
- 코드 작성/수정: Mac (Claude Code, Codex CLI) — 실제 빌드/실행 불가능한 환경
- 실제 빌드 및 실행: Linux 데스크톱의 Docker 컨테이너에서 git pull 후 재빌드
- Mac 환경에서는 "테스트 통과"라고 단정하지 말 것. 가능한 범위 내 정적 검증(문법/명백한 참조 오류)만 수행

## 빌드 로그 피드백 루프
- 리눅스 데스크톱에서 재빌드/실행 후 로그를 저장소 내 `logs/` 폴더에 커밋 (예: `logs/build_YYYYMMDD_HHMM.log`)
- Mac에서 작업 시작 전, `logs/` 폴더의 최신 로그부터 확인
- 로그에 에러가 있으면 최우선으로 해당 문제부터 수정

## 에이전트 간 인수인계 규칙
- 작업 시작 전 `git log --oneline -5`로 직전 작업 내역 확인
- `git diff HEAD~1`로 직전 커밋 변경사항만 검토 (전체 재스캔 금지)
- 리뷰 단계(Codex)에서는 새 기능 추가 금지, 기존 변경사항의 버그/스타일/예외처리만 점검
- 사용자가 리뷰 결과의 수정을 명시적으로 요청하면 Codex도 동일 브랜치에서 해당 버그/스타일/예외처리를 최소 범위로 수정 가능
- Codex가 코드를 수정한 경우 `git diff --check`와 가능한 범위의 정적 검증을 수행하되, Mac 환경에서 실제 빌드/테스트 통과를 단정하지 않음
- 마무리 단계(Claude)에서는 전체 diff 재확인 + PR 설명 초안 작성

## 금지 사항
- main 브랜치에 직접 push/commit 금지
- 인증정보(API 키, 토큰 등) 하드코딩 금지
- 리뷰 없이 임의 대규모 리팩토링 금지
- 코드 스타일 섹션이 이미 채워져 있는데 임의로 다시 전체 스캔하는 행위 금지
