<!--
Copyright (c) Stanford University, The Regents of the University of
California, and others.

All Rights Reserved.

See Copyright-SimVascular.txt for additional details.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject
to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-->

# SimVascular 저장소 구조 및 품질 스캔

## 1. 문서 목적

이 문서는 SimVascular 저장소를 공부하거나 변경 작업을 시작할 때 필요한 전체 구조,
빌드 흐름, 코드 계층, 실제 코드 스타일, 테스트 상태 및 주요 기술 부채를 한곳에
정리한다.

이 문서는 특정 기능의 구현 설명서가 아니다. 저장소 전체를 탐색하기 위한 지도이자,
코드를 읽을 때 주의해야 할 위험 요소를 기록한 기준 문서다.

스캔 기준은 다음과 같다.

| 항목 | 값 |
| --- | --- |
| 스캔 날짜 | 2026-07-19 |
| 브랜치 | `macOS_write` |
| 기준 커밋 | `1e16997b` (`[Claude] feat: 병합부위 곡률 기반 벽 두께 클램프 추가`) |
| 작업 트리 | clean |
| Linux 빌드 로그 | `logs/.gitkeep`만 존재하여 확인 불가 |
| 빌드·실행 | Mac 환경 정책에 따라 수행하지 않음 |

`Python/site-packages/svZeroDSolver`는 Git submodule 포인터만 있고 실제 내용이
초기화되지 않은 상태였다. 따라서 이 문서의 파일·줄 수 통계와 정적 검증에는 해당
서브모듈 내부 소스가 포함되지 않는다.

## 2. 저장소 규모

### 2.1 전체 인벤토리

| 항목 | 결과 |
| --- | ---: |
| Git 추적 항목 | 2,442개 |
| 추적 파일 총용량 | 약 94.26 MiB |
| `.cxx` | 514개, 295,491줄 |
| `.h` | 515개, 64,569줄 |
| `.hxx` | 8개, 27,309줄 |
| `.c` | 12개, 7,088줄 |
| `.cpp` | 3개, 124줄 |
| `.py` | 198개, 21,952줄 |
| `.sh` | 287개, 31,494줄 |
| `.cmake` | 114개, 13,561줄 |
| `CMakeLists.txt` | 73개 |
| `Makefile` | 63개 |

`Code/Source`만 보면 C/C++ 파일은 926개, 약 317,867줄이다. 자동 생성 또는
vendoring된 `sv4gui_json.hxx`를 포함한 값이다.

### 2.2 최상위 디렉터리 분포

| 디렉터리 | 추적 항목 | 용량 | 주요 역할 |
| --- | ---: | ---: | --- |
| `Code` | 1,567 | 20.83 MiB | 애플리케이션, 라이브러리, GUI, 내장 서드파티 |
| `Externals` | 509 | 2.05 MiB | 외부 라이브러리 빌드 정의와 패치 |
| `BuildWithMake` | 225 | 66.11 MiB | 레거시 빌드·배포 및 설치 파일 |
| `Python` | 78 | 0.52 MiB | Python 패키지와 submodule 포인터 |
| `Distribution` | 25 | 4.40 MiB | OS별 패키징 자료 |
| `Docker` | 16 | 0.04 MiB | 현재 Linux 빌드·실행 자동화 |
| `Documentation` | 3 | 0.21 MiB | 개발 문서와 Doxygen 설정 |

저장소 용량의 대부분은 소스가 아니라 레거시 배포 바이너리다. 추적된 `.exe`
14개가 약 52.79 MiB, `.msi` 1개가 약 14.46 MiB를 차지한다. Git LFS 규칙은
없다.

## 3. 전체 아키텍처

### 3.1 계층 관계

```text
외부 라이브러리
Qt / VTK / ITK / MITK / GDCM / HDF5 / OpenCASCADE / MMG
                         │
                         ▼
핵심 라이브러리
sv + sv2 + sv3 + vtkSV + 내장 ThirdParty
                         │
              ┌──────────┴──────────┐
              ▼                     ▼
        C++ Python API          sv4gui Modules
              │                     │
              ▼                     ▼
       Python sv 패키지       sv4gui Plugins
                                    │
                                    ▼
                           simvascular 실행 파일
```

하위 계층이 데이터 구조와 알고리즘을 제공하고, `sv4gui` 모듈과 플러그인이 이를
MITK/Qt 데이터 모델 및 사용자 인터페이스에 연결한다. `PythonAPI`는 같은 핵심
라이브러리를 CPython C API로 노출한다.

### 3.2 소스 계층별 역할

#### `Code/Source/sv`

가장 오래된 핵심 C++ 계층이다. `cvRepositoryData`, `cvMeshObject`처럼 `cv`
접두사를 가진 클래스와 C 스타일 API가 많다.

주요 하위 영역은 다음과 같다.

- `Utils`: 수학, VTK 유틸리티, 공용 함수
- `Repository`: 애플리케이션 객체 저장소
- `Geometry`: 형상 계산과 VTK 기반 처리
- `Model`: SolidModel 추상화, PolyData, OpenCASCADE 구현
- `Mesh`: MeshObject, TetGen, VMTK, MMG 및 적응형 메싱

#### `Code/Source/sv2`

레거시 계층과 최신 계층 사이의 과도기 코드다. 현재는 ImageProcessing과
PostProcessing에 집중되어 있다.

#### `Code/Source/sv3`

Path, Segmentation, Common, ITKSegmentation을 제공한다. 최신 `m_MemberName`
형식과 namespace를 사용하는 코드가 늘었지만, ITKSegmentation에는 이전 코드에서
이어진 탭과 다양한 헤더 가드가 많이 남아 있다.

#### `Code/Source/vtkSV`

VTK 알고리즘 라이브러리 모음이다.

- `Common`: 행렬, 수학, 렌더링, 일반 VTK 유틸리티
- `IO`: raw polydata/unstructured grid 입출력
- `Modules/Misc`: centerline, loft, 거리 및 데이터 처리
- `Modules/Geometry`: smoothing, subdivision, remeshing
- `Modules/Boolean`: 표면 교차 및 Boolean 연산
- `Modules/NURBS`: NURBS 표현과 계산
- `Modules/Parameterization`: 평면·구면·polycube parameterization
- `Modules/Segmentation`: centerline 및 vascular decomposition

#### `Code/Source/sv4gui`

현재 GUI의 중심이며 저장소에서 가장 큰 계층이다.

| 하위 계층 | 추적 파일 수 | C/C++ 파일·줄 수 |
| --- | ---: | ---: |
| `sv4gui` 전체 | 886 | 540개, 약 140,658줄 |
| `Modules` | 275 | 도메인 데이터와 서비스 |
| `Plugins` | 611 | Qt/MITK 사용자 인터페이스 |

모듈은 Path, Segmentation, Model, Mesh, Simulation, ROMSimulation,
MultiPhysics, ProjectManagement 등을 구현한다. 플러그인은 이 모듈을 화면,
메뉴, editor 및 Data Manager에 연결한다.

#### `Code/Source/PythonAPI`

CPython C API를 직접 사용해 geometry, imaging, math, meshing, modeling,
pathplanning, segmentation, simulation, ROM simulation 및 VMTK 기능을 노출한다.

파일 이름은 대체로 다음 구조를 따른다.

```text
ModuleName_PyModule.cxx
ModuleName_PyModule.h
ModuleNameClassName_PyClass.cxx
```

클래스 구현 `.cxx`를 모듈 구현에서 직접 `#include`하는 구조가 있으므로, 일반적인
C++ translation unit 구조와 다르게 읽어야 한다.

#### `Python/site-packages`

- `sv`: C++ 확장 모듈을 로드하는 기본 패키지
- `sv_vis`: VTK 기반 시각화 유틸리티
- `sv_rom_simulation`: 0D/1D reduced-order model 입력 생성
- `sv_rom_extract_results`: ROM 결과 추출과 후처리
- `sv_ml`: 머신러닝 2D 혈관 단면 분할
- `sv_auto_lv_modeling`: 심장 영상 기반 자동 LV 모델링
- `svZeroDSolver`: 별도 Git submodule

## 4. 빌드 시스템

### 4.1 상위 CMake 프로젝트

저장소 루트의 `CMakeLists.txt`는 `SV_TOP` 프로젝트다. 실제 소스를 직접
`add_subdirectory()`하지 않고 `ExternalProject_Add(SimVascular)`를 사용하여
`Code` 디렉터리를 별도 CMake 프로젝트로 구성한다.

외부 라이브러리를 자동으로 빌드하는 `ExternalProject_Add(Externals)` 블록은 현재
주석 처리되어 있다. 따라서 실제 빌드 환경은 외부 라이브러리 경로를 별도로 준비해
전달한다.

### 4.2 `Code` 프로젝트

`Code/CMakeLists.txt`는 실제 라이브러리와 실행 파일을 조립한다.

```text
SimVascularOptions
    ↓
SimVascularDependentOptions
    ↓
SimVascularInternals / Externals / ThirdParty / Licensed
    ↓
Core libraries
    ↓
sv4gui modules
    ↓
sv4gui plugins
    ↓
Application
```

기본 핵심 라이브러리에는 Repository, Geometry, SolidModel, MeshObject, Path,
Segmentation 및 vtkSV 모듈이 포함된다. TetGen, MMG, ITK, OpenCASCADE 등은
관련 옵션에 따라 추가된다.

MeshSim 및 Parasolid 관련 디렉터리는 CMake에서 조건부로 참조하지만 현재
체크아웃에는 없다. 이 기능들은 기본적으로 OFF이며, 활성화하려면 별도 라이선스
소스 또는 설치가 필요하다.

### 4.3 Docker 빌드

현재 Linux 전체 GUI 빌드 흐름은 `Docker/scripts/build-all.sh`에 가장 명확하게
정리되어 있다.

```text
10-qt6.sh
    ↓
20-python.sh
    ↓
30-small-libs.sh
    ↓
40-vtk.sh
    ↓
50-itk.sh
    ↓
60-opencascade.sh
    ↓
70-mitk.sh
    ↓
80-simvascular.sh
```

각 단계는 `/work/logs/<stage>.log`를 생성하고 첫 실패에서 중단한다. 저장소 정책상
Mac에서는 코드를 수정하고 Linux Docker 환경에서 이 흐름으로 빌드·실행한다.

### 4.4 외부 라이브러리 버전

현재 Code와 Docker가 주로 사용하는 버전은 다음과 같다.

| 라이브러리 | 버전 |
| --- | --- |
| Qt | 6.6.2 |
| HDF5 | 1.14.3 |
| TinyXML2 | 8.0.0 |
| Python | 3.11.0 |
| FreeType | 2.13.0 |
| MMG | 5.3.9 |
| GDCM | 3.0.10 |
| VTK | 9.3.0 |
| ITK | 5.4.0 |
| OpenCASCADE | 7.6.0 |
| MITK | Docker에서는 `v2024.06` |

버전 정보는 `Code/CMake/SimVascularExternalsVersions.cmake`,
`Externals/CMake/SvExtVersions.cmake`, `Docker/scripts/versions.env`에 분산되어
있으며 완전히 일치하지 않는다.

- Code와 Docker: TinyXML2 8.0.0
- Externals CMake: TinyXML2 6.2.0
- Code와 Externals CMake: MITK 2022.10
- Docker: MITK `v2024.06`

외부 버전을 변경할 때는 이 세 위치와 관련 빌드 스크립트를 함께 확인해야 한다.

## 5. 실제 코드 스타일

### 5.1 들여쓰기

“항상 스페이스, 탭 미사용”은 현재 저장소 전체에는 적용되지 않는다.

| 계층 | 탭이 있는 C/C++ 파일 | 탭 포함 줄 |
| --- | ---: | ---: |
| Application | 1 | 69 |
| PythonAPI | 1 | 2 |
| sv | 36 | 1,054 |
| sv2 | 7 | 200 |
| sv3 | 27 | 3,939 |
| sv4gui | 11 | 71 |
| vtkSV | 43 | 644 |
| 합계 | 126 | 5,979 |

특히 `sv3/ITKSegmentation`에 탭 사용이 집중되어 있다. 새 코드에서는 저장소
전체를 일괄 정리하지 말고 수정 대상 파일의 기존 들여쓰기를 유지하는 것이 안전하다.

### 5.2 중괄호

아래 수치는 함수와 제어문을 포함한 단순 패턴 집계이므로 절대적인 스타일 판정은
아니지만 계층별 경향은 잘 보여준다.

| 계층 | 같은 줄 `{` | 다음 줄 `{` | 경향 |
| --- | ---: | ---: | --- |
| Application | 31 | 32 | 혼합 |
| PythonAPI | 1,531 | 1,356 | 혼합 |
| sv | 1,127 | 2,819 | Allman 우세 |
| sv2 | 533 | 48 | 같은 줄 우세 |
| sv3 | 191 | 1,290 | Allman 우세 |
| sv4gui | 2,652 | 13,233 | Allman 우세 |
| vtkSV | 171 | 13,316 | Allman 강세 |

### 5.3 네이밍

- 레거시 클래스: `cvRepositoryData`, `cvMeshObject` 같은 `cv` 접두사
- 최신 클래스: PascalCase
- 멤버 함수: `Get...`, `Set...` 형태의 PascalCase
- 레거시 멤버: `type_`, `lockCnt_` 같은 trailing underscore
- 최신 멤버: `m_ControlPoints` 같은 `m_PascalCase`
- 지역 변수와 매개변수: 주로 camelCase
- 익스포트 매크로: `SV_EXPORT_*`

`SV_EXPORT_*` 매크로는 68개 헤더에서 387회 사용된다. DLL과 공유 라이브러리의
심볼 가시성을 제어하는 중요한 프로젝트 관행이다.

### 5.4 헤더 가드

`Code/Source` 헤더 457개의 가드는 다음과 같이 분포한다.

| 형식 | 개수 |
| --- | ---: |
| `__NAME...` 형태 | 53 |
| trailing underscore 형태 | 26 |
| 기타 형식 | 375 |
| 명시적인 가드 없음 | 3 |

따라서 `__MODULE_FILE_H__`를 저장소 전체의 단일 규칙으로 간주하면 안 된다.
새 헤더는 같은 모듈의 인접 헤더 형식을 따르는 것이 가장 일관적이다.

### 5.5 C++ 표준

저장소 안에서 `CMAKE_CXX_STANDARD 17`을 명시한 곳은
`sv4gui/Plugins/org.sv.gui.qt.datamanager/CMakeLists.txt` 한 곳뿐이다. 다른
타깃의 유효 C++ 표준은 컴파일러 기본값 또는 Qt/MITK target의 전이 compile
feature에 의존할 수 있다.

### 5.6 Python 스타일

- snake_case 함수 정의: 507개
- camelCase가 포함된 함수 정의: 84개
- Google식 `Args:`가 있는 파일: 19개
- Google식 `Returns:` 줄: 57개
- 탭이 있는 Python 파일: 2개, 35줄

Python은 대체로 Python 3 및 snake_case를 사용하지만 완전히 통일되어 있지는 않다.
`Python/site-packages/sv_vis/renfun.py`에는 Python 2의 `raw_input()` 호출도 남아
있다.

### 5.7 라이선스 헤더

자사 C/C++ 파일 922개 중 920개에서 표준 Stanford/SimVascular 라이선스
헤더를 확인했다. 예외 두 파일은 외부에서 유래한 것으로 보이는 JSON 헤더와 MITK
관련 구현이다.

Python 파일 71개 중 8개는 표준 헤더가 없으며, 대부분 `__init__.py` 또는 작은
옵션 정의 파일이다.

## 6. 테스트와 CI

### 6.1 저장소 내부 테스트

`Code/Testing`에는 CTest helper와 결과 비교 스크립트 7개만 있다. 실제 테스트
데이터는 저장소 밖의 `SV_TEST_DIR`을 요구한다.

`BUILD_TESTING`은 기본 OFF이며 상위 `ExternalProject_Add(SimVascular)`도
`BUILD_TESTING=OFF`를 전달한다. 테스트를 ON으로 설정하고 `SV_TEST_DIR`을
제공하지 않으면 configure 단계에서 중단된다.

현재 등록된 실질 테스트는 외부 Tcl 스크립트를 실행하는 `StartUpTest`뿐이다.

### 6.2 CI

`.github/workflows/build.yml`은 다음 세 플랫폼에서 빌드와 패키징을 정의한다.

- macOS 11
- Ubuntu 20.04, GCC 8, Qt 5 패키지
- Windows 2019, Visual Studio 2017

CI에는 `ctest`, Python 테스트, lint 또는 formatter 검사가 없다. Actions 버전도
`actions/checkout@v2`, `actions/upload-artifact@v3` 등에 고정되어 있다.
이 문서 작성 시점에 CI를 실제 실행하거나 원격 상태를 조회하지는 않았다.

### 6.3 정적 검증 결과

| 검사 | 결과 |
| --- | --- |
| Python 198개 구문 컴파일 | 오류 없음 |
| JSON 24개 파싱 | 오류 없음 |
| XML/UI/QRC 84개 파싱 | 오류 없음 |
| 직전 커밋 `git diff --check` | 오류 없음 |
| 주요 소스의 CRLF 검사 | 검출 없음 |
| 단순 비밀정보 패턴 검사 | 하드코딩된 키·토큰·개인키 검출 없음 |

`.sh` 287개 중 Bash 구문 검사에서 5개가 실패했지만, 하나는 Python shebang을
가진 Python 프로그램이고 나머지 네 개는 실행 스크립트가 아닌 sed 치환 규칙이다.
확장자와 실제 파일 형식이 일치하지 않는 레거시 사례다.

실제 C++ 빌드와 애플리케이션 실행은 수행하지 않았으므로 “빌드 통과” 또는
“테스트 통과”로 해석하면 안 된다.

## 7. 주요 위험 요소

### 7.1 테스트 자동화 부족

가장 큰 위험은 동작 변경을 자동 검증할 저장소 내부 테스트가 사실상 없다는 점이다.
현재 변경은 Linux Docker에서 수동 재빌드·실행한 뒤 로그를 저장소로 전달하는 흐름에
의존한다.

### 7.2 외부 버전 정의 분산

Code, Externals, Docker가 서로 다른 버전 값을 가질 수 있다. 특히 TinyXML2와
MITK가 실제로 다르다. 외부 라이브러리 업데이트 시 한쪽만 변경하면 configure,
compile 또는 ABI 문제가 발생할 수 있다.

### 7.3 CMake 오타 후보

정적 스캔에서 다음 참조가 확인됐다.

| 위치 | 현재 값 | 예상 값 또는 문제 |
| --- | --- | --- |
| `Code/CMakeLists.txt:212` | `SV_USE_MESHIM_ADAPTOR` | `SV_USE_MESHSIM_ADAPTOR` 오타 가능성 |
| `PolyDataSolidModel/CMakeLists.txt:50` | `VSV_LIB_VTKSVBOOLEAN_NAME` | `SV_LIB_VTKSVBOOLEAN_NAME` 오타 가능성 |
| `Testing/CMakeLists.txt` | `SV_TEST_SAVEOUT_DIR` 정의 | macro는 `SimVascular_TEST_SAVEOUT_DIR` 사용 |
| `Testing/CMakeLists.txt` | `SV_DEV_OUTPUT` 검사 | 실제 옵션은 `SV_DEVELOPER_OUTPUT` |

첫 번째 문제는 TetGen adaptor가 ON인 기본 구성에서는 가려질 수 있다. 두 번째는
빈 CMake 변수로 평가되어 필요한 직접 링크가 누락될 가능성이 있다. 테스트 관련
두 문제는 현재 테스트가 기본 OFF라 노출되지 않는다.

### 7.4 경고 억제

`SV_SUPPRESS_WARNINGS`의 기본값은 ON이며, 활성화되면 `-w`가 추가된다. GNU C++에는
`-fpermissive`도 적용된다. 오래된 코드의 빌드를 유지하는 데 도움이 되지만 새로운
타입·변환·API 문제를 조기에 발견하기 어렵게 만든다.

### 7.5 레거시 C 문자열 API

`Code/Source`에서 다음 호출을 확인했다.

- `strcpy`: 104회
- `sprintf`: 39회

호출이 존재한다는 사실만으로 취약점이라고 단정할 수는 없다. 다만 대상 버퍼 크기와
입력 길이를 검토해야 하는 우선 점검 위치다. `vtkSV`와 `sv` 계층에 집중되어 있다.

### 7.6 대형 파일과 복잡도

1,000줄 이상인 C/C++ 파일은 60개다. 주요 대형 파일은 다음과 같다.

| 파일 | 줄 수 |
| --- | ---: |
| `sv4gui_json.hxx` | 24,766, vendored/generated 헤더 |
| `vtkSVNURBSUtils.cxx` | 7,090 |
| `vtkSVSurfaceCuboidPatcher.cxx` | 5,963 |
| `sv4gui_ROMSimulationView.cxx` | 5,234 |
| `sv_sys_geom.cxx` | 4,360 |

`Code/Source`에는 TODO/FIXME/XXX/HACK 계열 표시도 371개 있다. 대규모 리팩터링보다
변경 대상 함수 주변을 제한적으로 이해하고 검증하는 접근이 안전하다.

### 7.7 저장소 위생

- `.orig`, `.modified`, `.old.cxx` 백업 파일 추적
- `.ui.1`, `.ui.2`, `.ui.3` 버전 파일 추적
- `CMakeFiles/cmake.check_cache` 같은 빌드 흔적 추적
- 실행 권한이 설정된 파일 502개
- 실행 비트가 있는 헤더 19개와 C++ 파일 3개
- `.gitattributes` 1,403줄 중 1,395줄이 과거 파일별 `-text` 규칙
- 대용량 설치 바이너리를 일반 Git object로 저장

이 항목들은 즉시 기능 오류를 의미하지 않지만 clone 크기, diff 가독성 및 플랫폼별
파일 처리의 유지보수 비용을 높인다.

### 7.8 초기화되지 않은 submodule

`.gitmodules`는 `Python/site-packages/svZeroDSolver`를
`https://github.com/SimVascular/svZeroDSolver.git`에 연결한다. 현재 체크아웃은
gitlink 커밋 `8a58d680...`만 기록하며 실제 소스는 내려받지 않았다.

0D solver를 공부하거나 빌드하려면 별도로 다음 작업이 필요하다.

```bash
git submodule update --init --recursive
```

## 8. 코드 변경 시 권장 접근

1. `git log --oneline -5`로 최근 작업 의도를 확인한다.
2. `git diff HEAD~1`로 직전 변경만 먼저 검토한다.
3. `logs/`의 최신 Linux 빌드 로그를 확인한다.
4. 대상 계층과 인접 파일 1~2개의 스타일을 기준으로 삼는다.
5. GUI 옵션을 추가할 때는 UI, command serialization, parser, core option,
   Python API가 모두 연결되는지 확인한다.
6. CMake 옵션을 추가할 때는 root superbuild와 `Code` 프로젝트 양쪽 전달 여부를
   확인한다.
7. Mac에서는 정적 검증까지만 수행하고 빌드 성공을 단정하지 않는다.
8. Linux Docker에서 증분 빌드·실행 후 로그를 `logs/`로 전달한다.

GUI에서 core 기능으로 이어지는 일반적인 옵션 전달 경로는 다음과 같다.

```text
Qt .ui 위젯
    ↓
Plugin command 생성
    ↓
sv4gui module parser
    ↓
legacy core SetMeshOptions 또는 대응 API
    ↓
실제 알고리즘
    ↑
PythonAPI option/type binding
```

한 계층만 변경하면 저장·복원, GUI, Python 또는 core 실행 경로 중 일부가 빠질 수
있으므로 이 흐름을 끝까지 추적해야 한다.

## 9. 학습 권장 순서

### 1단계: 프로젝트 목적과 빌드 진입점

1. `README.md`
2. 루트 `CMakeLists.txt`
3. `Code/CMakeLists.txt`
4. `Code/CMake/SimVascularOptions.cmake`
5. `Docker/WORKFLOW.md`

### 2단계: 핵심 데이터와 알고리즘

1. `sv/Repository`
2. `sv/Geometry`
3. `sv/Model`
4. `sv/Mesh`
5. `sv3/Path`
6. `sv3/Segmentation`
7. `vtkSV/Common`과 필요한 vtkSV module

### 3단계: GUI 아키텍처

1. `sv4gui/Modules/Common`
2. `sv4gui/Modules/ProjectManagement`
3. 관심 도메인의 `sv4gui/Modules/*`
4. 대응하는 `sv4gui/Plugins/org.sv.gui.qt.*`
5. `Source/Application`

### 4단계: Python 인터페이스

1. `Code/Source/PythonAPI/README.md`
2. `*_PyModule.cxx`
3. 대응하는 `*_PyClass.cxx`
4. `Python/site-packages/sv/__init__.py`
5. ROM, ML 또는 자동 모델링 패키지

### 5단계: 기능 하나를 수직으로 추적

한 기능을 선택하여 UI부터 core까지 수직으로 읽는 것이 가장 효과적이다. TetGen을
예로 들면 다음 순서가 적절하다.

1. `Documentation/TetGen_Mesh_Generation_Workflow.md`
2. `sv4gui_MeshEdit.ui`
3. `sv4gui_MeshEdit.cxx`
4. `sv4gui_MeshTetGen.cxx`
5. `sv_TetGenMeshObject.cxx`
6. `sv_tetgenmesh_utils.cxx`
7. `MeshingTetGenOptions_PyClass.cxx`

이 방식은 SimVascular의 레거시 core, 현대 GUI, Python binding이 어떻게 연결되는지
동시에 이해할 수 있게 해준다.

## 10. 결론

SimVascular는 여러 세대의 C++ 코드와 빌드 체계가 공존하는 대형 과학·의료 GUI
애플리케이션이다. 핵심 알고리즘, vtkSV 필터, MITK/Qt GUI, Python 바인딩이 잘
분리되어 있지만 실제 빌드 구성과 스타일은 완전히 통일되어 있지 않다.

작업 시 가장 중요한 원칙은 다음 세 가지다.

1. 대상 파일과 계층의 기존 관행을 우선한다.
2. GUI, module, core, Python 경로를 수직으로 함께 검토한다.
3. Mac 정적 검증과 Linux Docker 실빌드를 명확히 구분한다.

테스트 자동화, 외부 버전 단일화, CMake 오타 제거 및 경고 가시성 개선은 향후
저장소 안정성을 높이는 우선 과제다.
