# SimVascular 전체 GUI 빌드 (Docker / Ubuntu 22.04)

이 브랜치(메시 생성 신기능 포함)를 **GUI 포함**(Qt6 + MITK + sv4gui)으로
Docker에서 소스부터 빌드하기 위한 환경입니다. Windows Docker Desktop
(WSL2 백엔드) + D 드라이브 기준으로 작성했습니다.

## 왜 이런 구조인가 (배경 요약)

SimVascular의 공식 externals 배포 경로는 현재 전부 끊겨 있습니다:

- 소스/바이너리를 서빙하던 `simvascular.stanford.edu` 호스트 사망
- externals 빌드 스크립트 저장소 `SimVascular/svExternals` 삭제 (404)
- in-tree `Externals/CMake`는 관리자 개인 머신 절대경로
  (`/Users/parkerda/...`) 하드코딩 + 빌드 스텝 주석 처리 상태 (upstream
  master도 동일)
- GitHub 릴리스 자산은 소스 코드 스냅샷뿐 (externals 없음)

따라서 externals를 **각 upstream 공식 소스에서 직접 빌드**합니다.
빌드 플래그는 이 repo에 남아 있는 검증된 레시피
(`Externals/Make/2022.10/BuildHelpers/CompileScripts/*`)를 현재 버전에
맞게 이식한 것입니다.

## 버전 (repo가 요구하는 확정값)

`Code/CMake/SimVascularExternalsVersions.cmake`의 active(2024.05) 블록
기준 — `Docker/scripts/versions.env`가 단일 소스입니다:

| External | 버전 | 확보 방법 |
|---|---|---|
| Qt6 | 6.6.2 | aqtinstall (공식 Qt 바이너리 미러) |
| VTK | 9.3.0 | vtk.org 소스 |
| ITK | 5.4.0 | GitHub 릴리스 소스 |
| OpenCASCADE | 7.6.0 | GitHub(OCCT) 소스 |
| MITK | **v2023.12.2** ⚠️ | GitHub 소스 (아래 참고) |
| Python | 3.11.0 | python.org 소스 |
| HDF5 | 1.14.3 | GitHub 소스 |
| GDCM | 3.0.10 | GitHub 소스 |
| tinyxml2 | 8.0.0 | GitHub 소스 |
| freetype | 2.13.0 | savannah 소스 |
| MMG | 5.3.9 | GitHub 소스 (static) |

> ⚠️ **MITK 버전에 관한 판단**: repo의 버전 문자열은 "2022.10"이지만
> upstream MITK v2022.10은 Qt5 전용입니다. 이 트리는 Qt6 6.6.2
> (+Core5Compat)를 요구하고 "Fix code to use new MITK API" 커밋(862a1e8)이
> 있으므로, Qt6를 지원하는 첫 세대(v2023.04+) 이상이어야만 링크가
> 가능합니다. 기본값은 v2023.12.2이며, sv4gui 컴파일에서 MITK API
> 불일치가 나오면 **컴포넌트를 빼지 말고** 에러를 보고할 것
> (`versions.env`의 MITK_GIT_TAG로 v2024.06 / v2023.04 시도 가능).

## 요구 사양

- 디스크 **≥ 150 GB** (소스 + 빌드 트리 + 설치본)
- RAM **≥ 16 GB** (MITK/VTK 링킹; 부족하면 `NPROC`을 줄일 것)
- 전체 소요: 수 시간 (MITK가 가장 김; MITK superbuild는 빌드 중
  CTK/DCMTK/POCO 등을 추가 다운로드하므로 네트워크 필요)

## 사용법

### 0. D 드라이브 준비 (PowerShell)

```powershell
mkdir D:\sv
```

Docker Desktop → Settings → Resources: CPU 8+, Memory 16GB+, Disk 150GB+.

### 1. 이미지 빌드 (이 repo 체크아웃 위치에서)

```powershell
cd D:\sv
git clone https://github.com/halfmoon0109/SimVascular.git
cd SimVascular
git checkout claude/mesh-creation-code-review-2f7qol   # 또는 대상 브랜치
docker build -f Docker/Dockerfile -t simvascular-build:22.04 Docker
```

### 2. 컨테이너 실행 (D:\sv 를 /work 로 마운트)

```powershell
docker run -it --name sv -v D:\sv:/work simvascular-build:22.04 bash
```

### 3. 전체 빌드 (컨테이너 안에서)

```bash
bash /work/SimVascular/Docker/scripts/build-all.sh
```

단계별 로그는 `/work/logs/<단계>.log`. 특정 단계만 다시 돌리려면:

```bash
bash /work/SimVascular/Docker/scripts/40-vtk.sh   # 예: VTK만
```

**어느 단계든 실패하면 거기서 멈춥니다.** 해당 로그의 마지막 ~100줄을
가져오면 원인 분석 후 수정합니다. 컴포넌트를 빼고 우회하지 않습니다.

### 4. GUI 실행 (WSLg)

Docker Desktop(WSL2)은 WSLg X 서버를 쓸 수 있습니다. GUI 실행용으로는
컨테이너를 이렇게 띄웁니다:

```powershell
docker run -it --name sv-gui -v D:\sv:/work `
  -e DISPLAY=:0 `
  -v /run/desktop/mnt/host/wslg/.X11-unix:/tmp/.X11-unix `
  -v /run/desktop/mnt/host/wslg:/mnt/wslg `
  simvascular-build:22.04 bash
```

컨테이너 안에서:

```bash
export LD_LIBRARY_PATH=/work/externals/install/python/lib:$LD_LIBRARY_PATH
export QT_QPA_PLATFORM=xcb
/work/build/SimVascular-build/bin/simvascular   # 실제 launcher 경로는 빌드 후 확인
```

## 디렉토리 배치

```
D:\sv (= /work)
├── SimVascular/            이 repo (호스트에서 clone)
├── externals/
│   ├── src/                다운로드한 소스 tarball (재사용됨)
│   ├── build/              externals 빌드 트리 (완료 후 지워도 됨)
│   └── install/            ★ 설치본: qt6, vtk, itk, mitk, opencascade,
│                             python, gdcm, hdf5, freetype, tinyxml2, mmg, ml
├── build/                  SimVascular 빌드 (SimVascular-build/)
└── logs/                   단계별 로그
```

`install/<pkg>` 배치는 repo의 `run-cmake.sh`가 기대하는 레이아웃과
동일하며, `80-simvascular.sh`는 그 스크립트와 같은 `-DSV_*_DIR` 인자로
configure 합니다.

## 알려진 리스크 (막히면 여기부터 볼 것)

1. **MITK ↔ sv4gui API** — 위 박스 참고. 가장 가능성 높은 실패 지점.
2. **MITK superbuild vs 외부 VTK/ITK** — MITK가 자체 패치 VTK를 원할 수
   있음. 실패 시 fallback: `EXTERNAL_*_DIR` 없이 MITK가 VTK/ITK를 빌드하게
   하고, SimVascular의 `SV_VTK_DIR`/`SV_ITK_DIR`를 MITK `ep/` 트리로 변경.
3. **OpenCASCADE 7.6.0 + VTK 9.3** — in-repo 레시피는 VTK 9.2와 조합했음.
   IVtk 쪽 컴파일 에러 가능성 낮지만 존재.
4. **링커 OOM** — `NPROC=4 bash .../build-all.sh` 처럼 병렬도를 낮출 것.
