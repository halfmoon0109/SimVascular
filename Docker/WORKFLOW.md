# SimVascular Docker 빌드 — 전체 워크플로우 및 운영 가이드

이 문서는 이 저장소를 **인터넷 있는 개발 머신**에서 Docker(Ubuntu 22.04)로
**GUI 포함 전체 빌드**하고, 코드 수정 시 재빌드하고, 실행하는 전 과정을
정리한 것입니다. 실제 구축 과정에서 부딪힌 문제와 해결도 함께 기록합니다.

---

## 1. 목표

- 죽어버린 SimVascular 공식 externals 배포 경로(stanford 다운로드 호스트,
  `svExternals` 저장소)에 의존하지 않고,
- Qt6 / VTK / ITK / MITK / OpenCASCADE 등 **모든 external을 upstream 소스에서
  직접 빌드**하여,
- 이 브랜치의 **혈관 벽(solid wall) 메싱 신기능**(FSI용 boundary layer 유지 +
  벽 별도 region extrude, 전역/면별 두께, aspect-ratio 품질 리포트)을 포함한
  **GUI(Qt/MITK/sv4gui) 전체 빌드**를 만들고 실제 실행한다.

---

## 2. 구조 (2단계 빌드)

이 저장소는 "externals"와 "SimVascular 본체"를 **따로** 빌드한다.

1. **Stage 1 — externals**: `Docker/scripts/10..70-*.sh`
   각 external을 upstream 소스에서 받아 `/work/externals/install/<pkg>`에 설치.
2. **Stage 2 — SimVascular**: `Docker/scripts/80-simvascular.sh`
   `/work/SimVascular` 소스를, stage 1 결과물을 `-DSV_*_DIR`로 가리켜 빌드.
   (저장소 원래 `run-cmake.sh`와 동일한 인자 구성)

`build-all.sh`가 1→2를 순서대로 실행하며, **각 단계는 재개 가능**하다(이미
설치된 것은 자동 스킵).

### 확정 버전 (`scripts/versions.env`, `SimVascularExternalsVersions.cmake` 2024.05 블록)

| External | 버전 | 확보 |
|---|---|---|
| Qt6 | 6.6.2 | aqtinstall (공식 Qt 미러) |
| VTK | 9.3.0 | vtk.org |
| ITK | 5.4.0 | GitHub 릴리스 (+ GrowCut/IsotropicWavelets/OpenJPEG 모듈) |
| OpenCASCADE | 7.6.0 | Open-Cascade-SAS/OCCT |
| MITK | **v2024.06** | github.com/MITK/MITK (Qt6 지원 첫 릴리스) |
| Python | 3.11 | python.org |
| GDCM / HDF5 / tinyxml2 / freetype / mmg | 3.0.10 / 1.14.3 / 8.0.0 / 2.13.0 / 5.3.9 | upstream |

---

## 3. 최초 전체 빌드

### 3.1 준비 (호스트, Windows PowerShell)
```powershell
mkdir D:\sv; cd D:\sv
git clone https://github.com/halfmoon0109/SimVascular.git
cd SimVascular
git checkout <빌드할 브랜치>          # 예: master 또는 최신 기능 브랜치
docker build -f Docker/Dockerfile -t simvascular-build:22.04 Docker
```
> Docker Desktop → Settings → Resources: **CPU 8+, RAM 16GB+, Disk 150GB+**.
> `D:\sv` = 컨테이너의 `/work` 로 마운트된다.

### 3.2 빌드 (컨테이너)
```powershell
docker run -it --name sv -v D:\sv:/work simvascular-build:22.04 bash
```
```bash
bash /work/SimVascular/Docker/scripts/build-all.sh
```
- 오래 걸린다(MITK가 가장 김, 수 시간; 네트워크로 CTK/DCMTK 등 추가 다운로드).
- 창을 닫아도 계속 돌리려면:
  ```powershell
  docker exec -d sv bash -c "bash /work/SimVascular/Docker/scripts/build-all.sh > /work/logs/build.log 2>&1"
  docker exec -it sv tail -f /work/logs/build.log   # 진행 확인
  ```
- 어느 단계든 실패하면 그 자리에서 멈춘다. `/work/logs/<단계>.log` 확인.

---

## 4. 코드 수정 후 재빌드 ★ (가장 자주 쓰는 흐름)

수정 대상에 따라 방법이 다르다. **핵심 원칙: 소스는 호스트에서 pull → 컨테이너는
마운트로 자동 반영 → 하지만 재컴파일은 반드시 직접 실행해야 바이너리에 반영된다.**

### 4.1 SimVascular 소스(C++/sv4gui/.ui)만 바뀐 경우 — 대부분의 경우
```powershell
# 호스트: 최신 소스 받기
cd D:\sv\SimVascular
git fetch origin
git checkout <브랜치>        # 또는 git pull
```
```bash
# 컨테이너: SimVascular만 증분 재빌드 (externals는 그대로)
cd /work/build/SimVascular-build
export LD_LIBRARY_PATH="/work/externals/install/python/lib:$LD_LIBRARY_PATH"
make -j"$(nproc)"
```
> **Qt .ui 파일을 바꿨는데 UI가 안 바뀌면** (마운트 볼륨 타임스탬프 문제로 uic가
> 재생성을 건너뛴 경우), 해당 생성 헤더를 지우고 다시 make:
> ```bash
> rm -f /work/build/SimVascular-build/Source/sv4gui/Plugins/org.sv.gui.qt.meshing/ui_sv4gui_MeshEdit.h
> make -j"$(nproc)"
> # 검증: 새 위젯이 생성 헤더에 들어갔는지
> grep -c MeshingWallMeshWidgetPage \
>   /work/build/SimVascular-build/Source/sv4gui/Plugins/org.sv.gui.qt.meshing/ui_sv4gui_MeshEdit.h
> ```

`bash /work/SimVascular/Docker/scripts/80-simvascular.sh` 로도 되지만, 위처럼
`SimVascular-build`에서 직접 `make` 하는 게 빠르고 확실하다.

### 4.2 externals 설정(`versions.env`, `70-mitk.sh` 등)이 바뀐 경우
해당 external의 **build + install 디렉토리를 지워야** 새 설정이 반영된다
(재개 가드가 완성물이 있으면 스킵하기 때문). 예: MITK 옵션 변경 시
```bash
rm -rf /work/externals/build/mitk /work/externals/install/mitk
bash /work/SimVascular/Docker/scripts/build-all.sh   # MITK만 재빌드 후 SimVascular 재링크
```

---

## 5. 실행 (GUI)

디스플레이(WSLg) 연결이 있는 컨테이너에서:
```powershell
docker rm -f sv-gui
docker run -it --name sv-gui -v D:\sv:/work `
  -e DISPLAY=:0 `
  -v /run/desktop/mnt/host/wslg/.X11-unix:/tmp/.X11-unix `
  -v /run/desktop/mnt/host/wslg:/mnt/wslg `
  simvascular-build:22.04 bash
```
```bash
bash /work/SimVascular/Docker/scripts/run-gui.sh
```
`run-gui.sh`가 모든 런타임 환경(아래)을 자동 설정한다:
- `LD_LIBRARY_PATH` — 모든 externals의 lib
- `QT_PLUGIN_PATH` — Qt 플러그인(특히 CTK 플러그인 DB가 쓰는 sqlite 드라이버)
- `SV_PLUGIN_PATH` — SimVascular + MITK/BlueBerry 플러그인
- `LIBGL_ALWAYS_SOFTWARE=1` — WSLg 소프트웨어 OpenGL

### 빌드 반영 검증(실행 전 확인용)
```bash
# 벽 메싱 GUI가 컴파일됐는지
strings /work/build/SimVascular-build/lib/plugins/liborg_sv_gui_qt_meshing.so | grep -i "Wall Mesh"
# MITK Python 콘솔 플러그인이 빌드됐는지
find /work/externals/install/mitk -iname "*qt_python*"
```

---

## 6. 구축 중 부딪힌 문제와 해결 (트러블슈팅 이력)

| 단계 | 증상 | 원인 | 해결 |
|---|---|---|---|
| git | `build-all.sh`가 커밋 안 됨 | 저장소 `.gitignore`의 `build*` 규칙 | `git add -f` |
| git | Windows에서 `.sh` 실행/‑pull 충돌 반복 | `core.autocrlf`가 CRLF로 변환 | `.gitattributes`에 `Docker/** eol=lf` |
| mmg | 링크 시 `multiple definition` | mmg 5.3.9가 GCC≥10 `-fno-common` 이전 코드 | `-DCMAKE_C_FLAGS=-fcommon` |
| ITK | `find_package(Qt6)`/`Python3` 실패 | Qt/Python 래핑된 VTK config가 소비자에서 Qt6·Python3을 전이 탐색 | ITK/MITK에 `Qt6_DIR` + Python3 힌트 명시 |
| MITK | 태그 `v2023.12.2` 없음 | 존재하지 않는 버전 | 실제 태그 `v2024.06` (Qt6 첫 지원) |
| MITK | `Qt5` 요구 | v2023.x는 Qt5 전용 | v2024.06으로 상향 |
| MITK | `Qt6StateMachine` 없음 | aqt 기본 설치에 없음 | `qtscxml` 모듈 추가 |
| MITK | ACVD 서브빌드 Python3 실패 | 서브프로젝트가 시스템 3.10을 찾음 | 우리 Python 3.11을 `PATH` 최상단에 |
| MITK | `itkFastGrowCut.h`/wavelet 없음 | ITK remote 모듈 미활성 | `Module_GrowCut`/`IsotropicWavelets`/`OpenJPEG` ON |
| MITK | `fixup_bundle` READ_ELF 실패 | 외부 lib를 ldd로 못 찾음 | 모든 externals lib를 `LD_LIBRARY_PATH`에 |
| SV configure | `MITK_LIBS_MISSING`(DCMTK/Poco 등) | `cmake --install`은 MITK 자체 lib만 설치 | `MITK-build`+`ep` 병합 조립(원래 post-install 로직 이식) |
| SV configure | `mitkContourElement.h` 없음 | 헤더가 `Modules/*/include` 밖 하위폴더에 | 소스 트리 헤더 전체를 `include/mitk/`로 미러(FindMITK가 include를 글롭) |
| SV compile | `mitkITKImageImport.txx` 없음 | `.txx` 확장자 누락 | 헤더 미러에 `.txx` 포함 |
| 실행 | `QSQLITE driver not loaded` → abort | Qt sqldrivers 플러그인 경로 미지정 | `QT_PLUGIN_PATH` 설정 |
| 실행 | `libgdcmjpeg8.so` 못 찾음 | 런타임 lib 경로 누락 | 모든 externals lib를 `LD_LIBRARY_PATH`에 |
| 실행 | `org.mitk.gui.common` resolve 실패 | 플러그인 검색 경로 미지정 | `SV_PLUGIN_PATH`에 SV+MITK 플러그인 디렉토리 |
| GUI | Python Console 안 보임 | `org.mitk.gui.qt.python` 기본 OFF | `-DMITK_BUILD_org.mitk.gui.qt.python=ON` + MITK 재빌드 |
| GUI | Wall Mesh 체크박스 안 보임 | 재빌드가 옛 소스로 돎 | 최신 브랜치 pull 후 SimVascular 재빌드(필요 시 uic 강제 재생성) |

---

## 7. 요약 치트시트

```bash
# 최초 전체 빌드
bash /work/SimVascular/Docker/scripts/build-all.sh

# 코드만 수정 후 (빠른 증분)
cd /work/build/SimVascular-build && make -j"$(nproc)"

# externals 설정 수정 후 (예: MITK)
rm -rf /work/externals/build/mitk /work/externals/install/mitk
bash /work/SimVascular/Docker/scripts/build-all.sh

# 실행
bash /work/SimVascular/Docker/scripts/run-gui.sh
```
