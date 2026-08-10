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

# 고체 혈관벽 메시 생성 알고리즘

## 1. 목적과 범위

이 문서는 `thickness`를 이용해 고체 혈관벽(solid vessel wall) 메시를 만드는 현재
구현을 코드 호출 단위로 설명한다. 다음 내용을 포함한다.

- GUI와 Python API에서 코어 mesher까지 옵션이 전달되는 경로
- 원래 혈관 표면을 fluid/wall interface로 보존하는 방식
- point별 요청 두께와 최종 두께를 계산하는 알고리즘
- 분지 접합부에서 normal, 곡률, 두께, fold를 처리하는 방식
- 접합부 thinning의 원인이 형상인지 두께 패스인지 가르는 `t/R` 진단
- 바깥 면을 거리장 레벨셋으로 오프셋하고 벽을 TetGen으로 채우는 경로와, 그 경로가
  두께 축소 패스를 건너뛰는 이유
- VMTK가 표면을 여러 layer의 wedge로 압출하는 방식
- wedge를 tetrahedron으로 변환하는 방식
- TetGen 유체 코어, 유체 경계층, 고체벽을 병합하는 방식
- `ModelFaceID`, `ModelRegionID`, `GlobalNodeID`, `GlobalElementID` 처리
- 현재 알고리즘이 보장하지 않는 사항과 실행 후 검증 항목

전체 TetGen 워크플로의 요약은
[TetGen 메시 생성 로직](TetGen_Mesh_Generation_Workflow.md)을 참고한다.

## 2. 결과 메시의 공간적 배치

고체벽을 사용하는 전체 단면은 안쪽에서 바깥쪽으로 다음 순서다.

```text
혈관 중심
   │
   ├─ TetGen 유체 코어 tetrahedron
   │
   ├─ VMTK 유체 boundary-layer tetrahedron
   │
   ├─ 원래 모델 벽 표면
   │     = fluid/wall interface
   │     = 고체벽 압출의 시작 표면
   │
   ├─ 고체벽 요소
   │     = VMTK wedge 압출을 tet으로 변환(기본)
   │     = 또는 TetGen이 두 면 사이를 채운 tet(`WallMeshTetGenShell`, 12절)
   │
   └─ outer wall 표면
         = 안쪽 면을 normal 따라 압출한 면(기본)
         = 또는 거리장 오프셋 면(12절)
```

핵심 불변조건은 원래 모델 벽 표면의 point 좌표를 움직이지 않는 것이다. 유체
boundary layer는 이 표면에서 안쪽으로, 고체벽은 동일한 표면에서 바깥쪽으로
생성된다. 따라서 고체벽 두께를 스무딩하거나 줄여도 이동하는 것은 outer wall
정점이며 fluid/wall interface 정점은 그대로 남는다.

## 3. 코드 호출 경로

### 3.1 프런트엔드

GUI는 `sv4gui_MeshEdit.cxx`에서 다음 명령을 만든다.

```text
option GenerateWallMesh
option WallThickness <value>
option NumberOfWallLayers <value>
option WallThicknessSmoothingIterations <value>
option WallThicknessCurvatureFactor <value>
localWallThickness <face-name> <value>
```

`sv4gui_MeshTetGen.cxx`가 문자열 명령을 파싱해 코어의 `SetMeshOptions()`와
`SetWalls()` 호출로 변환한다.

Python API는 `MeshingTetGenOptions_PyClass.cxx`의 다음 속성을 같은 코어 옵션으로
변환한다.

| Python 옵션 | 코어 옵션 |
| --- | --- |
| `generate_wall_mesh` | `GenerateWallMesh` |
| `wall_thickness` | `WallThickness` |
| `number_of_wall_layers` | `NumberOfWallLayers` |
| `wall_thickness_smoothing_iterations` | `WallThicknessSmoothingIterations` |
| `wall_thickness_curvature_factor` | `WallThicknessCurvatureFactor` |
| `local_wall_thickness` | `LocalWallThickness` |
| `local_wall_thickness_on` | (Python 전용 boolean; `local_wall_thickness` 목록의 적용 여부) |

`MeshingTetGen_PyClass.cxx`의 `MesherTetGen_generate_mesh()`는 옵션을 설정한 뒤
wall face ID를 `SetWalls()`로 전달하고 `GenerateMesh()`를 호출한다.

### 3.2 코어 호출 그래프

```text
SetMeshOptions()
SetWalls()
SetBoundaryLayer()
        │
        ▼
cvTetGenMeshObject::GenerateMesh()
        │
        ├─ GenerateSurfaceRemesh()
        │
        ├─ GenerateBoundaryLayerMesh()
        │      ├─ normal 계산, surface clean
        │      ├─ 원래 interface surface 보관
        │      ├─ 유체 boundary layer 안쪽 압출
        │      └─ GenerateWallMesh()
        │             ├─ point별 WallThickness 계산
        │             ├─ 오목 영역 warp vector(normal) 스무딩
        │             ├─ t/R 진단(보고만, 수정 없음)
        │             ├─ 곡률 기반 두께 제한          ┐
        │             ├─ 두께 Laplacian 스무딩 + 재제한 │ wedge 경로 전용
        │             ├─ 두께 gradation 제한          │ (재제한/라운딩/fold)
        │             ├─ 바깥벽 라운딩(볼록 fillet)    │
        │             ├─ fold-over 방지               ┘
        │             │
        │             ├─ [wedge]  고체벽 바깥쪽 압출
        │             └─ [shell]  FillWallMeshWithTetGen()
        │                            ├─ 거리장 오프셋 바깥면
        │                            ├─ MMG 리메시
        │                            ├─ cap 평면 트리밍 + rim 스티칭
        │                            └─ TetGen shell 채움
        │
        ├─ GenerateAndMeshCaps()
        ├─ TetGen 내부 유체 tetrahedralization
        └─ AppendBoundaryLayerMesh()
               └─ VMTKUtils_AppendData()
                      ├─ 유체 코어 + 유체 BL 병합
                      ├─ 고체벽 wedge → tetrahedron
                      ├─ region/global ID 생성
                      └─ 유체 + 고체벽 conformal 병합
```

`GenerateWallMesh()`는 독립적인 최상위 실행 경로가 아니다.
`GenerateBoundaryLayerMesh()` 내부에서 호출되므로 고체벽 생성에는 surface mesh,
volume mesh, 유체 boundary layer가 모두 필요하다.

## 4. 옵션과 내부 변수

옵션은 `sv_TetGenMeshObject.h`의 `TGoptions`와
`sv_TetGenMeshObject.cxx`의 `localWallThickness_` 등에 저장된다.

| 외부 옵션 | 내부 변수 | 기본값 | 역할 |
| --- | --- | ---: | --- |
| `GenerateWallMesh` | `wallmeshflag` | `0` | 고체벽 생성 활성화 |
| `WallThickness` | `wallthickness` | `0.0` | global 요청 두께 |
| `NumberOfWallLayers` | `numwallsublayers` | `2` | 벽 두께 방향 layer 수 |
| `WallThicknessSmoothingIterations` | `wallthicknesssmoothingiterations` | `5` | normal 및 두께 스무딩 반복 수 |
| `WallThicknessCurvatureFactor` | `wallthicknesscurvaturefactor` | `0.8` | 오목 곡률 반경 대비 허용 두께 비율 |
| `LocalWallThickness` | `localWallThickness_` | 비어 있음 | `face ID → 두께` override |
| `WallMeshTetGenShell` | `walltetgenshell` | `0` | 벽을 wedge 압출 대신 TetGen tetrahedron으로 채움(12절) |
| `BoundaryLayerDirection` | `boundarylayerdirection` | `1` | 유체 boundary layer를 안쪽으로 압출 |

`WallMeshTetGenShell`은 벽의 **바깥 면을 만드는 방식**을 바꾼다. 꺼져 있으면
9~11절과 13절의 wedge 압출 경로, 켜져 있으면 12절의 거리장 오프셋 경로다. 두
경로 모두 안쪽 면(fluid/wall interface)의 좌표는 건드리지 않으므로 선택은
자유롭다.

`SetMeshOptions()`는 공통 진입점에서 값을 검증한다.

- global 두께는 유한한 0 이상의 값이어야 한다. 벽 생성을 켰다면
  `GenerateMesh()`에서 다시 `> 0`인지 확인한다.
- local 두께는 유한한 양수여야 하고 face ID는 정수여야 한다.
- wall layer 수는 1 이상의 정수여야 한다.
- smoothing 반복 수는 0~50 사이 정수다.
- curvature factor와 radius factor는 0~1 범위다.

`GenerateMesh()`는 벽 생성 전에 다음 조합도 검증한다.

1. surface mesh와 volume mesh가 모두 활성화되어야 한다.
2. boundary-layer mesh가 활성화되어야 한다.
3. `BoundaryLayerDirection == 1`이어야 한다.
4. `WallThickness > 0`이어야 한다.

## 5. 입력 표면과 식별 배열

### 5.1 wall face 선택

`SetWalls(numWalls, walls)`는 입력 `polydatasolid_`의 cell data
`ModelFaceID`를 읽는다. 지정된 face에는 `WallID=1`, cap을 포함한 나머지 face에는
`WallID=0`을 기록한다.

- VMTK remesher를 사용하는 경로는 `WallID=1`인 cell만 즉시 추출한다.
- MMG 경로는 remesh 중 face 배열을 유지한 다음 wall 영역을 사용한다.
- `wallFaceIDs_`에도 선택한 face ID를 저장한다.

고체벽 알고리즘은 서로 떨어진 wall face들을 나중에 boolean union하지 않는다.
분지를 포함하는 입력 혈관 모델이 이미 하나의 연결된 표면 topology를 가져야 한다.

### 5.2 interface 표면 준비

`GenerateBoundaryLayerMesh()`는 재메시된 `polydatasolid_`에 다음 처리를 한다.

1. `vtkPolyDataNormals`로 일관된 outward point normal을 계산한다.
   `Consistency=1`, `AutoOrientNormals=1`, `SplittingOff`이므로 face 경계에서도
   공유 point를 분리하지 않는다.
2. `vtkCleanPolyData`로 coincident point를 정리한다.
3. 결과를 `originalsurfpd`에 deep copy한다.
4. `TGenUtils_ReportSurfaceTriangleQuality(originalsurfpd, "fluid/wall interface surface")`로
   삼각형 품질을 보고한다. 이 표면은 유체 BL, TetGen 유체 체적, 고체벽 압출이 모두
   파생되는 공유 표면이므로 여기의 sliver는 셋 모두에 들어간다. 점 병합이 삼각형을
   바꿀 수 있어 clean **뒤에** 보고한다.
5. `MeshSizingFunction` point 배열을 만든다.
6. `SetCapBoundaryNormals()`로 열린 cap rim의 normal을 cap 평면 방향에 맞춘다.

이 `originalsurfpd`가 유체 boundary layer와 고체벽 양쪽의 공통 시작 표면이다.
유체 boundary layer를 압출하기 전에 별도 copy를 보관하므로, 유체 쪽 압출 결과가
고체벽의 시작 좌표를 바꾸지 않는다.

`SetCapBoundaryNormals()`는 각 열린 경계 loop의 중심을 구한 뒤 경계점에 대해 대략
다음 방향을 쓴다.

```text
capNormal_i = normalize(boundaryPoint_i - loopCenter)
```

따라서 끝단의 고체벽은 축 방향으로 기울기보다 cap 평면 안에서 바깥으로 확장된다.

### 5.3 주요 데이터 객체와 배열

| 이름 | 형식 | 의미 |
| --- | --- | --- |
| `originalsurfpd` | `vtkPolyData` | 정리된 원래 fluid/wall interface |
| `boundarylayermesh_` | `vtkUnstructuredGrid` | VMTK 유체 boundary layer |
| `innerblmesh_` | `vtkUnstructuredGrid` | 압출 전 interface copy |
| `polydatasolid_` | `vtkPolyData` | 유체 BL의 안쪽 표면 및 cap, TetGen 입력 |
| `wallmesh_` | `vtkUnstructuredGrid` | VMTK가 만든 고체벽 surface/volume 혼합 grid |
| `volumemesh_` | `vtkUnstructuredGrid` | 최종 유체+고체 체적 메시 |
| `WallID` | cell data | wall `1`, cap `0` |
| `ModelFaceID` | cell data | 모델 face 식별자 |
| `Normals` | point data | 압출 방향 |
| `MeshSizingFunction` | point data | 유체 BL 두께 또는 크기 기준 |
| `WallThickness` | point data | 고체벽 point별 최종 압출 거리 |
| `CellEntityIds` | cell data | BL volume/surface/sidewall 구분 |
| `ModelRegionID` | cell data | 유체와 고체 재료 영역 구분 |
| `GlobalNodeID` | point data | 병합된 전역 node 번호 |
| `GlobalElementID` | cell data | 병합된 전역 element 번호 |

## 6. 유체 boundary layer와 interface 보존

고체벽보다 먼저 유체 boundary layer를 생성한다. `Normals`를 반전해 원래 표면에서
혈관 안쪽으로 압출한다.

- constant thickness가 켜지면 전체 두께는
  `maxedgesize × blthicknessfactor`다.
- 꺼지면 point별 `MeshSizingFunction × blthicknessfactor`를 쓴다.
- `numsublayers`와 `sublayerratio`가 두께 방향 분할을 정한다.

VMTK가 반환한 가장 안쪽 표면은 local 변수 `innerSurface`에 있고, 이 표면이
`polydatasolid_`로 전달되어 cap 생성과 TetGen 유체 코어의 경계가 된다.
`innerblmesh_`는 압출 전 surface의 보관본이며 최종 inner surface 자체는 아니다.

그 다음 `GenerateWallMesh(originalsurfpd, markerListName)`를 호출한다. 즉 고체벽은
유체 boundary layer의 안쪽 표면이 아니라 원래 interface에서 시작한다.

## 7. point별 벽 두께 계산

`GenerateWallMesh()`는 입력 surface를 deep copy하고, point 수와 같은 길이의
`vtkDoubleArray WallThickness`를 만든다. 이후 두께는 다음 순서로 변한다.

```text
                       wedge 압출 경로              TetGen shell 경로
global/local 요청 두께        │                            │
  → 곡률 제한                 ●                            ─ (건너뜀)
  → Laplacian 두께 스무딩      ●                            ●
  → 곡률 재제한               ●                            ─ (건너뜀)
  → gradation 제한            ●                            ●
  → 바깥벽 라운딩             ●                            ─ (건너뜀)
  → fold 방지 축소            ●                            ─ (건너뜀)
  → 최종 gradation 제한       ●                            ─ (건너뜀)
                          = 최종 WallThickness          = 요청 두께 그대로
```

이와 별개로 오목 영역 warp vector(normal) 스무딩(8절)이 두께 계산 직후·곡률 제한
직전에 실행된다. 이 단계는 두께 값이 아니라 압출 방향(normal)만 바꾼다. 그 직후
곡률 제한 직전에 `t/R` 진단(9.2절)이 실행되지만 이것은 보고만 하며 위 흐름에
개입하지 않는다. 바깥벽
라운딩(11.1절)은 두께를 유지하면서 바깥면을 밀어내고 그 달성 거리를 `WallThickness`
값과 normal에 다시 인코딩하므로, 이후 fold 방지가 처리할 잔여 fold만 남긴다.

**건너뛰는 패스들의 공통점**은 하나다. 전부 1:1 바깥쪽 압출을 유효하게 만들려고
존재하며, 그 대가를 두께로 치른다. shell 경로는 압출하지 않으므로 막을 것이 없다.
남는 gradation 제한만이 압출의 유효성이 아니라 두께장 자체의 성질(급변 제거)을
다루기 때문에 양쪽 모두에서 실행된다. 자세한 근거는 12.2절.

`gradation 제한`(`TGenUtils_LimitThicknessGradation`)은 두께의 기울기를
`maxSlope = 0.5`(약 26.6도 taper)로 제한한다. 값을 낮추기만 하므로 곡률 상한과
local 두께를 천장으로 유지하며, 재클램프가 필요 없다. 이것이 필요한 이유는
fold 방지가 접힌 삼각형을 그 최솟값으로 평준화해서 절벽 하나가 다음 절벽을
만들기 때문이다 — 실측으로 감당 불가 접합부 15곳이 얇아진 구역 37곳이 되었다.

### 7.1 global/local 요청 두께

point `i`의 요청 두께를 `t_req(i)`라 하자. 먼저 모든 point를 global 두께
`t_global`로 초기화한다.

local 두께가 있으면 각 point가 속한 cell의 고유 `ModelFaceID` 집합 `F_i`를 만든다.
각 face의 두께는 local override가 있으면 그 값, 없으면 global 값을 사용한다.

```text
                 1
t_req(i) = ------------- × Σ t_face(f)
             |F_i|         f∈F_i
```

중요한 점은 cell 개수가 아니라 고유 face ID별로 한 번만 합산한다는 것이다.
따라서 한 face가 더 촘촘하게 삼각분할되어도 공유점의 두께 평균을 더 많이
가중하지 않는다.

예를 들어 두 face의 요청 두께가 `0.5`와 `1.0`이고 한 point에서 만난다면 그
point의 초기 요청 두께는 `0.75`다. 어느 face에 triangle이 더 많은지는 무관하다.

## 8. 오목 영역 warp vector(normal) 스무딩

두께 계산 직후, 곡률 제한 전에 `TGenUtils_SmoothWarpVectorsInConcaveRegions()`가
오목 영역의 `Normals`(압출 warp vector)를 스무딩한다. 표면을 normal 방향으로
바깥 압출하면 오목 영역(두 혈관이 만나는 crotch 등)의 warp vector가 서로 수렴하므로,
압출된 outer wall이 안쪽으로 파이고 그 부근의 벽 요소가 뒤틀린다. 두께가 이 뒤틀림의
원인은 아니며 두께를 줄여도 뒤틀림은 사라지지 않는다. 원인은 방향장이므로 각 오목
point의 normal을 one-ring 이웃 normal의 평균 쪽으로 완화해 수렴한 방향을 벌린 뒤
재정규화한다.

- 각 point는 **오목한 정도에 비례**해 완화한다. 오목도는 point의 접평면 위로 올라온
  이웃들의 rise 각 sine 평균(무차원)으로, 볼록·평평한 point에서는 0이다. 따라서
  볼록·평평한 영역과 직관(straight tube)의 normal은 그대로 유지된다.
- boundary edge에 놓인 point(cap rim, `SetCapBoundaryNormals`로 cap 평면 안에 놓인
  normal)는 pin 처리해 스무딩하지 않으므로 벽이 cap에서 평평하게 유지된다.
- 완화 계수는 `warpVectorRelaxation = 0.5`, 반복 수는
  `WallThicknessSmoothingIterations`다.

**두께 값과 표면 point 좌표(fluid/wall interface)는 이 단계에서 바뀌지 않는다.**
바뀌는 것은 normal 방향뿐이다. 로그에는
`Smoothed the wall extrusion warp vectors at N concave points (max direction change X degrees)`
형태로 스무딩한 point 수와 최대 방향 변화가 출력된다.

이 단계는 압출 방향의 뒤틀림/함몰을 다루고, 접합부의 **두께 감소**(완만한 thinning)는
곡률 제한(9.1절)이, 접합부 두께 **보존**은 바깥벽 라운딩(11.1절)이 담당한다. 각 단계의
역할이 다르므로 순서대로 적용된다.

## 9. 오목 곡률: 두께 제한과 t/R 진단

이 절의 두 단계는 같은 one-ring 곡률 추정을 공유하지만 역할이 다르다. 9.1의
클램프는 두께를 실제로 줄이고, 9.2의 진단은 아무것도 바꾸지 않고 요청 두께가
형상이 감당할 수 있는 범위인지만 보고한다. 실행 순서는 진단(9.2)이 먼저이고
클램프(9.1)가 나중이다. 진단이 두께 감소 전 값을 봐야 하기 때문이다.

### 9.1 곡률 기반 두께 제한

`TGenUtils_ClampThicknessToConcaveCurvature()`는 한 point의 one-ring 이웃에서
오목 곡률을 근사한다. 이웃까지의 제곱 거리 `d²`와 normal 방향 높이 `h`를 사용해
다음 값을 계산한다.

```text
κ_ij ≈ 2h_ij / d_ij², h_ij > 0
κ_i = max(κ_ij)
t_limit(i) = curvatureFactor / κ_i
t(i) = min(t(i), t_limit(i))
```

`κ_i = 0`인 볼록/평평한 point는 변경하지 않는다. factor가 0이면 이 단계 전체를
건너뛴다.

이 제한의 의미는 벽 두께가 국소 오목 곡률 반경의 일정 비율을 넘지 않게 하는
것이다. 제한은 현재 값보다 작은 경우에만 적용되므로 반복 호출해도 값을 다시
키우지 않는다.

### 9.2 요청 두께 대비 곡률 반경(t/R) 진단

`TGenUtils_ReportConcaveCurvatureVsThickness()`는 두께를 줄이는 어떤 패스보다
먼저 실행되어 요청 두께 `t`를 국소 오목 곡률 반경 `R`과 비교해 보고한다. 두께
배열도 surface도 수정하지 않는다.

접합부가 얇아졌을 때 원인이 **형상이 요청 두께를 감당할 수 없어서**인지 **두께
패스가 과하게 깎아서**인지를 가르는 값은 `t/R` 하나다. `t > R`이면 오목부를
normal 방향으로 offset한 바깥면은 반드시 자기교차하므로, 어떤 오프셋 기반
방법으로도 그 두께를 실을 수 없다. 반대로 `t < R`인데도 두께가 줄었다면 원인은
형상이 아니라 두께 패스 쪽에 있다.

기존 패스는 이 값을 간접적으로만 드러낸다. 곡률 클램프(9.1)는 factor가 `0`이면
함수 본문이 통째로 건너뛰어져 아무것도 계산하지 않고, fold 방지(11.2절)는 자기가
깎은 두께만 보고할 뿐 그렇게 강제한 비율은 보고하지 않는다. 게다가 fold 방지는
0.8배 계단식 축소와 세 정점 평준화를 쓰므로, 깎인 비율에서 `R`을 역산하면 실제보다
작게 나온다. 그래서 직접 측정이 필요하다.

#### 곡률의 두 요약

곡률 추정은 9.1과 같은 `2h/d²`를 쓰되 one-ring을 두 가지로 요약한다.

| 요약 | 정의 | 성질 |
| --- | --- | --- |
| `R_smallest` | `1 / max(κ_ij)` | 클램프가 쓰는 값과 동일. 분모가 `d²`이라 퇴화 삼각형 하나에 발산 |
| `R_typical` | `1 / median(κ_ij)` | 접평면 아래 이웃을 곡률 `0`으로 포함한 중앙값. 이웃 하나에 흔들리지 않음 |

`R_typical`이 중앙값이면서 비오목 이웃을 `0`으로 포함하는 이유는, 그래야 이 값이
**point 전체의 성질**이 되기 때문이다. 이웃 하나만 우연히 접평면 위로 올라온
평평한 지점은 중앙값이 `0`에 머물러 접합부로 오분류되지 않고, one-ring 대부분이
실제로 오목한 지점에서만 큰 값이 된다.

두 값을 함께 읽으면 형상 문제와 메시 문제가 분리된다. 둘 다 작으면 진짜로 날카로운
crotch이고, `R_typical`은 큰데 `R_smallest`만 작으면 형상이 아니라 퇴화 삼각형이
만든 아티팩트다. 합성 검증에서 오목 원호는 `R`을 오차 `5e-17`로 복원하고, 퇴화
이웃을 하나 추가하면 `R_smallest`는 `0.3 → 0.000505`로 붕괴하는 동안 `R_typical`은
변하지 않는다.

퇴화 삼각형(최소 altitude가 최장 edge의 5% 미만)을 쓰는 point는
`[near-degenerate triangle]`로 따로 표시되므로, sliver가 원인인 구역과 형상이
원인인 구역이 같은 목록에서 섞이지 않는다.

#### 호출 위치

warp vector 스무딩(8절) **직후**, 곡률 클램프(9.1) **직전**이다. 두 조건이 이
위치를 강제한다.

1. 높이 `h`를 재는 기준 normal이 실제 압출에 쓰이는 방향과 같아야 하므로 warp
   스무딩보다 뒤여야 한다.
2. 보고하는 `t`가 순수 요청값이어야 하므로 두께를 줄이는 모든 패스보다 앞이어야
   한다.

## 10. 두께 Laplacian 스무딩과 재클램프

local 두께, curvature factor, radius factor 중 하나라도 활성화되면
`TGenUtils_SmoothPointArray()`가 `WallThickness`를 스무딩한다. 기본 curvature
factor가 `0.8`이므로 기본 설정에서는 이 경로가 실행된다.

point `i`의 고유 one-ring 이웃 집합을 `N_i`라 하면 한 iteration은 다음과 같다.

```text
              t_old(i) + Σ t_old(j)
                         j∈N_i
t_new(i) = --------------------------
                    1 + |N_i|
```

자기 자신도 평균에 포함하며, 모든 point를 이전 iteration에서 동시에 계산하는
Jacobi 방식이다. 따라서 local face 경계의 계단과 곡률 제한 부근의 급격한 두께
차이가 완만해진다.

스무딩은 곡률로 낮춘 값을 이웃 값 쪽으로 다시 올릴 수 있다. 그래서 스무딩 뒤
`TGenUtils_ClampThicknessToConcaveCurvature()`를 한 번 더 호출해 곡률 상한을
복원한다.

## 11. 바깥벽 라운딩과 fold-over 방지

> 이 절 전체는 **wedge 압출 경로 전용**이다. `WallMeshTetGenShell`이 켜져 있으면
> 두 단계 모두 실행되지 않는다(12.2절).

두께가 최종 확정된 뒤, VMTK 압출 전에 두 단계가 순서대로 실행된다. 먼저 라운딩이
접합부 두께를 유지하면서 바깥면을 밀어 함몰을 메우고(11.1), 그다음 fold 방지가
어떤 방향으로도 벽을 실을 수 없는 퇴화 입력에 대해서만 안전망으로 두께를 깎는다
(11.2).

### 11.1 바깥벽 라운딩(볼록 fillet)

접합부 함몰의 근본 원인은, 안쪽면을 normal 따라 두께만큼 바깥으로 offset하는
방식에서 오목 crotch의 normal이 수렴해 순진한(naive) 바깥면이 자기교차하는 것이다.
모든 point가 배정 두께에 있어도 바깥면은 겹친다. fold 방지(11.2)는 이 자기교차를
두께를 깎아 해소했고, 그 결과 벽이 덜 뻗어 안으로 파였다.

`TGenUtils_RoundOuterWallToPreserveThickness()`는 대신 두께를 유지하면서 바깥면을
바깥으로 밀어 매끄러운 볼록 fillet을 만든다(두꺼운 용접 접합부의 바깥이 재료로
채워지는 것과 같음). 안쪽면(fluid/wall interface)은 고정, 바깥면 point만 이동한다.

- 각 바깥 point를 이웃 바깥 point 평균 쪽으로 relax한다. 함몰의 이웃은 더 바깥에
  있으므로 골을 메운다. 이동량은 point의 오목도에 비례하므로 볼록·평평·직관은
  불변이고, cap rim point는 pin한다.
- 이동 후 안쪽 point와의 normal 방향 거리가 배정 두께 아래로 내려가지 않게 바깥으로
  되민다.
- 아주 뾰족한 crotch가 무한정 튀지 않게 `maxFilletRatio = 3.0`(배정 두께의 3배)으로
  상한을 둔다.
- 완화 계수는 `outerRoundingRelaxation = 0.5`, 반복 수는
  `WallThicknessSmoothingIterations`다.

라운딩 결과는 `Normals`(압출 방향)와 두께 배열(압출 크기)로 다시 인코딩되므로 기존
VMTK 압출이 그대로 이 바깥면을 재현한다. 입력 두께 배열은 **달성한 outer 거리(최소
배정 두께 이상)**로 덮어써진다. VMTK의 tangle 기준(뒤집힘/면적 붕괴)과 동일한
기하이므로 fold-free한 라운딩 결과는 untangle이 건드리지 않는다.

퇴화 입력 sliver는 어떤 방향으로도 벽을 실을 수 없으므로 라운딩만으로는 그 fold를
제거하지 못한다. 그래서 뒤의 fold 방지(11.2)가 안전망으로 남는다. 로그에는
`Wall outer rounding: filled the junction depression by raising N concave points; largest fillet Xx the assigned thickness at (...)`
형태로 라운딩량과 최대 fillet 위치가 출력된다.

### 11.2 실제 outer triangle 기반 fold-over 방지

국소 곡률 근사는 거친 표면이나 sliver triangle에서 fold를 놓칠 수 있다.
`TGenUtils_LimitThicknessToPreventFold()`는 VMTK 호출 전에 예측 outer surface를
직접 구성해 검사한다.

point `i`의 예측 outer vertex는 다음과 같다.

```text
q_i = p_i + t_i × normalize(n_i)
```

각 inner triangle과 대응하는 outer triangle의 unit face normal을 각각
`N_inner`, `N_outer`라 한다. outer triangle이 퇴화했거나 다음 조건을 만족하면
fold로 판정한다.

```text
dot(N_inner, N_outer) <= 0.1
```

완전히 뒤집힌 경우뿐 아니라 거의 붕괴한 triangle도 검출하기 위해 임계값이 작은
양수다.

fold된 triangle의 세 point 두께 중 최소값을 `t_fold`라 하면 관련 point의 새
두께 후보는 다음과 같다.

```text
t_reduced(i) = 0.8 × min(t_current(i), t_fold(i))
t_new(i) = max(t_min(i), t_reduced(i))
```

먼저 세 두께를 최소값 쪽으로 평준화한 뒤 0.8배 줄인다. 두께 차이 자체가 fold를
만든 경우 세 값을 같은 비율로만 줄이면 차이의 비율이 계속 남으므로 평준화가
필요하다. 이 검사와 축소를 최대 30회 반복한다.

point별 하한은 fold 검사 시작 시 두께 `t_start(i)`와 그 point를 쓰는 삼각형의
최소 altitude `a_min(i)`로 계산한다. `t_start`는 요청 두께가 아니라 바깥벽
라운딩(11.1절)까지 끝난 뒤의 값이다. 라운딩은 두께 배열을 달성한 outer 거리로
덮어쓰므로, 라운딩이 크게 들어간 point에서는 이 하한도 그만큼 높아져 fold 방지가
덜 공격적으로 동작한다. 라운딩이 이미 fold를 해소했다는 전제에서는 정합적이지만,
sliver처럼 라운딩으로 풀리지 않는 지점이 라운딩을 받았다면 유의해야 한다.

```text
triangle altitude = 2 × area / longestEdge
t_min(i) = min(0.05 × t_start(i), 0.5 × a_min(i))
```

정상 triangle에서는 5% 하한이 주로 결정하고, altitude가 매우 작은 sliver에서는
기하 하한이 더 작아져 fold를 해소할 여지를 준다. 하한과 30회 반복 뒤에도 fold가
남으면 오류를 숨기지 않고 최대 10개 위치에 대해 edge 길이, area, altitude,
두께, 하한, point normal의 최소 dot product를 출력한다.

이 검사는 한 triangle의 winding만 보는 국소 검사다. 떨어진 두 표면 구간이
압출 후 서로 충돌하는 전역 자기 교차는 검출하지 못한다.

## 12. 거리장 오프셋 바깥면과 TetGen shell 채움

`WallMeshTetGenShell`이 켜지면 고체벽은 wedge 압출이 아니라 안쪽 면과 바깥 면
사이를 TetGen tetrahedron으로 채워 만든다. 이 경로는 8~10절의 두께 계산과
gradation 제한까지만 공유하고, 11절(라운딩, fold 방지)과 13절(VMTK 압출)을
대체한다.

### 12.1 왜 대응을 버리는가

wedge 압출은 바깥 절점 하나를 안쪽 절점 하나에 묶는다. 고체를 두께 `t`만큼
팽창시킨 경계(Minkowski dilation)를 생각하면 이 묶음이 왜 표현력이 부족한지가
나온다.

- 볼록한 곳에서 바깥 면은 반경 `t`로 둥글어진다.
- 오목한 crotch에서는 두 오프셋 시트가 서로를 파고들고, 그 **교차선(crease)**이
  경계다. 교차선 너머의 시트는 경계가 아니다.

즉 오프셋이 잘려나가는 쪽에 놓이는 안쪽 점들은 **대응하는 바깥 점이 아예 없는
것이 정답**이다. 1:1 대응을 유지한 채로는 이 배치가 존재하지 않는다. 라운딩에
클리어런스 구속을 넣었던 시도가 수렴하지 못한 것은 튜닝 문제가 아니라 고정된
삼각분할 위에서 존재하지 않는 답을 찾았기 때문이며, 그래서 되돌렸다.

`d(x) = t` 레벨셋에는 지켜야 할 대응이 없다. crease와 라운딩은 겨냥하는 것이
아니라 구성에서 나오고, 결과의 모든 점은 안쪽 면에서 최소 `t` 떨어져 있다.

### 12.2 두께 축소 패스를 건너뛴다

`GenerateWallMesh()`의 지역 플래그 `extrudeWedges = !walltetgenshell`가 다음
패스를 shell 경로에서 끈다.

| 패스 | shell 경로 | 이유 |
| --- | --- | --- |
| 곡률 클램프(9.1)와 재클램프 | 건너뜀 | 오프셋은 `t > R`에서도 자기교차하지 않는다 |
| 바깥벽 라운딩(11.1) | 건너뜀 | 오프셋이 이미 crease/fillet을 만든다 |
| fold 방지(11.2) | 건너뜀 | 막을 fold가 없다 |
| 최종 gradation 제한 | 건너뜀 | 위 패스가 만드는 절벽이 없다 |
| `TGenUtils_ReportAchievedWallThickness` | 건너뜀 | 만들어지지 않는 면을 측정한다 |

유지하는 것은 요청 두께 계산(7절), warp vector 스무딩(8절), `t/R` 진단(9.2절),
두께 Laplacian 스무딩(10절), 그리고 요청 두께에 대한 gradation 제한이다.
gradation 제한만 남는 이유는 그것이 압출의 유효성이 아니라 두께장 자체의
성질이기 때문이다.

### 12.3 오프셋 바깥면 생성

`TGenUtils_BuildOffsetOuterSurface()`가 다음 순서로 만든다.

1. **cap을 닫는다.** 부호 있는 거리를 정의하려면 닫힌 면이 필요하다.
   `TGenUtils_ExtractBoundaryLoops()`가 경계 edge(셀 하나만 쓰는 edge)를 소유 셀의
   winding 순서로 걸어 rim 루프를 만들고, 각 루프를 중심점으로의 삼각형 팬으로
   닫는다. 팬은 rim edge를 벽 삼각형과 반대 방향으로 traverse하므로 닫힌 면의
   방향이 일관된다.
2. **격자.** 간격은 `min(최소 두께, maxedgesize) / 2`이고, voxel 예산
   (`FillWallMeshWithTetGen`이 6400만을 넘긴다)을 넘으면 맞을 때까지 거칠어진다.
   거칠어졌으면 그 사실과 실제 해상도를 로그에 출력한다. 격자는 모델 경계상자를
   `1.25 × 최대두께 + 5 × 간격`만큼 넓힌 범위를 덮는다.
3. **밴드.** 각 삼각형의 경계상자를 `1.25 × 국소두께 + 3 × 간격`만큼 부풀린
   범위의 voxel만 표시한다. 거리는 이 밴드에서만 평가한다.
4. **값.** 밴드 voxel마다 `vtkImplicitPolyDataDistance`의 부호 거리 `d(x)`와,
   `vtkStaticPointLocator`로 찾은 **안쪽 면**(팬 중심점 제외) 최근접 점의 두께
   `t`를 써서 `d(x) - t`를 넣는다. 부호 규약은 단정하지 않고 격자 모서리에서
   실측해 정한다. 모서리 거리의 크기까지 확인해 "뒤집힌 것"과 "격자 가정이
   깨진 것"을 구분한다.
5. **밴드 밖.** 격자 모서리에서 flood fill한 voxel은 `+LARGE`, 도달하지 못한
   voxel은 lumen 안쪽이므로 `-LARGE`.
6. **검증.** 밴드 voxel이 반대 부호의 채워진 voxel과 맞닿아 있으면 레벨셋이
   밴드를 벗어난 것이므로 오류로 보고한다. 그대로 컨투어하면 오프셋이 아니라
   밴드 가장자리를 따라가는 잘못된 면이 나온다.
7. **컨투어.** `vtkFlyingEdges3D`로 `0`을 컨투어한다. 연결된 shell 개수를
   보고한다(2개 이상이면 lumen을 막았거나 모델이 분리되어 있다는 뜻).

### 12.4 리메시, 트리밍, 스티칭

컨투어는 격자 해상도이고 marching cubes 특유의 sliver가 많으므로, `MMG`로
`maxedgesize`에 맞춰 리메시한다(각도 임계 45도가 접합부 crease를 ridge로
보존한다). MMG 없이 빌드된 경우 리메시를 건너뛰고 그 사실을 로그에 남긴다.

**리메시가 트리밍보다 먼저다.** 반대로 하면 리메시가 trimmed rim을 cap 평면 밖으로
옮겨 안쪽 rim과의 연결이 깨진다.

`TGenUtils_TrimOffsetSurfaceAtCaps()`가 cap 위의 돔을 잘라낸다.

- 자르는 면은 cap 평면이되 **그 cap 근처에서만** 자른다. 무한 평면은 반대편에
  우연히 놓인 다른 부위까지 자르며, 되꺾이는 혈관에서는 그것이 실제 벽이다.
  유지 조건을 점 스칼라 `max(평면 아래, rim 반경 2배 밖)`로 만들어
  `vtkClipPolyData`에 넘긴다.
- 평면 방향은 point normal이 아니라 rim 자체에서 구한다. `SetCapBoundaryNormals()`가
  이미 normal을 cap 평면 안에 눕혀 놓았기 때문이다. rim은 벽 삼각형의 winding으로
  걸리고 그 순서는 바깥 방향에 대해 시계방향이므로, **바깥 방향은 rim의 Newell
  법선의 반대**다.
- 트리밍 후 rim 개수가 cap 개수와 다르면 오류다. 잘림이 돔 이상을 가져간
  경우이며, 조용히 넘어가면 벽에 구멍이 남는다.

`TGenUtils_StitchCapAnnulus()`가 각 혈관 끝에서 두 rim 사이를 삼각분할한다. 두
rim은 점도 점 개수도 다르므로 1:1 스트립을 쓸 수 없다. 두 rim 모두 같은 끝을
감으므로 cap 축에 대한 각도가 둘 다를 정렬하고, 각 단계에서 뒤처진 rim을
전진시키며 병합한다. 각 rim이 축을 정확히 한 바퀴 감는지 검사하고 아니면 오류다
(되꺾이는 rim에는 이 순서가 없어 겹친 facet이 된다). 삼각형 방향은 경우별로
따지지 않고 바깥 방향과의 내적으로 측정해 맞춘다. 면적 0 삼각형은 버리지 않고
세어서 보고한다 — 끝면에 구멍을 남기는 쪽이 더 나쁘다.

### 12.5 shell 구성과 TetGen 채움

`TGenUtils_BuildWallShellSurface()`가 다음을 합쳐 닫힌 면을 만든다.

- 안쪽 면 삼각형(뒤집어 벽 바깥을 향하게), point index `0 .. numPts-1`
- 트리밍된 오프셋 면 삼각형, point index `numPts ..`
- 각 cap의 annulus

컨투어 삼각형의 방향은 필터의 성질이지 이 벽의 성질이 아니므로, 표본 삼각형의
법선이 최근접 안쪽 점에서 바깥으로 향하는지 측정해 정한다. 표본의 90% 미만이
다수와 어긋나면 오류다(일관되게 감기지 않은 면은 어느 쪽으로 채워도 틀린다).

`FillWallMeshWithTetGen()`이 이 shell을 `plc=1`, `nobisect=1`, `quality=1`,
`minratio=1.414`, `mindihedral=10`으로 채운다. `nobisect=1`이 입력 facet에 절점을
추가하지 못하게 하므로 fluid/wall interface 절점이 그대로 보존된다. 안쪽 면이
닫혀 있어 cap이 없는 입력에서는 shell이 lumen까지 감싸므로
`TGenUtils_FindLumenHolePoint()`로 hole을 지정한다.

TetGen 출력에는 체적 tetrahedron만 있으므로, 하위
`VMTKUtils_CreateBoundaryLayerSurfaceAndCaps()`가 기대하는 모양에 맞춰 shell
삼각형을 `CellEntityIds`/`ModelFaceID`와 함께 다시 넣는다. 삼각형의 세 점이 모두
`numPts` 미만이면 interface(`1`), 모두 이상이면 outer wall(`2`), 섞여 있으면
혈관 끝(`9999`)이다.

### 12.6 이 경로가 보장하는 것과 실패하는 방식

레벨셋은 자기교차하지 않는다. 두 벽이 만나면 교차가 아니라 **융합**한다. 따라서
간격이 두께의 2배 미만인 두 혈관은 오류 없이 하나의 고체로 합쳐지며, 이는
오프셋 진단(18절)에서 "요청보다 훨씬 두꺼운 구역"으로 읽힌다. 물리적으로 원치
않으면 해당 face의 local wall thickness를 간격 절반 아래로 내리거나 모델을
블렌드해야 한다.

TetGen이 shell에서 자기교차를 보고한다면 원인은 접합부가 아니다. 오프셋이
만들지 않은 유일한 부분인 annulus이거나, 이미 자기교차한 안쪽 면이다.

## 13. VMTK 고체벽 압출

> 이 절과 14절은 **wedge 압출 경로 전용**이다. `WallMeshTetGenShell`이 켜져 있으면
> 12절이 이를 대체한다.

최종 `WallThickness` 배열을 surface에 추가하고
`VMTKUtils_BoundaryLayerMesh()`를 호출한다.

| 인자 | 고체벽 값 | 의미 |
| --- | ---: | --- |
| `edgeSize` | global `wallthickness` | wrapper 인자이며 variable thickness 모드에서는 직접 두께로 사용되지 않음 |
| `blThicknessFactor` | `1.0` | `WallThickness`에 곱하는 비율 |
| `numSublayers` | `numwallsublayers` | 벽 layer 수 |
| `sublayerRatio` | `1.0` | 벽 layer를 균등 분할 |
| `sidewallCellEntityId` | `9999` | 열린 끝단 sidewall 표식 |
| `innerSurfaceCellEntityId` | `1` | 시작 interface 표식 |
| `negateWarpVectors` | `0` | outward normal 유지 |
| `useConstantThickness` | `0` | point별 배열 사용 |
| `layerThicknessArrayName` | `WallThickness` | point별 전체 두께 |

wrapper는 VMTK generator에 다음 값도 고정한다.

- `WarpVectorsArrayName = "Normals"`
- `UseWarpVectorMagnitudeAsThickness = 0`
- `SurfaceCellIdsArrayName = "ModelFaceID"`
- surface cell과 sidewall cell 모두 포함
- untangle을 위한 substep 수 `100`

wrapper가 설정하지 않는 값 중 하나는 주의해서 읽어야 한다. VMTK generator에는
point별 두께를 잘라내는 `MaximumLayerThickness`가 있지만 생성자 기본값이
`VTK_VMTK_LARGE_DOUBLE`이고 wrapper가 이를 덮어쓰지 않으므로 실질적으로 걸리지
않는다. 바깥벽 라운딩(11.1절)이 두께 배열에 배정 두께의 최대 3배까지 써 넣기
때문에 이 상한이 라운딩 결과를 자를 수 있는지 의심할 수 있으나, 현재 구성에서는
자르지 않는다.

### 13.1 warp vector 구성

`vtkvmtkBoundaryLayerGenerator::BuildWarpVectors()`는 각 point normal을 정규화하고
point별 두께를 곱한다.

```text
w_i = normalize(n_i) × WallThickness(i) × LayerThicknessRatio
```

고체벽의 `LayerThicknessRatio`는 `1.0`이므로 `|w_i|`는 최종
`WallThickness(i)`다.

### 13.2 VMTK 내부 tangle 검사와 방향 보정

VMTK는 warp를 한 번에 적용하지 않고 총 100 substep으로 점진적으로 적용한다.
현재 구성은 initial 1회, intermediate 10회, final 89회로 나뉜다.

`CheckTangle()`은 base와 warped triangle에 대해 다음 중 하나면 tangle로 본다.

```text
dot(baseNormal, warpedNormal) < 0
warpedArea / baseArea <= 0.1
```

tangle이 있으면 `LocalUntangle()`이 주변 triangle tangent 성분을 사용해 warp
방향을 보정한 다음 다시 검사한다. 보정 후에도 warp vector의 크기는 원래
`WallThickness`로 유지하고 방향만 바꾼다.

따라서 코어의 fold 검사는 코어에서 스무딩된 normal 방향을 검사하지만, VMTK의
untangle이 최종 방향을 추가로 기울일 수 있다. 이 방향 변경이 크면 요청한 normal
방향보다 바깥으로 덜 나가 접합부가 움푹 들어가 보일 수 있다.

### 13.3 layer 위치

전체 layer 수를 `L`, sublayer ratio를 `r`, 0부터 시작하는 layer index를 `k`라
하면 layer 가중치는 다음과 같다.

```text
weight(k) = r^(L-k-1) / Σ r^(L-j-1), j=0..L-1
offset(k) = Σ weight(j), j=0..k
p(i,k) = p_i + offset(k) × w_i
```

고체벽은 `r=1`이므로 모든 layer가 같은 두께다.

```text
p(i,k) = p_i + ((k+1)/L) × w_i
```

따라서 `NumberOfWallLayers`는 총 벽 두께를 늘리지 않는다. 동일한 총 두께를 몇
개의 요소 layer로 나눌지만 결정한다.

## 14. wedge 생성과 tetrahedron 변환

VMTK는 시작 surface cell과 각 layer의 대응 cell을 연결한다.

- 입력 triangle 하나와 다음 layer triangle을 연결해 `VTK_WEDGE` 하나를 만든다.
- 입력 quad라면 `VTK_HEXAHEDRON`을 만든다.
- generator에는 `VTK_QUADRATIC_TRIANGLE`을 `VTK_QUADRATIC_WEDGE`로 압출하는
  경로도 있으나, SimVascular의 표면 재메시와 TetGen은 선형 삼각형만 만들므로 이
  경로는 사용되지 않는다.
- 입력 surface의 열린 boundary edge마다 layer 사이를 잇는 `VTK_QUAD`
  sidewall을 만든다.
- 설정에 따라 원래 inner surface와 최종 outer surface cell도 결과에 포함한다.
- `ModelFaceID`를 surface 계열 cell에 전달한다.

즉 triangle surface와 벽 layer 수 `L`에 대해 기본 체적 셀은 triangle당 wedge
`L`개다.

최종 결합 때 `vtkvmtkUnstructuredGridTetraFilter`가
`vtkOrderedTriangulator`의 topology template을 사용해 wedge/hex를 tetrahedron으로
분할한다. surface quad도 triangle으로 단순화된다.

그 뒤 `VMTKUtils_ReorderTetElements()`가 각 tetrahedron의 signed volume을
검사한다. 음수인 요소는 첫 번째와 두 번째 node를 교환해 orientation을 바꾼다.
이는 node 순서만 바로잡는 과정이다. 기하적으로 뒤틀리거나 자기 교차한 요소를
복구하는 알고리즘은 아니다.

## 15. TetGen 유체 코어 생성

유체 boundary layer의 가장 안쪽 표면은 열린 입구/출구를 갖는다.
`GenerateAndMeshCaps()`가 VMTK capper로 이를 닫고 원본 `ModelFaceID`를 복원한
뒤 cap을 재메시한다. 이 닫힌 PLC가 `polydatasolid_`에 저장된다.

TetGen은 이 표면 내부를 tetrahedron으로 채운다. boundary layer를 사용하는
경로에서는 `quality=3`, `metric=1`, `mindihedral=10`, `nobisect=1`을 강제로
설정한다. `nobisect=1`은 interface facet을 다시 나누지 않아 유체 코어와
boundary layer의 node 대응을 유지하는 데 필요하다. sizing function과 refinement
설정도 TetGen 입력에 반영된다.

TetGen은 고체벽 체적을 채우지 않는다. 고체벽은 VMTK의 surface extrusion과
wedge-to-tet 변환으로 만들어진다.

## 16. 유체와 고체벽 병합

`GenerateMesh()`의 TetGen 단계가 끝나면 `AppendBoundaryLayerMesh()`가
`VMTKUtils_AppendData()`를 호출한다.

### 16.1 유체 코어와 유체 boundary layer

1. 유체 boundary-layer wedge/hex를 tetrahedron으로 변환한다.
2. 음수 signed-volume tetrahedron의 node 순서를 보정한다.
3. TetGen 코어가 가진 `ModelRegionID` 범위를 읽는다. 현재 append 경로는 단일
   TetGen region을 전제로 한다.
4. boundary-layer volume에 같은 유체 `ModelRegionID`를 준다.
5. `vtkvmtkAppendFilter`로 TetGen 코어와 유체 boundary layer를 합친다.

`vtkvmtkAppendFilter`는 기본적으로 `MergeDuplicatePoints=1`이고
`vtkMergePoints`를 사용한다. 따라서 같은 좌표의 interface point가 하나의
mesh point로 병합된다.

고체벽이 있으면 `NewRegionBoundaryLayer` 옵션과 무관하게 유체 boundary layer는
TetGen 코어와 같은 유체 region에 남는다. 새 region ID는 고체벽에 사용된다.

### 16.2 고체벽 region

고체벽 `wallmesh_`도 tetrahedron으로 변환하고 orientation을 보정한다.
유체 region ID가 `modelId`라면 고체벽의 모든 volume cell에 다음 값을 준다.

```text
wall ModelRegionID = modelId + 1
```

`VMTKUtils_CreateBoundaryLayerSurfaceAndCaps()`는 cell type을 보고
`isSurface` 배열을 만든다.

- triangle/quad: `isSurface=1`
- tetrahedron 등 volume cell: `isSurface=0`

이 배열로 wall surface와 wall volume을 분리한다. `WallID=0`인 원본 cap
surface를 따로 얻고, `CellEntityIds=9999`인 wall 끝단 sidewall의
`ModelFaceID`를 대응 cap에서 복원한다.

### 16.3 fluid/wall interface의 node 병합

`VMTKUtils_CreateNewBoundaryLayerRegion()`이 합쳐진 유체 volume과 wall volume을
받는다.

1. 유체 point 좌표를 `vtkCoincidentPoints`에 먼저 등록한다.
2. 유체 point에 1부터 순서대로 `GlobalNodeID`를 준다.
3. wall point 좌표를 조회한다.
4. 같은 좌표의 유체 point가 있으면 그 유체 `GlobalNodeID`를 재사용한다.
5. 일치하지 않는 wall point에는 새 ID를 준다.
6. 유체 요소 뒤에 wall 요소가 이어지도록 `GlobalElementID`를 부여한다.
7. `vtkAppendFilter(SetMergePoints(true))`로 두 volume을 합친다.

원래 interface 좌표를 양쪽에서 보존했기 때문에 interface의 wall 시작점은 유체
외부 경계점과 일치한다. 병합 후 유체와 고체는 서로 다른 `ModelRegionID`를
가지지만 같은 mesh point와 `GlobalNodeID`를 공유하는 conformal interface가 된다.

좌표 해시를 만들 때 wall point도 미리 등록하므로 새 ID 시퀀스에 빈 번호가 생길
가능성은 있다. 중요한 보장은 ID가 조밀하다는 것이 아니라 같은 interface 좌표가
같은 ID를 사용한다는 것이다.

### 16.4 최종 surface 병합

volume 병합과 surface 병합은 별도로 수행된다.

1. wall volume의 surface를 추출하고 wall surface source로 `ModelFaceID`를
   복원한다.
2. fluid volume의 surface를 추출하고 fluid boundary-layer shell과 내부 cap을
   합친 source로 `ModelFaceID`를 복원한다.
3. 두 surface를 `vtkAppendPolyData`로 이어 최종 `surfacemesh_`를 만든다.

최종 surface에는 외부 유체 경계, cap, outer wall, 벽 끝단과 함께 재료 경계인
fluid/wall interface가 포함될 수 있다. volume은 interface point를 병합하지만,
surface 표현은 두 region의 경계를 표현하기 위해 양쪽 source에서 온 coincident
face를 가질 수 있다.

outer wall은 원래 모델보다 바깥에 있으므로 전체 surface에 대해 원본 face로 단순
최근접 매핑하지 않는다. fluid와 wall source를 분리해 각자의 `ModelFaceID`를
복원한다.

## 17. “접합부 병합”의 두 의미

혈관 분지에서 말하는 병합은 다음 두 개를 구분해야 한다.

### 17.1 혈관 branch의 기하 병합

mesher 안에서 독립적인 원통 branch들을 boolean union하는 단계는 없다. 입력
모델링 단계에서 이미 branch가 하나의 연결된 surface로 병합되어 있어야 한다.
mesher는 그 접합부를 다음 방식으로 연속적으로 처리한다.

- 공유 point topology 유지
- 공유점의 local 두께를 고유 face별 산술 평균
- 오목 영역 normal(warp vector) 스무딩
- Laplacian 두께 스무딩과 gradation 제한

여기까지는 두 경로가 같고, 접합부의 바깥면을 만드는 방식이 갈린다.

| | wedge 압출 경로 | TetGen shell 경로 |
| --- | --- | --- |
| 두께 | 곡률 기반 제한으로 축소 | 요청값 유지 |
| 바깥면 | 라운딩(볼록 fillet)으로 함몰을 메움 | 오프셋의 crease가 그 자리 |
| 잔여 fold | outer triangle 검사 후 국소 두께 축소 | 발생하지 않음 |

### 17.2 생성된 유체/고체 mesh의 병합

이는 boolean 연산이 아니라 같은 좌표의 point를 합치고 ID를 일관되게 만드는
데이터 병합이다.

- 유체 코어와 유체 BL: duplicate point merge, 같은 fluid region
- 유체와 고체벽: duplicate point merge, interface의 `GlobalNodeID` 공유
- 고체벽: 별도 `ModelRegionID`
- 최종 surface: fluid/wall source별 face ID 복원 후 append

## 18. 진단 출력의 의미

현재 구현은 다음 정보를 로그로 출력한다.

- 실제 적용된 벽 메시 옵션값(`Wall mesh options in effect: ...`)
- interface surface의 triangle 품질
- 오목 영역에서 스무딩한 warp vector point 수와 최대 normal 방향 변화
  (`Smoothed the wall extrusion warp vectors at ...`)
- 요청 두께가 국소 오목 곡률 반경을 넘는 point 수(`t/R >= 1 / 2 / 4`)와 해당
  구역별 `t`, `R_typical`, `R_smallest`, 퇴화 삼각형 인접 여부
  (`Concave curvature vs requested wall thickness (t/R, ...)`)
- 바깥벽 라운딩으로 끌어올린 오목 point 수와 최대 fillet 배율, 위치
  (`Wall outer rounding: ...`, wedge 경로 전용)
- 곡률 제한으로 매우 얇아진 point 수(wedge 경로 전용)
- fold 방지 반복 후에도 남은 triangle과 위치/형상/normal 정보(wedge 경로 전용)
- shell 경로에서는 대신 다음이 출력된다.
  - 오프셋 격자 크기/간격, 예산으로 거칠어졌다면 그 사실과 실제 해상도
  - 밴드에서 평가한 voxel 수, 부호로 채운 안팎 voxel 수
  - 오프셋 면의 점/삼각형 수와 연결된 shell 개수
  - 리메시 전후 삼각형 품질(`offset outer wall, before/after remesh`)
  - cap별 안쪽 rim과 트리밍된 rim의 점 수
  - **양방향 오프셋 두께**(`Offset wall thickness [...]`). 바깥 방향(오프셋 면 →
    안쪽 면)은 구성 오차를 재고 — 레벨셋은 모든 오프셋 점을 요청 거리에 놓으므로
    미달은 격자 간격이나 리메시 탓이다 — 안쪽 방향(안쪽 면 → 오프셋 면)은
    인터페이스 위에 실제로 서 있는 벽을 잰다. 접합부에서 1을 넘는 것은 crease가
    채운 것이므로 의도된 결과이고, 1 미만이 벽이 없는 것이다. 안쪽 방향 비율은
    `OffsetThicknessRatio` 필드로 `wall_offset_diagnostics.vtp`에 기록된다
  - 면적 0인 혈관 끝 삼각형 수
  - shell의 점/삼각형 수와 닫은 혈관 끝 수, TetGen이 채운 tetrahedron/절점 수
  - interface/outer/side wall 삼각형 수
- 최종 두께가 반경 보정 후 기준 두께의 90%/50%/25% 미만인 point 수
  (`Wall thickness reduction (final vs requested): ...`)
- 얇아진 point를 **공간적으로 분리된 구역(thinned regions) 단위**로, 최악 구역부터
  각 구역의 대표 비율/좌표를 보고(단일 최악 point가 아니라 "한 곳만 나쁜지 vs
  모든 접합부가 얇은지"를 구분하기 위함)

구역 목록은 두 곳에서 나온다. 두께 감소 진단의 `thinned regions`와 t/R 진단의
`concave regions`이며, 둘 다 `TGenUtils_ClusterPointsIntoRegions()`로 같은 그리디
방식과 같은 반경(모델 대각선의 2%)을 쓴다. 다만 **정렬 키가 서로 달라(두께는 비율
오름차순, t/R은 비율 내림차순) 두 목록의 구역 번호는 대응하지 않는다.** 같은
접합부를 대조할 때는 번호가 아니라 출력된 좌표를 기준으로 맞춰야 한다.
- VMTK untangle 전후 warp tilt가 15/30/45도를 넘는 point 수와 최악 위치
- VMTK `CheckTangle()`이 찾은 triangle 수
- 최종 tetrahedron aspect ratio와 최악 요소 위치

두께 감소 진단의 코드 변수 `baseThickness`는 local/global 요청값에 선택적 반경
보정까지 적용한 값이다. 로그 문구의 “requested”는 순수 `t_req`가 아니라 이
후처리 전 기준값을 뜻한다.

warp tilt는 VMTK wrapper가 유체 boundary layer와 고체벽 양쪽에서 사용되므로
호출 문맥과 압출 방향을 함께 보고 해석해야 한다.

## 19. 보장 범위와 한계

현재 알고리즘은 다음을 보장하거나 의도한다.

- fluid/wall interface의 시작 좌표를 두께 처리 중 이동하지 않는다.
- local 두께 공유점 평균이 triangle 밀도에 의존하지 않는다.
- normal/두께 스무딩이 point 순회 순서에 의존하지 않는다.
- 고체벽 layer 전체 두께가 point별 `WallThickness`를 따른다(wedge 경로).
- 유체와 벽이 interface node와 `GlobalNodeID`를 공유한다.
- 유체와 벽을 서로 다른 `ModelRegionID`로 구분한다.

shell 경로(12절)는 여기에 더해 다음을 보장한다.

- 바깥 면의 모든 점이 안쪽 면에서 최소 요청 두께만큼 떨어져 있다(격자 해상도와
  리메시 오차 범위 내. 그 오차 자체를 18절의 양방향 진단이 측정한다).
- 바깥 면은 자기교차하지 않는다. 레벨셋이기 때문이다.
- `nobisect=1`이므로 TetGen이 입력 facet에 절점을 추가하지 않는다.

다음은 보장하지 않는다.

- 떨어진 surface 구간 사이의 전역 outer-wall 자기 교차(wedge 경로).
  shell 경로에서는 대신 두 벽이 **융합**하며, 이는 오류로 검출되지 않는다(12.6절)
- 입력 surface의 0-area triangle 또는 0-length normal 복구
- VMTK untangle 뒤 모든 tetrahedron의 양의 Jacobian과 충분한 품질
- node 순서 보정만으로 기하학적으로 접힌 요소 복구
- 여러 TetGen `ModelRegionID`를 가진 입력의 일반적인 append
- 최종 surface에 coincident interface face가 전혀 없다는 보장
- shell 경로에서 벽 요소가 두께 방향으로 층을 이룬다는 보장. 채움은 비구조
  tetrahedron이며 `NumberOfWallLayers`는 이 경로에서 의미가 없다

## 20. 실행 후 검증 체크리스트

Linux Docker 실행 환경에서 최소한 다음을 확인한다.

1. `ModelRegionID`가 유체와 고체벽 두 영역으로 올바르게 나뉘는가.
2. fluid/wall interface에서 같은 좌표와 `GlobalNodeID`를 공유하는가.
3. interface에 crack 또는 중복 volume point가 없는가.
4. 모든 tetrahedron의 signed volume/Jacobian이 양수인가.
5. aspect ratio가 나쁜 상위 요소가 분지 crotch와 sliver 주변에 집중되는가.
6. outer wall에 국소 fold 또는 떨어진 branch 사이 전역 자기 교차가 없는가.
7. local 두께 face 경계에서 outer wall이 계단처럼 급변하지 않는가.
8. VMTK warp tilt가 큰 위치와 시각적으로 움푹 들어간 위치가 일치하는가.
9. `t/R >= 1`인 point 수가 접합부 규모인가(형상이 요청 두께를 감당 못 함), 아니면
    거의 없는데도 `thinned regions`가 여러 개인가(두께 패스가 원인). 구역별
    `R_typical`과 `R_smallest`가 함께 작은지(날카로운 형상), `R_smallest`만
    작은지(퇴화 삼각형 아티팩트)도 함께 본다.
10. wall layer 수를 바꿔도 총 두께가 유지되고 두께 방향 해상도만 변하는가
    (wedge 경로 한정. shell 경로에서는 이 옵션이 아무것도 하지 않는다).
11. cap 및 sidewall의 `ModelFaceID`가 입력 face와 일치하는가.
12. mesh-complete 출력 후에도 region, face, node, element ID가 보존되는가.

shell 경로(12절)에서는 추가로 다음을 확인한다.

13. 오프셋 격자 간격이 예산 때문에 거칠어졌는가. 거칠어졌다면 그 해상도가 벽
    두께 대비 충분한가.
14. 양방향 오프셋 두께 진단에서 **바깥 방향**이 거의 1인가. 1보다 작다면 격자나
    리메시가 두께를 갉아먹은 것이므로, 이 접근 자체의 이득이 사라진다.
15. **안쪽 방향**에서 90% 미만 point 수가 0에 가까운가. 접합부에서 1을 넘는 것은
    의도된 결과(crease가 채움)이므로 문제가 아니다.
16. 요청보다 훨씬 두꺼운 구역이 있다면 그것이 두 혈관의 융합인가(12.6절).
    융합이라면 해당 face의 local wall thickness를 간격 절반 아래로 내린다.
17. 면적 0인 혈관 끝 삼각형 수가 0인가. 0이 아니라면 cap rim이 퇴화에 가깝다.
18. 오프셋 면의 연결된 shell 개수가 예상과 맞는가(정상은 1개, 모델이 분리되어
    있으면 그 개수).

Mac 개발 환경에서는 실제 SimVascular/VMTK/TetGen 빌드와 실행 성공을 단정하지
않는다. 실행 검증 결과는 저장소의 `logs/`에 남긴다.

## 21. 코드 위치

라인 번호는 변경에 따라 이동할 수 있으므로 함수명을 기준으로 찾는다.

| 파일 | 핵심 함수/구조 |
| --- | --- |
| `Code/Source/sv/Mesh/TetGenMeshObject/sv_TetGenMeshObject.h` | `TGoptions`(`walltetgenshell` 포함), `GenerateBoundaryLayerMesh`, `GenerateWallMesh`, `FillWallMeshWithTetGen`, `AppendBoundaryLayerMesh` |
| `Code/Source/sv/Mesh/TetGenMeshObject/sv_TetGenMeshObject.cxx` | constructor 기본값, `SetMeshOptions`, `SetBoundaryLayer`, `SetWalls`, `GenerateMesh`, `GenerateSurfaceRemesh`, `SetCapBoundaryNormals`, `GenerateBoundaryLayerMesh`, `GenerateWallMesh`(`extrudeWedges` 게이팅), `FillWallMeshWithTetGen`, `GenerateAndMeshCaps`, `AppendBoundaryLayerMesh` |
| `Code/Source/sv/Mesh/TetGenMeshObject/sv_tetgenmesh_utils.cxx` (wedge 경로) | `TGenUtils_SmoothWarpVectorsInConcaveRegions`, `TGenUtils_SmoothPointArray`, `TGenUtils_ClampThicknessToConcaveCurvature`, `TGenUtils_ReportConcaveCurvatureVsThickness`, `TGenUtils_LimitThicknessGradation`, `TGenUtils_ClusterPointsIntoRegions`, `TGenUtils_RoundOuterWallToPreserveThickness`, `TGenUtils_LimitThicknessToPreventFold`, `TGenUtils_ReportAchievedWallThickness`, 품질 보고 함수 |
| `Code/Source/sv/Mesh/TetGenMeshObject/sv_tetgenmesh_utils.cxx` (shell 경로) | `TGenUtils_ExtractBoundaryLoops`, `TGenUtils_BuildOffsetOuterSurface`, `TGenUtils_TrimOffsetSurfaceAtCaps`, `TGenUtils_StitchCapAnnulus`, `TGenUtils_BuildWallShellSurface`, `TGenUtils_FindLumenHolePoint`, `TGenUtils_ReportOffsetWallThickness`, 구조체 `TGenUtilsCapRim` |
| `Code/Source/sv/Mesh/MMGMeshUtils/sv_mmg_mesh_utils.cxx` | `MMGUtils_SurfaceRemeshing` (오프셋 면 리메시) |
| `Code/Source/sv/Mesh/VMTKUtils/sv_vmtk_utils.cxx` | `VMTKUtils_BoundaryLayerMesh`, `VMTKUtils_AppendData`, `VMTKUtils_CreateNewBoundaryLayerRegion`, `VMTKUtils_CreateBoundaryLayerSurfaceAndCaps`, `VMTKUtils_ReorderTetElements` |
| `Code/ThirdParty/vmtk/simvascular_vmtk/vtkvmtkBoundaryLayerGenerator.cxx` | `RequestData`, `BuildWarpVectors`, `CheckTangle`, `LocalUntangle`, `WarpPoints` |
| `Code/ThirdParty/vmtk/simvascular_vmtk/vtkvmtkUnstructuredGridTetraFilter.cxx` | `Execute` |
| `Code/ThirdParty/vmtk/simvascular_vmtk/vtkvmtkAppendFilter.cxx` | duplicate point 병합 |
| `Code/Source/sv4gui/Plugins/org.sv.gui.qt.meshing/sv4gui_MeshEdit.cxx` | GUI 명령 생성 |
| `Code/Source/sv4gui/Modules/Mesh/Common/sv4gui_MeshTetGen.cxx` | GUI/.msh 명령 파싱 |
| `Code/Source/PythonAPI/MeshingTetGen_PyClass.cxx` | Python `generate_mesh()` 진입 |
| `Code/Source/PythonAPI/MeshingTetGenOptions_PyClass.cxx` | Python 옵션 정의와 코어 옵션 변환 |
