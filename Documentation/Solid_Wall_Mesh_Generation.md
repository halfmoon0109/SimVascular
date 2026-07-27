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
   ├─ VMTK 고체벽 tetrahedron
   │
   └─ 압출된 outer wall 표면
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
option WallThicknessRadiusFactor <value>
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
| `wall_thickness_radius_factor` | `WallThicknessRadiusFactor` |
| `local_wall_thickness` | `LocalWallThickness` |

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
        │             ├─ 곡률 기반 두께 제한
        │             ├─ 두께 Laplacian 스무딩 + 곡률 재제한
        │             ├─ 바깥벽 라운딩(볼록 fillet)
        │             ├─ fold-over 방지
        │             └─ 고체벽 바깥쪽 압출
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
| `WallThicknessRadiusFactor` | `wallthicknessradiusfactor` | `0.0` | 국소 반경 대비 두께 비율 |
| `LocalWallThickness` | `localWallThickness_` | 비어 있음 | `face ID → 두께` override |
| `BoundaryLayerDirection` | `boundarylayerdirection` | `1` | 유체 boundary layer를 안쪽으로 압출 |

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
4. `MeshSizingFunction` point 배열을 만든다.
5. `SetCapBoundaryNormals()`로 열린 cap rim의 normal을 cap 평면 방향에 맞춘다.

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
global/local 요청 두께
  → 선택적 반경 적응 두께
  → 곡률 제한
  → Laplacian 두께 스무딩
  → 곡률 재제한
  → 바깥벽 라운딩(두께 유지, 값은 달성 outer 거리로 갱신)
  → fold 방지 축소
  = 최종 WallThickness
```

이와 별개로 오목 영역 warp vector(normal) 스무딩(8절)이 두께 계산 직후·곡률 제한
직전에 실행된다. 이 단계는 두께 값이 아니라 압출 방향(normal)만 바꾼다. 바깥벽
라운딩(11.1절)은 두께를 유지하면서 바깥면을 밀어내고 그 달성 거리를 `WallThickness`
값과 normal에 다시 인코딩하므로, 이후 fold 방지가 처리할 잔여 fold만 남긴다.

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

### 7.2 반경 적응 두께

`wallthicknessradiusfactor = α > 0`이면 centerline까지의 거리로 국소 반경
`R_i`를 근사한다.

```text
R_i = DistanceToCenterlines(i)
t_radius(i) = α × R_i
t_base(i) = clamp(t_radius(i), 0.05 × t_req(i), t_req(i))
```

`sys_geom_distancetocenterlines()`가 `DistanceToCenterlines`를 계산한다.
centerline은 `GenerateWallMesh()`가 만들지 않으며, radius meshing 경로가
계산한 것을 `SetCenterlines()`로 미리 전달해야 한다. radius factor가 켜졌는데
centerline이 없거나 point 수가 맞지 않으면 벽 생성을 중단한다.

상한을 global 값이 아닌 해당 point의 `t_req(i)`로 두므로 local override를
존중한다. 5% 하한은 잘못된 반경 추정이나 작은 혈관에서 0에 가까운 두께가 생기는
것을 막는다. radius factor가 0이면 `t_base(i) = t_req(i)`다.

분지 crotch에서는 가장 가까운 centerline이 옆 branch의 centerline일 수 있으므로
`R_i`는 근사값이다. 이 단계는 혈관 구경에 따른 두께 분배를 위한 것이며 접합부
fold를 직접 해결하는 단계는 아니다.

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
곡률 제한(9절)이, 접합부 두께 **보존**은 바깥벽 라운딩(11.1절)이 담당한다. 각 단계의
역할이 다르므로 순서대로 적용된다.

## 9. 오목 곡률 기반 두께 제한

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
최소 altitude `a_min(i)`로 계산한다.

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

## 12. VMTK 고체벽 압출

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

### 12.1 warp vector 구성

`vtkvmtkBoundaryLayerGenerator::BuildWarpVectors()`는 각 point normal을 정규화하고
point별 두께를 곱한다.

```text
w_i = normalize(n_i) × WallThickness(i) × LayerThicknessRatio
```

고체벽의 `LayerThicknessRatio`는 `1.0`이므로 `|w_i|`는 최종
`WallThickness(i)`다.

### 12.2 VMTK 내부 tangle 검사와 방향 보정

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

### 12.3 layer 위치

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

## 13. wedge 생성과 tetrahedron 변환

VMTK는 시작 surface cell과 각 layer의 대응 cell을 연결한다.

- 입력 triangle 하나와 다음 layer triangle을 연결해 `VTK_WEDGE` 하나를 만든다.
- 입력 quad라면 `VTK_HEXAHEDRON`을 만든다.
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

## 14. TetGen 유체 코어 생성

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

## 15. 유체와 고체벽 병합

`GenerateMesh()`의 TetGen 단계가 끝나면 `AppendBoundaryLayerMesh()`가
`VMTKUtils_AppendData()`를 호출한다.

### 15.1 유체 코어와 유체 boundary layer

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

### 15.2 고체벽 region

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

### 15.3 fluid/wall interface의 node 병합

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

### 15.4 최종 surface 병합

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

## 16. “접합부 병합”의 두 의미

혈관 분지에서 말하는 병합은 다음 두 개를 구분해야 한다.

### 16.1 혈관 branch의 기하 병합

mesher 안에서 독립적인 원통 branch들을 boolean union하는 단계는 없다. 입력
모델링 단계에서 이미 branch가 하나의 연결된 surface로 병합되어 있어야 한다.
mesher는 그 접합부를 다음 방식으로 연속적으로 처리한다.

- 공유 point topology 유지
- 공유점의 local 두께를 고유 face별 산술 평균
- centerline 반경에 따른 두께 변화
- 오목 영역 normal(warp vector) 스무딩
- 곡률 기반 두께 제한
- Laplacian 두께 스무딩
- 바깥벽 라운딩(두께 유지 볼록 fillet)
- outer triangle fold 검사와 국소 두께 축소

### 16.2 생성된 유체/고체 mesh의 병합

이는 boolean 연산이 아니라 같은 좌표의 point를 합치고 ID를 일관되게 만드는
데이터 병합이다.

- 유체 코어와 유체 BL: duplicate point merge, 같은 fluid region
- 유체와 고체벽: duplicate point merge, interface의 `GlobalNodeID` 공유
- 고체벽: 별도 `ModelRegionID`
- 최종 surface: fluid/wall source별 face ID 복원 후 append

## 17. 진단 출력의 의미

현재 구현은 다음 정보를 로그로 출력한다.

- 실제 적용된 벽 메시 옵션값(`Wall mesh options in effect: ...`)
- interface surface의 triangle 품질
- 오목 영역에서 스무딩한 warp vector point 수와 최대 normal 방향 변화
  (`Smoothed the wall extrusion warp vectors at ...`)
- 바깥벽 라운딩으로 끌어올린 오목 point 수와 최대 fillet 배율, 위치
  (`Wall outer rounding: ...`)
- 곡률 제한으로 매우 얇아진 point 수
- fold 방지 반복 후에도 남은 triangle과 위치/형상/normal 정보
- 최종 두께가 반경 보정 후 기준 두께의 90%/50%/25% 미만인 point 수
  (`Wall thickness reduction (final vs requested): ...`)
- 얇아진 point를 **공간적으로 분리된 구역(thinned regions) 단위**로, 최악 구역부터
  각 구역의 대표 비율/좌표를 보고(단일 최악 point가 아니라 "한 곳만 나쁜지 vs
  모든 접합부가 얇은지"를 구분하기 위함)
- VMTK untangle 전후 warp tilt가 15/30/45도를 넘는 point 수와 최악 위치
- VMTK `CheckTangle()`이 찾은 triangle 수
- 최종 tetrahedron aspect ratio와 최악 요소 위치

두께 감소 진단의 코드 변수 `baseThickness`는 local/global 요청값에 선택적 반경
보정까지 적용한 값이다. 로그 문구의 “requested”는 순수 `t_req`가 아니라 이
후처리 전 기준값을 뜻한다.

warp tilt는 VMTK wrapper가 유체 boundary layer와 고체벽 양쪽에서 사용되므로
호출 문맥과 압출 방향을 함께 보고 해석해야 한다.

## 18. 보장 범위와 한계

현재 알고리즘은 다음을 보장하거나 의도한다.

- fluid/wall interface의 시작 좌표를 두께 처리 중 이동하지 않는다.
- local 두께 공유점 평균이 triangle 밀도에 의존하지 않는다.
- normal/두께 스무딩이 point 순회 순서에 의존하지 않는다.
- 고체벽 layer 전체 두께가 point별 `WallThickness`를 따른다.
- 유체와 벽이 interface node와 `GlobalNodeID`를 공유한다.
- 유체와 벽을 서로 다른 `ModelRegionID`로 구분한다.

다음은 보장하지 않는다.

- 떨어진 surface 구간 사이의 전역 outer-wall 자기 교차
- 입력 surface의 0-area triangle 또는 0-length normal 복구
- VMTK untangle 뒤 모든 tetrahedron의 양의 Jacobian과 충분한 품질
- node 순서 보정만으로 기하학적으로 접힌 요소 복구
- 여러 TetGen `ModelRegionID`를 가진 입력의 일반적인 append
- 최종 surface에 coincident interface face가 전혀 없다는 보장

## 19. 실행 후 검증 체크리스트

Linux Docker 실행 환경에서 최소한 다음을 확인한다.

1. `ModelRegionID`가 유체와 고체벽 두 영역으로 올바르게 나뉘는가.
2. fluid/wall interface에서 같은 좌표와 `GlobalNodeID`를 공유하는가.
3. interface에 crack 또는 중복 volume point가 없는가.
4. 모든 tetrahedron의 signed volume/Jacobian이 양수인가.
5. aspect ratio가 나쁜 상위 요소가 분지 crotch와 sliver 주변에 집중되는가.
6. outer wall에 국소 fold 또는 떨어진 branch 사이 전역 자기 교차가 없는가.
7. local 두께 face 경계에서 outer wall이 계단처럼 급변하지 않는가.
8. 반경 적응 사용 시 `DistanceToCenterlines`가 branch별 실제 반경을 합리적으로
   나타내는가.
9. VMTK warp tilt가 큰 위치와 시각적으로 움푹 들어간 위치가 일치하는가.
10. wall layer 수를 바꿔도 총 두께가 유지되고 두께 방향 해상도만 변하는가.
11. cap 및 sidewall의 `ModelFaceID`가 입력 face와 일치하는가.
12. mesh-complete 출력 후에도 region, face, node, element ID가 보존되는가.

Mac 개발 환경에서는 실제 SimVascular/VMTK/TetGen 빌드와 실행 성공을 단정하지
않는다. 실행 검증 결과는 저장소의 `logs/`에 남긴다.

## 20. 코드 위치

라인 번호는 변경에 따라 이동할 수 있으므로 함수명을 기준으로 찾는다.

| 파일 | 핵심 함수/구조 |
| --- | --- |
| `Code/Source/sv/Mesh/TetGenMeshObject/sv_TetGenMeshObject.h` | `TGoptions`, `GenerateBoundaryLayerMesh`, `GenerateWallMesh`, `AppendBoundaryLayerMesh` |
| `Code/Source/sv/Mesh/TetGenMeshObject/sv_TetGenMeshObject.cxx` | constructor 기본값, `SetMeshOptions`, `SetBoundaryLayer`, `SetWalls`, `GenerateMesh`, `GenerateSurfaceRemesh`, `SetCapBoundaryNormals`, `GenerateBoundaryLayerMesh`, `GenerateWallMesh`, `GenerateAndMeshCaps`, `AppendBoundaryLayerMesh` |
| `Code/Source/sv/Mesh/TetGenMeshObject/sv_tetgenmesh_utils.cxx` | `TGenUtils_SmoothWarpVectorsInConcaveRegions`, `TGenUtils_SmoothPointArray`, `TGenUtils_ClampThicknessToConcaveCurvature`, `TGenUtils_RoundOuterWallToPreserveThickness`, `TGenUtils_LimitThicknessToPreventFold`, 품질 보고 함수 |
| `Code/Source/sv/Mesh/VMTKUtils/sv_vmtk_utils.cxx` | `VMTKUtils_BoundaryLayerMesh`, `VMTKUtils_AppendData`, `VMTKUtils_CreateNewBoundaryLayerRegion`, `VMTKUtils_CreateBoundaryLayerSurfaceAndCaps`, `VMTKUtils_ReorderTetElements` |
| `Code/ThirdParty/vmtk/simvascular_vmtk/vtkvmtkBoundaryLayerGenerator.cxx` | `RequestData`, `BuildWarpVectors`, `CheckTangle`, `LocalUntangle`, `WarpPoints` |
| `Code/ThirdParty/vmtk/simvascular_vmtk/vtkvmtkUnstructuredGridTetraFilter.cxx` | `Execute` |
| `Code/ThirdParty/vmtk/simvascular_vmtk/vtkvmtkAppendFilter.cxx` | duplicate point 병합 |
| `Code/Source/sv4gui/Plugins/org.sv.gui.qt.meshing/sv4gui_MeshEdit.cxx` | GUI 명령 생성 |
| `Code/Source/sv4gui/Modules/Mesh/Common/sv4gui_MeshTetGen.cxx` | GUI/.msh 명령 파싱 |
| `Code/Source/PythonAPI/MeshingTetGen_PyClass.cxx` | Python `generate_mesh()` 진입 |
| `Code/Source/PythonAPI/MeshingTetGenOptions_PyClass.cxx` | Python 옵션 정의와 코어 옵션 변환 |
