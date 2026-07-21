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

# TetGen 메시 생성 로직

## 1. 문서 목적

이 문서는 SimVascular의 TetGen 메시 생성 과정이 GUI와 Python API에서 시작되어
VMTK/MMG 표면 처리, 유체 경계층 생성, TetGen 체적 메시 생성, 선택적인 고체 혈관벽
메시 생성 및 최종 메시 결합으로 이어지는 흐름을 설명한다.

주요 대상 코드는 다음과 같다.

- `Code/Source/sv/Mesh/TetGenMeshObject/sv_TetGenMeshObject.cxx`
- `Code/Source/sv/Mesh/TetGenMeshObject/sv_tetgenmesh_utils.cxx`
- `Code/Source/sv/Mesh/VMTKUtils/sv_vmtk_utils.cxx`
- `Code/Source/sv4gui/Modules/Mesh/Common/sv4gui_MeshTetGen.cxx`
- `Code/Source/sv4gui/Plugins/org.sv.gui.qt.meshing/sv4gui_MeshEdit.cxx`
- `Code/Source/PythonAPI/MeshingTetGen_PyClass.cxx`
- `Code/Source/PythonAPI/MeshingTetGenOptions_PyClass.cxx`

## 2. 전체 처리 흐름

```text
PolyData 모델 입력
        │
        ├─ GUI 명령 생성 또는 Python TetGenOptions 설정
        ├─ wall face 지정
        ├─ global/local 크기 및 refinement 설정
        └─ boundary layer / solid wall 옵션 설정
        │
        ▼
cvTetGenMeshObject::GenerateMesh()
        │
        ├─ 옵션 조합 검증
        ├─ 표면 재메싱
        │
        ├─ [경계층 사용]
        │      ├─ 유체 경계층을 모델 안쪽으로 압출
        │      ├─ [solid wall 사용] 모델 벽을 바깥쪽으로 압출
        │      └─ 내부 유체 표면의 입구/출구 cap 생성
        │
        ├─ TetGen 입력 및 sizing function 구성
        ├─ TetGen 내부 체적 tetrahedralization
        │
        ├─ [경계층 사용]
        │      └─ 내부 유체 + 유체 경계층 + 선택적 solid wall 결합
        │
        ├─ surface/volume 데이터 배열 정리
        └─ 메시 품질 계산
```

모든 프런트엔드는 최종적으로 `cvTetGenMeshObject::GenerateMesh()`에 합류한다.
GUI와 Python API의 차이는 옵션과 입력 데이터를 코어 mesher에 전달하는 방식에 있다.

## 3. 주요 객체와 데이터

### 3.1 핵심 메시 객체

| 객체 | 형식 | 역할 |
| --- | --- | --- |
| `originalpolydata_` | `vtkPolyData` | 최초 입력 모델 보존. 원래 face ID 복원에 사용한다. |
| `polydatasolid_` | `vtkPolyData` | 표면 재메싱, 경계층 생성, cap 생성 과정에서 변경되는 작업 표면이다. |
| `inmesh_` | `tetgenio` | TetGen에 전달되는 입력 표면, marker, sizing 및 region 정보다. |
| `outmesh_` | `tetgenio` | TetGen이 생성한 내부 tetrahedral mesh다. |
| `boundarylayermesh_` | `vtkUnstructuredGrid` | VMTK가 생성한 유체 경계층 메시다. |
| `innerblmesh_` | `vtkUnstructuredGrid` | 경계층 생성 전의 기준 표면 또는 경계층 내부 표면 처리에 사용된다. |
| `wallmesh_` | `vtkUnstructuredGrid` | 선택적으로 생성되는 고체 혈관벽 메시다. |
| `surfacemesh_` | `vtkPolyData` | 외부 경계와 region interface를 포함한 최종 표면 메시다. |
| `volumemesh_` | `vtkUnstructuredGrid` | 내부 유체, 유체 경계층, 선택적 고체벽을 포함한 최종 체적 메시다. |

### 3.2 주요 VTK 데이터 배열

| 배열 | 위치 | 역할 |
| --- | --- | --- |
| `ModelFaceID` | cell data | 원본 모델의 wall, inlet, outlet 등 경계면을 식별한다. |
| `WallID` | cell data | wall cell과 cap cell을 구분한다. |
| `CellEntityIds` | cell data | VMTK 경계층 처리 중 surface와 sidewall entity를 구분한다. |
| `MeshSizingFunction` | point data | TetGen 또는 VMTK에 전달할 위치별 메시 크기를 저장한다. |
| `WallThickness` | point data | 고체 혈관벽을 압출할 위치별 두께를 저장한다. |
| `ModelRegionID` | cell data | 유체와 고체벽 같은 재료 영역을 구분한다. |
| `GlobalNodeID` | point data | 결합된 최종 메시의 전역 노드 번호다. |
| `GlobalElementID` | cell data | 결합된 최종 메시의 전역 요소 번호다. |
| `AspectRatio` | cell data | 생성된 tetrahedral element의 품질 지표다. |

## 4. 프런트엔드 진입 경로

### 4.1 GUI 경로

`sv4guiMeshEdit::RunCommands()`는 다음 순서로 메시 생성을 시작한다.

1. PolyData 모델인지 확인한다.
2. `sv4guiMeshFactory::CreateMesh()`로 TetGen GUI mesher를 생성한다.
3. `sv4guiMeshEdit::CreateCmdsT()`에서 GUI 값을 문자열 명령으로 변환한다.
4. `sv4guiMesh::ExecuteCommands()`가 명령을 순서대로 실행한다.
5. 각 명령은 `sv4guiMeshTetGen::ParseCommand()`와
   `sv4guiMeshTetGen::Execute()`를 통과한다.
6. 마지막 `generateMesh` 명령이 `cvTetGenMeshObject::GenerateMesh()`를 호출한다.

대표적인 GUI 명령 기록은 다음과 같다.

```text
option surface 1
option volume 1
option GlobalEdgeSize 0.5
setWalls
boundaryLayer 3 0.5 0.8 0
option BoundaryLayerDirection 1
option GenerateWallMesh
option WallThickness 0.5
option NumberOfWallLayers 2
option WallThicknessSmoothingIterations 5
option WallThicknessCurvatureFactor 0.8
option WallThicknessRadiusFactor 0.0
localWallThickness wall_aorta 0.8
generateMesh
writeMesh
```

명령 기록은 GUI 프로젝트에 저장되며, 기존 메시를 다시 생성할 때 동일한 명령을
재실행한다.

주의: `WallThicknessCurvatureFactor`의 코어 기본값은 `0.8`이므로, 이 명령이
없는 기존 명령 기록을 재실행해도 곡률 기반 두께 제한이 활성화된다. 이는 병합
부위의 자기 교차를 줄이기 위한 의도된 동작 변경이며, 이전과 동일한 결과가
필요하면 `option WallThicknessCurvatureFactor 0.0`을 명시해 비활성화한다.

### 4.2 Python API 경로

Python API에서는 `meshing.TetGenOptions`에 옵션을 저장한 후
`meshing.TetGen.generate_mesh(options)`를 호출한다.

```python
mesher = sv.meshing.create_mesher(sv.meshing.Kernel.TETGEN)
mesher.set_model(model)
mesher.set_walls(wall_face_ids)

options = sv.meshing.TetGenOptions()
options.surface_mesh_flag = True
options.volume_mesh_flag = True
options.global_edge_size = 0.5
options.generate_wall_mesh = True
options.wall_thickness = 0.5
options.number_of_wall_layers = 2
options.wall_thickness_smoothing_iterations = 5
options.wall_thickness_curvature_factor = 0.8
# 반경 적응 두께(선택): centerline이 계산되어 있어야 한다.
# radius_meshing_compute_centerlines 기본값은 False이므로 명시적으로 켜거나
# radius_meshing_centerlines로 centerline을 직접 제공해야 한다.
# options.radius_meshing_on = True
# options.radius_meshing_compute_centerlines = True
# options.wall_thickness_radius_factor = 0.15

mesher.generate_mesh(options)
```

`MesherTetGen_generate_mesh()`의 처리 순서는 다음과 같다.

1. `SetOptions()`로 일반 옵션을 `SetMeshOptions()`에 전달한다.
2. Python mesher에 저장된 wall face ID를 `SetWalls()`에 전달한다.
3. local edge size, centerline radius, sphere refinement용 배열을 생성한다.
4. `cvTetGenMeshObject::GenerateMesh()`를 호출한다.

## 5. 모델 입력과 wall face 설정

### 5.1 모델 로드

`cvTetGenMeshObject::LoadModel()`은 입력 `vtkPolyData`를 다음 두 객체에 복사한다.

- `originalpolydata_`: 변경하지 않는 원본
- `polydatasolid_`: 실제 메시 생성에 사용하는 작업본

입력 모델에는 최소한 각 cell의 경계면을 식별하는 `ModelFaceID` 배열이 있어야 한다.

### 5.2 wall face 설정

`cvTetGenMeshObject::SetWalls()`는 지정된 `ModelFaceID`를 기준으로 `WallID` 배열을
생성하고 `meshwallfirst` 옵션을 활성화한다.

VMTK 경로에서는 wall cell을 직접 추출할 수 있다. MMG 경로에서는 재메싱 후
전달된 `WallID`를 이용하여 wall surface를 다시 추출한다. 결과적으로 경계층과
solid wall 생성에 사용되는 표면은 wall 영역이며, 입구와 출구 cap은 이후 별도
단계에서 다시 생성된다.

## 6. 옵션 검증

`GenerateMesh()`는 실제 메시 생성을 시작하기 전에 옵션 조합을 검증한다.

solid wall mesh를 생성하려면 다음 조건을 모두 만족해야 한다.

- 유체 boundary layer가 활성화되어 있어야 한다.
- `BoundaryLayerDirection`이 `1`이어야 한다.
- 즉, 유체 boundary layer가 원래 모델 표면에서 안쪽으로 생성되어야 한다.
- global wall thickness가 0보다 커야 한다.

이 제약은 다음과 같은 영역 배치를 만들기 위한 것이다.

```text
혈관 중심
   │
   ├─ TetGen 내부 유체 영역
   ├─ VMTK 유체 boundary layer
   ├─ 원래 모델 표면 = fluid/wall interface
   └─ VMTK solid wall mesh
바깥쪽
```

## 7. 표면 재메싱

`cvTetGenMeshObject::GenerateSurfaceRemesh()`는 `polydatasolid_`를 재메싱한다.

주요 동작은 다음과 같다.

1. `meshwallfirst` 여부에 따라 edge 보존과 triangle split 설정을 결정한다.
2. local size 또는 refinement가 있으면 `MeshSizingFunction`을 사용한다.
3. 빌드 설정과 `UseMMG` 옵션에 따라 MMG 또는 VMTK surface remesher를 호출한다.
4. 재메싱 후 원본의 `ModelFaceID`를 가능한 범위에서 복원한다.
5. 연결 영역, free edge, non-manifold/bad edge를 검사한다.

surface mesh만 요청된 경우 이 단계의 `polydatasolid_`가 최종 `surfacemesh_`가 되고
체적 메시 생성은 수행하지 않는다.

## 8. 유체 boundary layer 생성

`cvTetGenMeshObject::GenerateBoundaryLayerMesh()`는 VMTK를 사용하여 wall surface에서
유체 방향의 boundary layer를 생성한다.

### 8.1 표면 준비

1. `vtkPolyDataNormals`로 point normal을 계산하고 방향을 일관되게 만든다.
2. `vtkCleanPolyData`로 중복 point를 정리한다.
3. `VMTKUtils_ComputeSizingFunction()`으로 `MeshSizingFunction`을 생성한다.
4. `SetCapBoundaryNormals()`로 열린 입구/출구 경계의 normal을 보정한다.

cap 경계 normal을 보정하지 않으면 경계의 평균 normal이 cap 평면 밖을 향할 수 있고,
이 경우 boundary layer 끝단이 평평하게 생성되지 않을 수 있다.

### 8.2 경계층 압출

표면을 `vtkUnstructuredGrid`로 변환한 후
`VMTKUtils_BoundaryLayerMesh()`를 호출한다.

VMTK boundary layer generator는 다음 값을 사용한다.

- point normal: 압출 방향
- `MeshSizingFunction`: 위치별 layer 두께
- layer 수
- 전체 두께 비율
- sublayer ratio
- constant thickness 사용 여부
- `CellEntityIds`: surface와 sidewall 구분

`BoundaryLayerDirection == 1`이면 normal을 반전하여 경계층을 모델 안쪽으로 생성한다.
생성된 가장 안쪽 표면은 TetGen이 채울 내부 유체 영역의 새로운 경계가 되며,
`polydatasolid_`에 저장된다.

## 9. 고체 혈관벽 메시 생성

solid wall 옵션이 활성화된 경우 `GenerateBoundaryLayerMesh()` 안에서
`cvTetGenMeshObject::GenerateWallMesh()`가 호출된다.

유체 boundary layer는 원래 표면에서 안쪽으로 생성되지만, solid wall은 같은 원래
표면에서 바깥쪽으로 생성된다. 따라서 원래 모델 표면의 좌표가 fluid/wall interface로
유지된다.

### 9.1 point별 두께 배열 생성

`WallThickness` point data 배열은 다음 원칙으로 만든다.

1. 각 point의 **요청 두께**를 정한다. 기본은 global wall thickness이고,
   `LocalWallThickness`가 지정된 face는 local 값을 사용한다.
2. 여러 face가 만나는 point는 연결된 고유 `ModelFaceID`별 두께를 평균한다.
   local override가 없는 face는 global 두께를 쓴다.
3. 같은 face의 triangle 수는 평균 가중치에 영향을 주지 않아야 한다.
4. `WallThicknessRadiusFactor`가 활성화되면 반경 기반 값으로 요청 두께를
   **조정**한다(9.2 참고). 요청 두께를 대체하지 않고 그 범위 안에서만
   낮추므로, local 두께를 준 face는 그 값이 상한으로 유지된다.

예를 들어 두 face의 두께가 각각 `0.5`, `1.0`이면 공유점의 초기 두께는
두 face의 triangle 밀도와 관계없이 `0.75`가 되어야 한다.

### 9.2 반경 적응 두께

`WallThicknessRadiusFactor`가 `0`보다 크면 각 벽 node의 두께를 단일
상수가 아니라 **국소 혈관 반경 × factor**로 설정한다. 소구경 혈관(예: PCoA)은
자동으로 얇아지고 대구경 혈관(예: ICA)은 두꺼워져, 서로 다른 구경의 혈관이
만나는 접합부 두께가 기하에 따라 전이된다. 하나의 face 안에서도 혈관이
가늘어지는 만큼 두께가 변하므로, face당 상수인 local 두께로는 표현할 수 없는
테이퍼를 처리한다.

- 국소 반경은 `sys_geom_distancetocenterlines()`로 계산한 각 벽 표면 node의
  **centerline까지의 거리**(`DistanceToCenterlines`)를 사용한다.
- centerline은 mesher가 직접 계산하지 않는다. radius(centerline) 메싱 경로가
  계산한 centerline을 `cvTetGenMeshObject::SetCenterlines()`로 전달받아
  재사용한다. 따라서 이 옵션을 쓰려면 centerline이 계산되어 있어야 한다.
  Python에서는 `radius_meshing_on`과 `radius_meshing_compute_centerlines`를
  함께 켜거나(`radius_meshing_compute_centerlines` 기본값은 False),
  `radius_meshing_centerlines`로 centerline을 직접 제공한다. GUI에서는
  centerline(radius) 메싱을 켠다. centerline이 없으면 명확한 에러로 종료한다.
- 두께는 **그 point의 요청 두께**(local 두께가 있으면 그 값, 없으면 global
  `WallThickness`)를 상한, 그 5%를 하한으로 클램프해 반경 추정이 튀거나
  극소혈관에서 두께가 0에 수렴하는 것을 막는다. 상한이 global이 아니라 요청
  두께이므로 local 두께 지정과 반경 적응이 서로를 덮어쓰지 않는다. 하한이
  5%인 것은 대동맥~뇌동맥처럼 혈관 구경이 20배 이상 차이 나는 모델에서
  하한 자체가 반경 적응을 무력화하지 않게 하기 위한 것이다.
- 반경 추정(`DistanceToCenterlines`)은 분지 crotch에서 가장 가까운 centerline이
  옆 가지의 것일 수 있어 접합부에서 가장 부정확하다. 접합부 겹침(fold-over)의
  처방으로는 쓰지 말고(9.5 참고), 혈관 구경차가 큰 모델의 두께 배분 용도로만
  사용한다.
- 유효 범위는 0.0~1.0이며 `0.0`이면 상수 `WallThickness`를 그대로 쓴다.
- 표면 point 좌표(fluid/wall interface)는 이동하지 않는다.

### 9.3 곡률 기반 두께 제한

`WallThicknessCurvatureFactor`가 `0`보다 크면
`TGenUtils_ClampThicknessToConcaveCurvature()`가 스무딩 전에 `WallThickness`
배열을 제한한다.

- 각 point의 오목 곡률을 one-ring 이웃으로부터 추정한다(이웃까지의 거리 `d`,
  법선 방향 높이 `h`에 대해 곡률 근사 `2h/d^2`의 최댓값).
- 두께를 `factor / 곡률`(오목 곡률 반경의 `factor`배)로 제한한다. 제한은 항상
  적용되어 클램프는 멱등적이며, 요청 두께의 10% 미만으로 줄어든 point 개수는
  경고로 보고된다.
- 볼록하거나 평평한 영역은 변경하지 않으며, 표면 point 좌표(fluid/wall
  interface)는 이동하지 않는다.

추정이 국소적이므로 이 제한은 압출된 outer wall의 자기 교차 가능성을 줄일 뿐
전역적인 자기 교차 부재를 보장하지 않는다. 생성된 메시의 품질(음수 체적,
aspect ratio, 표면 교차)은 별도로 확인해야 한다.

### 9.4 두께 스무딩

local wall thickness가 존재하거나 곡률 기반 두께 제한 또는 반경 적응 두께가
활성화된 경우 `TGenUtils_SmoothPointArray()`가 `WallThickness` 배열을 반복
Laplacian averaging으로 스무딩한다.

각 반복에서 point의 새 값은 이전 반복의 다음 값들을 평균하여 계산한다.

- point 자신의 두께
- 같은 cell을 공유하는 one-ring 이웃 point의 두께

모든 새 값은 이전 반복 값에서 동시에 계산하는 Jacobi 방식이므로 point 순회 순서에
결과가 의존하지 않는다.

스무딩은 두께 scalar만 변경한다. 원래 표면 point 좌표와 fluid/wall interface는
이동하지 않는다. 반복 횟수가 `0`이면 스무딩을 수행하지 않는다.

스무딩이 클램프된 값을 이웃 방향으로 다시 끌어올리므로, 곡률 기반 두께 제한이
활성화된 경우 스무딩 후 클램프를 한 번 더 적용해 제한을 복원한다.

### 9.5 fold-over 방지 (압출 기하 검출)

곡률 클램프(9.3)는 one-ring 국소 추정이라 거칠게 메싱된 접합부에서 fold를
놓칠 수 있고, 압출 후처리(`VMTKUtils_ReorderTetElements`)는 음수 체적 tet를
노드 재정렬로 뒤집어 **접힌 요소를 감춘다**(관례상 음수와 실제 fold를 구분
못 함). 그래서 압출 직전에 `TGenUtils_LimitThicknessToPreventFold()`가 실제
압출 기하를 검사한다.

- 각 표면 point의 바깥 벽 정점은 `p + t·n`(두께 × 법선). 각 삼각형에 대해
  바깥 삼각형의 winding을 안쪽과 비교해, 오목부에서 두께가 커 바깥 삼각형이
  반전(또는 거의 붕괴)하면 fold로 판정한다.
- fold로 판정된 삼각형의 정점 두께를 반복적으로 축소(반복당 0.8배, 원래
  두께의 5% 하한)하고, 반전이 사라질 때까지 최대 20회 재검사한다.
- 하한에서도 fold가 남으면 마스킹하지 않고 경고로 보고한다(표면 리메시 또는
  접합부 두께 축소 필요).
- fold 발생 지점의 두께만 줄이므로 항상 실행해도 안전하며, 표면 point 좌표
  (fluid/wall interface)는 이동하지 않는다.
- 국소 검사이므로 서로 떨어진 두 표면 구간이 압출 후 충돌하는 **전역 자기
  교차는 검출하지 못한다.**

### 9.6 바깥쪽 압출

`WallThickness` 배열을 surface에 추가한 후
`VMTKUtils_BoundaryLayerMesh()`를 다시 호출한다.

solid wall 생성 시 주요 설정은 다음과 같다.

- normal 방향을 반전하지 않음: 바깥쪽 압출
- constant thickness 비활성화
- 위치별 `WallThickness` 배열 사용
- wall layer 수는 `NumberOfWallLayers` 사용
- sublayer ratio는 `1.0`: 두께 방향의 layer를 균등하게 분할

결과는 `wallmesh_`에 저장된다.

## 10. 내부 유체 영역 cap 생성

유체 boundary layer의 안쪽 표면은 입구와 출구가 열린 상태다.
`cvTetGenMeshObject::GenerateAndMeshCaps()`는 이 표면을 TetGen 입력으로 사용할 수
있도록 닫는다.

처리 순서는 다음과 같다.

1. VMTK capper로 열린 경계를 닫는다.
2. 원본 모델에서 cap의 `ModelFaceID`를 복원한다.
3. wall surface를 제외하고 새 cap만 재메싱한다.
4. `WallID` 배열을 다시 생성한다.

이 단계가 끝나면 `polydatasolid_`는 TetGen이 내부 유체 체적을 생성할 수 있는
닫힌 PLC(Piecewise Linear Complex) 표면이 된다.

## 11. TetGen sizing function과 입력 생성

### 11.1 sizing function 생성

`cvTetGenMeshObject::GenerateMeshSizingFunction()`은 현재 닫힌 내부 유체 표면에
`MeshSizingFunction` point data 배열을 생성한다.

boundary layer와 TetGen 내부 체적 메시가 만나는 interface에서 node가 일치하려면
TetGen이 interface의 표면 크기를 따라야 하므로 이 단계가 필요하다.

### 11.2 `tetgenio` 변환

`cvTetGenMeshObject::NewMesh()`는 다음 작업을 수행한다.

1. 이전 `inmesh_`와 `outmesh_`를 정리한다.
2. 새 `tetgenio` 입력과 출력을 생성한다.
3. `polydatasolid_`를 clean한다.
4. surface point와 triangle을 TetGen 입력 형식으로 변환한다.
5. 필요하면 `MeshSizingFunction`을 point metric으로 추가한다.
6. `ModelFaceID` 또는 `CellEntityIds` facet marker를 추가한다.
7. hole과 subdomain/region 정보를 추가한다.

boundary layer를 사용하는 경우 TetGen이 interface triangle을 다시 분할하지 않도록
sizing metric과 `nobisect` 설정을 함께 사용한다.

## 12. TetGen 체적 메시 생성

`GenerateMesh()`는 `tetgenbehavior`를 생성하고 옵션에 따라 다음 값을 설정한다.

- `plc`: 닫힌 surface PLC 사용
- `neighout`: tetrahedron adjacency 출력
- global edge size 기반 최대 tetrahedron volume
- quality ratio
- minimum dihedral angle
- optimization level
- epsilon
- no-bisect
- variable region volume
- sizing metric

global edge size `a`에 대응하는 최대 tetrahedron volume은 다음 근사식을 사용한다.

```text
max_volume = a³ / (6 × √2)
```

boundary layer, function-based sizing 또는 refinement를 사용하는 경우 TetGen의
metric 기반 품질 옵션이 활성화된다.

최종적으로 `tetrahedralize()`가 호출되어 `outmesh_`에 내부 유체 tetrahedral mesh를
생성한다.

## 13. 최종 메시 결합

boundary layer를 사용하지 않는 경우 TetGen 결과를 바로 VTK surface/volume mesh로
변환한다.

boundary layer를 사용하는 경우 `cvTetGenMeshObject::AppendBoundaryLayerMesh()`가
`VMTKUtils_AppendData()`를 호출한다.

### 13.1 유체 메시 결합

다음 두 메시를 먼저 결합한다.

- TetGen 내부 유체 tetrahedral mesh
- VMTK 유체 boundary layer mesh

VMTK boundary layer가 wedge, triangle, quad 요소를 포함할 수 있으므로 결합 전에
tetrahedron으로 변환한다. 변환 과정에서 음수 Jacobian이 생기지 않도록 tetrahedron의
node 순서를 보정한다.

### 13.2 solid wall 결합

`wallmesh_`가 존재하면 다음 과정을 추가로 수행한다.

1. wall mesh를 tetrahedron으로 변환한다.
2. tetrahedron node 순서를 보정한다.
3. 유체와 다른 `ModelRegionID`를 부여한다.
4. 유체와 wall의 공유 interface node를 병합한다.
5. 최종 `GlobalNodeID`와 `GlobalElementID`를 생성한다.

일반적인 단일 유체 영역에서는 다음 region 구성이 만들어진다.

| 영역 | 일반적인 `ModelRegionID` |
| --- | --- |
| TetGen 내부 유체 + 유체 boundary layer | `1` |
| solid vessel wall | `2` |

solid wall을 생성하면 `NewRegionBoundaryLayer` 옵션은 무시된다. 이 경우 boundary
layer는 유체 영역에 포함되고, 새 region은 solid wall에 사용된다.

## 14. 표면 ID와 전역 ID 처리

최종 결합 과정은 다음 ID를 정리한다.

- 유체 외부 경계와 cap의 `ModelFaceID`
- wall 내부/외부 표면과 끝단의 `ModelFaceID`
- 유체와 wall의 `ModelRegionID`
- 공유 interface를 포함한 `GlobalNodeID`
- 전체 요소의 `GlobalElementID`

solid wall이 없는 안쪽 boundary layer 메시에서는 최종 surface의 `ModelFaceID`를
원본 모델에서 다시 복원한다.

solid wall이 있는 경우 wall 외부 표면이 원본 모델보다 바깥에 있으므로 최종 전체
surface에 대해 단순한 최근접 원본 face ID 복원을 수행하지 않는다. 대신 결합 과정에서
유체와 wall surface source를 분리하여 face ID를 설정한다.

## 15. 메시 품질 계산

최종 `volumemesh_`가 생성되면 `TGenUtils_ReportMeshQuality()`가 tetrahedral element의
aspect ratio를 계산한다.

결과는 다음 두 방식으로 제공된다.

- 각 cell의 `AspectRatio` 배열
- 표준 출력의 최소/평균/최대값과 품질 평가

현재 평가는 왜곡 요소와 품질이 낮은 요소의 비율에 따라
`GOOD`, `ACCEPTABLE`, `POOR`로 분류한다.

## 16. 주요 실행 분기

| Surface | Volume | Boundary layer | Solid wall | 결과 |
| --- | --- | --- | --- | --- |
| On | Off | Off | Off | 재메싱된 surface만 생성 |
| On/Off | On | Off | Off | TetGen 단일 체적 메시 생성 |
| On | On | On | Off | TetGen 내부 메시와 유체 boundary layer 결합 |
| On | On | On | On | 내부 유체, 유체 boundary layer, solid wall을 결합한 FSI 메시 생성 |

solid wall은 boundary layer 없이 단독으로 생성할 수 없다.

## 17. 오류 발생 지점

메시 생성은 다음 조건에서 중단될 수 있다.

- 모델 또는 `ModelFaceID`가 없음
- wall face가 지정되지 않음
- 잘못된 surface topology
- 여러 연결 영역이 있지만 multiple region이 허용되지 않음
- free edge 또는 bad edge 존재
- wall mesh 옵션과 boundary layer 방향이 호환되지 않음
- wall thickness가 설정되지 않음
- VMTK/MMG surface remeshing 실패
- boundary layer 또는 wall 압출 실패
- cap 생성 또는 cap 재메싱 실패
- TetGen 입력 변환 실패
- TetGen tetrahedralization 예외
- 메시 결합 또는 region/global ID 생성 실패

## 18. 구현 불변조건

향후 로직을 수정할 때 다음 조건을 유지해야 한다.

1. fluid/wall interface point 좌표는 wall thickness 스무딩으로 이동하지 않는다.
2. local wall thickness 접합부의 초기값은 triangle 수가 아니라 고유 face를 기준으로
   계산한다.
3. wall thickness 스무딩은 이전 iteration 값만 읽는 순서 독립 방식이어야 한다.
4. TetGen 내부 유체와 boundary layer interface의 node는 일치해야 한다.
5. fluid/wall interface의 전역 노드 ID는 양쪽 영역에서 일관되어야 한다.
6. solid wall 생성 시 boundary layer는 유체 region에 남아야 한다.
7. 최종 surface에는 의미 있는 `ModelFaceID`, volume에는 `ModelRegionID`가 있어야 한다.
8. VMTK 요소를 tetrahedron으로 변환한 후 음수 Jacobian이 생기지 않도록 node 순서를
   보정해야 한다.

## 19. 검증 권장 항목

Mac 개발 환경에서는 실제 SimVascular 빌드와 VMTK/TetGen 실행을 검증할 수 없으므로,
Linux Docker 빌드 및 실행 환경에서 다음 항목을 확인한다.

1. 일반 TetGen surface/volume 메시 생성
2. local edge size가 있는 메시 생성
3. 안쪽 fluid boundary layer 생성
4. uniform solid wall 생성
5. 서로 다른 local wall thickness를 가진 두 face의 접합부
6. 세 개 이상의 face가 만나는 wall thickness junction
7. smoothing iteration `0`, 기본값, 최대 허용값
8. 최종 유체와 wall의 `ModelRegionID`
9. fluid/wall interface의 중복 point 및 `GlobalNodeID`
10. tetrahedron Jacobian과 `AspectRatio`
11. mesh-complete 파일 출력 후 surface ID 보존

빌드 및 실행 로그는 저장소의 `logs/` 디렉터리에 저장하고 커밋한다.
