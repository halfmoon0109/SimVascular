/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/** @file sv_tetgenmesh_utils.cxx
 *  @brief The implementations of functions in cv_tetgenmesh_utils
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "SimVascular.h"

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkMeshQuality.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkImplicitPolyDataDistance.h"
#include "vtkFlyingEdges3D.h"
#include "vtkStaticPointLocator.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkClipPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkMath.h"

#include "simvascular_tetgen.h"

#include "sv_polydatasolid_utils.h"
#include "sv_misc_utils.h"
#include "sv_vtk_utils.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

#define MAXPATHLEN 1024

#ifdef SV_USE_ZLIB
  #ifdef SV_USE_SYSTEM_ZLIB
    #include <zlib.h>
  #else
    #include "simvascular_zlib.h"
  #endif
#else
  #include <stdlib.h>
  #define gzopen fopen
  #define gzprintf fprintf
  #define gzFile FILE*
  #define gzclose fclose
#endif

#include "sv_tetgenmesh_utils.h"

// -----------------------------
// cvTetGenMeshObjectUtils_Init()
// -----------------------------
int TGenUtils_Init()
{
  return SV_OK;
}

// -----------------------------
// cvTGenUtils_ConvertToTetGen()
// -----------------------------
/**
 * @brief Takes vtkPolyData and turns it into TetGen data structures
 * @param *inmesh tetgenio structure in which to hold the input surface
 * @param *polydatasolid the solid in which to take the discrete points
 * from in order to form the mesh to put into TetGen
 * @return SV_OK if function completes properly
 */

int TGenUtils_ConvertSurfaceToTetGen(tetgenio *inmesh,vtkPolyData *polydatasolid)
{
  tetgenio::facet *f;
  tetgenio::polygon *p;

  //All input numbers start from zero, all outmesh_put number start from zero
  inmesh->firstnumber = 0;
  inmesh->numberofpoints = polydatasolid->GetNumberOfPoints();
  inmesh->pointlist = new REAL[inmesh->numberofpoints*3];

  //Do Point transition from polydatasolid into pointlist
  //fprintf(stderr,"Converting Points...\n");
  vtkSmartPointer<vtkPoints> inPts = vtkSmartPointer<vtkPoints>::New();
  inPts = polydatasolid->GetPoints();
  for (int i=0; i<inmesh->numberofpoints;i++)
  {
    double polyPt[3];
    inPts->GetPoint(i,polyPt);
    inmesh->pointlist[i*3] = polyPt[0];
    inmesh->pointlist[i*3+1] = polyPt[1];
    inmesh->pointlist[i*3+2] = polyPt[2];
  }

  // Convert faces
  inmesh->numberoffacets = (int) polydatasolid->GetNumberOfPolys();
  inmesh->facetlist = new tetgenio::facet[inmesh->numberoffacets];
  inmesh->facetmarkerlist = new int[inmesh->numberoffacets];

  //fprintf(stderr,"Converting Faces...\n");
  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds->SetNumberOfIds(3);
  for (int i=0;i<inmesh->numberoffacets;i++)
  {
    polydatasolid->GetCellPoints(i,ptIds);

    f = &inmesh->facetlist[i];
    f->numberofpolygons=1;

    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

    f->numberofholes = 0;
    f->holelist = nullptr;

    p = &f->polygonlist[0];
    p->numberofvertices=3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] =  (int) ptIds->GetId(0);
    p->vertexlist[1] =  (int) ptIds->GetId(1);
    p->vertexlist[2] =  (int) ptIds->GetId(2);

  }

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_AddPointSizingFunction
// -----------------------------
/**
 */
int TGenUtils_AddPointSizingFunction(tetgenio *inmesh,vtkPolyData *polydatasolid,
    std::string meshSizingFunctionName, double maxEdgeSize)
{
  // check number of points
  if (polydatasolid->GetNumberOfPoints() != inmesh->numberofpoints)
  {
    fprintf(stderr,"surface and tetgen object must match and must have already been converted to add sizing function\n");
    return SV_ERROR;
  }

  //Check if applying a mesh sizing function and initiate point metrics list
  //if ok
  if (VtkUtils_PDCheckArrayName(polydatasolid, 0, meshSizingFunctionName) != SV_OK)
  {
    fprintf(stderr,"Array name 'MeshSizingFunction' does not exist. \
        Something may have gone wrong when setting up BL\n");
    return SV_ERROR;
  }
  vtkDoubleArray *meshSizingFunction = vtkDoubleArray::SafeDownCast(
    polydatasolid->GetPointData()->GetArray(meshSizingFunctionName.c_str()));

  inmesh->numberofpointmtrs = 1;
  inmesh->pointmtrlist = new REAL[inmesh->numberofpoints];

  //Do Point transition from polydatasolid into pointlist
  fprintf(stderr,"Adding mesh sizing metric...\n");
  for (int i=0; i<inmesh->numberofpoints;i++)
  {
    inmesh->pointmtrlist[i] = meshSizingFunction->GetComponent(i,0);
    if (inmesh->pointmtrlist[i] == 0.0)
      inmesh->pointmtrlist[i] = maxEdgeSize;
  }

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_AddFacetMarkers
// -----------------------------
/**
 */
int TGenUtils_AddFacetMarkers(tetgenio *inmesh,vtkPolyData *polydatasolid,
    std::string markerListArrayName)
{
  // Check to make sure number of facets matches
  if (polydatasolid->GetNumberOfPolys() != inmesh->numberoffacets)
  {
    fprintf(stderr,"surface and tetgen object must match and must have already been converted to add facet markers\n");
    return SV_ERROR;
  }

  //Do Poly transition from polydatasolid into facetlist
  if (VtkUtils_PDCheckArrayName(polydatasolid,1,markerListArrayName) != SV_OK)
  {
    fprintf(stderr,"Array name does not exist in polydata. Regions must be identified \
      and named prior to this function call\n");
    return SV_ERROR;
  }
  vtkIntArray *boundaryScalars = static_cast<vtkIntArray*>(polydatasolid->GetCellData()->GetScalars(markerListArrayName.c_str()));

  fprintf(stderr,"Adding Facet Markers...\n");
  for (int i=0;i<inmesh->numberoffacets;i++)
  {
    double boundarymarker = (int) boundaryScalars->GetValue(i);
    inmesh->facetmarkerlist[i]=boundarymarker;
  }

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_AddHoles
// -----------------------------
/**
 */
int TGenUtils_AddHoles(tetgenio *inmesh, vtkPoints *holeList)
{
  inmesh->numberofholes = holeList->GetNumberOfPoints();
  inmesh->holelist = new REAL[inmesh->numberofholes * 3];

  for (int i=0; i<inmesh->numberofholes; i++)
  {
    double pt[3];
    holeList->GetPoint(i, pt);
    inmesh->holelist[3*i] =   pt[0];
    inmesh->holelist[3*i+1] = pt[1];
    inmesh->holelist[3*i+2] = pt[2];
  }

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_AddRegions
// -----------------------------
/**
 */
int TGenUtils_AddRegions(tetgenio *inmesh, vtkPoints *regionList, vtkDoubleArray *regionSizeList)
{
  inmesh->numberofregions = regionList->GetNumberOfPoints();
  inmesh->regionlist = new REAL[inmesh->numberofregions * 5];

  for (int i=0; i<inmesh->numberofregions; i++)
  {
    double pt[3];
    regionList->GetPoint(i, pt);
    inmesh->regionlist[5*i] =   pt[0];
    inmesh->regionlist[5*i+1] = pt[1];
    inmesh->regionlist[5*i+2] = pt[2];
    inmesh->regionlist[5*i+3] = i;

    double mES = regionSizeList->GetTuple1(i);
    double maxvol = (mES*mES*mES)/(6*sqrt(2.));
    inmesh->regionlist[5*i+4] = maxvol;
  }

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_ConvertVolumeToTetGen()
// -----------------------------
/**
 * @brief Function to convert the current mesh to a tetgen mesh object to be
 * able to remesh
 * @param mesh This is the full mesh to be remeshed
 * @param surfaceMesh This is the intial mesh; If we don't need the final
 * mesh regions, then we don't have to actually use this
 * @param inmesh This is the tegen mesh object to be transferred to
 */

int TGenUtils_ConvertVolumeToTetGen(vtkUnstructuredGrid *mesh,vtkPolyData *surfaceMesh,
    tetgenio *inmesh)
{
  int numTets,numPolys;
  int numPoints,numSurfacePoints;
  double tetPts[3];
  tetgenio::facet *f;
  tetgenio::polygon *p;
  vtkIdType i,j;
  vtkIdType npts = 0;
  const vtkIdType *pts;
  vtkIdType cellId;
  vtkSmartPointer<vtkPoints> uPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> pPolys = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> uTets = vtkSmartPointer<vtkCellArray>::New();
  vtkIntArray *boundaryScalars;
  vtkDoubleArray *errorMetricArray;

  mesh->BuildLinks();
  numTets = mesh->GetNumberOfCells();
  numPoints = mesh->GetNumberOfPoints();
  uPoints = mesh->GetPoints();
  uTets = mesh->GetCells();

  numSurfacePoints = surfaceMesh->GetNumberOfPoints();
  numPolys = surfaceMesh->GetNumberOfPolys();
  pPolys = surfaceMesh->GetPolys();
  boundaryScalars = vtkIntArray::SafeDownCast(surfaceMesh->GetCellData()->GetArray("ModelFaceID"));
  errorMetricArray = vtkDoubleArray::SafeDownCast(mesh->GetPointData()->GetArray("errormetric"));

  cout<<"Num Cells "<<numTets<<endl;
  cout<<"Num Points "<<numPoints<<endl;
  inmesh->firstnumber = 0;
  inmesh->numberofcorners = 4;
  inmesh->numberoftetrahedra = numTets;
  inmesh->numberofpoints = numPoints;
  inmesh->pointlist = new double[numPoints*3];
  inmesh->tetrahedronlist = new int[numTets*4];
  inmesh->numberofpointmtrs = 1;
  inmesh->pointmtrlist = new REAL[numPoints*inmesh->numberofpointmtrs];

  cout<<"Converting to Adapt Points..."<<endl;
  for (i = 0; i < numPoints; i++)
  {
    uPoints->GetPoint(i,tetPts);
    inmesh->pointlist[i*3] = tetPts[0];
    inmesh->pointlist[i*3+1] = tetPts[1];
    inmesh->pointlist[i*3+2] = tetPts[2];
    inmesh->pointmtrlist[i] = errorMetricArray->GetValue(i);
  }

  cout<<"Converting to Adapt Tets..."<<endl;
  for (i=0,uTets->InitTraversal();uTets->GetNextCell(npts,pts);i++)
  {
    for (j = 0;j < npts;j++)
    {
      inmesh->tetrahedronlist[i*npts+j] = pts[j];
    }
  }

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_ConvertToVTK()
// -----------------------------
//
/**
 * @brief Takes tetgenio and turns it into an output vtkPolyData and
 * vtkUnstructuredGrid
 * @param *outmesh tetgen structure for which the mesh is output
 * @param *volumemesh vtkPolyData on which to save the surface mesh
 * @param *surfacemesh vtkUnstructuredGrid on which to save the volume mesh
 * @return SV_OK if function completes properly
 */

int TGenUtils_ConvertToVTK(tetgenio *outmesh,vtkUnstructuredGrid *volumemesh,vtkPolyData *surfacemesh,int *modelRegions,int getBoundary)
{
  int modelId = 1;
  int globalId = 1;
  int count=0;
  int totRegions=0;
  double tmp;
  vtkIdType i, j;
  vtkIdType vtkId;
  vtkIdType npts = 0;
  const vtkIdType *pts;

  vtkIdType numPts,numPolys,numFaces;

  vtkSmartPointer<vtkUnstructuredGrid> fullUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPolyData> fullPolyData = vtkSmartPointer<vtkPolyData>::New();

  //Create pointers to vtk scalar lists, point lists, and element lists
  vtkSmartPointer<vtkIdList> polyPointIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> facePointIds = vtkSmartPointer<vtkIdList>::New();

  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> vtpPoints = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkIntArray> modelRegionIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> globalNodeIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> globalElementIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> vtpNodeIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> vtpFaceIds = vtkSmartPointer<vtkIntArray>::New();

  vtkSmartPointer<vtkIntArray> boundaryScalars = vtkSmartPointer<vtkIntArray>::New();

  //Get number of points, polys, and faces
  numPts = outmesh->numberofpoints;
  numPolys = outmesh->numberoftetrahedra;
  numFaces = outmesh->numberoftrifaces;

  bool *pointOnSurface = new bool[numPts];
  int *pointMapping = new int[numPts];

  //Save all point information in a vtkPoints list
  //fprintf(stderr,"Converting Points to VTK Structures...\n");
  points->SetNumberOfPoints(numPts);
  for (i=0;i< outmesh->numberofpoints; i++)
  {
    points->SetPoint(i,outmesh->pointlist[i*3],outmesh->pointlist[i*3+1],outmesh->pointlist[i*3+2]);
    globalNodeIds->InsertValue(i,globalId);
    pointOnSurface[i] = false;
    globalId++;
  }

  //Save all point information in a vtkPoints list
  for (i=0;i<numFaces;i++)
  {
    for (j=0;j<3;j++)
    {
      if (pointOnSurface[outmesh->trifacelist[3*i+j]] == false)
      {
        pointOnSurface[outmesh->trifacelist[3*i+j]] = true;
        pointMapping[outmesh->trifacelist[3*i+j]] = count++;
      }
    }
  }

  vtpPoints->SetNumberOfPoints(count);
  //Create face point list
  for (i=0;i<numPts;i++)
  {
    if (pointOnSurface[i] == true)
    {
      vtpPoints->SetPoint(pointMapping[i],outmesh->pointlist[i*3],outmesh->pointlist[i*3+1],outmesh->pointlist[i*3+2]);
      vtpNodeIds->InsertValue(pointMapping[i],i+1);
    }
  }

  //Save all element information in a vtkCellArray list
  //fprintf(stderr,"Converting Elements to VTK Structures...\n");
  polyPointIds->SetNumberOfIds(4);
  globalId=1;
  for (i=0;i< numPolys;i++)
  {
    for (j=0; j< outmesh->numberofcorners;j++)
    {
      vtkId = outmesh->tetrahedronlist[i*outmesh->numberofcorners+j];
      polyPointIds->SetId(j,vtkId);
    }

    if (outmesh->numberoftetrahedronattributes > 0)
      modelRegionIds->InsertValue(i, outmesh->tetrahedronattributelist[i] + 1);
    else
      modelRegionIds->InsertValue(i,modelId);

    globalElementIds->InsertValue(i,globalId);
    globalId++;
    polys->InsertNextCell(polyPointIds);
  }

  //Create an unstructured grid and link scalar information to nodes and
  //elements
  fullUGrid->SetPoints(points);
  fullUGrid->SetCells(VTK_TETRA, polys);

  modelRegionIds->SetName("ModelRegionID");
  fullUGrid->GetCellData()->AddArray(modelRegionIds);
  fullUGrid->GetCellData()->SetActiveScalars("ModelRegionID");

  globalNodeIds->SetName("GlobalNodeID");
  fullUGrid->GetPointData()->AddArray(globalNodeIds);

  globalElementIds->SetName("GlobalElementID");
  fullUGrid->GetCellData()->AddArray(globalElementIds);

  //Save all external faces to a vtkCellArray list
  //fprintf(stdout,"Converting Faces to VTK Structures...\n");
  facePointIds->SetNumberOfIds(3);

  for (i=0;i< numFaces;i++)
  {
    for (j=0; j<3;j++)
    {
      facePointIds->SetId(j,pointMapping[outmesh->trifacelist[i*3+j]]);
    }

    faces->InsertNextCell(facePointIds);

    if (!(outmesh->adjtetlist[2*i] >= numPolys || outmesh->adjtetlist[2*i] < 0))
    {
      vtpFaceIds->InsertValue(i,globalElementIds->GetValue(outmesh->adjtetlist[2*i]));
    }
    else if (!(outmesh->adjtetlist[2*i+1] >= numPolys || outmesh->adjtetlist[2*i+1] < 0))
    {
      vtpFaceIds->InsertValue(i,globalElementIds->GetValue(outmesh->adjtetlist[2*i+1]));
    }
    else
    {
      fprintf(stderr,"WARNING: TetGen says face has no adjacent tetrahedron\n");
      vtpFaceIds->InsertValue(i,globalElementIds->GetValue(outmesh->adjtetlist[2*i+1]));
    }

    if (getBoundary)
    {
      if (outmesh->trifacemarkerlist != nullptr)
      {
        boundaryScalars->InsertValue(i,outmesh->trifacemarkerlist[i]);
      }
      if (boundaryScalars->GetValue(i)>totRegions)
      {
        totRegions = outmesh->trifacemarkerlist[i];
      }
    }
  }

  //Create a polydata grid and link scalar information to nodes and elements
  fullPolyData->SetPoints(vtpPoints);
  fullPolyData->SetPolys(faces);
  fullPolyData->BuildCells();

  vtpNodeIds->SetName("GlobalNodeID");
  fullPolyData->GetPointData()->AddArray(vtpNodeIds);
  fullPolyData->GetPointData()->SetActiveScalars("GlobalNodeID");

  vtpFaceIds->SetName("GlobalElementID");
  fullPolyData->GetCellData()->AddArray(vtpFaceIds);
  fullPolyData->GetCellData()->SetActiveScalars("GlobalElementID");

  if (getBoundary)
  {
    boundaryScalars->SetName("ModelFaceID");
    fullPolyData->GetCellData()->AddArray(boundaryScalars);
    fullPolyData->GetCellData()->SetActiveScalars("ModelFaceID");

    *modelRegions = totRegions;
  }

//  //Flip the cells on the polydata for presolver
  for (i=0;i<fullPolyData->GetNumberOfCells();i++)
  {
    fullPolyData->GetCellPoints(i,npts,pts);

    // // this could fail if pts has more than 2 components?
    const vtkIdType pts_[3] = {pts[1], pts[0], pts[2]};
    // tmp = pts[0];
    // pts[0] = pts[1];
    // pts[1] = tmp;
    fullPolyData->ReplaceCell(i,npts,pts);
  }

  delete [] pointMapping;
  delete [] pointOnSurface;

  surfacemesh->DeepCopy(fullPolyData);
  volumemesh->DeepCopy(fullUGrid);

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_WriteVTU()
// -----------------------------
/**
 * @brief Writes a vtu file file
 * @param *filename Name of desired file location
 * @param *UGrid vtkUnstructuredGrid to be written
 * @return SV_OK if function completes properly
 */

int TGenUtils_WriteVTU(char *filename,vtkUnstructuredGrid *UGrid)
{
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename);

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(UGrid);
#else
  writer->SetInputData(UGrid);
#endif

  writer->Write();
  return SV_OK;
}

// -----------------------------
// cvTGenUtils_WriteVTP()
// -----------------------------
/**
 * @brief Writes a vtp file file
 * @param *filename Name of desired file location
 * @param *UGrid vtkPolyData to be written
 * @return SV_OK if function completes properly
 */

int TGenUtils_WriteVTP(char *filename,vtkPolyData *PData)
{
  vtkSmartPointer<vtkXMLPolyDataWriter> writer  = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  std::string fn = "out.vtp";
  writer->SetFileName(fn.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(PData);
#else
  writer->SetInputData(PData);
#endif
  //writer->SetDataModeToAscii();

  writer->Write();

  return SV_OK;
}


// -----------------------------
// cvTGenUtils_GetFacePolyData()
// -----------------------------
/**
 * @brief Based on Scalars Defined by the GetBoundaryFaces filter,
 * separate into face VTKs
 * @param *mesh vtkPolyData on which to extract the face from
 * @param *face vtkPolyData on which to set the face PolyData
 * @param angle double that specifies the extraction angle. Any faces
 * @param id int that specifies the face id to extract
 * @return SV_OK if function completes properly
 * @note There is another method to do this that does not retain id
 * information. It may be faster, but doesn't reatain info
 */
//

int TGenUtils_GetFacePolyData(int id,vtkPolyData *mesh, vtkPolyData *face)
{
  //Initiate variable used by function
  int i,j;
  int count=0;
  vtkIdType cellId;
  vtkIdType npts = 0;
 const vtkIdType *pts;
  vtkIdType globalElement2=-1;
  double ptCmps[3];

  vtkSmartPointer<vtkPolyData> tempFace = vtkSmartPointer<vtkPolyData>::New();

  vtkSmartPointer<vtkIdList> facePointIds = vtkSmartPointer<vtkIdList>::New();

  vtkSmartPointer<vtkCellArray> meshFaces = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> selectFaces = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> meshPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> selectPoints = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkIntArray> globalNodeIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> globalElementIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> boundaryScalars = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> lessNodeIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> lessElementIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> globalElement2Ids = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> modelFaceIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> modelFaceRegionIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> modelRegionIds = vtkSmartPointer<vtkIntArray>::New();

  if (VtkUtils_PDCheckArrayName(mesh,0,"GlobalNodeID") != SV_OK)
  {
    fprintf(stderr,"Array name 'GlobalNodeID' does not exist.");
    fprintf(stderr," IDs on mesh may not have been assigned properly\n");
    return SV_ERROR;
  }

  if (VtkUtils_PDCheckArrayName(mesh,1,"GlobalElementID") != SV_OK)
  {
    fprintf(stderr,"Array name 'GlobalElementID' does not exist.");
    fprintf(stderr," IDs on mesh may not have been assigned properly\n");
    return SV_ERROR;
  }

  if (VtkUtils_PDCheckArrayName(mesh,1,"ModelFaceID") != SV_OK)
  {
    fprintf(stderr,"Array name 'ModelFaceID' does not exist. Regions must be identified");
		fprintf(stderr," and named 'ModelFaceID' prior to this function call\n");
    return SV_ERROR;
  }

  bool has_ModelRegionIDs = false;

  if (VtkUtils_PDCheckArrayName(mesh,1,"ModelRegionID") == SV_OK)
  {
    has_ModelRegionIDs = true;
    modelRegionIds = vtkIntArray::SafeDownCast(mesh->GetCellData()->GetScalars("ModelRegionID"));
  }

  globalNodeIds = vtkIntArray::SafeDownCast(mesh->GetPointData()->GetScalars("GlobalNodeID"));
  globalElementIds = vtkIntArray::SafeDownCast(mesh->GetCellData()->GetScalars("GlobalElementID"));
  boundaryScalars = vtkIntArray::SafeDownCast(mesh->GetCellData()->GetScalars("ModelFaceID"));

  meshFaces = mesh->GetPolys();
  meshPoints = mesh->GetPoints();
  int numPts = mesh->GetNumberOfPoints();
  int numFaces = mesh->GetNumberOfPolys();

  bool *cellOnFace = new bool[numFaces];
  bool *pointOnFace = new bool[numPts];
  int *pointMapping = new int[numPts];

  for (i = 0;i<numPts;i++)
  {
    pointOnFace[i] = false;
    pointMapping[i] = -1;
  }

  //Set up point mapping and boolean whether point is on face
  for (cellId = 0,meshFaces->InitTraversal();meshFaces->GetNextCell(npts,pts);cellId++)
  {
    if (boundaryScalars->GetValue(cellId) == id)
    {
      cellOnFace[cellId] = true;
      for(j=0;j<npts;j++)
      {
	if (pointOnFace[pts[j]] == false)
	{
	  pointOnFace[pts[j]] = true;
	  pointMapping[pts[j]] = count++;
	}
      }
    }
    else
    {
      cellOnFace[cellId] = false;
    }
  }

  selectPoints->SetNumberOfPoints(count);

  for (i=0;i<numPts;i++)
  {
    if (pointOnFace[i] == true)
    {
      meshPoints->GetPoint(i,ptCmps);
      selectPoints->SetPoint(pointMapping[i],ptCmps[0],ptCmps[1],ptCmps[2]);
      lessNodeIds->InsertValue(pointMapping[i],globalNodeIds->GetValue(i));
    }
  }

  facePointIds->SetNumberOfIds(3);
  //Get node and element information for the current boundary on the full
  //polydata and save to a smaller polydata
  for(cellId = 0,meshFaces->InitTraversal();meshFaces->GetNextCell(npts,pts); cellId++)
  {
    if (cellOnFace[cellId] == true)
    {
      for (j=0; j<npts; j++)
      {
        facePointIds->SetId(j,pointMapping[pts[j]]);
      }
      selectFaces->InsertNextCell(facePointIds);
      globalElement2Ids->InsertNextValue(globalElement2);
      lessElementIds->InsertNextValue(globalElementIds->GetValue(cellId));
      modelFaceIds->InsertNextValue(id);
      if (has_ModelRegionIDs) { 
        modelFaceRegionIds->InsertNextValue(modelRegionIds->GetValue(cellId));
      }
    }
  }

  //Create links between points and faces and respective global node and
  //element information
  tempFace->SetPoints(selectPoints);
  tempFace->SetPolys(selectFaces);

  lessNodeIds->SetName("GlobalNodeID");
  tempFace->GetPointData()->AddArray(lessNodeIds);
  tempFace->GetPointData()->SetActiveScalars("GlobalNodeID");

  globalElement2Ids->SetName("GlobalElementID2");
  tempFace->GetCellData()->AddArray(globalElement2Ids);
  tempFace->GetCellData()->SetActiveScalars("GlobalElementID2");

  lessElementIds->SetName("GlobalElementID");
  tempFace->GetCellData()->AddArray(lessElementIds);
  tempFace->GetCellData()->SetActiveScalars("GlobalElementID");

  modelFaceIds->SetName("ModelFaceID");
  tempFace->GetCellData()->AddArray(modelFaceIds);
  tempFace->GetCellData()->SetActiveScalars("ModelFaceID");

  if (has_ModelRegionIDs) { 
    modelFaceRegionIds->SetName("ModelRegionID");
    tempFace->GetCellData()->AddArray(modelFaceRegionIds);
    tempFace->GetCellData()->SetActiveScalars("ModelRegionID");
  }

  // Add cell normals.
  //
  auto normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SplittingOff();
  normals->ConsistencyOn();
  normals->AutoOrientNormalsOn();
  normals->ComputeCellNormalsOn();
  normals->ComputePointNormalsOff();
  normals->SetInputData(tempFace);
  normals->Update();

  tempFace->DeepCopy(normals->GetOutput());
  tempFace->GetCellData()->GetNormals()->SetName("Normals");

  delete [] pointOnFace;
  delete [] pointMapping;
  delete [] cellOnFace;

  face->DeepCopy(tempFace);

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_writeDiffAdj()
// -----------------------------
/**
 * @brief This is the new way to write an adjacency file based on the mesh
 * @note now implemented in the presolver as new command
 */

int TGenUtils_writeDiffAdj(vtkUnstructuredGrid *volumemesh)
{
  gzFile myfile = nullptr;

  std::string filename("compareAdjacency.xadj");

  #ifdef SV_USE_ZLIB
  char filenamegz[MAXPATHLEN];
  filenamegz[0]='\0';
  sprintf (filenamegz, "%s.gz", filename.c_str());
  myfile = gzopen (filenamegz, "wb");
  if (myfile == nullptr) {
      fprintf(stderr,"Error: Could not open output file %s.\n",filenamegz);
      return SV_ERROR;
  }
  #else
  myfile = gzopen (filename.c_str(), "wb");
  if (myfile == nullptr) {
      fprintf(stderr,"Error: Could not open output file %s.\n",filename.c_str());
      return SV_ERROR;
  }
  #endif

  int i;
  int numCells;
  int *xadj;
  int *adjacency;
  vtkIdType cellId;
  vtkIdType meshCellId;
  vtkIdType p1,p2,p3;
  vtkIdType ns = 0;
  vtkIdType npts = 0;
  const vtkIdType *pts;
  vtkSmartPointer<vtkCellArray> volCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIntArray> globalIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  volumemesh->BuildLinks();

  if (VtkUtils_UGCheckArrayName(volumemesh,1,"GlobalElementID") != SV_OK)
  {
    fprintf(stderr,"Array name 'GlobalElementID' does not exist. IDs on mesh may not have been assigned properly\n");
    return SV_ERROR;
  }
  globalIds = vtkIntArray::SafeDownCast(volumemesh->GetCellData()->GetScalars("GlobalElementID"));
  numCells = volumemesh->GetNumberOfCells();
  volCells = volumemesh->GetCells();

  xadj = new int[numCells];
  adjacency = new int[4*numCells];
  int adj = 0;
  int xcheck = 0;
  xadj[xcheck] = 0;

  ptIds->SetNumberOfIds(3);
  for (cellId = 0;cellId<numCells;cellId++)
  {
    meshCellId = globalIds->LookupValue(cellId+1);
    volumemesh->GetCellPoints(meshCellId,npts,pts);
    for (i=0;i < npts; i++)
    {
      p1 = pts[i];
      p2 = pts[(i+1)%(npts)];
      p3 = pts[(i+2)%(npts)];

      ptIds->InsertId(0,p1);
      ptIds->InsertId(1,p2);
      ptIds->InsertId(2,p3);

      volumemesh->GetCellNeighbors(meshCellId,ptIds,cellIds);

      //If it is zero, it is a face on the exterior. Otherwise, it has
      //neighbors
      if (cellIds->GetNumberOfIds() != 0)
      {
	adjacency[adj++] = (int) globalIds->GetValue(cellIds->GetId(0)-1);
      }

    }
    xadj[++xcheck] = adj;
  }

  gzprintf(myfile,"xadj: %i\n",numCells+1);
  gzprintf(myfile,"adjncy: %i\n",adj);

  for (i=0;i < numCells+1; i++)
  {
      gzprintf(myfile,"%i\n",xadj[i]);
  }
  for (i=0;i < adj; i++)
  {
      gzprintf(myfile,"%i\n",adjacency[i]);
  }

  delete xadj;
  delete adjacency;

  gzclose(myfile);
  return SV_OK;
}

// -----------------------------
// cvTGenUtils_SetRefinementCylinder()
// -----------------------------
/**
 * @brief computes the distance between each point on surface and center
 * @brief of cylinder. Then, if inside radius, the meshsizing function at the
 * @brief is set to the reduced size,
 * @param size This is the smaller refined of the edges within cylinder region.
 * @param radius This is the radius of the refinement cylinder.
 * @param center This is the center of the refinement cylinder.
 * @param length This is the length of the cylinder. Center is half the length.
 * @param normal This is the normal direction of the length of the cylinder.
 * It is normalized before being used for compuation.
 * @return SV_OK if function completes properly
 */

int TGenUtils_SetRefinementCylinder(vtkPolyData *polydatasolid,
    std::string sizingFunctionArrayName,double size,double radius, double *center,
    double length, double *normal, int secondarray,double maxedgesize,
    std::string refineIDArrayName, int refinecount)
{
  int numPts;
  double disttopoint;
  double distalonglength;
  double pts[3];
  double norm[3];
  for (int i=0;i < 3;i++)
    norm[i] = normal[i];
  vtkIdType pointId;
  vtkSmartPointer<vtkDoubleArray> meshSizeArray = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> refineIDArray = vtkSmartPointer<vtkIntArray>::New();

  //Set sizing function params
  numPts = polydatasolid->GetNumberOfPoints();
  if (secondarray)
  {
    if (VtkUtils_PDCheckArrayName(polydatasolid,0,sizingFunctionArrayName) != SV_OK)
    {
      fprintf(stderr,"Solid does not contain a double array of name %s. Regions must be identified \
		      Reset or remake the array and try again\n",sizingFunctionArrayName.c_str());
      return SV_ERROR;
    }
    meshSizeArray = vtkDoubleArray::SafeDownCast(polydatasolid->GetPointData()->GetArray(sizingFunctionArrayName.c_str()));
    if (VtkUtils_PDCheckArrayName(polydatasolid,0,refineIDArrayName) != SV_OK)
    {
      fprintf(stderr,"Solid does not contain an int array of name %s. Regions must be identified \
		      Reset or remake the array and try again\n",refineIDArrayName.c_str());
      return SV_ERROR;
    }
    refineIDArray = vtkIntArray::SafeDownCast(polydatasolid->GetPointData()->GetArray(refineIDArrayName.c_str()));
  }
  else
  {
    meshSizeArray->SetNumberOfComponents(1);
    meshSizeArray->Allocate(numPts,1000);
    meshSizeArray->SetNumberOfTuples(numPts);
    meshSizeArray->SetName(sizingFunctionArrayName.c_str());
    refineIDArray->SetNumberOfComponents(1);
    refineIDArray->Allocate(numPts,1000);
    refineIDArray->SetNumberOfTuples(numPts);
    refineIDArray->SetName(refineIDArrayName.c_str());
    for (pointId = 0;pointId<numPts;pointId++)
    {
      meshSizeArray->SetValue(pointId,0.0);
      refineIDArray->SetValue(pointId,0);
    }
  }

  for (pointId = 0;pointId<numPts;pointId++)
  {
    polydatasolid->GetPoint(pointId,pts);
    //compute distance
    double pvec[3];
    double scale;
    vtkMath::Norm(norm);
    vtkMath::Subtract(pts,center,pvec);
    scale = vtkMath::Dot(pvec,norm);
    vtkMath::MultiplyScalar(norm,scale);
    disttopoint = sqrt(pow(pts[0]-norm[0],2)+
	pow(pts[1]-norm[1],2)+
	pow(pts[2]-norm[2],2));

    distalonglength = sqrt(pow(norm[0]-center[0],2)+
	pow(norm[1]-center[1],2)+
	pow(norm[2]-center[2],2));

    //set value to new size
    if (disttopoint <= radius && distalonglength <= length/2)
    {
      meshSizeArray->SetValue(pointId,size);
      refineIDArray->SetValue(pointId,refinecount+1);
    }
    else
    {
      if (meshSizeArray->GetValue(pointId) == 0)
        meshSizeArray->SetValue(pointId,maxedgesize);
    }
  }

  if (secondarray)
  {
    polydatasolid->GetPointData()->RemoveArray(sizingFunctionArrayName.c_str());
  }
  polydatasolid->GetPointData()->AddArray(meshSizeArray);
  polydatasolid->GetPointData()->SetActiveScalars(sizingFunctionArrayName.c_str());
  polydatasolid->GetPointData()->AddArray(refineIDArray);

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_SetRefinementSphere()
// -----------------------------
/**
 * @brief computes the distance between each point on surface and center
 * @brief of sphere. Then, if inside radius, the meshsizing function at the
 * @brief is set to the reduced size,
 * @param size This is the smaller refined of the edges within sphere region.
 * @param radius This is the radius of the refinement sphere.
 * @param center This is the center of the refinement sphere.
 * @return SV_OK if function completes properly
 */

int TGenUtils_SetRefinementSphere(vtkPolyData *polydatasolid,
    std::string sizingFunctionArrayName,double size,double radius, double *center,
    int secondarray,double maxedgesize, std::string refineIDArrayName, int refinecount)
{
  int numPts;
  double dist;
  double pts[3];
  vtkIdType pointId, cellId;
  vtkSmartPointer<vtkDoubleArray> meshSizeArray = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> refineIDArray = vtkSmartPointer<vtkIntArray>::New();

  //Set sizing function params
  numPts = polydatasolid->GetNumberOfPoints();
  if (secondarray)
  {
    if (VtkUtils_PDCheckArrayName(polydatasolid,0,sizingFunctionArrayName) != SV_OK)
    {
      fprintf(stderr,"Solid does not contain a double array of name %s. Regions must be identified \
		      Reset or remake the array and try again\n",sizingFunctionArrayName.c_str());
      return SV_ERROR;
    }
    meshSizeArray = vtkDoubleArray::SafeDownCast(polydatasolid->GetPointData()->GetArray(sizingFunctionArrayName.c_str()));
  }
  else
  {
    meshSizeArray->SetNumberOfComponents(1);
    meshSizeArray->Allocate(numPts,1000);
    meshSizeArray->SetNumberOfTuples(numPts);
    meshSizeArray->SetName(sizingFunctionArrayName.c_str());
    for (pointId = 0;pointId<numPts;pointId++)
    {
      meshSizeArray->SetValue(pointId,0.0);
    }
  }
  if (refinecount != 0)
  {
    if (VtkUtils_PDCheckArrayName(polydatasolid,0,refineIDArrayName) != SV_OK)
    {
      fprintf(stderr,"Solid does not contain an int array of name %s. Regions must be identified \
		      Reset or remake the array and try again\n",refineIDArrayName.c_str());
      return SV_ERROR;
    }
    refineIDArray = vtkIntArray::SafeDownCast(polydatasolid->GetPointData()->GetArray(refineIDArrayName.c_str()));
  }
  else
  {
    refineIDArray->SetNumberOfComponents(1);
    refineIDArray->Allocate(numPts,1000);
    refineIDArray->SetNumberOfTuples(numPts);
    refineIDArray->SetName(refineIDArrayName.c_str());
    for (pointId = 0;pointId<numPts;pointId++)
    {
      refineIDArray->SetValue(pointId,0);
    }
  }

  for (pointId = 0;pointId<numPts;pointId++)
  {
    polydatasolid->GetPoint(pointId,pts);
    //compute distance
    dist = sqrt(pow(pts[0]-center[0],2)+
	pow(pts[1]-center[1],2)+
	pow(pts[2]-center[2],2));

    //set value to new size
    if (dist <= radius)
    {
      meshSizeArray->SetValue(pointId,size);
      refineIDArray->SetValue(pointId,refinecount+1);
    }
    else
    {
      if (meshSizeArray->GetValue(pointId) == 0)
        meshSizeArray->SetValue(pointId,maxedgesize);
    }
  }

  if (secondarray)
  {
    polydatasolid->GetPointData()->RemoveArray(sizingFunctionArrayName.c_str());
    polydatasolid->GetPointData()->RemoveArray(refineIDArrayName.c_str());
  }
  polydatasolid->GetPointData()->AddArray(meshSizeArray);
  polydatasolid->GetPointData()->SetActiveScalars(sizingFunctionArrayName.c_str());
  polydatasolid->GetPointData()->AddArray(refineIDArray);

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_SetSizeFunctionArray()
// -----------------------------
/**
 * @brief set a mesh size function based on given array.
 * @brief Values of given array are normalized based on minimum value. Then
 * @brief normalized values are multiplied by size in order to give the mesh
 * @brief size function for the mesher
 * @param size This is the smaller refined of the edges within sphere region.
 * @param sizingFunctionArrayName Name for which to pull values from
 * @param functionname This is the desired function name to be sent to the
 * mesher
 * @param secondarray This designates whether a previous function is already
 * applied.
 * @return SV_OK if function completes properly
 */

int TGenUtils_SetSizeFunctionArray(vtkPolyData *polydatasolid,
    std::string sizingFunctionArrayName,double size,char *functionname,
    int secondarray)
{
  int numPts,numCells;
  double dist;
  double value;
  double factor;
  vtkIdType npts;
  const vtkIdType *pts;
  vtkIdType pointId,cellId;
  double min = 0;
  double range[2];
  vtkSmartPointer<vtkDoubleArray> arrayonmesh = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> regionarray = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkDoubleArray> meshSizeArray = vtkSmartPointer<vtkDoubleArray>::New();

  //Set sizing function params
  numPts = polydatasolid->GetNumberOfPoints();
  numCells = polydatasolid->GetNumberOfCells();
  if (secondarray)
  {
    if (VtkUtils_PDCheckArrayName(polydatasolid,0,sizingFunctionArrayName) != SV_OK)
    {
      fprintf(stderr,"Solid does not contain a double array of name %s. Regions must be identified \
		      Reset or remake the array and try again\n",sizingFunctionArrayName.c_str());
      return SV_ERROR;
    }
    meshSizeArray = vtkDoubleArray::SafeDownCast(polydatasolid->GetPointData()->GetArray(sizingFunctionArrayName.c_str()));
  }
  else
  {
    meshSizeArray->SetNumberOfComponents(1);
    meshSizeArray->Allocate(numPts,1000);
    meshSizeArray->SetNumberOfTuples(numPts);
    meshSizeArray->SetName(sizingFunctionArrayName.c_str());
    for (pointId = 0;pointId<numPts;pointId++)
    {
      meshSizeArray->SetValue(pointId,0.0);
    }
  }

  if (VtkUtils_PDCheckArrayName(polydatasolid,0,functionname) != SV_OK)
  {
    fprintf(stderr,"Solid does not contain a double array of name %s.",
		    functionname);
    return SV_ERROR;
  }

  arrayonmesh = vtkDoubleArray::SafeDownCast(polydatasolid->GetPointData()->GetArray(functionname));

  if (!strncmp(functionname,"DistanceToCenterlines",21))
  {
    arrayonmesh->GetRange(range,0);
    min = range[0];
    fprintf(stderr,"Size Function minimum is: %.4f\n",min);
    fprintf(stderr,"Size Function maximum is: %.4f\n",range[1]);
    if (min <= 0)
    {
      fprintf(stderr,"Min is Zero or negative. This will not work!!!\n",min);
      return SV_ERROR;
    }
    if (min < size)
    {
      std::cout<<"Given mesh size is smaller than minimum radius!!"<<endl;
      std::cout<<"Setting new mesh size to minimum radius :)"<<endl;
      size = min;
    }

    for (pointId = 0;pointId<numPts;pointId++)
    {
      value = arrayonmesh->GetValue(pointId);
      factor = value/min;
  //    fprintf(stderr,"Value is : %.4f\n",factor);
      //compute distance
      //set value to reduced size
      meshSizeArray->SetValue(pointId,factor*size);
    }
    polydatasolid->GetPointData()->RemoveArray(functionname);
  }
  else
  {
  }

  polydatasolid->GetPointData()->AddArray(meshSizeArray);
  polydatasolid->GetPointData()->SetActiveScalars(sizingFunctionArrayName.c_str());


  fprintf(stderr,"Sizing function set\n");
  return SV_OK;
}

// -----------------------------
// cvTGenUtils_LoadMesh()
// -----------------------------
/**
 * @brief Function to load in a vtkUnstructuredGrid
 * @note This is only used by LoadMesh in vtkTetGenMeshObject
 */
//

int TGenUtils_LoadMesh(char *filename,vtkUnstructuredGrid *result)
{
  const char *extension = strrchr(filename,'.');
  extension = extension +1;

  if (!strncmp(extension,"vtu",3)) {
    vtkSmartPointer<vtkXMLUnstructuredGridReader> ugreader =
      vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    ugreader->SetFileName(filename);
    ugreader->Update();

    result->DeepCopy(ugreader->GetOutput());
    result->BuildLinks();
  }
  else {
    fprintf(stderr,"Cannot load the mesh. \
	It must be of type vtkUnstructuredGrid\n");
    return SV_ERROR;
  }

  return SV_OK;
}

int TGenUtils_ResetOriginalRegions(vtkPolyData *newgeom,
    vtkPolyData *originalgeom,
    std::string regionName)
{
  int i,j,k;
  int subId;
  int maxIndex;
  int temp;
  int flag = 1;
  int count;
  int bigcount;
  vtkIdType npts;
  const vtkIdType *pts;
  double distance;
  double closestPt[3];
  double tolerance = 1.0;
  double centroid[3];
  int range;
  vtkIdType closestCell;
  vtkIdType cellId;
  vtkIdType currentValue;
  vtkIdType realValue;
  vtkSmartPointer<vtkCellLocator> locator =
    vtkSmartPointer<vtkCellLocator>::New();
  vtkSmartPointer<vtkGenericCell> genericCell =
    vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkIntArray> currentRegions =
    vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> realRegions =
    vtkSmartPointer<vtkIntArray>::New();

  newgeom->BuildLinks();
  originalgeom->BuildLinks();
  locator->SetDataSet(originalgeom);
  locator->BuildLocator();

  if (VtkUtils_PDCheckArrayName(originalgeom,1,regionName) != SV_OK)
  {
    fprintf(stderr,"Array name 'ModelFaceID' does not exist. Regions must be identified \
		    and named 'ModelFaceID' prior to this function call\n");
    return SV_ERROR;
  }

  realRegions = static_cast<vtkIntArray*>(originalgeom->GetCellData()->GetScalars(regionName.c_str()));


  for (cellId=0;cellId<newgeom->GetNumberOfCells();cellId++)
  {
      newgeom->GetCellPoints(cellId,npts,pts);
      //int eachValue[npts];
      vtkSmartPointer<vtkPoints> polyPts = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkIdTypeArray> polyPtIds = vtkSmartPointer<vtkIdTypeArray>::New();
      for (i=0;i<npts;i++)
      {
	polyPtIds->InsertValue(i,i);
	polyPts->InsertNextPoint(newgeom->GetPoint(pts[i]));
      }
      vtkPolygon::ComputeCentroid(polyPtIds,polyPts,centroid);

      locator->FindClosestPoint(centroid,closestPt,genericCell,closestCell,
	  subId,distance);
      currentRegions->InsertValue(cellId,realRegions->GetValue(closestCell));
  }

  newgeom->GetCellData()->RemoveArray(regionName.c_str());
  currentRegions->SetName(regionName.c_str());
  newgeom->GetCellData()->AddArray(currentRegions);

  newgeom->GetCellData()->SetActiveScalars(regionName.c_str());

  return SV_OK;
}

int TGenUtils_ResetOriginalRegions(vtkPolyData *newgeom,
    vtkPolyData *originalgeom,
    std::string regionName,
    vtkIdList *excludeList)
{
  int i,j,k;
  int subId;
  int region;
  int temp;
  int flag = 1;
  int count;
  int bigcount;
  vtkIdType npts;
  const vtkIdType *pts;
  double distance;
  double closestPt[3];
  double tolerance = 1.0;
  double centroid[3];
  int range;
  vtkIdType closestCell;
  vtkIdType cellId;
  vtkIdType currentValue;
  vtkIdType realValue;
  vtkSmartPointer<vtkCellLocator> locator =
    vtkSmartPointer<vtkCellLocator>::New();
  vtkSmartPointer<vtkGenericCell> genericCell =
    vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkPolyData> originalCopy =
    vtkSmartPointer<vtkPolyData>::New();

  if (excludeList == nullptr)
  {
    fprintf(stderr,"Cannot give nullptr excludeList. Use other reset function without exclude list\n");
    return SV_ERROR;
  }

  newgeom->BuildLinks();
  originalgeom->BuildLinks();
  originalCopy->DeepCopy(originalgeom);

  if (VtkUtils_PDCheckArrayName(originalCopy,1, regionName) != SV_OK)
  {
    fprintf(stderr,"Array name %s does not exist. Regions must be identified \
		    and named 'ModelFaceID' prior to this function call\n",  regionName.c_str());
    return SV_ERROR;
  }

  vtkDataArray *testRegions = originalCopy->GetCellData()->GetScalars( regionName.c_str());

    if (VtkUtils_PDCheckArrayName(newgeom,1, regionName.c_str()) != SV_OK)
    {
      fprintf(stderr,"Array name %s does not exist. Regions must be identified \
          and named 'ModelFaceID' prior to this function call\n", regionName.c_str());
      return SV_ERROR;
    }

    vtkDataArray *currentRegions = newgeom->GetCellData()->GetArray(regionName.c_str());

    for (int i=0; i<originalCopy->GetNumberOfCells(); i++)
    {
      region = testRegions->GetTuple1(i);
      if (excludeList->IsId(region) != -1)
      {
        originalCopy->DeleteCell(i);
      }
    }

    originalCopy->RemoveDeletedCells();

    vtkSmartPointer<vtkCleanPolyData> cleaner =
      vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(originalCopy);
    cleaner->Update();

    originalCopy->DeepCopy(cleaner->GetOutput());
    originalCopy->BuildLinks();

  locator->SetDataSet(originalCopy);
  locator->BuildLocator();
  vtkDataArray *realRegions = originalCopy->GetCellData()->GetScalars( regionName.c_str());

  for (cellId=0;cellId<newgeom->GetNumberOfCells();cellId++)
  {
    currentValue = currentRegions->GetTuple1(cellId);
    if (excludeList->IsId(currentValue) != -1)
    {
      continue;
    }

    newgeom->GetCellPoints(cellId,npts,pts);
    vtkSmartPointer<vtkPoints> polyPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdTypeArray> polyPtIds = vtkSmartPointer<vtkIdTypeArray>::New();
    for (i=0;i<npts;i++)
    {
      polyPtIds->InsertValue(i,i);
      polyPts->InsertNextPoint(newgeom->GetPoint(pts[i]));
    }
    vtkPolygon::ComputeCentroid(polyPtIds,polyPts,centroid);

    locator->FindClosestPoint(centroid,closestPt,genericCell,closestCell,
	subId,distance);
    currentRegions->SetTuple1(cellId,realRegions->GetTuple1(closestCell));
  }

  newgeom->GetCellData()->SetActiveScalars(regionName.c_str());

  return SV_OK;
}

int TGenUtils_ResetOriginalRegions(vtkPolyData *newgeom,
    vtkPolyData *originalgeom,
    std::string regionName,
    vtkIdList *onlyList,
    int dummy)
{
  int i,j,k;
  int subId;
  int region;
  int temp;
  int flag = 1;
  int count;
  int bigcount;
  vtkIdType npts;
  const vtkIdType *pts;
  double distance;
  double closestPt[3];
  double tolerance = 1.0;
  double centroid[3];
  int range;
  vtkIdType closestCell;
  vtkIdType cellId;
  vtkIdType currentValue;
  vtkIdType realValue;
  vtkSmartPointer<vtkCellLocator> locator =
    vtkSmartPointer<vtkCellLocator>::New();
  vtkSmartPointer<vtkGenericCell> genericCell =
    vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkPolyData> originalCopy =
    vtkSmartPointer<vtkPolyData>::New();

  if (onlyList == nullptr)
  {
    fprintf(stderr,"Cannot give nullptr onlyList. Use other reset function without only list\n");
    return SV_ERROR;
  }

  newgeom->BuildLinks();
  originalgeom->BuildLinks();
  originalCopy->DeepCopy(originalgeom);

  if (VtkUtils_PDCheckArrayName(originalCopy,1, regionName) != SV_OK)
  {
    fprintf(stderr,"Array name %s does not exist. Regions must be identified \
		    and named 'ModelFaceID' prior to this function call\n",  regionName.c_str());
    return SV_ERROR;
  }

  vtkDataArray *testRegions = originalCopy->GetCellData()->GetScalars( regionName.c_str());

  if (VtkUtils_PDCheckArrayName(newgeom,1, regionName.c_str()) != SV_OK)
  {
    fprintf(stderr,"Array name %s does not exist. Regions must be identified \
        and named 'ModelFaceID' prior to this function call\n", regionName.c_str());
    return SV_ERROR;
  }

  vtkDataArray *currentRegions = newgeom->GetCellData()->GetArray(regionName.c_str());

  locator->SetDataSet(originalCopy);
  locator->BuildLocator();
  vtkDataArray *realRegions = originalCopy->GetCellData()->GetScalars( regionName.c_str());

  for (cellId=0;cellId<newgeom->GetNumberOfCells();cellId++)
  {
    currentValue = currentRegions->GetTuple1(cellId);
    if (onlyList->IsId(currentValue) == -1)
    {
      continue;
    }

    newgeom->GetCellPoints(cellId,npts,pts);
    vtkSmartPointer<vtkPoints> polyPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdTypeArray> polyPtIds = vtkSmartPointer<vtkIdTypeArray>::New();
    for (i=0;i<npts;i++)
    {
      polyPtIds->InsertValue(i,i);
      polyPts->InsertNextPoint(newgeom->GetPoint(pts[i]));
    }
    vtkPolygon::ComputeCentroid(polyPtIds,polyPts,centroid);

    locator->FindClosestPoint(centroid,closestPt,genericCell,closestCell,
	subId,distance);
    currentRegions->SetTuple1(cellId,realRegions->GetTuple1(closestCell));
  }

  newgeom->GetCellData()->SetActiveScalars(regionName.c_str());

  return SV_OK;
}

// -----------------------------
// cvTGenUtils_CheckSurfaceMesh()
// -----------------------------
/**
 * @brief Function to load in a vtkUnstructuredGrid
 * @note This is only used by LoadMesh in vtkTetGenMeshObject
 */
//

int TGenUtils_CheckSurfaceMesh(vtkPolyData *pd, int meshInfo[3])
{
  fprintf(stdout,"Checking surface mesh\n");
  vtkIdType npts,p0,p1;
  const vtkIdType *pts;
  int NonManifoldEdges = 0,FreeEdges = 0;
  int Regions=0;
  vtkSmartPointer<vtkCleanPolyData> cleaner =
    vtkSmartPointer<vtkCleanPolyData>::New();
  vtkSmartPointer<vtkIdList> edgeneigh =
    vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkConnectivityFilter> connector =
    vtkSmartPointer<vtkConnectivityFilter>::New();
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfacer =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();

  //Clean the input surface
  cleaner->SetInputData(pd);
  cleaner->Update();
  pd->DeepCopy(cleaner->GetOutput());
  pd->BuildLinks();

  surfacer->SetInputData(cleaner->GetOutput());
  surfacer->Update();

  connector->SetInputData(surfacer->GetOutput());
  connector->ColorRegionsOn();
  connector->Update();

  vtkDataArray *regionarray = connector->GetOutput()->
      GetCellData()->GetScalars("RegionId");

  //Loop through the surface and find edges with cells that have either more
  //than one neighbor or no neighbors. No neighbors can be okay,as this can
  //indicate a free edge. Multiple neighbors indicates a
  //non-manifold edge. This can cause issues as well in certain cases.
  for (int i = 0;i<pd->GetNumberOfCells();i++)
  {
    pd->GetCellPoints(i,npts,pts);
    for (int j=0;j<npts;j++)
    {
      p0 = pts[j];
      p1 = pts[(j+1)%npts];

      pd->GetCellEdgeNeighbors(i,p0,p1,edgeneigh);
      if (edgeneigh->GetNumberOfIds() > 1)
        NonManifoldEdges++;
      else if (edgeneigh->GetNumberOfIds() < 1)
        FreeEdges++;
    }
    int val = regionarray->GetTuple1(i);
    if (val > Regions)
      Regions = val;
  }

  fprintf(stdout,"Regions: %d\n", Regions + 1);
  fprintf(stdout,"Number of Free Edges on Surface: %d\n", FreeEdges);
  fprintf(stdout,"Number of Non-Manifold Edges on Surface: %d\n", NonManifoldEdges);

  meshInfo[0] = Regions+1;
  meshInfo[1] = FreeEdges;
  meshInfo[2] = NonManifoldEdges;

  return SV_OK;
}

int TGenUtils_SetLocalMeshSize(vtkPolyData *pd,int regionId,double size)
{
  vtkIdType pointId, cellId;
  vtkIdType npts;
  const vtkIdType *pts;
  vtkSmartPointer<vtkIntArray> regionarray =
    vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkDoubleArray> meshSizeArray =
    vtkSmartPointer<vtkDoubleArray>::New();

  int numPts = pd->GetNumberOfPoints();
  int numCells = pd->GetNumberOfCells();
  regionarray = vtkIntArray::SafeDownCast(pd->GetCellData()->GetArray("ModelFaceID"));
  if (VtkUtils_PDCheckArrayName(pd,0,"MeshSizingFunction") != SV_OK)
  {
    meshSizeArray->SetNumberOfComponents(1);
    meshSizeArray->Allocate(numPts,1000);
    meshSizeArray->SetNumberOfTuples(numPts);
    meshSizeArray->SetName("MeshSizingFunction");
    for (pointId = 0;pointId<numPts;pointId++)
    {
      meshSizeArray->SetValue(pointId,0.0);
    }
  }
  else
  {
    meshSizeArray = vtkDoubleArray::SafeDownCast(pd->GetPointData()->GetArray("MeshSizingFunction"));
  }
  pd->BuildLinks();
  for (cellId = 0;cellId<numCells;cellId++)
  {
    if (regionarray->GetValue(cellId) == regionId)
    {
      pd->GetCellPoints(cellId,npts,pts);
      for (int j=0;j<npts;j++)
      {
	meshSizeArray->SetValue(pts[j],size);
      }
    }
  }

  pd->GetPointData()->RemoveArray("MeshSizingFunction");
  meshSizeArray->SetName("MeshSizingFunction");
  pd->GetPointData()->AddArray(meshSizeArray);

  return SV_OK;
}

// -----------------------------
// TGenUtils_ReportMeshQuality
// -----------------------------
/**
 * @brief Computes the quality of a volume mesh and prints a summary.
 * @note The aspect ratio (the ratio of the longest edge length to the
 * shortest tetrahedron height, normalized so that 1.0 is an equilateral
 * tetrahedron) is computed for each tetrahedral element and stored in
 * a cell data array named 'AspectRatio' so it can be visualized.
 * @note The commonly used rule of thumb thresholds are used to assess
 * the mesh: elements with an aspect ratio above 3 are distorted and
 * elements above 10 are considered poor.
 * @param mesh The volume mesh to check.
 * @return SV_OK if the quality is computed.
 */

int TGenUtils_ReportMeshQuality(vtkUnstructuredGrid *mesh)
{
  if (mesh == nullptr || mesh->GetNumberOfCells() == 0)
  {
    fprintf(stderr,"Cannot compute the quality of an empty mesh\n");
    return SV_ERROR;
  }

  auto qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
  qualityFilter->SetInputData(mesh);
  qualityFilter->SetTetQualityMeasureToAspectRatio();
  qualityFilter->Update();

  auto qualityArray = vtkDoubleArray::SafeDownCast(
    qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
  if (qualityArray == nullptr)
  {
    fprintf(stderr,"Could not compute the mesh quality\n");
    return SV_ERROR;
  }

  // Compute the aspect ratio statistics over the tetrahedral elements
  // and store the aspect ratio in an 'AspectRatio' cell data array.
  //
  int numTets = 0;
  int numDistorted = 0;
  int numPoor = 0;
  double minRatio = 0.0;
  double maxRatio = 0.0;
  double sumRatio = 0.0;

  // Track the worst few tetrahedra (highest aspect ratio) so their location
  // is reported; this tells whether the worst element is at a junction or at
  // a known input-surface sliver.
  const int numWorstToReport = 5;
  std::vector<std::pair<double,vtkIdType>> worstTets;

  auto aspectRatio = vtkSmartPointer<vtkDoubleArray>::New();
  aspectRatio->SetNumberOfComponents(1);
  aspectRatio->SetNumberOfTuples(mesh->GetNumberOfCells());
  aspectRatio->FillComponent(0, 0.0);
  aspectRatio->SetName("AspectRatio");

  for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); cellId++)
  {
    if (mesh->GetCellType(cellId) != VTK_TETRA)
    {
      continue;
    }
    double ratio = qualityArray->GetValue(cellId);
    aspectRatio->SetValue(cellId, ratio);

    if (numTets == 0 || ratio < minRatio)
    {
      minRatio = ratio;
    }
    if (numTets == 0 || ratio > maxRatio)
    {
      maxRatio = ratio;
    }
    sumRatio += ratio;
    numTets++;

    if (ratio > 10.0)
    {
      numPoor++;
    }
    else if (ratio > 3.0)
    {
      numDistorted++;
    }

    // Keep the top 'numWorstToReport' tetrahedra sorted by descending aspect
    // ratio.
    if ((int)worstTets.size() < numWorstToReport || ratio > worstTets.back().first)
    {
      auto pos = std::lower_bound(worstTets.begin(), worstTets.end(), ratio,
          [](const std::pair<double,vtkIdType>& entry, double value)
          { return entry.first > value; });
      worstTets.insert(pos, std::make_pair(ratio, cellId));
      if ((int)worstTets.size() > numWorstToReport)
      {
        worstTets.pop_back();
      }
    }
  }

  mesh->GetCellData()->RemoveArray("AspectRatio");
  mesh->GetCellData()->AddArray(aspectRatio);

  if (numTets == 0)
  {
    fprintf(stderr,"No tetrahedral elements found to compute the mesh quality\n");
    return SV_ERROR;
  }

  double avgRatio = sumRatio / numTets;
  double pctDistorted = 100.0 * numDistorted / numTets;
  double pctPoor = 100.0 * numPoor / numTets;

  fprintf(stdout,"Mesh quality (aspect ratio, 1.0 is an equilateral tetrahedron):\n");
  fprintf(stdout,"  Number of elements: %d\n", numTets);
  fprintf(stdout,"  Min / Avg / Max: %.3f / %.3f / %.3f\n", minRatio, avgRatio, maxRatio);
  fprintf(stdout,"  Elements with aspect ratio > 3 (distorted): %d (%.2f%%)\n", numDistorted, pctDistorted);
  fprintf(stdout,"  Elements with aspect ratio > 10 (poor): %d (%.2f%%)\n", numPoor, pctPoor);

  if (numPoor == 0 && pctDistorted < 5.0)
  {
    fprintf(stdout,"  Mesh quality assessment: GOOD\n");
  }
  else if (pctPoor < 1.0)
  {
    fprintf(stdout,"  Mesh quality assessment: ACCEPTABLE (some distorted elements)\n");
  }
  else
  {
    fprintf(stdout,"  Mesh quality assessment: POOR (consider adjusting the mesh size options or remeshing)\n");
  }

  // Report the location (centroid) of the worst elements so a high aspect
  // ratio can be traced to a junction or to a known input-surface sliver.
  auto worstPtIds = vtkSmartPointer<vtkIdList>::New();
  for (auto& worst : worstTets)
  {
    mesh->GetCellPoints(worst.second, worstPtIds);
    double centroid[3] = {0.0, 0.0, 0.0};
    vtkIdType npts = worstPtIds->GetNumberOfIds();
    for (vtkIdType k = 0; k < npts; k++)
    {
      double point[3];
      mesh->GetPoint(worstPtIds->GetId(k), point);
      centroid[0] += point[0];
      centroid[1] += point[1];
      centroid[2] += point[2];
    }
    if (npts > 0)
    {
      centroid[0] /= npts;
      centroid[1] /= npts;
      centroid[2] /= npts;
    }
    fprintf(stdout,"  aspect ratio %.3f at element centroid (%.5g, %.5g, %.5g)\n",
        worst.first, centroid[0], centroid[1], centroid[2]);
  }

  return SV_OK;
}

// -----------------------------
// TGenUtils_SmoothPointArray
// -----------------------------
/**
 * @brief Smooths a point data array over a surface by iterative Laplacian
 * averaging.
 * @note Each iteration replaces every point value with the average of its
 * own value and the values of its one-ring neighbors (the points sharing a
 * cell with it). All new values are computed from the previous iteration's
 * values (Jacobi iteration) so the result does not depend on the order in
 * which the points are visited. The surface geometry is not modified; only
 * the array values change.
 * @param surface The surface whose cells define the point neighbors.
 * @param array The point data array to smooth; must have exactly one
 * component and one tuple per surface point.
 * @param iterations The number of averaging iterations; a value less than
 * one leaves the array unchanged.
 * @return SV_OK if the array is smoothed.
 */

int TGenUtils_SmoothPointArray(vtkPolyData *surface, vtkDoubleArray *array, int iterations)
{
  if (iterations <= 0)
  {
    return SV_OK;
  }

  if (surface == nullptr || array == nullptr)
  {
    fprintf(stderr,"Cannot smooth a point array without a surface and an array\n");
    return SV_ERROR;
  }

  if (array->GetNumberOfComponents() != 1)
  {
    fprintf(stderr,"The point array to smooth must have exactly one component\n");
    return SV_ERROR;
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The point array size does not match the number of surface points\n");
    return SV_ERROR;
  }

  // Build the one-ring point neighbors of each point once. The points
  // sharing a cell with a point are its neighbors.
  //
  surface->BuildLinks();
  std::vector<std::vector<vtkIdType>> neighbors(numPts);
  auto cellIds = vtkSmartPointer<vtkIdList>::New();
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    surface->GetPointCells(ptId, cellIds);
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellIds->GetId(i), npts, pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (pts[j] == ptId)
        {
          continue;
        }
        auto& ptNeighbors = neighbors[ptId];
        if (std::find(ptNeighbors.begin(), ptNeighbors.end(), pts[j]) == ptNeighbors.end())
        {
          ptNeighbors.push_back(pts[j]);
        }
      }
    }
  }

  std::vector<double> values(numPts);
  std::vector<double> smoothed(numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    values[ptId] = array->GetValue(ptId);
  }

  for (int iter = 0; iter < iterations; iter++)
  {
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      double sum = values[ptId];
      int count = 1;
      for (auto neighborId : neighbors[ptId])
      {
        sum += values[neighborId];
        count++;
      }
      smoothed[ptId] = sum / count;
    }
    values.swap(smoothed);
  }

  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    array->SetValue(ptId, values[ptId]);
  }

  return SV_OK;
}

// -------------------------------------------------
// TGenUtils_SmoothWarpVectorsInConcaveRegions
// -------------------------------------------------
/**
 * @brief Smooths the extrusion warp vectors (the point normals) of a surface
 * in its concave regions so the outward wall extrusion does not dip inward
 * and skew where the normals converge.
 * @note When a surface is extruded outward along its point normals, the warp
 * vectors of a concave region (such as the crotch where two vessels merge)
 * converge toward each other, so the extruded outer wall dips inward and the
 * wall elements there are twisted even when the wall does not fully fold over.
 * The wall thickness is not the cause and reducing it does not fix the twist;
 * the direction field is. This relaxes each concave point's normal toward the
 * average of its one-ring neighbors' normals so the converging directions
 * spread apart, and renormalizes it. Only the normal direction changes: the
 * wall thickness (taken from a separate array) and the surface points (the
 * fluid/wall interface) never move.
 *
 * Each point is relaxed in proportion to how concave it is. The concavity is
 * the average, over the one-ring neighbors that rise above the point's
 * tangent plane, of the sine of their rise angle (a dimensionless value that
 * is zero on convex and flat regions), so convex and flat regions keep their
 * normals and a straight tube is left unchanged. Points on a boundary edge
 * (the cap rims, whose normals are set to lie in the cap plane) are pinned
 * and never smoothed, so the wall stays flat at the caps.
 * @param surface The surface being extruded; must have a 3-component normals
 * point data array with the outward point normals.
 * @param normalsArrayName The name of the normals point data array to smooth.
 * @param iterations The number of relaxation iterations; a value less than
 * one leaves the array unchanged.
 * @param maxRelaxation The largest fraction of the neighbor average blended
 * into a fully concave point's normal per iteration (between 0 and 1); a
 * value of zero or less leaves the array unchanged.
 * @return SV_OK if the normals are smoothed.
 */

int TGenUtils_SmoothWarpVectorsInConcaveRegions(vtkPolyData *surface, const char *normalsArrayName,
    int iterations, double maxRelaxation)
{
  if (iterations <= 0 || maxRelaxation <= 0.0)
  {
    return SV_OK;
  }

  if (surface == nullptr || normalsArrayName == nullptr)
  {
    fprintf(stderr,"Cannot smooth the warp vectors without a surface and an array name\n");
    return SV_ERROR;
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  auto normals = surface->GetPointData()->GetArray(normalsArrayName);
  if (normals == nullptr || normals->GetNumberOfComponents() != 3 ||
      normals->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The surface must have a 3-component '%s' point array to smooth the warp vectors\n",
        normalsArrayName);
    return SV_ERROR;
  }

  // Read the current normals into a flat working buffer (component k of point
  // ptId is at index 3*ptId+k) so the array's own storage type does not
  // matter.
  std::vector<double> vectors(3*numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    normals->GetTuple(ptId, &vectors[3*ptId]);
  }

  // Build the one-ring point neighbors of each point and, at the same time,
  // detect the points on a boundary edge (an edge used by a single cell, such
  // as the cap rims). Boundary points are pinned so the cap normals set to
  // lie in the cap plane are kept and the wall stays flat at the caps.
  surface->BuildLinks();
  std::vector<std::vector<vtkIdType>> neighbors(numPts);
  std::vector<char> pinned(numPts, 0);
  auto cellIds = vtkSmartPointer<vtkIdList>::New();
  auto edgeNeighbors = vtkSmartPointer<vtkIdList>::New();
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    surface->GetPointCells(ptId, cellIds);
    auto& ptNeighbors = neighbors[ptId];
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
      vtkIdType cellId = cellIds->GetId(i);
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId, npts, pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (pts[j] == ptId)
        {
          continue;
        }
        if (std::find(ptNeighbors.begin(), ptNeighbors.end(), pts[j]) == ptNeighbors.end())
        {
          ptNeighbors.push_back(pts[j]);
        }
        // The edge (ptId, pts[j]) is a boundary edge when no other cell
        // shares it, which makes both its endpoints boundary points.
        surface->GetCellEdgeNeighbors(cellId, ptId, pts[j], edgeNeighbors);
        if (edgeNeighbors->GetNumberOfIds() == 0)
        {
          pinned[ptId] = 1;
        }
      }
    }
  }

  // Precompute a per-point relaxation weight from the initial geometry and
  // normals. The concavity is the average sine of the rise angle over the
  // neighbors above the tangent plane and is zero on convex and flat points,
  // so only concave points are relaxed.
  std::vector<double> weight(numPts, 0.0);
  vtkIdType numConcavePts = 0;
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    if (pinned[ptId] || neighbors[ptId].empty())
    {
      continue;
    }
    const double *normal = &vectors[3*ptId];
    double normalLength = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
        normal[2]*normal[2]);
    if (normalLength <= 0.0)
    {
      continue;
    }
    double point[3];
    surface->GetPoint(ptId, point);
    double concavitySum = 0.0;
    int concaveCount = 0;
    for (auto neighborId : neighbors[ptId])
    {
      double neighbor[3];
      surface->GetPoint(neighborId, neighbor);
      double offset[3] = {neighbor[0]-point[0], neighbor[1]-point[1], neighbor[2]-point[2]};
      double distance = std::sqrt(offset[0]*offset[0] + offset[1]*offset[1] +
          offset[2]*offset[2]);
      if (distance <= 0.0)
      {
        continue;
      }
      double height = (offset[0]*normal[0] + offset[1]*normal[1] +
          offset[2]*normal[2]) / normalLength;
      if (height <= 0.0)
      {
        continue;
      }
      concavitySum += height/distance;   // sine of the rise angle, in [0,1)
      concaveCount++;
    }
    if (concaveCount > 0)
    {
      weight[ptId] = maxRelaxation * (concavitySum/concaveCount);
      numConcavePts++;
    }
  }

  // Relax each concave point's normal toward its neighbor average (Jacobi
  // iteration, so the result is independent of the point visiting order) and
  // renormalize; pinned and convex/flat points keep their normals. Track the
  // largest direction change so the effect is observable in the log.
  double maxAngleChange = 0.0;
  std::vector<double> smoothed(3*numPts);
  for (int iter = 0; iter < iterations; iter++)
  {
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      smoothed[3*ptId] = vectors[3*ptId];
      smoothed[3*ptId+1] = vectors[3*ptId+1];
      smoothed[3*ptId+2] = vectors[3*ptId+2];
      double w = weight[ptId];
      if (w <= 0.0 || neighbors[ptId].empty())
      {
        continue;
      }
      double average[3] = {0.0, 0.0, 0.0};
      for (auto neighborId : neighbors[ptId])
      {
        average[0] += vectors[3*neighborId];
        average[1] += vectors[3*neighborId+1];
        average[2] += vectors[3*neighborId+2];
      }
      double count = (double)neighbors[ptId].size();
      average[0] /= count; average[1] /= count; average[2] /= count;
      double blended[3] = {
        (1.0-w)*vectors[3*ptId]   + w*average[0],
        (1.0-w)*vectors[3*ptId+1] + w*average[1],
        (1.0-w)*vectors[3*ptId+2] + w*average[2]};
      double length = std::sqrt(blended[0]*blended[0] + blended[1]*blended[1] +
          blended[2]*blended[2]);
      if (length <= 0.0)
      {
        continue;   // degenerate average; keep the current normal
      }
      smoothed[3*ptId]   = blended[0]/length;
      smoothed[3*ptId+1] = blended[1]/length;
      smoothed[3*ptId+2] = blended[2]/length;
    }
    vectors.swap(smoothed);
  }

  // Write the smoothed normals back and report the largest direction change,
  // measured against the original normals.
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double original[3];
    normals->GetTuple(ptId, original);
    double dot = original[0]*vectors[3*ptId] + original[1]*vectors[3*ptId+1] +
        original[2]*vectors[3*ptId+2];
    if (dot > 1.0) { dot = 1.0; }
    if (dot < -1.0) { dot = -1.0; }
    double angle = std::acos(dot);
    if (angle > maxAngleChange)
    {
      maxAngleChange = angle;
    }
    normals->SetTuple(ptId, &vectors[3*ptId]);
  }

  std::cout << "Smoothed the wall extrusion warp vectors at " << numConcavePts
      << " concave points (max direction change " << maxAngleChange*180.0/M_PI
      << " degrees)" << std::endl;

  return SV_OK;
}

// --------------------------------------------
// TGenUtils_ClampThicknessToConcaveCurvature
// --------------------------------------------
/**
 * @brief Limits a wall thickness point array so the thickness nowhere
 * exceeds a fraction of the local concave radius of curvature of the
 * surface.
 * @note When a surface is extruded outward along its point normals, the
 * warp vectors of a concave region (such as the crotch where two vessels
 * merge) converge; if the thickness is larger than the concave radius of
 * curvature the extruded outer wall folds over and self-intersects. For
 * each point the concave curvature is estimated from its one-ring
 * neighbors: a neighbor at distance d that rises a height h above the
 * point's tangent plane (measured along the outward normal) implies a
 * curvature of about 2*h/d^2, and the largest such value over the
 * neighbors is used. The thickness at the point is then clamped to
 * factor divided by that curvature; the limit is always enforced so the
 * clamp is idempotent, and a warning is printed when it reduces a
 * thickness below a small fraction of its requested value because the
 * resulting wall elements may be very thin there. The one-ring estimate
 * is local, so the clamp reduces but does not guarantee the absence of
 * global self-intersections of the extruded outer wall.
 * Convex and flat regions (h <= 0 for all neighbors) are left unchanged.
 * Only the thickness values change; the surface points (the fluid/wall
 * interface) never move.
 * @param surface The surface being extruded; must have a 3-component
 * 'Normals' point data array with the outward point normals.
 * @param array The wall thickness point array to clamp; must have exactly
 * one component and one tuple per surface point.
 * @param factor The maximum allowed thickness as a fraction of the concave
 * radius of curvature; a value of zero or less leaves the array unchanged.
 * @return SV_OK if the array is clamped.
 */

int TGenUtils_ClampThicknessToConcaveCurvature(vtkPolyData *surface, vtkDoubleArray *array, double factor)
{
  if (factor <= 0.0)
  {
    return SV_OK;
  }

  if (surface == nullptr || array == nullptr)
  {
    fprintf(stderr,"Cannot clamp a thickness array without a surface and an array\n");
    return SV_ERROR;
  }

  if (array->GetNumberOfComponents() != 1)
  {
    fprintf(stderr,"The thickness array to clamp must have exactly one component\n");
    return SV_ERROR;
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array size does not match the number of surface points\n");
    return SV_ERROR;
  }

  auto normals = surface->GetPointData()->GetArray("Normals");
  if (normals == nullptr || normals->GetNumberOfComponents() != 3 ||
      normals->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The surface must have a 3-component 'Normals' point array to clamp the thickness\n");
    return SV_ERROR;
  }

  // The curvature limit is always enforced; when it reduces a thickness
  // below this fraction of its requested value the wall elements there
  // may be very thin, so the reduction is reported as a warning.
  const double thinThicknessRatio = 0.1;
  vtkIdType numThinPts = 0;

  surface->BuildLinks();
  auto cellIds = vtkSmartPointer<vtkIdList>::New();
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double normal[3];
    normals->GetTuple(ptId, normal);
    double normalLength = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
        normal[2]*normal[2]);
    if (normalLength <= 0.0)
    {
      continue;
    }

    double point[3];
    surface->GetPoint(ptId, point);

    // Estimate the largest concave curvature at the point from its
    // one-ring neighbors. Visiting a neighbor more than once is harmless
    // because only the maximum is kept.
    double maxCurvature = 0.0;
    surface->GetPointCells(ptId, cellIds);
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellIds->GetId(i), npts, pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (pts[j] == ptId)
        {
          continue;
        }
        double neighbor[3];
        surface->GetPoint(pts[j], neighbor);
        double offset[3] = {neighbor[0]-point[0], neighbor[1]-point[1],
            neighbor[2]-point[2]};
        double distanceSquared = offset[0]*offset[0] + offset[1]*offset[1] +
            offset[2]*offset[2];
        if (distanceSquared <= 0.0)
        {
          continue;
        }
        double height = (offset[0]*normal[0] + offset[1]*normal[1] +
            offset[2]*normal[2]) / normalLength;
        if (height <= 0.0)
        {
          continue;
        }
        double curvature = 2.0*height/distanceSquared;
        if (curvature > maxCurvature)
        {
          maxCurvature = curvature;
        }
      }
    }

    if (maxCurvature <= 0.0)
    {
      continue;
    }

    double thickness = array->GetValue(ptId);
    double maxThickness = factor/maxCurvature;
    if (thickness > maxThickness)
    {
      array->SetValue(ptId, maxThickness);
      if (maxThickness < thinThicknessRatio*thickness)
      {
        numThinPts++;
      }
    }
  }

  if (numThinPts > 0)
  {
    fprintf(stderr,"Warning: the curvature limit reduced the wall thickness below %g%% of its\
 requested value at %lld points; the wall elements there may be very thin\n",
        100.0*thinThicknessRatio, (long long)numThinPts);
  }

  return SV_OK;
}

// -----------------------------------------
// TGenUtils_ClusterPointsIntoRegions
// -----------------------------------------
/**
 * @brief Groups flagged surface points into spatially separated regions,
 * worst first.
 * @note A diagnostic that reports only its worst point is always pulled to the
 * same single severe defect, which hides every other junction, so the log
 * cannot distinguish "one bad spot" from "every junction is affected" - the
 * question that decides whether a fix belongs in the input surface or in the
 * thickness passes. The worst remaining point seeds a region, every flagged
 * point within a radius of it is absorbed, and the process repeats. Only the
 * first few seeds are wanted, so the cost is maxRegions passes over the flagged
 * points rather than a full clustering. Nothing on the surface is modified.
 * @param surface The surface the flagged points belong to.
 * @param points The flagged points as (sort key, point id), sorted in place.
 * The worst point is the one with the smallest key, so a caller whose worst
 * value is its largest passes the negated value as the key.
 * @param maxRegions The maximum number of regions to return.
 * @param radiusFraction The region radius as a fraction of the diagonal of the
 * surface bounding box.
 * @param regions The seed and size of each region found, worst first.
 * @param radius The absolute region radius used, for the caller to report.
 * @param numOutside The flagged points left outside the returned regions, which
 * is non-zero only when maxRegions regions were filled.
 * @param numRegionsTotal Every region the flagged points form, not only the
 * maxRegions reported. Reporting a truncated list cannot distinguish "these
 * few junctions are affected" from "every junction is affected and only the
 * worst few are shown", which is the question the caller is asking, so the
 * clustering always runs to exhaustion and only the reported list is capped.
 * The extra passes cost nothing because they walk the flagged points, not the
 * surface.
 * @return SV_OK if the points are clustered.
 */

int TGenUtils_ClusterPointsIntoRegions(vtkPolyData *surface,
    std::vector<std::pair<double,vtkIdType> > &points, int maxRegions,
    double radiusFraction, std::vector<TGenUtilsPointRegion> &regions,
    double &radius, int &numOutside, int &numRegionsTotal)
{
  regions.clear();
  radius = 0.0;
  numOutside = 0;
  numRegionsTotal = 0;

  if (surface == nullptr)
  {
    fprintf(stderr,"Cannot cluster surface points into regions without a surface\n");
    return SV_ERROR;
  }

  if (points.empty() || maxRegions <= 0)
  {
    return SV_OK;
  }

  std::sort(points.begin(), points.end());

  double bounds[6];
  surface->GetBounds(bounds);
  double dx = bounds[1]-bounds[0], dy = bounds[3]-bounds[2], dz = bounds[5]-bounds[4];
  radius = radiusFraction*std::sqrt(dx*dx + dy*dy + dz*dz);
  double radius2 = radius*radius;

  // The clustering runs over every flagged point so numRegionsTotal counts all
  // of them; only the first maxRegions are kept for reporting, and the points
  // absorbed by the regions past that cap are counted as "outside" so the
  // reported list plus numOutside still accounts for every flagged point.
  std::vector<bool> absorbed(points.size(), false);
  for (size_t i = 0; i < points.size(); i++)
  {
    if (absorbed[i])
    {
      continue;
    }
    vtkIdType seedId = points[i].second;
    double seed[3];
    surface->GetPoint(seedId, seed);
    absorbed[i] = true;
    int members = 1;
    for (size_t j = i+1; j < points.size(); j++)
    {
      if (absorbed[j])
      {
        continue;
      }
      double p[3];
      surface->GetPoint(points[j].second, p);
      double d2 = (p[0]-seed[0])*(p[0]-seed[0]) + (p[1]-seed[1])*(p[1]-seed[1])
                + (p[2]-seed[2])*(p[2]-seed[2]);
      if (d2 < radius2)
      {
        absorbed[j] = true;
        members++;
      }
    }
    numRegionsTotal++;
    if ((int)regions.size() < maxRegions)
    {
      TGenUtilsPointRegion region;
      region.seedId = seedId;
      region.numPoints = members;
      regions.push_back(region);
    }
    else
    {
      numOutside += members;
    }
  }

  return SV_OK;
}

// -------------------------------------------------
// TGenUtils_ReportConcaveCurvatureVsThickness
// -------------------------------------------------
/**
 * @brief Reports the requested wall thickness against the local concave radius
 * of curvature of the surface it is extruded from.
 * @note Extruding a surface outward along its point normals makes the warp
 * vectors of a concave region converge, so the outer wall self-intersects once
 * the thickness exceeds the concave radius of curvature there. Whether that
 * happens is decided by the single ratio t/R, and the thickness passes only
 * ever observe it indirectly: the curvature clamp is skipped entirely when its
 * factor is zero, and the fold-prevention pass reports the thickness it had to
 * remove rather than the ratio that forced it. Reporting t/R directly separates
 * a junction whose shape cannot carry the requested thickness (t/R > 1, where
 * no offset-based method can avoid the self-intersection) from a thinning the
 * thickness passes produce on a junction that could have carried it.
 *
 * The curvature at a point is estimated from its one-ring as in the curvature
 * clamp: a neighbor at distance d rising a height h above the point's tangent
 * plane implies a curvature of about 2*h/d^2. Two summaries of that one-ring
 * are reported. The maximum matches the clamp and gives the smallest radius
 * present, but it has d^2 in the denominator, so a single near-degenerate
 * triangle makes it diverge and it cannot tell a sharp junction from a sliver.
 * The median over the whole one-ring, counting the neighbors at or below the
 * tangent plane as zero curvature, is insensitive to one bad neighbor and
 * describes the shape of the junction: it is only large where most of the
 * one-ring is genuinely concave. A region whose median radius is large while
 * its smallest radius is tiny is a mesh artifact rather than a sharp shape, so
 * the points using a near-degenerate triangle are flagged as well.
 *
 * This is a report; neither the thickness array nor the surface is modified. It
 * is meant to be called after the warp vectors are smoothed, because the
 * smoothing changes the normals the heights are measured along, and before any
 * thickness reduction, so the thickness it reports is the requested one.
 * @param surface The surface being extruded; must have a 3-component 'Normals'
 * point data array with the outward point normals. The per-point ratio and
 * radii are left on it as the 'ThicknessOverRadius', 'ConcaveRadiusTypical'
 * and 'ConcaveRadiusSmallest' point arrays so the field can be viewed directly
 * rather than only through the region summary in the log.
 * @param array The wall thickness point array; one component, one tuple per
 * surface point.
 * @param inward Zero when the extrusion goes outward along the normals (the
 * solid wall), non-zero when it goes inward (the fluid boundary layer). Only
 * the side the extrusion converges on is at risk, and it is the opposite side
 * for the two, so this flips which neighbors count as concave. Without it the
 * fluid boundary layer cannot be measured on the same footing as the wall, and
 * the two being extruded from the same surface in opposite directions is
 * exactly what has to be compared.
 * @param label Names the extrusion the report belongs to, so the wall and the
 * fluid boundary layer are distinguishable in one log.
 * @return SV_OK if the surface is reported on.
 */

int TGenUtils_ReportConcaveCurvatureVsThickness(vtkPolyData *surface, vtkDoubleArray *array,
    int inward, const char *label)
{
  if (surface == nullptr || array == nullptr)
  {
    fprintf(stderr,"Cannot report the concave curvature without a surface and a thickness array\n");
    return SV_ERROR;
  }

  if (label == nullptr)
  {
    label = "wall";
  }
  // The extrusion converges on the side it moves toward, so the heights that
  // decide concavity are measured along the extrusion direction, not along the
  // stored outward normal.
  const double normalSign = inward ? -1.0 : 1.0;

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array must have one component and one tuple per surface point\n");
    return SV_ERROR;
  }

  auto normals = surface->GetPointData()->GetArray("Normals");
  if (normals == nullptr || normals->GetNumberOfComponents() != 3 ||
      normals->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The surface must have a 3-component 'Normals' point array to report the concave curvature\n");
    return SV_ERROR;
  }

  // A triangle whose smallest altitude has collapsed against its longest edge
  // cannot carry a wall at all, and it also drives the one-ring curvature
  // estimate to a meaningless value, so the points using one are flagged and
  // read separately from the junctions.
  const double sliverAltitudeRatio = 0.05;

  surface->BuildLinks();

  std::vector<double> ratioTypical(numPts, 0.0);
  std::vector<double> radiusTypical(numPts, 0.0);
  std::vector<double> radiusSmallest(numPts, 0.0);
  std::vector<char> sliverAdjacent(numPts, 0);

  std::vector<std::pair<double,vtkIdType> > flagged;
  int numAtLeast1 = 0, numAtLeast2 = 0, numAtLeast4 = 0, numConcave = 0;
  int numSliverAdjacent = 0;

  auto cellIds = vtkSmartPointer<vtkIdList>::New();
  std::vector<vtkIdType> neighbors;
  std::vector<double> curvatures;

  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double normal[3];
    normals->GetTuple(ptId, normal);
    double normalLength = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
        normal[2]*normal[2]);
    if (normalLength <= 0.0)
    {
      continue;
    }

    double thickness = array->GetValue(ptId);
    if (thickness <= 0.0)
    {
      continue;
    }

    double point[3];
    surface->GetPoint(ptId, point);

    // Unique one-ring neighbors; unlike the clamp, which keeps only the
    // maximum and is therefore indifferent to visiting a neighbor once per
    // incident cell, the median would be biased by the duplicates. The
    // smallest altitude relative to the longest edge over the incident
    // triangles is collected in the same pass.
    neighbors.clear();
    double minAltitudeRatio = -1.0;
    surface->GetPointCells(ptId, cellIds);
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellIds->GetId(i), npts, pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (pts[j] == ptId)
        {
          continue;
        }
        if (std::find(neighbors.begin(), neighbors.end(), pts[j]) == neighbors.end())
        {
          neighbors.push_back(pts[j]);
        }
      }

      if (npts != 3)
      {
        continue;
      }
      double a[3], b[3], c[3];
      surface->GetPoint(pts[0], a);
      surface->GetPoint(pts[1], b);
      surface->GetPoint(pts[2], c);
      double ab[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
      double ac[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
      double bc[3] = {c[0]-b[0], c[1]-b[1], c[2]-b[2]};
      double cross[3] = {ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2],
          ab[0]*ac[1]-ab[1]*ac[0]};
      double area = 0.5*std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
      double lenAB = std::sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      double lenAC = std::sqrt(ac[0]*ac[0] + ac[1]*ac[1] + ac[2]*ac[2]);
      double lenBC = std::sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      double longest = lenAB;
      if (lenAC > longest) { longest = lenAC; }
      if (lenBC > longest) { longest = lenBC; }
      if (longest <= 0.0)
      {
        continue;
      }
      // The smallest altitude is 2*area/longestEdge, so this ratio is the
      // altitude measured in units of the longest edge.
      double altitudeRatio = 2.0*area/(longest*longest);
      if (minAltitudeRatio < 0.0 || altitudeRatio < minAltitudeRatio)
      {
        minAltitudeRatio = altitudeRatio;
      }
    }

    if (neighbors.empty())
    {
      continue;
    }

    if (minAltitudeRatio >= 0.0 && minAltitudeRatio < sliverAltitudeRatio)
    {
      sliverAdjacent[ptId] = 1;
      numSliverAdjacent++;
    }

    // Curvature of every unique neighbor, with the neighbors at or below the
    // tangent plane contributing zero. Keeping them makes the median describe
    // the point as a whole: it stays near zero where only one neighbor happens
    // to rise, and is only large where most of the one-ring is concave.
    curvatures.clear();
    double maxCurvature = 0.0;
    for (auto neighborId : neighbors)
    {
      double neighbor[3];
      surface->GetPoint(neighborId, neighbor);
      double offset[3] = {neighbor[0]-point[0], neighbor[1]-point[1],
          neighbor[2]-point[2]};
      double distanceSquared = offset[0]*offset[0] + offset[1]*offset[1] +
          offset[2]*offset[2];
      if (distanceSquared <= 0.0)
      {
        continue;
      }
      double height = normalSign*(offset[0]*normal[0] + offset[1]*normal[1] +
          offset[2]*normal[2]) / normalLength;
      double curvature = (height <= 0.0) ? 0.0 : 2.0*height/distanceSquared;
      curvatures.push_back(curvature);
      if (curvature > maxCurvature)
      {
        maxCurvature = curvature;
      }
    }

    if (curvatures.empty() || maxCurvature <= 0.0)
    {
      continue;   // convex or flat: the outer wall spreads instead of folding
    }

    std::sort(curvatures.begin(), curvatures.end());
    size_t middle = curvatures.size()/2;
    double medianCurvature = (curvatures.size() % 2 == 1) ? curvatures[middle] :
        0.5*(curvatures[middle-1] + curvatures[middle]);

    numConcave++;
    radiusSmallest[ptId] = 1.0/maxCurvature;
    radiusTypical[ptId] = (medianCurvature > 0.0) ? 1.0/medianCurvature : 0.0;
    double ratio = thickness*medianCurvature;   // t / R_typical
    ratioTypical[ptId] = ratio;

    if (ratio >= 1.0)
    {
      numAtLeast1++;
      // The clustering seeds from the smallest key, so the worst (largest)
      // ratio is passed negated.
      flagged.push_back(std::make_pair(-ratio, ptId));
    }
    if (ratio >= 2.0) { numAtLeast2++; }
    if (ratio >= 4.0) { numAtLeast4++; }
  }

  // Leave the field on the surface so it can be viewed directly. The region
  // summary below is a list of seeds, which answers "how bad is the worst" but
  // not "which junctions light up"; that one is answered by looking at the
  // field, so it has to leave the log behind and reach the written mesh.
  {
    auto ratioArray = vtkSmartPointer<vtkDoubleArray>::New();
    ratioArray->SetName("ThicknessOverRadius");
    ratioArray->SetNumberOfComponents(1);
    ratioArray->SetNumberOfTuples(numPts);
    auto typicalArray = vtkSmartPointer<vtkDoubleArray>::New();
    typicalArray->SetName("ConcaveRadiusTypical");
    typicalArray->SetNumberOfComponents(1);
    typicalArray->SetNumberOfTuples(numPts);
    auto smallestArray = vtkSmartPointer<vtkDoubleArray>::New();
    smallestArray->SetName("ConcaveRadiusSmallest");
    smallestArray->SetNumberOfComponents(1);
    smallestArray->SetNumberOfTuples(numPts);
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      ratioArray->SetValue(ptId, ratioTypical[ptId]);
      typicalArray->SetValue(ptId, radiusTypical[ptId]);
      smallestArray->SetValue(ptId, radiusSmallest[ptId]);
    }
    surface->GetPointData()->RemoveArray("ThicknessOverRadius");
    surface->GetPointData()->RemoveArray("ConcaveRadiusTypical");
    surface->GetPointData()->RemoveArray("ConcaveRadiusSmallest");
    surface->GetPointData()->AddArray(ratioArray);
    surface->GetPointData()->AddArray(typicalArray);
    surface->GetPointData()->AddArray(smallestArray);
  }

  fprintf(stdout,"Concave curvature vs requested thickness (t/R, before any thickness reduction) [%s]:\n", label);
  fprintf(stdout,"  concave points: %d of %lld; points with t/R >= 1 / 2 / 4: %d/%d/%d\n",
      numConcave, (long long)numPts, numAtLeast1, numAtLeast2, numAtLeast4);
  fprintf(stdout,"  t/R > 1 means the extruded outer surface must self-intersect there whatever the extrusion does;\n");
  fprintf(stdout,"  points using a near-degenerate triangle (altitude < %g of the longest edge): %d\n",
      sliverAltitudeRatio, numSliverAdjacent);

  if (flagged.empty())
  {
    fprintf(stdout,"  no point requests a thickness above its concave radius of curvature\n");
    return SV_OK;
  }

  const int maxRegions = 8;
  const double radiusFraction = 0.02;
  std::vector<TGenUtilsPointRegion> regions;
  double regionRadius = 0.0;
  int numOutside = 0;
  int numRegionsTotal = 0;
  if (TGenUtils_ClusterPointsIntoRegions(surface, flagged, maxRegions, radiusFraction,
        regions, regionRadius, numOutside, numRegionsTotal) != SV_OK)
  {
    fprintf(stderr,"Problem clustering the points whose thickness exceeds the concave radius of curvature\n");
    return SV_ERROR;
  }

  // The total is the number that answers whether every junction is affected or
  // only a few, so it is reported before the truncated list of the worst ones.
  fprintf(stdout,"  concave regions with t/R >= 1: %d in total (separated by %.4g), worst %d shown:\n",
      numRegionsTotal, regionRadius, (int)regions.size());
  for (size_t i = 0; i < regions.size(); i++)
  {
    vtkIdType seedId = regions[i].seedId;
    double seed[3];
    surface->GetPoint(seedId, seed);
    fprintf(stdout,"    [%d] t/R %.3g (t %.5g, R_typical %.5g, R_smallest %.5g) at (%.5g, %.5g, %.5g), %d points%s\n",
        (int)(i+1), ratioTypical[seedId], array->GetValue(seedId), radiusTypical[seedId],
        radiusSmallest[seedId], seed[0], seed[1], seed[2], regions[i].numPoints,
        sliverAdjacent[seedId] ? "  [near-degenerate triangle]" : "");
  }
  if (numOutside > 0)
  {
    fprintf(stdout,"    ... %d further points with t/R >= 1 in the remaining %d regions\n",
        numOutside, numRegionsTotal - (int)regions.size());
  }

  return SV_OK;
}

// ------------------------------------
// TGenUtils_LimitThicknessGradation
// ------------------------------------
/**
 * @brief Limits how fast the wall thickness may change from point to point, by
 * lowering the thickness wherever it stands too far above a neighbor.
 * @note Every pass that reduces the thickness does so at the points that need
 * it and leaves their neighbors alone, which turns a local reduction into a
 * cliff in the thickness field. That cliff is itself destructive: the fold
 * prevention pass notes that an imbalance between the three points of a
 * triangle folds it on its own, so a reduction meant to remove a fold creates
 * the conditions for new ones. The curvature clamp dropping a point from its
 * requested thickness to the local concave radius in one step, and the fold
 * prevention pass levelling a triangle to its smallest thickness, are both this
 * same act at different doses; measured, the clamp alone took the worst element
 * aspect ratio from 1603 to 150779 while gaining nothing.
 *
 * The cure is to bound the gradient of the thickness field, the way a mesh
 * sizing field is gradation-limited before it is used (MMG reports its own
 * GRADATION for exactly this reason). This enforces
 *
 *     t_i <= t_j + maxSlope * |P_i - P_j|
 *
 * over every edge, by lowering t_i. Two properties matter and both come from
 * only ever lowering: the pass cannot raise a thickness back above a ceiling
 * another pass established, so it composes with the clamp and the fold
 * prevention pass without a re-clamp afterwards; and since the values decrease
 * monotonically and are bounded below, the relaxation terminates. This replaces
 * the role the Laplacian thickness smoothing was meant to play - that one
 * averages, so it raises a clamped point back above its ceiling and needs the
 * clamp applied again, which restores the very cliff it was there to remove.
 *
 * A bounded gradient does not by itself guarantee a fold-free extrusion: a fold
 * is driven by the thickness against the concave radius of curvature, which
 * this pass does not change. It removes the cliffs, not the infeasibility.
 * @param surface The surface the thickness belongs to.
 * @param array The thickness point array, modified in place. Values are only
 * ever lowered.
 * @param maxSlope The largest allowed change in thickness per unit distance.
 * Zero or negative disables the pass. A value of one lets the thickness change
 * by the edge length over one edge, which is a 45 degree taper of the outer
 * surface relative to the inner one.
 * @param label Names the field being limited, for the log.
 * @return SV_OK if the thickness is limited.
 */

int TGenUtils_LimitThicknessGradation(vtkPolyData *surface, vtkDoubleArray *array,
    double maxSlope, const char *label)
{
  if (maxSlope <= 0.0)
  {
    return SV_OK;
  }

  if (surface == nullptr || array == nullptr)
  {
    fprintf(stderr,"Cannot limit the thickness gradation without a surface and an array\n");
    return SV_ERROR;
  }

  if (label == nullptr)
  {
    label = "wall thickness";
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array must have one component and one tuple per surface point\n");
    return SV_ERROR;
  }

  surface->BuildLinks();

  // One-ring neighbors with the edge lengths, built once. The lengths are kept
  // alongside so the relaxation below does not re-read the coordinates.
  std::vector<std::vector<std::pair<vtkIdType,double> > > neighbors(numPts);
  auto cellIds = vtkSmartPointer<vtkIdList>::New();
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double p[3];
    surface->GetPoint(ptId, p);
    auto& ptNeighbors = neighbors[ptId];
    surface->GetPointCells(ptId, cellIds);
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellIds->GetId(i), npts, pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (pts[j] == ptId)
        {
          continue;
        }
        bool seen = false;
        for (auto& existing : ptNeighbors)
        {
          if (existing.first == pts[j])
          {
            seen = true;
            break;
          }
        }
        if (seen)
        {
          continue;
        }
        double q[3];
        surface->GetPoint(pts[j], q);
        double offset[3] = {q[0]-p[0], q[1]-p[1], q[2]-p[2]};
        double distance = std::sqrt(offset[0]*offset[0] + offset[1]*offset[1] +
            offset[2]*offset[2]);
        ptNeighbors.push_back(std::make_pair(pts[j], distance));
      }
    }
  }

  std::vector<double> original(numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    original[ptId] = array->GetValue(ptId);
  }

  // A point whose thickness has just been lowered can put its own neighbors
  // over the limit, so the lowered points are revisited rather than the whole
  // surface being swept a fixed number of times. Every point starts queued
  // because any of them may be the one others have to come down to.
  std::deque<vtkIdType> queue;
  std::vector<char> queued(numPts, 1);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    queue.push_back(ptId);
  }

  // Relative, so the guard is meaningful whatever the model units are. Without
  // it a chain of reductions each smaller than the floating point resolution
  // could keep requeueing points.
  const double relativeTolerance = 1.0e-9;

  while (!queue.empty())
  {
    vtkIdType ptId = queue.front();
    queue.pop_front();
    queued[ptId] = 0;

    double thickness = array->GetValue(ptId);
    for (auto& neighbor : neighbors[ptId])
    {
      double limit = thickness + maxSlope*neighbor.second;
      double neighborThickness = array->GetValue(neighbor.first);
      if (neighborThickness <= limit + relativeTolerance*std::fabs(limit))
      {
        continue;
      }
      array->SetValue(neighbor.first, limit);
      if (!queued[neighbor.first])
      {
        queued[neighbor.first] = 1;
        queue.push_back(neighbor.first);
      }
    }
  }

  int numLowered = 0;
  double maxReduction = 0.0;
  vtkIdType maxReductionId = -1;
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double reduction = original[ptId] - array->GetValue(ptId);
    if (reduction <= 0.0)
    {
      continue;
    }
    numLowered++;
    if (reduction > maxReduction)
    {
      maxReduction = reduction;
      maxReductionId = ptId;
    }
  }

  fprintf(stdout,"Thickness gradation limit (max change %g per unit distance) [%s]: "
      "lowered %d of %lld points; largest reduction %.5g",
      maxSlope, label, numLowered, (long long)numPts, maxReduction);
  if (maxReductionId >= 0)
  {
    double p[3];
    surface->GetPoint(maxReductionId, p);
    fprintf(stdout," at (%.5g, %.5g, %.5g) (%.5g -> %.5g)", p[0], p[1], p[2],
        original[maxReductionId], array->GetValue(maxReductionId));
  }
  fprintf(stdout,"\n");

  return SV_OK;
}

// -------------------------------------
// TGenUtils_ExtractBoundaryLoops
// -------------------------------------
/**
 * @brief Collects the boundary edges of a surface into ordered closed loops.
 * @note The passes that only need to know whether a point sits on the boundary
 * test one edge at a time. Closing the surface needs more than that: the cap
 * has to be filled against the rim in the rim's own order, so the loop has to
 * be walked. A boundary edge is an edge used by a single cell, the same test
 * used elsewhere, and it is recorded in the winding order of that cell. Each
 * boundary point then starts exactly one boundary edge, so following that map
 * from any point returns to it and traverses its loop once.
 *
 * The loop order is the order of the cells that own it, so the edge
 * loop[m] -> loop[m+1] is directed as the wall triangle traverses it. A facet
 * closing the loop has to traverse the same edge the other way round, which is
 * what keeps the closed surface consistently oriented.
 * @param surface The surface; its cells must be triangles.
 * @param loops Set to the boundary loops, each a list of point ids in order.
 * @return SV_OK if every boundary edge belongs to a simple closed loop.
 */

int TGenUtils_ExtractBoundaryLoops(vtkPolyData *surface,
    std::vector<std::vector<vtkIdType> > &loops)
{
  loops.clear();

  if (surface == nullptr)
  {
    fprintf(stderr,"Cannot extract the boundary loops without a surface\n");
    return SV_ERROR;
  }

  surface->BuildLinks();
  auto edgeNeighbors = vtkSmartPointer<vtkIdList>::New();

  std::map<vtkIdType,vtkIdType> nextPoint;
  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    surface->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }

    for (vtkIdType j = 0; j < npts; j++)
    {
      vtkIdType a = pts[j];
      vtkIdType b = pts[(j+1)%npts];
      surface->GetCellEdgeNeighbors(cellId, a, b, edgeNeighbors);
      if (edgeNeighbors->GetNumberOfIds() != 0)
      {
        continue;
      }
      if (nextPoint.find(a) != nextPoint.end())
      {
        fprintf(stderr,"Point %lld starts more than one boundary edge, so the surface boundary is not a set of simple loops and cannot be capped\n",
            (long long)a);
        return SV_ERROR;
      }
      nextPoint[a] = b;
    }
  }

  std::set<vtkIdType> visited;
  for (std::map<vtkIdType,vtkIdType>::const_iterator it = nextPoint.begin();
       it != nextPoint.end(); ++it)
  {
    vtkIdType start = it->first;
    if (visited.find(start) != visited.end())
    {
      continue;
    }

    std::vector<vtkIdType> loop;
    vtkIdType current = start;
    while (visited.find(current) == visited.end())
    {
      visited.insert(current);
      loop.push_back(current);

      std::map<vtkIdType,vtkIdType>::const_iterator step = nextPoint.find(current);
      if (step == nextPoint.end())
      {
        fprintf(stderr,"The boundary edge chain from point %lld ends at point %lld instead of closing\n",
            (long long)start, (long long)current);
        return SV_ERROR;
      }
      current = step->second;
    }

    if (current != start)
    {
      fprintf(stderr,"A boundary edge chain closed onto point %lld rather than onto its start %lld\n",
          (long long)current, (long long)start);
      return SV_ERROR;
    }

    if (loop.size() >= 3)
    {
      loops.push_back(loop);
    }
  }

  return SV_OK;
}

// -------------------------------------
// TGenUtils_BuildOffsetOuterSurface
// -------------------------------------
/**
 * @brief Builds the outer wall surface as the true offset of the inner surface
 * at the requested thickness, by contouring a signed distance field.
 * @note Moving each point out along its own normal gives every outer point one
 * inner point, and that correspondence is what cannot represent the answer at a
 * junction. Dilating a solid by t rounds its convex features to radius t, but
 * at a concave crotch the two offset sheets run into each other: the correct
 * outer surface is their intersection curve, a crease, and the parts of both
 * sheets beyond it are not on the boundary at all. The inner points whose
 * offset lands in the discarded part simply have no outer point. Every pass
 * that kept the correspondence had to pay for that with thickness, which is why
 * the wall thinned exactly where the geometry is concave.
 *
 * Contouring the distance field has no correspondence to keep. The level set
 * d(x) = t is the dilated boundary by construction, so the crease and the
 * rounding come out of it rather than being aimed at, and every point of the
 * result is at least t from the inner surface - the invariant the clearance
 * constraint tried and failed to impose on a fixed triangulation.
 *
 * The distance has to be signed to tell the wall side from the lumen, which
 * needs a closed surface, so the cap rims are filled with a fan first. The
 * result therefore also covers the caps with a t-thick dome, which the caller
 * trims back to the cap planes.
 *
 * The field is only evaluated in a band around the surface, since that is where
 * the level set is; the rest of the grid is flooded from outside so that the
 * lumen keeps a negative sign and no spurious sheet is contoured in it.
 * @param surface The inner surface, open at the caps.
 * @param array The requested thickness per point; one tuple per surface point.
 * @param targetEdgeSize The mesh edge size, or a non-positive value if unknown.
 * Used with the thickness to choose the grid spacing.
 * @param maxVoxels The largest grid to build; the spacing is coarsened until
 * the grid fits and the outcome is reported.
 * @param outer Set to the contoured offset surface.
 * @return SV_OK if the offset surface is built.
 */

int TGenUtils_BuildOffsetOuterSurface(vtkPolyData *surface, vtkDoubleArray *array,
    double targetEdgeSize, vtkIdType maxVoxels, vtkPolyData *outer)
{
  if (surface == nullptr || array == nullptr || outer == nullptr)
  {
    fprintf(stderr,"Cannot build the offset outer surface without a surface, a thickness array and an output\n");
    return SV_ERROR;
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array must have one component and one tuple per surface point\n");
    return SV_ERROR;
  }

  double thicknessMin = std::numeric_limits<double>::max();
  double thicknessMax = 0.0;
  int numPositive = 0;
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double t = array->GetValue(ptId);
    if (t <= 0.0)
    {
      continue;
    }
    thicknessMin = std::min(thicknessMin, t);
    thicknessMax = std::max(thicknessMax, t);
    numPositive++;
  }
  if (numPositive == 0)
  {
    fprintf(stderr,"Every wall thickness is zero or negative, so there is no wall to offset\n");
    return SV_ERROR;
  }

  // Close the cap openings so the distance can be signed. The fan traverses
  // each rim edge opposite to the wall triangle that owns it, which is what
  // makes the closed surface consistently oriented and therefore what decides
  // the sign; it is checked against a grid corner below rather than assumed.
  std::vector<std::vector<vtkIdType> > loops;
  if (TGenUtils_ExtractBoundaryLoops(surface, loops) != SV_OK)
  {
    fprintf(stderr,"Problem extracting the cap rims of the wall surface\n");
    return SV_ERROR;
  }

  auto closedPoints = vtkSmartPointer<vtkPoints>::New();
  closedPoints->DeepCopy(surface->GetPoints());

  std::vector<double> closedThickness((size_t)numPts, 0.0);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    closedThickness[(size_t)ptId] = std::max(array->GetValue(ptId), 0.0);
  }

  auto closedCells = vtkSmartPointer<vtkCellArray>::New();
  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    surface->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }
    vtkIdType triangle[3] = {pts[0], pts[1], pts[2]};
    closedCells->InsertNextCell(3, triangle);
  }

  for (size_t i = 0; i < loops.size(); i++)
  {
    const std::vector<vtkIdType> &loop = loops[i];
    double centroid[3] = {0.0, 0.0, 0.0};
    double meanThickness = 0.0;
    for (size_t m = 0; m < loop.size(); m++)
    {
      double p[3];
      closedPoints->GetPoint(loop[m], p);
      for (int k = 0; k < 3; k++)
      {
        centroid[k] += p[k];
      }
      meanThickness += closedThickness[(size_t)loop[m]];
    }
    for (int k = 0; k < 3; k++)
    {
      centroid[k] /= (double)loop.size();
    }
    meanThickness /= (double)loop.size();

    vtkIdType centroidId = closedPoints->InsertNextPoint(centroid);
    closedThickness.push_back(meanThickness);

    for (size_t m = 0; m < loop.size(); m++)
    {
      vtkIdType a = loop[m];
      vtkIdType b = loop[(m+1)%loop.size()];
      vtkIdType fan[3] = {b, a, centroidId};
      closedCells->InsertNextCell(3, fan);
    }
  }

  auto closed = vtkSmartPointer<vtkPolyData>::New();
  closed->SetPoints(closedPoints);
  closed->SetPolys(closedCells);
  closed->BuildLinks();

  // The grid has to hold the whole offset surface, so it covers the model plus
  // the largest thickness, and a few cells beyond that so the band around the
  // surface never touches the grid face and the outside flood always has a
  // corner to start from.
  double bounds[6];
  closed->GetBounds(bounds);

  double resolution = thicknessMin;
  if (targetEdgeSize > 0.0 && targetEdgeSize < resolution)
  {
    resolution = targetEdgeSize;
  }
  double requestedSpacing = 0.5*resolution;
  double spacing = requestedSpacing;

  int dims[3] = {0, 0, 0};
  double margin = 0.0;
  for (int attempt = 0; attempt < 64; attempt++)
  {
    // Two cells wider than the widest band a triangle can mark, so the grid
    // corners the outside flood starts from are never part of the band.
    margin = 1.25*thicknessMax + 5.0*spacing;
    double total = 1.0;
    for (int k = 0; k < 3; k++)
    {
      double extent = bounds[2*k+1] - bounds[2*k] + 2.0*margin;
      dims[k] = (int)std::ceil(extent/spacing) + 1;
      if (dims[k] < 2)
      {
        dims[k] = 2;
      }
      total *= (double)dims[k];
    }
    if (total <= (double)maxVoxels)
    {
      break;
    }
    spacing *= 1.25;
  }

  double origin[3];
  for (int k = 0; k < 3; k++)
  {
    origin[k] = bounds[2*k] - margin;
  }

  size_t numVoxels = (size_t)dims[0]*(size_t)dims[1]*(size_t)dims[2];

  fprintf(stdout,"Wall outer surface by distance field offset:\n");
  fprintf(stdout,"  thickness %.5g to %.5g, %zu cap rims closed with a fan\n",
      thicknessMin, thicknessMax, loops.size());
  fprintf(stdout,"  grid %d x %d x %d = %zu voxels at spacing %.5g, which is %.2f cells per smallest wall thickness\n",
      dims[0], dims[1], dims[2], numVoxels, spacing, thicknessMin/spacing);
  if (spacing > requestedSpacing*1.001)
  {
    fprintf(stdout,"  the spacing was coarsened from %.5g to fit the %lld voxel budget\n",
        requestedSpacing, (long long)maxVoxels);
  }
  if (thicknessMin/spacing < 1.5)
  {
    fprintf(stdout,"  WARNING: fewer than 1.5 cells per wall thickness. The smooth part of the offset survives this, but the crease at a junction is rounded over about half a cell, so the wall there will read short in the offset thickness report below. Contouring in slabs rather than over one grid is what buys resolution here; raising the voxel budget alone will not reach far on a model this long.\n");
  }

  // Mark the band the level set can pass through: every voxel within the local
  // thickness of the surface, plus two cells so the contour has values on both
  // sides of it everywhere. Marking is per triangle and the marks overlap, but
  // they are idempotent writes, and what the band buys is that the distance is
  // only evaluated where it can matter.
  std::vector<unsigned char> state(numVoxels, 0);
  const unsigned char kUnknown = 0, kBand = 1, kOutside = 2;

  for (vtkIdType cellId = 0; cellId < closed->GetNumberOfCells(); cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    closed->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }

    double lo[3] = {0.0, 0.0, 0.0}, hi[3] = {0.0, 0.0, 0.0};
    double radius = 0.0;
    for (vtkIdType j = 0; j < npts; j++)
    {
      double p[3];
      closedPoints->GetPoint(pts[j], p);
      for (int k = 0; k < 3; k++)
      {
        if (j == 0 || p[k] < lo[k]) { lo[k] = p[k]; }
        if (j == 0 || p[k] > hi[k]) { hi[k] = p[k]; }
      }
      radius = std::max(radius, closedThickness[(size_t)pts[j]]);
    }

    // The band has to reach past the level set on both sides, and the value at
    // a voxel uses the thickness of its nearest point, which need not be one of
    // this triangle's. The quarter of slack absorbs that difference where the
    // thickness varies; the fill below checks that it was in fact enough.
    radius = 1.25*radius + 3.0*spacing;

    int begin[3], end[3];
    for (int k = 0; k < 3; k++)
    {
      begin[k] = (int)std::floor((lo[k] - radius - origin[k])/spacing);
      end[k] = (int)std::ceil((hi[k] + radius - origin[k])/spacing);
      begin[k] = std::max(begin[k], 0);
      end[k] = std::min(end[k], dims[k]-1);
    }

    for (int k = begin[2]; k <= end[2]; k++)
    {
      for (int j = begin[1]; j <= end[1]; j++)
      {
        size_t row = (size_t)begin[0] + (size_t)dims[0]*((size_t)j + (size_t)dims[1]*(size_t)k);
        for (int i = begin[0]; i <= end[0]; i++, row++)
        {
          state[row] = kBand;
        }
      }
    }
  }

  auto implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  implicit->SetInput(closed);

  // Which side of the surface counts as negative is a convention, and the
  // offset has to be built outward whichever way round it is. The grid corner
  // sits outside the model by at least the margin in every direction, so the
  // sign there is the sign of the outside and everything else follows from it.
  // Its magnitude is checked as well: a corner that is not that far from the
  // surface means the grid or the margin is not what this assumes, which would
  // make the reading meaningless rather than merely backwards.
  double corner[3] = {origin[0], origin[1], origin[2]};
  double cornerDistance = implicit->EvaluateFunction(corner);
  if (std::abs(cornerDistance) < 0.9*margin)
  {
    fprintf(stderr,"The signed distance at the grid corner is %.5g but the corner lies at least %.5g from the model, so the distance field is not measuring what the offset needs\n",
        cornerDistance, margin);
    return SV_ERROR;
  }
  double outwardSign = (cornerDistance > 0.0) ? 1.0 : -1.0;
  if (outwardSign < 0.0)
  {
    fprintf(stdout,"  the distance field is negative outside the model, so its sign is used inverted\n");
  }

  // The thickness is looked up on the inner surface rather than on the closed
  // one, so the fan centres never contribute a thickness of their own; a voxel
  // over a cap takes the thickness of the nearest rim point.
  auto thicknessLocator = vtkSmartPointer<vtkStaticPointLocator>::New();
  thicknessLocator->SetDataSet(surface);
  thicknessLocator->BuildLocator();

  auto values = vtkSmartPointer<vtkFloatArray>::New();
  values->SetName("WallOffsetLevel");
  values->SetNumberOfComponents(1);
  values->SetNumberOfTuples((vtkIdType)numVoxels);

  double farValue = thicknessMax + (bounds[1]-bounds[0]) + (bounds[3]-bounds[2]) + (bounds[5]-bounds[4]);
  size_t numBandVoxels = 0;

  for (int k = 0; k < dims[2]; k++)
  {
    for (int j = 0; j < dims[1]; j++)
    {
      size_t index = (size_t)dims[0]*((size_t)j + (size_t)dims[1]*(size_t)k);
      for (int i = 0; i < dims[0]; i++, index++)
      {
        if (state[index] != kBand)
        {
          continue;
        }
        numBandVoxels++;

        double x[3] = {origin[0] + spacing*i, origin[1] + spacing*j, origin[2] + spacing*k};
        double distance = outwardSign*implicit->EvaluateFunction(x);
        vtkIdType nearest = thicknessLocator->FindClosestPoint(x);
        double thickness = (nearest >= 0) ? std::max(array->GetValue(nearest), 0.0) : 0.0;
        values->SetValue((vtkIdType)index, (float)(distance - thickness));
      }
    }
  }

  // Everything the band does not cover is either well outside the offset or
  // inside the lumen, and the two need opposite signs or the contour would find
  // a sheet between them. The band is at least two cells thick and closed
  // around the model, so a flood that starts outside cannot leak through it.
  std::deque<size_t> queue;
  const int cornerIndex[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
                                 {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
  for (int c = 0; c < 8; c++)
  {
    size_t i = (size_t)(cornerIndex[c][0] ? dims[0]-1 : 0);
    size_t j = (size_t)(cornerIndex[c][1] ? dims[1]-1 : 0);
    size_t k = (size_t)(cornerIndex[c][2] ? dims[2]-1 : 0);
    size_t index = i + (size_t)dims[0]*(j + (size_t)dims[1]*k);
    if (state[index] == kUnknown)
    {
      state[index] = kOutside;
      queue.push_back(index);
    }
  }

  if (queue.empty())
  {
    fprintf(stderr,"Every corner of the distance grid is inside the evaluated band, so there is nowhere to start the outside flood from and the lumen cannot be told from the exterior\n");
    return SV_ERROR;
  }

  while (!queue.empty())
  {
    size_t index = queue.front();
    queue.pop_front();

    int i = (int)(index % (size_t)dims[0]);
    int j = (int)((index / (size_t)dims[0]) % (size_t)dims[1]);
    int k = (int)(index / ((size_t)dims[0]*(size_t)dims[1]));

    const int step[6][3] = {{-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};
    for (int s = 0; s < 6; s++)
    {
      int ni = i + step[s][0], nj = j + step[s][1], nk = k + step[s][2];
      if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1] || nk < 0 || nk >= dims[2])
      {
        continue;
      }
      size_t neighbor = (size_t)ni + (size_t)dims[0]*((size_t)nj + (size_t)dims[1]*(size_t)nk);
      if (state[neighbor] != kUnknown)
      {
        continue;
      }
      state[neighbor] = kOutside;
      queue.push_back(neighbor);
    }
  }

  size_t numInside = 0;
  for (size_t index = 0; index < numVoxels; index++)
  {
    if (state[index] == kBand)
    {
      continue;
    }
    if (state[index] == kOutside)
    {
      values->SetValue((vtkIdType)index, (float)farValue);
    }
    else
    {
      values->SetValue((vtkIdType)index, (float)(-farValue));
      numInside++;
    }
  }

  fprintf(stdout,"  distance evaluated at %zu band voxels; %zu enclosed and %zu outside voxels filled by sign\n",
      numBandVoxels, numInside, numVoxels - numBandVoxels - numInside);

  // The signs the fill wrote are only right if the level set stayed inside the
  // band. Where it did not, a band voxel sits against a filled one of the
  // opposite sign and the contour would follow the edge of the band instead of
  // the offset. That is a wrong surface rather than a rough one, so it is
  // reported as an error rather than contoured.
  size_t numStraddling = 0;
  for (int k = 0; k < dims[2]; k++)
  {
    for (int j = 0; j < dims[1]; j++)
    {
      size_t index = (size_t)dims[0]*((size_t)j + (size_t)dims[1]*(size_t)k);
      for (int i = 0; i < dims[0]; i++, index++)
      {
        if (state[index] != kBand)
        {
          continue;
        }
        double value = values->GetValue((vtkIdType)index);

        const int step[6][3] = {{-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};
        for (int s = 0; s < 6; s++)
        {
          int ni = i + step[s][0], nj = j + step[s][1], nk = k + step[s][2];
          if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1] || nk < 0 || nk >= dims[2])
          {
            continue;
          }
          size_t neighbor = (size_t)ni + (size_t)dims[0]*((size_t)nj + (size_t)dims[1]*(size_t)nk);
          if (state[neighbor] == kBand)
          {
            continue;
          }
          double filled = (state[neighbor] == kOutside) ? farValue : -farValue;
          if ((value < 0.0) != (filled < 0.0))
          {
            numStraddling++;
            break;
          }
        }
      }
    }
  }

  if (numStraddling > 0)
  {
    fprintf(stderr,"The offset level set leaves the evaluated band at %zu voxels, so the band around the surface is too thin to hold it. This happens when the wall thickness varies faster than the band allows for; contouring it would follow the edge of the band rather than the offset.\n",
        numStraddling);
    return SV_ERROR;
  }

  auto image = vtkSmartPointer<vtkImageData>::New();
  image->SetDimensions(dims[0], dims[1], dims[2]);
  image->SetOrigin(origin[0], origin[1], origin[2]);
  image->SetSpacing(spacing, spacing, spacing);
  image->GetPointData()->SetScalars(values);

  auto contour = vtkSmartPointer<vtkFlyingEdges3D>::New();
  contour->SetInputData(image);
  contour->SetValue(0, 0.0);
  contour->ComputeNormalsOff();
  contour->ComputeGradientsOff();
  contour->ComputeScalarsOff();
  contour->Update();

  outer->Initialize();
  outer->DeepCopy(contour->GetOutput());

  if (outer->GetNumberOfCells() == 0)
  {
    fprintf(stderr,"The offset level set is empty, so no outer wall surface was produced\n");
    return SV_ERROR;
  }

  // More than one shell means the offset is not simply the outside of the wall:
  // a thickness large enough to close a lumen leaves a sheet inside it, and a
  // model with detached vessels offsets into one shell per vessel. Neither is
  // an error here, but both change what the fill will be asked to do.
  auto connectivity = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  connectivity->SetInputData(outer);
  connectivity->SetExtractionModeToAllRegions();
  connectivity->Update();

  fprintf(stdout,"  offset surface has %lld points and %lld triangles in %d connected shells\n",
      (long long)outer->GetNumberOfPoints(), (long long)outer->GetNumberOfCells(),
      connectivity->GetNumberOfExtractedRegions());

  return SV_OK;
}

// -------------------------------------
// TGenUtils_TrimOffsetSurfaceAtCaps
// -------------------------------------
/**
 * @brief Trims the offset outer surface back to the cap planes of the inner
 * surface, and reports the rim it was trimmed to alongside the inner rim.
 * @note The offset is built from a surface whose cap openings were filled, so
 * it covers each vessel end with a dome a thickness deep. The wall does not
 * extend over the end - the lumen opens there - so the dome comes off and the
 * two rims left behind are closed to each other instead.
 *
 * The cut is the cap plane, but only near the cap. An infinite plane would also
 * cut whatever else of the model happens to lie on its far side, which for a
 * vessel that curves back on itself is real wall; keeping the cut inside a
 * sphere around the rim bounds what it can reach to the end it belongs to.
 *
 * The plane is oriented from the rim itself rather than from the point normals,
 * which the capping pass has already laid into the cap plane. The rim is
 * traversed in the winding of the wall triangles that own it, and that winding
 * runs clockwise about the outward direction, so the outward direction is the
 * reverse of the rim's own normal.
 * @param surface The inner surface, open at the caps.
 * @param outer The offset surface; trimmed in place.
 * @param maxThickness The largest wall thickness, which bounds how far past a
 * rim the dome to be cut off can reach.
 * @param caps Set to one entry per vessel end, holding both rims and the plane.
 * @return SV_OK if every cap was trimmed and its two rims paired.
 */

int TGenUtils_TrimOffsetSurfaceAtCaps(vtkPolyData *surface, vtkPolyData *outer,
    double maxThickness, std::vector<TGenUtilsCapRim> &caps)
{
  caps.clear();

  if (surface == nullptr || outer == nullptr)
  {
    fprintf(stderr,"Cannot trim the offset surface without an inner surface and an offset surface\n");
    return SV_ERROR;
  }

  std::vector<std::vector<vtkIdType> > innerLoops;
  if (TGenUtils_ExtractBoundaryLoops(surface, innerLoops) != SV_OK)
  {
    fprintf(stderr,"Problem extracting the cap rims of the wall surface\n");
    return SV_ERROR;
  }

  if (innerLoops.empty())
  {
    fprintf(stdout,"  the inner surface is closed, so the offset surface needs no trimming\n");
    return SV_OK;
  }

  std::vector<double> capRadius(innerLoops.size(), 0.0);
  caps.resize(innerLoops.size());

  for (size_t c = 0; c < innerLoops.size(); c++)
  {
    const std::vector<vtkIdType> &loop = innerLoops[c];
    TGenUtilsCapRim &cap = caps[c];
    cap.innerLoop = loop;

    for (int k = 0; k < 3; k++)
    {
      cap.origin[k] = 0.0;
    }
    for (size_t m = 0; m < loop.size(); m++)
    {
      double p[3];
      surface->GetPoint(loop[m], p);
      for (int k = 0; k < 3; k++)
      {
        cap.origin[k] += p[k];
      }
    }
    for (int k = 0; k < 3; k++)
    {
      cap.origin[k] /= (double)loop.size();
    }

    // Newell's normal, which is the rim's own normal for the order it is
    // stored in; the outward direction is its reverse.
    double normal[3] = {0.0, 0.0, 0.0};
    for (size_t m = 0; m < loop.size(); m++)
    {
      double p[3], q[3];
      surface->GetPoint(loop[m], p);
      surface->GetPoint(loop[(m+1)%loop.size()], q);
      normal[0] += (p[1]-q[1])*(p[2]+q[2]);
      normal[1] += (p[2]-q[2])*(p[0]+q[0]);
      normal[2] += (p[0]-q[0])*(p[1]+q[1]);
    }
    if (vtkMath::Normalize(normal) <= 0.0)
    {
      fprintf(stderr,"A cap rim of %zu points at (%.5g, %.5g, %.5g) encloses no area, so its plane cannot be found\n",
          loop.size(), cap.origin[0], cap.origin[1], cap.origin[2]);
      return SV_ERROR;
    }
    for (int k = 0; k < 3; k++)
    {
      cap.outward[k] = -normal[k];
    }

    for (size_t m = 0; m < loop.size(); m++)
    {
      double p[3];
      surface->GetPoint(loop[m], p);
      capRadius[c] = std::max(capRadius[c], std::sqrt(vtkMath::Distance2BetweenPoints(p, cap.origin)));
    }
    if (capRadius[c] <= 0.0)
    {
      fprintf(stderr,"A cap rim at (%.5g, %.5g, %.5g) has zero radius\n",
          cap.origin[0], cap.origin[1], cap.origin[2]);
      return SV_ERROR;
    }
  }

  // Clip once per cap. The scalar is positive on everything that is kept: below
  // the plane, or far enough from this rim that the plane has no business
  // reaching it. The dome is the only place both are negative.
  //
  // The dome of a rim of radius R under a wall of thickness t meets the cap
  // plane at R + t, so the window has to hold that and no more than it needs
  // to. Sizing it off the thickness rather than off R keeps it tight where the
  // vessel is wide, and keeps it valid where the wall is thick relative to the
  // vessel - a window of a fixed multiple of R would fall inside the rim it is
  // meant to cut once t approached R.
  for (size_t c = 0; c < caps.size(); c++)
  {
    const TGenUtilsCapRim &cap = caps[c];
    double window = capRadius[c] + 2.5*maxThickness;

    auto level = vtkSmartPointer<vtkDoubleArray>::New();
    level->SetName("CapTrimLevel");
    level->SetNumberOfComponents(1);
    level->SetNumberOfTuples(outer->GetNumberOfPoints());
    for (vtkIdType ptId = 0; ptId < outer->GetNumberOfPoints(); ptId++)
    {
      double x[3];
      outer->GetPoint(ptId, x);
      double offset[3];
      vtkMath::Subtract(cap.origin, x, offset);
      double below = vtkMath::Dot(offset, cap.outward);
      double away = std::sqrt(vtkMath::Distance2BetweenPoints(x, cap.origin)) - window;
      level->SetValue(ptId, std::max(below, away));
    }
    outer->GetPointData()->SetScalars(level);

    auto clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputData(outer);
    clipper->GenerateClipScalarsOff();
    clipper->GenerateClippedOutputOff();
    clipper->InsideOutOff();
    clipper->SetValue(0.0);

    auto triangles = vtkSmartPointer<vtkTriangleFilter>::New();
    triangles->SetInputConnection(clipper->GetOutputPort());
    triangles->PassLinesOff();
    triangles->PassVertsOff();

    // Clipping leaves a pair of coincident points on every cut edge, and the
    // rim cannot be walked until they are one point.
    auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputConnection(triangles->GetOutputPort());
    cleaner->Update();

    // Only the triangles are carried on. Cleaning can turn a collapsed one into
    // a line, and a line sharing an edge with a triangle would make that edge
    // look interior when the rim is walked.
    auto trimmed = vtkSmartPointer<vtkPolyData>::New();
    trimmed->SetPoints(cleaner->GetOutput()->GetPoints());
    trimmed->SetPolys(cleaner->GetOutput()->GetPolys());
    outer->DeepCopy(trimmed);

    if (outer->GetNumberOfCells() == 0)
    {
      fprintf(stderr,"Trimming the offset surface at the cap at (%.5g, %.5g, %.5g) removed all of it\n",
          cap.origin[0], cap.origin[1], cap.origin[2]);
      return SV_ERROR;
    }
  }

  std::vector<std::vector<vtkIdType> > outerLoops;
  if (TGenUtils_ExtractBoundaryLoops(outer, outerLoops) != SV_OK)
  {
    fprintf(stderr,"Problem extracting the trimmed rims of the offset surface\n");
    return SV_ERROR;
  }

  if (outerLoops.size() != caps.size())
  {
    fprintf(stderr,"The trimmed offset surface has %zu rims but the wall has %zu cap openings. A cut has taken more than the dome off its end, which happens when a vessel curves back within a rim radius of another vessel's cap.\n",
        outerLoops.size(), caps.size());
    return SV_ERROR;
  }

  // Pair each trimmed rim with the cap whose plane it lies on. Being on the
  // plane is the test that matters, since two caps can be near each other but
  // only one cut produced this rim.
  std::vector<bool> used(outerLoops.size(), false);
  for (size_t c = 0; c < caps.size(); c++)
  {
    TGenUtilsCapRim &cap = caps[c];
    size_t best = outerLoops.size();
    double bestDeviation = 0.0;

    for (size_t l = 0; l < outerLoops.size(); l++)
    {
      if (used[l])
      {
        continue;
      }
      double deviation = 0.0;
      for (size_t m = 0; m < outerLoops[l].size(); m++)
      {
        double x[3];
        outer->GetPoint(outerLoops[l][m], x);
        double offset[3];
        vtkMath::Subtract(x, cap.origin, offset);
        deviation = std::max(deviation, std::abs(vtkMath::Dot(offset, cap.outward)));
      }
      if (best == outerLoops.size() || deviation < bestDeviation)
      {
        best = l;
        bestDeviation = deviation;
      }
    }

    if (best == outerLoops.size() || bestDeviation > 0.05*capRadius[c])
    {
      fprintf(stderr,"No trimmed rim lies on the cap plane at (%.5g, %.5g, %.5g); the closest is off it by %.5g against a rim radius of %.5g\n",
          cap.origin[0], cap.origin[1], cap.origin[2], bestDeviation, capRadius[c]);
      return SV_ERROR;
    }

    used[best] = true;
    cap.outerLoop = outerLoops[best];
  }

  fprintf(stdout,"  trimmed the offset surface at %zu cap planes, leaving %lld points and %lld triangles\n",
      caps.size(), (long long)outer->GetNumberOfPoints(), (long long)outer->GetNumberOfCells());
  for (size_t c = 0; c < caps.size(); c++)
  {
    fprintf(stdout,"    cap at (%.5g, %.5g, %.5g): inner rim %zu points, trimmed rim %zu points\n",
        caps[c].origin[0], caps[c].origin[1], caps[c].origin[2],
        caps[c].innerLoop.size(), caps[c].outerLoop.size());
  }

  return SV_OK;
}

// -------------------------------------
// TGenUtils_StitchCapAnnulus
// -------------------------------------
/**
 * @brief Closes the wall at a vessel end by triangulating between the inner cap
 * rim and the trimmed outer rim.
 * @note The extrusion could close its ends with one quad per rim edge, because
 * every outer point was one inner point offset. The offset surface is contoured
 * from a grid, so its rim has neither the same points nor the same number of
 * them, and the two have to be triangulated against each other instead.
 *
 * Both rims run around the same vessel end, so the angle about the cap axis
 * orders them both and the two can be merged on it, taking whichever rim is
 * behind at each step. This holds while each rim winds once around the axis,
 * which is checked rather than assumed - a rim that doubles back has no such
 * order and would be triangulated into overlapping facets.
 *
 * Neither rim is moved. The inner rim in particular is part of the fluid/wall
 * interface, and the solver matches it against the fluid mesh.
 * @param points The points both rims index into.
 * @param innerLoop The inner cap rim, in order.
 * @param outerLoop The trimmed outer rim, in order.
 * @param outward The direction out of the vessel end.
 * @param cells The annulus triangles are appended here, facing outward.
 * @param numDegenerate Set to the number of zero-area triangles produced.
 * @return SV_OK if the annulus is built.
 */

int TGenUtils_StitchCapAnnulus(vtkPoints *points,
    const std::vector<vtkIdType> &innerLoop,
    const std::vector<vtkIdType> &outerLoop,
    const double outward[3],
    vtkCellArray *cells,
    int &numDegenerate)
{
  numDegenerate = 0;

  if (points == nullptr || cells == nullptr)
  {
    fprintf(stderr,"Cannot stitch a cap annulus without points and an output cell array\n");
    return SV_ERROR;
  }
  if (innerLoop.size() < 3 || outerLoop.size() < 3)
  {
    fprintf(stderr,"Cannot stitch a cap annulus between rims of %zu and %zu points\n",
        innerLoop.size(), outerLoop.size());
    return SV_ERROR;
  }

  double center[3] = {0.0, 0.0, 0.0};
  for (size_t m = 0; m < innerLoop.size(); m++)
  {
    double p[3];
    points->GetPoint(innerLoop[m], p);
    for (int k = 0; k < 3; k++)
    {
      center[k] += p[k];
    }
  }
  for (int k = 0; k < 3; k++)
  {
    center[k] /= (double)innerLoop.size();
  }

  // A frame on the cap plane in which the angle increases counterclockwise
  // about the outward direction. Crossing with the axis the outward direction
  // leans on least keeps the first vector well away from degenerate.
  double axis[3] = {0.0, 0.0, 0.0};
  int smallest = 0;
  for (int k = 1; k < 3; k++)
  {
    if (std::abs(outward[k]) < std::abs(outward[smallest]))
    {
      smallest = k;
    }
  }
  axis[smallest] = 1.0;

  double u[3], v[3];
  vtkMath::Cross(axis, outward, u);
  if (vtkMath::Normalize(u) <= 0.0)
  {
    fprintf(stderr,"The cap outward direction is not a usable axis for stitching the annulus\n");
    return SV_ERROR;
  }
  vtkMath::Cross(outward, u, v);
  if (vtkMath::Normalize(v) <= 0.0)
  {
    fprintf(stderr,"The cap outward direction is not a usable axis for stitching the annulus\n");
    return SV_ERROR;
  }

  auto angleOf = [&](vtkIdType ptId)
  {
    double p[3], offset[3];
    points->GetPoint(ptId, p);
    vtkMath::Subtract(p, center, offset);
    return std::atan2(vtkMath::Dot(offset, v), vtkMath::Dot(offset, u));
  };

  auto wrap = [](double d)
  {
    while (d > vtkMath::Pi()) { d -= 2.0*vtkMath::Pi(); }
    while (d <= -vtkMath::Pi()) { d += 2.0*vtkMath::Pi(); }
    return d;
  };

  std::vector<vtkIdType> inner = innerLoop;
  std::vector<vtkIdType> outerRim = outerLoop;

  // Both rims have to run the same way round before they can be merged, and
  // each has to run round exactly once for the angle to order it at all.
  for (int side = 0; side < 2; side++)
  {
    std::vector<vtkIdType> &loop = (side == 0) ? inner : outerRim;
    double turning = 0.0;
    for (size_t m = 0; m < loop.size(); m++)
    {
      turning += wrap(angleOf(loop[(m+1)%loop.size()]) - angleOf(loop[m]));
    }
    if (std::abs(std::abs(turning) - 2.0*vtkMath::Pi()) > 0.5)
    {
      fprintf(stderr,"A cap rim of %zu points turns %.4g radians about the cap axis instead of one full turn, so it does not wind once around the vessel end and cannot be stitched by angle\n",
          loop.size(), turning);
      return SV_ERROR;
    }
    if (turning < 0.0)
    {
      std::reverse(loop.begin(), loop.end());
    }
  }

  size_t n = inner.size();
  size_t m = outerRim.size();

  // Start the outer rim at the point nearest in angle to where the inner rim
  // starts, so the first triangle is not a sliver spanning most of the cap.
  size_t startOuter = 0;
  double startAngle = angleOf(inner[0]);
  double bestGap = 0.0;
  for (size_t j = 0; j < m; j++)
  {
    double gap = std::abs(wrap(angleOf(outerRim[j]) - startAngle));
    if (j == 0 || gap < bestGap)
    {
      bestGap = gap;
      startOuter = j;
    }
  }

  // The angle swept from each rim's start, rescaled so both end at a full turn.
  // Rescaling matters because the two rims start a little apart in angle;
  // without it the merge would run one rim out before the other.
  std::vector<double> innerSweep(n+1, 0.0), outerSweep(m+1, 0.0);
  for (size_t k = 0; k < n; k++)
  {
    innerSweep[k+1] = innerSweep[k] + wrap(angleOf(inner[(k+1)%n]) - angleOf(inner[k]));
  }
  for (size_t k = 0; k < m; k++)
  {
    outerSweep[k+1] = outerSweep[k] +
        wrap(angleOf(outerRim[(startOuter+k+1)%m]) - angleOf(outerRim[(startOuter+k)%m]));
  }
  for (size_t k = 0; k <= n; k++)
  {
    innerSweep[k] /= innerSweep[n];
  }
  for (size_t k = 0; k <= m; k++)
  {
    outerSweep[k] /= outerSweep[m];
  }

  size_t i = 0, j = 0;
  while (i < n || j < m)
  {
    bool advanceInner;
    if (i >= n)
    {
      advanceInner = false;
    }
    else if (j >= m)
    {
      advanceInner = true;
    }
    else
    {
      advanceInner = (innerSweep[i+1] <= outerSweep[j+1]);
    }

    vtkIdType triangle[3];
    if (advanceInner)
    {
      triangle[0] = inner[i%n];
      triangle[1] = inner[(i+1)%n];
      triangle[2] = outerRim[(startOuter+j)%m];
      i++;
    }
    else
    {
      triangle[0] = outerRim[(startOuter+j)%m];
      triangle[1] = outerRim[(startOuter+j+1)%m];
      triangle[2] = inner[i%n];
      j++;
    }

    // The annulus is the end face of the wall, so it faces out of the vessel
    // end. Which of the two orders gives that depends on which rim was
    // advanced, so it is measured rather than worked out per case.
    double p0[3], p1[3], p2[3], e1[3], e2[3], normal[3];
    points->GetPoint(triangle[0], p0);
    points->GetPoint(triangle[1], p1);
    points->GetPoint(triangle[2], p2);
    vtkMath::Subtract(p1, p0, e1);
    vtkMath::Subtract(p2, p0, e2);
    vtkMath::Cross(e1, e2, normal);

    if (vtkMath::Norm(normal) <= 0.0)
    {
      // Dropping it would leave a hole in the wall's end face, which is worse
      // than a facet the volume mesher will complain about, so it is kept and
      // counted for the caller to report.
      numDegenerate++;
    }
    else if (vtkMath::Dot(normal, outward) < 0.0)
    {
      std::swap(triangle[1], triangle[2]);
    }

    cells->InsertNextCell(3, triangle);
  }

  return SV_OK;
}

// -------------------------------------
// TGenUtils_BuildWallShellSurface
// -------------------------------------
/**
 * @brief Builds the closed surface bounding the solid wall from the inner
 * surface and the trimmed offset surface, so the wall can be filled with
 * tetrahedra.
 * @note The wedge extrusion ties each outer node to exactly one inner node.
 * That tie is what forces the wall to be thinned at a junction: the outer nodes
 * of a concave crotch converge on each other and the only way to keep the
 * extrusion valid with a fixed node correspondence is to shorten it. Filling
 * the volume between two surfaces has no such tie, so this builds that volume's
 * boundary and leaves the filling to a volume mesher.
 *
 * The outer surface is no longer the inner one pushed along its normals. It is
 * the offset surface, contoured from a distance field and trimmed at the caps,
 * which shares no point with the inner surface and need not even have the same
 * number of points on a cap rim. The two are joined by an annulus at each
 * vessel end instead of by a strip per rim edge.
 *
 * The inner points are the first numPts points of the result, in the input
 * order and at the input coordinates, so the fluid/wall interface nodes are
 * carried through unchanged, which is what the solver requires of them. The
 * caller relies on that split as well: a facet of the result is interface,
 * outer wall or vessel end according to how many of its points fall below
 * numPts.
 * @param surface The inner surface.
 * @param outer The offset surface, already trimmed at the cap planes.
 * @param caps The rim pairs to close the wall between, one per vessel end;
 * empty when the inner surface is closed, in which case the result encloses the
 * lumen as well and the caller must mark it as a hole.
 * @param shell Set to the closed boundary of the wall.
 * @param numDegenerate Set to the number of zero-area annulus triangles.
 * @return SV_OK if the shell surface is built.
 */

int TGenUtils_BuildWallShellSurface(vtkPolyData *surface, vtkPolyData *outer,
    const std::vector<TGenUtilsCapRim> &caps, vtkPolyData *shell, int &numDegenerate)
{
  numDegenerate = 0;

  if (surface == nullptr || outer == nullptr || shell == nullptr)
  {
    fprintf(stderr,"Cannot build the wall shell without an inner surface, an offset surface and an output\n");
    return SV_ERROR;
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  vtkIdType numOuterPts = outer->GetNumberOfPoints();
  if (numPts == 0 || numOuterPts == 0)
  {
    fprintf(stderr,"Cannot build the wall shell from an inner surface of %lld points and an offset surface of %lld points\n",
        (long long)numPts, (long long)numOuterPts);
    return SV_ERROR;
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numPts + numOuterPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double p[3];
    surface->GetPoint(ptId, p);
    points->SetPoint(ptId, p);
  }
  for (vtkIdType ptId = 0; ptId < numOuterPts; ptId++)
  {
    double p[3];
    outer->GetPoint(ptId, p);
    points->SetPoint(numPts + ptId, p);
  }

  auto cells = vtkSmartPointer<vtkCellArray>::New();

  // The input normals point out of the lumen, so an inner triangle in its
  // input winding faces into the wall; reversing it makes it face out.
  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    surface->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }
    vtkIdType innerTriangle[3] = {pts[2], pts[1], pts[0]};
    cells->InsertNextCell(3, innerTriangle);
  }

  // Which way the contoured triangles face is a property of the contouring
  // filter rather than of this wall, and a shell wound inside out is filled
  // inside out. Measure it: an offset point lies away from the inner surface,
  // so the outward direction there is the direction from the nearest inner
  // point to it. A sample settles it, because the contour is wound
  // consistently; a sample that does not agree with itself means it is not,
  // which is worse than either answer and is reported rather than voted on.
  auto locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(surface);
  locator->BuildLocator();
  auto genericCell = vtkSmartPointer<vtkGenericCell>::New();

  vtkIdType numOuterCells = outer->GetNumberOfCells();
  const vtkIdType maxSamples = 5000;
  vtkIdType sampleStride = (numOuterCells > maxSamples) ? (numOuterCells/maxSamples) : 1;
  int numAgree = 0, numDisagree = 0;

  for (vtkIdType cellId = 0; cellId < numOuterCells; cellId += sampleStride)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    outer->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }

    double p0[3], p1[3], p2[3];
    outer->GetPoint(pts[0], p0);
    outer->GetPoint(pts[1], p1);
    outer->GetPoint(pts[2], p2);

    double e1[3], e2[3], normal[3];
    vtkMath::Subtract(p1, p0, e1);
    vtkMath::Subtract(p2, p0, e2);
    vtkMath::Cross(e1, e2, normal);
    if (vtkMath::Norm(normal) <= 0.0)
    {
      continue;
    }

    double centroid[3];
    for (int k = 0; k < 3; k++)
    {
      centroid[k] = (p0[k] + p1[k] + p2[k])/3.0;
    }

    double closest[3];
    vtkIdType closestCell = -1;
    int subId = 0;
    double distanceSquared = 0.0;
    locator->FindClosestPoint(centroid, closest, genericCell, closestCell, subId, distanceSquared);

    double away[3];
    vtkMath::Subtract(centroid, closest, away);
    if (vtkMath::Dot(normal, away) >= 0.0)
    {
      numAgree++;
    }
    else
    {
      numDisagree++;
    }
  }

  if (numAgree + numDisagree == 0)
  {
    fprintf(stderr,"The offset surface has no triangle with an area, so which way it faces cannot be measured\n");
    return SV_ERROR;
  }

  int numMajority = std::max(numAgree, numDisagree);
  if (numMajority < 0.9*(numAgree + numDisagree))
  {
    fprintf(stderr,"The offset surface faces outward on %d of its sampled triangles and inward on %d, so it is not wound consistently and cannot bound a volume\n",
        numAgree, numDisagree);
    return SV_ERROR;
  }

  bool reverseOuter = (numDisagree > numAgree);
  if (reverseOuter)
  {
    fprintf(stdout,"  the offset surface is wound facing the wall, so its triangles are reversed into the shell\n");
  }

  for (vtkIdType cellId = 0; cellId < numOuterCells; cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    outer->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }
    vtkIdType outerTriangle[3];
    if (reverseOuter)
    {
      outerTriangle[0] = pts[2] + numPts;
      outerTriangle[1] = pts[1] + numPts;
      outerTriangle[2] = pts[0] + numPts;
    }
    else
    {
      outerTriangle[0] = pts[0] + numPts;
      outerTriangle[1] = pts[1] + numPts;
      outerTriangle[2] = pts[2] + numPts;
    }
    cells->InsertNextCell(3, outerTriangle);
  }

  // Close each vessel end between its two rims. The offset rim ids are those of
  // the offset surface, so they move with it into the shell's numbering.
  for (size_t c = 0; c < caps.size(); c++)
  {
    std::vector<vtkIdType> outerLoop(caps[c].outerLoop.size());
    for (size_t m = 0; m < caps[c].outerLoop.size(); m++)
    {
      outerLoop[m] = caps[c].outerLoop[m] + numPts;
    }

    int numCapDegenerate = 0;
    if (TGenUtils_StitchCapAnnulus(points, caps[c].innerLoop, outerLoop,
          caps[c].outward, cells, numCapDegenerate) != SV_OK)
    {
      fprintf(stderr,"Problem closing the wall at the cap at (%.5g, %.5g, %.5g)\n",
          caps[c].origin[0], caps[c].origin[1], caps[c].origin[2]);
      return SV_ERROR;
    }
    numDegenerate += numCapDegenerate;
  }

  shell->Initialize();
  shell->SetPoints(points);
  shell->SetPolys(cells);
  shell->BuildLinks();

  return SV_OK;
}

// -------------------------------
// TGenUtils_FindLumenHolePoint
// -------------------------------
/**
 * @brief Finds a point strictly inside the region a closed surface encloses.
 * @note When the wall shell's inner surface is closed it encloses the lumen as
 * well as the wall, and the volume mesher has to be told that the lumen is not
 * part of the wall. That is done with a point inside it, which has to be found
 * rather than assumed: stepping a fixed distance inward from a surface point
 * leaves the lumen wherever the vessel is thinner than the step. Instead a ray
 * is cast inward along the normal and the midpoint of the first chord it cuts
 * is taken, which is inside the region for any vessel width. The point whose
 * chord is longest is used, so the result sits in the widest part of the model
 * and is the least sensitive to a ray that grazes the surface.
 * @param surface The closed surface; must have a 3-component 'Normals' point
 * array with the outward point normals.
 * @param holePoint Set to a point inside the enclosed region.
 * @return SV_OK if a point is found.
 */

int TGenUtils_FindLumenHolePoint(vtkPolyData *surface, double holePoint[3])
{
  if (surface == nullptr)
  {
    fprintf(stderr,"Cannot find a hole point without a surface\n");
    return SV_ERROR;
  }

  auto normals = surface->GetPointData()->GetArray("Normals");
  if (normals == nullptr || normals->GetNumberOfComponents() != 3)
  {
    fprintf(stderr,"The surface must have a 3-component 'Normals' point array to find a hole point\n");
    return SV_ERROR;
  }

  double bounds[6];
  surface->GetBounds(bounds);
  double dx = bounds[1]-bounds[0], dy = bounds[3]-bounds[2], dz = bounds[5]-bounds[4];
  double diagonal = std::sqrt(dx*dx + dy*dy + dz*dz);
  if (diagonal <= 0.0)
  {
    fprintf(stderr,"The surface is degenerate; cannot find a hole point\n");
    return SV_ERROR;
  }

  auto locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(surface);
  locator->BuildLocator();

  // Stepping off the surface before casting keeps the ray from immediately
  // hitting the triangles at its own origin.
  const double startOffset = 1.0e-6*diagonal;

  vtkIdType numPts = surface->GetNumberOfPoints();
  double bestChord = 0.0;
  bool found = false;

  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double n[3];
    normals->GetTuple(ptId, n);
    double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (len <= 0.0)
    {
      continue;
    }

    double p[3];
    surface->GetPoint(ptId, p);

    double start[3], end[3];
    for (int k = 0; k < 3; k++)
    {
      double inward = -n[k]/len;
      start[k] = p[k] + startOffset*inward;
      end[k] = p[k] + diagonal*inward;
    }

    double t = 0.0, hit[3], pcoords[3];
    int subId = 0;
    if (locator->IntersectWithLine(start, end, 0.0, t, hit, pcoords, subId) == 0)
    {
      continue;
    }

    double chord[3] = {hit[0]-start[0], hit[1]-start[1], hit[2]-start[2]};
    double chordLength = std::sqrt(chord[0]*chord[0] + chord[1]*chord[1] + chord[2]*chord[2]);
    if (chordLength <= bestChord)
    {
      continue;
    }
    bestChord = chordLength;
    for (int k = 0; k < 3; k++)
    {
      holePoint[k] = 0.5*(start[k] + hit[k]);
    }
    found = true;
  }

  if (!found)
  {
    fprintf(stderr,"No inward ray from the surface hit it again; cannot find a hole point\n");
    return SV_ERROR;
  }

  return SV_OK;
}

// -------------------------------------------
// TGenUtils_ReportAchievedWallThickness
// -------------------------------------------
/**
 * @brief Reports the wall thickness the extrusion actually achieves, measured
 * as the clearance from each outer point to the whole inner surface.
 * @note Every existing pass reasons about the extrusion length, not the wall
 * thickness. The rounding pass holds the outward distance from a point to its
 * own outer point at or above the assigned thickness; the fold prevention pass
 * tests whether an outer triangle inverts; the thickness reduction report
 * divides the final extrusion length by the requested one. None of these is
 * the thickness of the wall. The wall is only as thick as the closest approach
 * of the outer surface to *any* part of the inner surface, and at a concave
 * junction the outer point of one vessel moves toward the inner surface of the
 * other vessel, which no distance along a point's own normal and no triangle
 * orientation test can see. The existing report is therefore an upper bound:
 * it says how much of the requested extrusion length survived, not how much
 * wall was produced, so a junction can read as thinned to 90% while the wall
 * there is far thinner, or read as full thickness while the outer surface has
 * come within a fraction of it of the opposite side.
 *
 * The clearance is measured against the surface itself rather than against its
 * points, so the value does not depend on where the vertices happen to fall.
 * It also catches two distinct vessels whose walls interpenetrate without any
 * junction being involved, which is the same defect and is otherwise invisible.
 *
 * This is a report; neither the thickness array nor the geometry is modified,
 * although the achieved ratio is left on the surface as the
 * 'AchievedThicknessRatio' point array so the field can be viewed. It is meant
 * to be called after every thickness pass, on the final extrusion inputs.
 * @param surface The surface being extruded; must have a 3-component 'Normals'
 * point data array holding the final extrusion directions.
 * @param array The final extrusion length per point, one tuple per point.
 * @param requested The originally requested thickness per point, which the
 * achieved clearance is reported against.
 * @param label Names the extrusion the report belongs to.
 * @return SV_OK if the surface is reported on.
 */

int TGenUtils_ReportAchievedWallThickness(vtkPolyData *surface, vtkDoubleArray *array,
    const std::vector<double> &requested, const char *label)
{
  if (surface == nullptr || array == nullptr)
  {
    fprintf(stderr,"Cannot report the achieved wall thickness without a surface and a thickness array\n");
    return SV_ERROR;
  }

  if (label == nullptr)
  {
    label = "wall";
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array must have one component and one tuple per surface point\n");
    return SV_ERROR;
  }

  if ((vtkIdType)requested.size() != numPts)
  {
    fprintf(stderr,"The requested thickness must have one value per surface point\n");
    return SV_ERROR;
  }

  auto normals = surface->GetPointData()->GetArray("Normals");
  if (normals == nullptr || normals->GetNumberOfComponents() != 3 ||
      normals->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The surface must have a 3-component 'Normals' point array to report the achieved thickness\n");
    return SV_ERROR;
  }

  auto locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(surface);
  locator->BuildLocator();

  auto genericCell = vtkSmartPointer<vtkGenericCell>::New();

  auto ratioArray = vtkSmartPointer<vtkDoubleArray>::New();
  ratioArray->SetName("AchievedThicknessRatio");
  ratioArray->SetNumberOfComponents(1);
  ratioArray->SetNumberOfTuples(numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    ratioArray->SetValue(ptId, 1.0);
  }

  std::vector<double> achieved(numPts, 0.0);
  std::vector<std::pair<double,vtkIdType> > flagged;
  int numBelow90 = 0, numBelow50 = 0, numBelow25 = 0, numMeasured = 0;

  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double want = requested[ptId];
    if (want <= 0.0)
    {
      continue;
    }

    double n[3];
    normals->GetTuple(ptId, n);
    double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (len <= 0.0)
    {
      continue;
    }

    double p[3];
    surface->GetPoint(ptId, p);
    double extrusion = array->GetValue(ptId);
    double outer[3];
    for (int k = 0; k < 3; k++)
    {
      outer[k] = p[k] + extrusion*n[k]/len;
    }

    // The closest point of the inner surface to this outer point; its distance
    // is the wall thickness there, because the wall occupies the space between
    // the two surfaces and nothing is thicker than that closest approach.
    double closest[3];
    vtkIdType cellId = -1;
    int subId = 0;
    double distanceSquared = 0.0;
    locator->FindClosestPoint(outer, closest, genericCell, cellId, subId, distanceSquared);
    double clearance = std::sqrt(distanceSquared);

    achieved[ptId] = clearance;
    numMeasured++;
    double ratio = clearance/want;
    ratioArray->SetValue(ptId, ratio);
    if (ratio < 0.90)
    {
      numBelow90++;
      flagged.push_back(std::make_pair(ratio, ptId));
    }
    if (ratio < 0.50) { numBelow50++; }
    if (ratio < 0.25) { numBelow25++; }
  }

  surface->GetPointData()->RemoveArray("AchievedThicknessRatio");
  surface->GetPointData()->AddArray(ratioArray);

  fprintf(stdout,"Achieved thickness (outer point to the whole inner surface, vs requested) [%s]:\n", label);
  fprintf(stdout,"  points below 90%%/50%%/25%%: %d/%d/%d of %d measured\n",
      numBelow90, numBelow50, numBelow25, numMeasured);
  fprintf(stdout,"  this is the wall actually produced; the extrusion length reported separately is an upper bound on it\n");

  if (flagged.empty())
  {
    fprintf(stdout,"  every point achieves at least 90%% of its requested thickness\n");
    return SV_OK;
  }

  const int maxRegions = 8;
  const double radiusFraction = 0.02;
  std::vector<TGenUtilsPointRegion> regions;
  double regionRadius = 0.0;
  int numOutside = 0;
  int numRegionsTotal = 0;
  if (TGenUtils_ClusterPointsIntoRegions(surface, flagged, maxRegions, radiusFraction,
        regions, regionRadius, numOutside, numRegionsTotal) != SV_OK)
  {
    fprintf(stderr,"Problem clustering the points whose achieved wall thickness is below the requested one\n");
    return SV_ERROR;
  }

  fprintf(stdout,"  regions below 90%%: %d in total (separated by %.4g), worst %d shown:\n",
      numRegionsTotal, regionRadius, (int)regions.size());
  for (size_t i = 0; i < regions.size(); i++)
  {
    vtkIdType seedId = regions[i].seedId;
    double seed[3];
    surface->GetPoint(seedId, seed);
    fprintf(stdout,"    [%d] ratio %.3f (achieved %.5g / requested %.5g, extruded %.5g) at (%.5g, %.5g, %.5g), %d points\n",
        (int)(i+1), achieved[seedId]/requested[seedId], achieved[seedId], requested[seedId],
        array->GetValue(seedId), seed[0], seed[1], seed[2], regions[i].numPoints);
  }
  if (numOutside > 0)
  {
    fprintf(stdout,"    ... %d further points in the remaining %d regions\n",
        numOutside, numRegionsTotal - (int)regions.size());
  }

  return SV_OK;
}

// -----------------------------------------
// TGenUtils_ReportOffsetWallThickness
// -----------------------------------------
/**
 * @brief Reports the wall the offset surface actually produces, from both
 * sides.
 * @note The extrusion's report divides the clearance of each extruded point by
 * what was asked for, which needs every outer point to belong to one inner
 * point. The offset surface has no such relation, so the wall it makes has to
 * be measured as the distance between two surfaces, and the two directions of
 * that distance answer different questions.
 *
 * Outward, every point of the offset is at the requested distance from the
 * inner surface by construction, so measuring it checks the construction rather
 * than the wall: the grid the level set was contoured on has a spacing, and the
 * remesh that followed moved the points again. A shortfall here is that error,
 * and it is the one thing that can quietly eat the thickness this whole
 * approach exists to keep.
 *
 * Inward, the nearest offset point to an inner point is the wall over it. This
 * is the thickness in the sense that matters, and it is not the same number: at
 * a concave junction the offset creases outward, so the wall there comes out
 * thicker than requested, the way the outside of a welded joint fills with
 * material. A ratio below one on this side means wall is missing.
 * @param surface The inner surface.
 * @param array The requested thickness per inner surface point.
 * @param outer The trimmed offset surface.
 * @param label Names the wall in the report.
 * @return SV_OK if both directions were measured.
 */

int TGenUtils_ReportOffsetWallThickness(vtkPolyData *surface, vtkDoubleArray *array,
    vtkPolyData *outer, const char *label)
{
  if (surface == nullptr || array == nullptr || outer == nullptr)
  {
    fprintf(stderr,"Cannot report the offset wall thickness without both surfaces and a thickness array\n");
    return SV_ERROR;
  }

  if (label == nullptr)
  {
    label = "wall";
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array must have one component and one tuple per inner surface point\n");
    return SV_ERROR;
  }

  auto innerCells = vtkSmartPointer<vtkCellLocator>::New();
  innerCells->SetDataSet(surface);
  innerCells->BuildLocator();

  auto outerCells = vtkSmartPointer<vtkCellLocator>::New();
  outerCells->SetDataSet(outer);
  outerCells->BuildLocator();

  auto innerPoints = vtkSmartPointer<vtkStaticPointLocator>::New();
  innerPoints->SetDataSet(surface);
  innerPoints->BuildLocator();

  auto genericCell = vtkSmartPointer<vtkGenericCell>::New();

  fprintf(stdout,"Offset wall thickness [%s]:\n", label);

  // Outward: how far each offset point ended up from the inner surface against
  // the thickness asked for where it sits.
  {
    vtkIdType numOuterPts = outer->GetNumberOfPoints();
    std::vector<double> ratio((size_t)numOuterPts, 1.0);
    std::vector<std::pair<double,vtkIdType> > flagged;
    int numBelow90 = 0, numBelow50 = 0, numBelow25 = 0, numMeasured = 0;
    double worst = 0.0;

    for (vtkIdType ptId = 0; ptId < numOuterPts; ptId++)
    {
      double x[3];
      outer->GetPoint(ptId, x);

      vtkIdType nearest = innerPoints->FindClosestPoint(x);
      if (nearest < 0)
      {
        continue;
      }
      double want = array->GetValue(nearest);
      if (want <= 0.0)
      {
        continue;
      }

      double closest[3];
      vtkIdType cellId = -1;
      int subId = 0;
      double distanceSquared = 0.0;
      innerCells->FindClosestPoint(x, closest, genericCell, cellId, subId, distanceSquared);

      double value = std::sqrt(distanceSquared)/want;
      ratio[(size_t)ptId] = value;
      numMeasured++;
      if (numMeasured == 1 || value < worst)
      {
        worst = value;
      }
      if (value < 0.90) { numBelow90++; flagged.push_back(std::make_pair(value, ptId)); }
      if (value < 0.50) { numBelow50++; }
      if (value < 0.25) { numBelow25++; }
    }

    fprintf(stdout,"  offset surface to inner surface, over %d of its points: below 90%%/50%%/25%% of the requested thickness at %d/%d/%d, worst %.3f\n",
        numMeasured, numBelow90, numBelow50, numBelow25, worst);
    fprintf(stdout,"    this is the construction, not the shape: the level set puts every one of these points at the requested distance, so a shortfall is the grid spacing or the remesh giving it back\n");

    if (!flagged.empty())
    {
      const int maxRegions = 8;
      const double radiusFraction = 0.02;
      std::vector<TGenUtilsPointRegion> regions;
      double regionRadius = 0.0;
      int numOutside = 0, numRegionsTotal = 0;
      if (TGenUtils_ClusterPointsIntoRegions(outer, flagged, maxRegions, radiusFraction,
            regions, regionRadius, numOutside, numRegionsTotal) != SV_OK)
      {
        fprintf(stderr,"Problem clustering the offset points that fell short of the requested thickness\n");
        return SV_ERROR;
      }
      fprintf(stdout,"    regions below 90%%: %d in total (separated by %.4g), worst %d shown:\n",
          numRegionsTotal, regionRadius, (int)regions.size());
      for (size_t i = 0; i < regions.size(); i++)
      {
        double seed[3];
        outer->GetPoint(regions[i].seedId, seed);
        fprintf(stdout,"      [%d] ratio %.3f at (%.5g, %.5g, %.5g), %d points\n",
            (int)(i+1), ratio[(size_t)regions[i].seedId], seed[0], seed[1], seed[2],
            regions[i].numPoints);
      }
      if (numOutside > 0)
      {
        fprintf(stdout,"      ... %d further points in the remaining %d regions\n",
            numOutside, numRegionsTotal - (int)regions.size());
      }
    }
  }

  // Inward: the wall standing over each point of the fluid/wall interface.
  {
    std::vector<double> ratio((size_t)numPts, 1.0);
    std::vector<std::pair<double,vtkIdType> > flagged;
    std::vector<std::pair<double,vtkIdType> > flaggedThick;
    int numBelow90 = 0, numBelow50 = 0, numBelow25 = 0, numMeasured = 0;
    int numAbove2 = 0;
    double worst = 0.0, thickest = 0.0;

    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      double want = array->GetValue(ptId);
      if (want <= 0.0)
      {
        continue;
      }

      double p[3];
      surface->GetPoint(ptId, p);

      double closest[3];
      vtkIdType cellId = -1;
      int subId = 0;
      double distanceSquared = 0.0;
      outerCells->FindClosestPoint(p, closest, genericCell, cellId, subId, distanceSquared);

      double value = std::sqrt(distanceSquared)/want;
      ratio[(size_t)ptId] = value;
      numMeasured++;
      if (numMeasured == 1 || value < worst) { worst = value; }
      if (value > thickest) { thickest = value; }
      if (value < 0.90) { numBelow90++; flagged.push_back(std::make_pair(value, ptId)); }
      if (value < 0.50) { numBelow50++; }
      if (value < 0.25) { numBelow25++; }

      // The high tail is the only signal there is for two walls having merged.
      // A level set does not cross itself, so vessels closer together than
      // twice the wall never produce an error - the space between them simply
      // fills in and the mesh that comes out is valid and wrong. It reads here
      // as an interface carrying several times the wall it asked for, which a
      // junction crease also does, so the two are reported together by
      // location for the eye to separate.
      if (value > 2.0) { numAbove2++; flaggedThick.push_back(std::make_pair(-value, ptId)); }
    }

    fprintf(stdout,"  inner surface to offset surface, over %d of its points: below 90%%/50%%/25%% at %d/%d/%d, worst %.3f, thickest %.2fx requested\n",
        numMeasured, numBelow90, numBelow50, numBelow25, worst, thickest);
    fprintf(stdout,"    this is the wall over the interface. Above one at a junction is the crease filling it, which is intended; below one is wall that is missing\n");

    if (!flaggedThick.empty())
    {
      const int maxThickRegions = 8;
      const double thickRadiusFraction = 0.02;
      std::vector<TGenUtilsPointRegion> thickRegions;
      double thickRadius = 0.0;
      int numThickOutside = 0, numThickTotal = 0;
      if (TGenUtils_ClusterPointsIntoRegions(surface, flaggedThick, maxThickRegions,
            thickRadiusFraction, thickRegions, thickRadius, numThickOutside, numThickTotal) != SV_OK)
      {
        fprintf(stderr,"Problem clustering the interface points carrying more wall than requested\n");
        return SV_ERROR;
      }
      fprintf(stdout,"    %d points carry more than twice the wall asked for, in %d regions (separated by %.4g), thickest %d shown:\n",
          numAbove2, numThickTotal, thickRadius, (int)thickRegions.size());
      for (size_t i = 0; i < thickRegions.size(); i++)
      {
        vtkIdType seedId = thickRegions[i].seedId;
        double seed[3];
        surface->GetPoint(seedId, seed);
        fprintf(stdout,"      [%d] ratio %.2f (requested %.5g) at (%.5g, %.5g, %.5g), %d points\n",
            (int)(i+1), ratio[(size_t)seedId], array->GetValue(seedId),
            seed[0], seed[1], seed[2], thickRegions[i].numPoints);
      }
      if (numThickOutside > 0)
      {
        fprintf(stdout,"      ... %d further points in the remaining %d regions\n",
            numThickOutside, numThickTotal - (int)thickRegions.size());
      }
      fprintf(stdout,"      a region at a branch crotch is the crease and is expected; a region on the plain side of a vessel means the wall has reached another vessel passing close by, and the local wall thickness on those faces has to come below half the gap\n");
    }

    auto ratioArray = vtkSmartPointer<vtkDoubleArray>::New();
    ratioArray->SetName("OffsetThicknessRatio");
    ratioArray->SetNumberOfComponents(1);
    ratioArray->SetNumberOfTuples(numPts);
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      ratioArray->SetValue(ptId, ratio[(size_t)ptId]);
    }
    surface->GetPointData()->RemoveArray("OffsetThicknessRatio");
    surface->GetPointData()->AddArray(ratioArray);

    if (flagged.empty())
    {
      fprintf(stdout,"    every interface point carries at least 90%% of its requested wall\n");
      return SV_OK;
    }

    const int maxRegions = 8;
    const double radiusFraction = 0.02;
    std::vector<TGenUtilsPointRegion> regions;
    double regionRadius = 0.0;
    int numOutside = 0, numRegionsTotal = 0;
    if (TGenUtils_ClusterPointsIntoRegions(surface, flagged, maxRegions, radiusFraction,
          regions, regionRadius, numOutside, numRegionsTotal) != SV_OK)
    {
      fprintf(stderr,"Problem clustering the interface points whose wall fell short\n");
      return SV_ERROR;
    }
    fprintf(stdout,"    regions below 90%%: %d in total (separated by %.4g), worst %d shown:\n",
        numRegionsTotal, regionRadius, (int)regions.size());
    for (size_t i = 0; i < regions.size(); i++)
    {
      vtkIdType seedId = regions[i].seedId;
      double seed[3];
      surface->GetPoint(seedId, seed);
      fprintf(stdout,"      [%d] ratio %.3f (requested %.5g) at (%.5g, %.5g, %.5g), %d points\n",
          (int)(i+1), ratio[(size_t)seedId], array->GetValue(seedId),
          seed[0], seed[1], seed[2], regions[i].numPoints);
    }
    if (numOutside > 0)
    {
      fprintf(stdout,"      ... %d further points in the remaining %d regions\n",
          numOutside, numRegionsTotal - (int)regions.size());
    }
  }

  return SV_OK;
}

// -----------------------------------------
// TGenUtils_LimitThicknessToPreventFold
// -----------------------------------------
/**
 * @brief Reduces a wall thickness point array where extruding the surface
 * outward by the thickness would fold the outer wall over itself, so the
 * generated wall mesh does not self-intersect.
 * @note The curvature clamp (TGenUtils_ClampThicknessToConcaveCurvature) is
 * a local, one-ring estimate and can under-predict the fold at coarsely
 * meshed concave junctions, so this pass checks the actual extruded outer
 * geometry. The outer wall vertex of surface point p is p + t*n (thickness
 * along the unit point normal, matching how the extrusion normalizes the
 * warp vectors before scaling them). For each surface triangle the winding of the
 * outer triangle (from the extruded vertices) is compared with the winding
 * of the inner triangle: when the thickness is too large in a concave region
 * the outer triangle collapses and inverts, flipping the winding. The
 * thickness at the vertices of every inverted (or near-collapsed) triangle is
 * levelled to the smallest of the three (an imbalance between them folds the
 * triangle on its own and a proportional reduction would never remove it),
 * then reduced, and the check repeated until no triangle folds or the
 * iteration limit is reached. Reduction is bounded by a small fraction of each
 * point's original thickness, or by a fraction of the smallest altitude of the
 * triangles using the point where that is smaller: a sliver in the input
 * surface cannot carry a thickness of the order of its altitude no matter what
 * fraction of the requested value that is. Triangles still folded at that
 * bound are reported so the fold is surfaced rather than silently produced.
 * Only the thickness
 * values change; the surface points (the fluid/wall interface) never move.
 * This is a local test and does not detect a global collision of two
 * separate surface regions.
 * @param surface The surface being extruded; must have a 3-component
 * 'Normals' point data array with the outward point normals.
 * @param array The wall thickness point array; one component, one tuple per
 * surface point.
 * @param maxIterations The maximum number of reduce-and-recheck iterations.
 * @return SV_OK if the array is processed.
 */

int TGenUtils_LimitThicknessToPreventFold(vtkPolyData *surface, vtkDoubleArray *array, int maxIterations)
{
  if (surface == nullptr || array == nullptr)
  {
    fprintf(stderr,"Cannot limit a thickness array without a surface and an array\n");
    return SV_ERROR;
  }

  if (maxIterations <= 0)
  {
    return SV_OK;
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array must have one component and one tuple per surface point\n");
    return SV_ERROR;
  }

  auto normals = surface->GetPointData()->GetArray("Normals");
  if (normals == nullptr || normals->GetNumberOfComponents() != 3 ||
      normals->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The surface must have a 3-component 'Normals' point array to limit the thickness\n");
    return SV_ERROR;
  }

  // The cells are walked below with GetCellPoints, which needs the cell array
  // built. The caller is not required to have built it, and relying on an
  // earlier pass having done so would break if the passes are reordered.
  if (surface->NeedToBuildCells())
  {
    surface->BuildCells();
  }

  // A triangle is treated as folded when the outer winding has turned by
  // more than this much from the inner winding (a dot product of the unit
  // face normals at or below the threshold). A small positive value also
  // catches nearly collapsed outer triangles, not only fully inverted ones.
  const double foldThreshold = 0.1;
  // Each folded point's thickness is scaled by this factor per iteration.
  const double reductionFactor = 0.8;
  // The thickness is never reduced below this fraction of its original value,
  // unless the geometric bound below is smaller.
  const double minThicknessRatio = 0.05;
  // A triangle can only be extruded without folding while the outer vertices
  // move apart by less than the triangle's smallest altitude, so the thickness
  // is also allowed down to this fraction of that altitude. On a well shaped
  // triangle the smallest altitude is comparable with the edge lengths and
  // this bound sits above the ratio floor, which then decides; on a sliver
  // (a nearly degenerate triangle in the input surface) the altitude collapses
  // and the fold cannot be removed at the ratio floor at all, so the bound
  // lets the thickness go further down there and only there.
  const double minAltitudeRatio = 0.5;

  std::vector<double> originalThickness(numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    originalThickness[ptId] = array->GetValue(ptId);
  }

  // The smallest altitude of a triangle is twice its area over its longest
  // edge; the bound for a point is the smallest one over the triangles using
  // it. Points not used by any triangle keep the ratio floor.
  std::vector<double> altitudeBound(numPts, std::numeric_limits<double>::max());
  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    surface->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }

    double corner[3][3];
    for (int i = 0; i < 3; i++)
    {
      surface->GetPoint(pts[i], corner[i]);
    }

    double maxEdge = 0.0;
    for (int i = 0; i < 3; i++)
    {
      const double *a = corner[i];
      const double *b = corner[(i+1)%3];
      double edge = std::sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) +
          (b[2]-a[2])*(b[2]-a[2]));
      if (edge > maxEdge)
      {
        maxEdge = edge;
      }
    }
    if (maxEdge <= 0.0)
    {
      continue;
    }

    double ab[3], ac[3], cross[3];
    for (int k = 0; k < 3; k++)
    {
      ab[k] = corner[1][k] - corner[0][k];
      ac[k] = corner[2][k] - corner[0][k];
    }
    cross[0] = ab[1]*ac[2] - ab[2]*ac[1];
    cross[1] = ab[2]*ac[0] - ab[0]*ac[2];
    cross[2] = ab[0]*ac[1] - ab[1]*ac[0];
    double area = 0.5*std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
    double bound = minAltitudeRatio*2.0*area/maxEdge;

    for (int i = 0; i < 3; i++)
    {
      if (bound < altitudeBound[pts[i]])
      {
        altitudeBound[pts[i]] = bound;
      }
    }
  }

  std::vector<double> minThickness(numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    minThickness[ptId] = minThicknessRatio*originalThickness[ptId];
    if (altitudeBound[ptId] < minThickness[ptId])
    {
      minThickness[ptId] = altitudeBound[ptId];
    }
  }

  // The extrusion moves each point by thickness*unit normal (the vmtk
  // boundary layer generator normalizes the warp vectors before scaling them
  // by the thickness), so the point normals are normalized here to compute
  // the same outer vertex. The stored normals are not always unit vectors:
  // the surface point normals are averaged where coincident points are
  // merged, which shortens them. Points with a degenerate normal cannot be
  // extruded in any direction and their thickness cannot fix a fold, so the
  // triangles using them are left out of the check.
  std::vector<double> unitNormals(3*numPts, 0.0);
  std::vector<bool> validNormal(numPts, false);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double normal[3];
    normals->GetTuple(ptId, normal);
    double length = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
        normal[2]*normal[2]);
    if (length <= 0.0)
    {
      continue;
    }
    unitNormals[3*ptId] = normal[0]/length;
    unitNormals[3*ptId+1] = normal[1]/length;
    unitNormals[3*ptId+2] = normal[2]/length;
    validNormal[ptId] = true;
  }

  // Returns the unit normal of a triangle from three points, or false if
  // the triangle is degenerate.
  auto triangleNormal = [](const double a[3], const double b[3],
      const double c[3], double normal[3]) -> bool
  {
    double ab[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    double ac[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    normal[0] = ab[1]*ac[2] - ab[2]*ac[1];
    normal[1] = ab[2]*ac[0] - ab[0]*ac[2];
    normal[2] = ab[0]*ac[1] - ab[1]*ac[0];
    double length = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
        normal[2]*normal[2]);
    if (length <= 0.0)
    {
      return false;
    }
    normal[0] /= length;
    normal[1] /= length;
    normal[2] /= length;
    return true;
  };

  std::vector<vtkIdType> foldedCells;
  int iter = 0;
  for (; iter < maxIterations; iter++)
  {
    std::vector<bool> foldedPoint(numPts, false);
    // The thickness a folded point is pulled down to before it is scaled: the
    // smallest thickness on the folded triangles using it. The outer vertices
    // of a triangle move apart by the difference of the thicknesses as well as
    // by the spread of the normals, so an imbalance between the three points
    // folds the triangle on its own. Scaling all three by the same factor keeps
    // that imbalance forever, which is why a fold driven by it survives every
    // iteration; levelling the three first removes it in one step.
    std::vector<double> foldTarget(numPts, std::numeric_limits<double>::max());
    foldedCells.clear();

    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId, npts, pts);
      if (npts != 3)
      {
        continue;
      }

      if (!validNormal[pts[0]] || !validNormal[pts[1]] || !validNormal[pts[2]])
      {
        continue;
      }

      double inner[3][3];
      double outer[3][3];
      for (int i = 0; i < 3; i++)
      {
        surface->GetPoint(pts[i], inner[i]);
        const double *normal = &unitNormals[3*pts[i]];
        double thickness = array->GetValue(pts[i]);
        outer[i][0] = inner[i][0] + thickness*normal[0];
        outer[i][1] = inner[i][1] + thickness*normal[1];
        outer[i][2] = inner[i][2] + thickness*normal[2];
      }

      double innerNormal[3];
      double outerNormal[3];
      if (!triangleNormal(inner[0], inner[1], inner[2], innerNormal))
      {
        continue;
      }
      // A degenerate outer triangle means the thickness has collapsed the
      // face, which is itself a fold.
      bool outerOk = triangleNormal(outer[0], outer[1], outer[2], outerNormal);
      double dot = outerOk ? (innerNormal[0]*outerNormal[0] +
          innerNormal[1]*outerNormal[1] + innerNormal[2]*outerNormal[2]) : -1.0;

      if (dot <= foldThreshold)
      {
        foldedCells.push_back(cellId);
        double smallest = array->GetValue(pts[0]);
        for (int i = 1; i < 3; i++)
        {
          smallest = std::min(smallest, array->GetValue(pts[i]));
        }
        for (int i = 0; i < 3; i++)
        {
          foldedPoint[pts[i]] = true;
          foldTarget[pts[i]] = std::min(foldTarget[pts[i]], smallest);
        }
      }
    }

    if (foldedCells.empty())
    {
      break;
    }

    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      if (!foldedPoint[ptId])
      {
        continue;
      }
      double reduced = reductionFactor*std::min(array->GetValue(ptId), foldTarget[ptId]);
      if (reduced < minThickness[ptId])
      {
        reduced = minThickness[ptId];
      }
      array->SetValue(ptId, reduced);
    }
  }

  if (!foldedCells.empty())
  {
    fprintf(stderr,"Warning: the extruded outer wall still folds over at %lld triangles after %d\
 thickness reduction iterations; the wall mesh may self-intersect there. Refine the surface mesh or\
 reduce the wall thickness at the junction\n", (long long)foldedCells.size(), iter);

    // Report where the triangles that could not be fixed are, what shape they
    // have and how far apart their point normals are, so a fold the thickness
    // cannot fix can be told apart from one the thickness caused. Reducing the
    // thickness moves the outer vertices by at most twice the thickness, which
    // can only turn the outer face normal enough to register as folded when
    // that displacement is comparable with the size of the triangle. So a
    // triangle still folded at the reduction bound is either nearly
    // degenerate, which shows up as a minimum altitude (twice the area over
    // the longest edge) far below its edge lengths, or it has point normals
    // pointing in very different directions, which shows up as a small
    // smallest normal dot product. Neither is fixed by the thickness: the
    // first needs the surface remeshed there, the second needs the normals
    // at that junction fixed.
    const size_t maxReported = 10;
    for (size_t i = 0; i < foldedCells.size() && i < maxReported; i++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(foldedCells[i], npts, pts);

      double center[3] = {0.0, 0.0, 0.0};
      double corner[3][3];
      for (int j = 0; j < 3; j++)
      {
        surface->GetPoint(pts[j], corner[j]);
        for (int k = 0; k < 3; k++)
        {
          center[k] += corner[j][k]/3.0;
        }
      }

      double edges[3];
      double maxEdge = 0.0;
      for (int j = 0; j < 3; j++)
      {
        const double *a = corner[j];
        const double *b = corner[(j+1)%3];
        edges[j] = std::sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) +
            (b[2]-a[2])*(b[2]-a[2]));
        if (edges[j] > maxEdge)
        {
          maxEdge = edges[j];
        }
      }

      // The area from the cross product of two edges, and the altitude on the
      // longest edge; a nearly degenerate triangle has an altitude orders of
      // magnitude below its edges.
      double ab[3], ac[3], cross[3];
      for (int k = 0; k < 3; k++)
      {
        ab[k] = corner[1][k] - corner[0][k];
        ac[k] = corner[2][k] - corner[0][k];
      }
      cross[0] = ab[1]*ac[2] - ab[2]*ac[1];
      cross[1] = ab[2]*ac[0] - ab[0]*ac[2];
      cross[2] = ab[0]*ac[1] - ab[1]*ac[0];
      double area = 0.5*std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
      double altitude = (maxEdge > 0.0) ? 2.0*area/maxEdge : 0.0;

      // The smallest dot product between the point normals; a value near one
      // means the three points are extruded in nearly the same direction.
      double minNormalDot = 1.0;
      for (int j = 0; j < 3; j++)
      {
        const double *a = &unitNormals[3*pts[j]];
        const double *b = &unitNormals[3*pts[(j+1)%3]];
        double normalDot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
        if (normalDot < minNormalDot)
        {
          minNormalDot = normalDot;
        }
      }

      fprintf(stderr,"  folded triangle at (%.6g, %.6g, %.6g):\n", center[0], center[1], center[2]);
      fprintf(stderr,"    edges %.6g %.6g %.6g, area %.6g, altitude %.6g\n", edges[0], edges[1],
          edges[2], area, altitude);
      fprintf(stderr,"    thickness %.6g %.6g %.6g (requested %.6g %.6g %.6g)\n",
          array->GetValue(pts[0]), array->GetValue(pts[1]), array->GetValue(pts[2]),
          originalThickness[pts[0]], originalThickness[pts[1]], originalThickness[pts[2]]);
      fprintf(stderr,"    reduction floor %.6g %.6g %.6g\n", minThickness[pts[0]],
          minThickness[pts[1]], minThickness[pts[2]]);
      fprintf(stderr,"    smallest point normal dot %.6g\n", minNormalDot);
    }
  }

  return SV_OK;
}

// -----------------------------------------
// TGenUtils_RoundOuterWallToPreserveThickness
// -----------------------------------------
/**
 * @brief Rounds the outer wall surface outward at concave junctions so the
 * assigned wall thickness is preserved there instead of the wall being thinned
 * (which shows up as an inward depression).
 * @note The wall is built by extruding the inner surface outward by the
 * thickness along the point normals. At a concave junction (the crotch where
 * two vessels merge) the outward normals converge, so the naive outer surface
 * (each point at its thickness along its normal) self-intersects even though
 * every point sits at the full thickness. The fold prevention pass removes
 * that self-intersection by thinning the wall, which then reaches less far
 * outward and caves in. This instead keeps the thickness and moves the outer
 * surface outward into a smooth convex fillet, the way the outer side of a
 * thick welded junction fills with material rather than denting inward.
 *
 * The inner surface (the fluid/wall interface) is fixed and never moves; only
 * the outer surface points move. Each outer point is relaxed toward the
 * average of its one-ring neighbors' outer points (which fills a dip because a
 * dip's neighbors sit further out), in proportion to how concave the point is,
 * so convex and flat points and a straight tube are left unchanged. It is then
 * pushed back out so its outward (normal) distance from the inner point is
 * never below the assigned thickness, and capped so a very sharp crotch cannot
 * spike outward without bound.
 *
 * This distance is measured from each point's own inner point, so it is the
 * extrusion length and not the thickness of the wall: at a junction the outer
 * point of one vessel can close up against the inner surface of the vessel it
 * is merging with while its own normal distance still reads full. Enforcing
 * the thickness against the whole inner surface was tried and removed. It is
 * the right invariant but it cannot be reached from this representation: with
 * one outer node tied to each inner node, no placement satisfies it wherever
 * the thickness exceeds the concave radius of curvature, so the pass pushes
 * points off their normals every iteration without ever converging and skews
 * the elements instead. Measured over four variants it bought about 35% of the
 * thickness deficit and cost an order of magnitude on the worst element aspect
 * ratio. TGenUtils_ReportAchievedWallThickness still measures the real
 * thickness, so the deficit stays visible; closing it needs the node
 * correspondence dropped (a shell filled with tetrahedra, or layers terminated
 * locally rather than thinned), not a stronger constraint here.
 * Boundary (cap rim) points are pinned so the
 * wall stays flat at the caps. The rounded outer surface is encoded back into
 * the normals (the extrusion direction) and the thickness array (the extrusion
 * magnitude) so the existing extrusion reproduces exactly this surface; the
 * tangle test in the extrusion uses the same inverted/collapsed-triangle
 * criterion, so a fold-free rounded surface leaves it nothing to undo.
 *
 * A degenerate junction triangle (an input sliver) still cannot carry a wall
 * in any direction, so this only fills the junction depression; the following
 * fold prevention pass remains the safety net that thins a sliver fold.
 * @param surface The surface being extruded; must have a 3-component 'Normals'
 * point data array with the outward point normals.
 * @param array The wall thickness point array (the assigned thickness on
 * entry); overwritten with the achieved outer distance, at least the assigned
 * thickness.
 * @param iterations The number of relaxation iterations.
 * @param relaxation The fraction of the neighbor-average move applied to a
 * fully concave point per iteration (between 0 and 1).
 * @param maxFilletRatio The largest multiple of the assigned thickness the
 * outer surface may bulge out to.
 * @return SV_OK if the outer surface is rounded.
 */

int TGenUtils_RoundOuterWallToPreserveThickness(vtkPolyData *surface, vtkDoubleArray *array,
    int iterations, double relaxation, double maxFilletRatio)
{
  if (iterations <= 0 || relaxation <= 0.0)
  {
    return SV_OK;
  }

  if (surface == nullptr || array == nullptr)
  {
    fprintf(stderr,"Cannot round the outer wall without a surface and a thickness array\n");
    return SV_ERROR;
  }

  vtkIdType numPts = surface->GetNumberOfPoints();
  if (array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The thickness array must have one component and one tuple per surface point\n");
    return SV_ERROR;
  }

  auto normals = surface->GetPointData()->GetArray("Normals");
  if (normals == nullptr || normals->GetNumberOfComponents() != 3 ||
      normals->GetNumberOfTuples() != numPts)
  {
    fprintf(stderr,"The surface must have a 3-component 'Normals' point array to round the outer wall\n");
    return SV_ERROR;
  }

  if (maxFilletRatio < 1.0)
  {
    maxFilletRatio = 1.0;
  }

  // Unit point normals; the stored normals are averaged where coincident points
  // merge and are not always unit length. A degenerate normal cannot be
  // extruded, so those points are left untouched.
  std::vector<double> unitNormals(3*numPts, 0.0);
  std::vector<bool> validNormal(numPts, false);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double n[3];
    normals->GetTuple(ptId, n);
    double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (len <= 0.0)
    {
      continue;
    }
    unitNormals[3*ptId]   = n[0]/len;
    unitNormals[3*ptId+1] = n[1]/len;
    unitNormals[3*ptId+2] = n[2]/len;
    validNormal[ptId] = true;
  }

  // One-ring neighbors and the boundary (cap rim) points, built as in the
  // warp-vector smoothing: a boundary edge is used by a single cell, and its
  // endpoints are pinned so the wall stays flat at the caps.
  surface->BuildLinks();
  std::vector<std::vector<vtkIdType>> neighbors(numPts);
  std::vector<char> pinned(numPts, 0);
  auto cellIds = vtkSmartPointer<vtkIdList>::New();
  auto edgeNeighbors = vtkSmartPointer<vtkIdList>::New();
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    surface->GetPointCells(ptId, cellIds);
    auto& ptNeighbors = neighbors[ptId];
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
      vtkIdType cellId = cellIds->GetId(i);
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId, npts, pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (pts[j] == ptId)
        {
          continue;
        }
        if (std::find(ptNeighbors.begin(), ptNeighbors.end(), pts[j]) == ptNeighbors.end())
        {
          ptNeighbors.push_back(pts[j]);
        }
        surface->GetCellEdgeNeighbors(cellId, ptId, pts[j], edgeNeighbors);
        if (edgeNeighbors->GetNumberOfIds() == 0)
        {
          pinned[ptId] = 1;
        }
      }
    }
  }

  // Per-point concavity weight (the average sine of the rise angle over the
  // neighbors above the tangent plane), zero on convex and flat points, so only
  // concave junctions are rounded. The assigned thickness is captured now
  // because the array is overwritten with the achieved distance at the end.
  std::vector<double> weight(numPts, 0.0);
  std::vector<double> assignedThickness(numPts, 0.0);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    assignedThickness[ptId] = array->GetValue(ptId);
    if (pinned[ptId] || !validNormal[ptId] || neighbors[ptId].empty())
    {
      continue;
    }
    const double *n = &unitNormals[3*ptId];
    double p[3];
    surface->GetPoint(ptId, p);
    double concavitySum = 0.0;
    int concaveCount = 0;
    for (auto neighborId : neighbors[ptId])
    {
      double q[3];
      surface->GetPoint(neighborId, q);
      double offset[3] = {q[0]-p[0], q[1]-p[1], q[2]-p[2]};
      double distance = std::sqrt(offset[0]*offset[0] + offset[1]*offset[1] + offset[2]*offset[2]);
      if (distance <= 0.0)
      {
        continue;
      }
      double height = offset[0]*n[0] + offset[1]*n[1] + offset[2]*n[2];
      if (height <= 0.0)
      {
        continue;
      }
      concavitySum += height/distance;   // sine of the rise angle, in [0,1)
      concaveCount++;
    }
    if (concaveCount > 0)
    {
      weight[ptId] = concavitySum/concaveCount;
    }
  }

  // Outer surface positions, initialized to the naive extrusion (each point at
  // its assigned thickness along its normal); Jacobi updates use a second
  // buffer so the result does not depend on the point visiting order.
  std::vector<double> outer(3*numPts);
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double p[3];
    surface->GetPoint(ptId, p);
    const double *n = &unitNormals[3*ptId];
    double t = assignedThickness[ptId];
    outer[3*ptId]   = p[0] + t*n[0];
    outer[3*ptId+1] = p[1] + t*n[1];
    outer[3*ptId+2] = p[2] + t*n[2];
  }

  std::vector<double> nextOuter(outer);
  for (int iter = 0; iter < iterations; iter++)
  {
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      // Pinned, convex/flat, degenerate, or isolated points keep the naive
      // outer position, so only concave junctions move.
      if (pinned[ptId] || !validNormal[ptId] || weight[ptId] <= 0.0 || neighbors[ptId].empty())
      {
        nextOuter[3*ptId]   = outer[3*ptId];
        nextOuter[3*ptId+1] = outer[3*ptId+1];
        nextOuter[3*ptId+2] = outer[3*ptId+2];
        continue;
      }

      // Move toward the average of the neighbors' outer points (a dip's
      // neighbors sit further out, so this fills the dip), scaled by concavity.
      double centroid[3] = {0.0, 0.0, 0.0};
      for (auto neighborId : neighbors[ptId])
      {
        centroid[0] += outer[3*neighborId];
        centroid[1] += outer[3*neighborId+1];
        centroid[2] += outer[3*neighborId+2];
      }
      double inv = 1.0/(double)neighbors[ptId].size();
      centroid[0] *= inv; centroid[1] *= inv; centroid[2] *= inv;

      double blend = relaxation*weight[ptId];
      double moved[3];
      for (int k = 0; k < 3; k++)
      {
        moved[k] = outer[3*ptId+k] + blend*(centroid[k] - outer[3*ptId+k]);
      }

      // Keep the thickness: the outward (normal) distance from the inner point
      // must not drop below the assigned thickness, and the fillet is capped so
      // a very sharp crotch cannot spike outward without bound.
      double p[3];
      surface->GetPoint(ptId, p);
      const double *n = &unitNormals[3*ptId];
      double disp[3] = {moved[0]-p[0], moved[1]-p[1], moved[2]-p[2]};
      double h = disp[0]*n[0] + disp[1]*n[1] + disp[2]*n[2];
      double t = assignedThickness[ptId];
      if (h < t)
      {
        double push = t - h;
        moved[0] += push*n[0];
        moved[1] += push*n[1];
        moved[2] += push*n[2];
        h = t;
      }

      // The fillet is capped so a very sharp crotch cannot spike outward
      // without bound.
      double maxH = maxFilletRatio*t;
      if (h > maxH)
      {
        double pull = h - maxH;
        moved[0] -= pull*n[0];
        moved[1] -= pull*n[1];
        moved[2] -= pull*n[2];
      }

      nextOuter[3*ptId]   = moved[0];
      nextOuter[3*ptId+1] = moved[1];
      nextOuter[3*ptId+2] = moved[2];
    }
    outer.swap(nextOuter);
  }

  // Encode the rounded outer surface back into the extrusion inputs: the normal
  // is the unit direction to the outer point and the thickness is the distance
  // to it, so the existing extrusion places the outer node exactly here. Points
  // that were not moved reproduce their original normal and thickness. The
  // distance is euclidean, matching the wedge edge length the extrusion builds.
  int numRaised = 0;
  double maxRatio = 1.0;
  vtkIdType maxRatioId = -1;
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    if (!validNormal[ptId])
    {
      continue;
    }
    double p[3];
    surface->GetPoint(ptId, p);
    double disp[3] = {outer[3*ptId]-p[0], outer[3*ptId+1]-p[1], outer[3*ptId+2]-p[2]};
    double dist = std::sqrt(disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2]);
    if (dist <= 0.0)
    {
      continue;
    }
    double unit[3] = {disp[0]/dist, disp[1]/dist, disp[2]/dist};
    normals->SetTuple(ptId, unit);
    array->SetValue(ptId, dist);

    double t = assignedThickness[ptId];
    if (t > 0.0)
    {
      double ratio = dist/t;
      if (ratio > 1.001)
      {
        numRaised++;
      }
      if (ratio > maxRatio)
      {
        maxRatio = ratio;
        maxRatioId = ptId;
      }
    }
  }

  fprintf(stdout,"Wall outer rounding: filled the junction depression by raising %d concave points; "
      "largest fillet %.3gx the assigned thickness", numRaised, maxRatio);
  if (maxRatioId >= 0)
  {
    double p[3];
    surface->GetPoint(maxRatioId, p);
    fprintf(stdout," at (%.5g, %.5g, %.5g)", p[0], p[1], p[2]);
  }
  fprintf(stdout,"\n");

  return SV_OK;
}

// -----------------------------------------
// TGenUtils_ReportSurfaceTriangleQuality
// -----------------------------------------
/**
 * @brief Reports the shape quality of the triangles of a surface.
 * @note A wall thickness cannot be extruded off a triangle without folding it
 * unless the outer vertices move apart by less than the triangle's smallest
 * altitude, so a sliver (a triangle whose altitude has collapsed against its
 * edge lengths) cannot carry a wall at all and no thickness value fixes it.
 * The fold prevention pass (TGenUtils_LimitThicknessToPreventFold) can only
 * thin the wall down to such a triangle; the triangle itself has to be gone
 * before the extrusion. This report is called on both sides of the surface
 * remeshing so a sliver reaching the extrusion can be traced to the model
 * surface or to the remeshing. Nothing is modified.
 * @param surface The surface to report on.
 * @param label A short name for the pipeline stage, printed with the report.
 * @return SV_OK if the surface is reported on.
 */

int TGenUtils_ReportSurfaceTriangleQuality(vtkPolyData *surface, const char *label)
{
  if (surface == nullptr)
  {
    fprintf(stderr,"Cannot report the triangle quality of a null surface\n");
    return SV_ERROR;
  }

  if (surface->NeedToBuildCells())
  {
    surface->BuildCells();
  }

  // The aspect ratio of a triangle is its longest edge over its smallest
  // altitude, scaled so an equilateral triangle is 1.0, matching the
  // convention of the tetrahedron quality report.
  const double equilateralScale = 0.5*std::sqrt(3.0);
  // Triangles at or above this aspect ratio are listed individually.
  const double reportAspect = 10.0;
  const size_t maxReported = 5;

  vtkIdType numTris = 0;
  double minAspect = 0.0;
  double maxAspect = 0.0;
  double sumAspect = 0.0;
  vtkIdType numAbove10 = 0;
  vtkIdType numAbove30 = 0;
  vtkIdType numAbove100 = 0;
  vtkIdType numDegenerate = 0;

  // The worst triangles, kept sorted by decreasing aspect ratio.
  std::vector<std::pair<double,vtkIdType> > worst;

  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    surface->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }

    double corner[3][3];
    for (int i = 0; i < 3; i++)
    {
      surface->GetPoint(pts[i], corner[i]);
    }

    double maxEdge = 0.0;
    for (int i = 0; i < 3; i++)
    {
      const double *a = corner[i];
      const double *b = corner[(i+1)%3];
      double edge = std::sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) +
          (b[2]-a[2])*(b[2]-a[2]));
      if (edge > maxEdge)
      {
        maxEdge = edge;
      }
    }

    double ab[3], ac[3], cross[3];
    for (int k = 0; k < 3; k++)
    {
      ab[k] = corner[1][k] - corner[0][k];
      ac[k] = corner[2][k] - corner[0][k];
    }
    cross[0] = ab[1]*ac[2] - ab[2]*ac[1];
    cross[1] = ab[2]*ac[0] - ab[0]*ac[2];
    cross[2] = ab[0]*ac[1] - ab[1]*ac[0];
    double area = 0.5*std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

    // A triangle with no area or no extent has no aspect ratio to report; it
    // is counted separately because it is a defect on its own.
    if (maxEdge <= 0.0 || area <= 0.0)
    {
      numDegenerate++;
      continue;
    }

    double altitude = 2.0*area/maxEdge;
    double aspect = equilateralScale*maxEdge/altitude;

    if (numTris == 0 || aspect < minAspect)
    {
      minAspect = aspect;
    }
    if (aspect > maxAspect)
    {
      maxAspect = aspect;
    }
    sumAspect += aspect;
    numTris++;

    if (aspect > 10.0)
    {
      numAbove10++;
    }
    if (aspect > 30.0)
    {
      numAbove30++;
    }
    if (aspect > 100.0)
    {
      numAbove100++;
    }

    if (aspect >= reportAspect)
    {
      worst.push_back(std::make_pair(aspect, cellId));
      std::sort(worst.begin(), worst.end(),
          std::greater<std::pair<double,vtkIdType> >());
      if (worst.size() > maxReported)
      {
        worst.resize(maxReported);
      }
    }
  }

  fprintf(stdout,"Surface triangle quality (%s), aspect ratio (1.0 is an equilateral triangle):\n",
      (label == nullptr) ? "surface" : label);
  if (numTris == 0)
  {
    fprintf(stdout,"  No triangles\n");
    return SV_OK;
  }

  fprintf(stdout,"  Number of triangles: %lld\n", (long long)numTris);
  fprintf(stdout,"  Min / Avg / Max: %.3f / %.3f / %.3f\n", minAspect,
      sumAspect/(double)numTris, maxAspect);
  fprintf(stdout,"  Aspect ratio > 10: %lld, > 30: %lld, > 100: %lld\n",
      (long long)numAbove10, (long long)numAbove30, (long long)numAbove100);
  if (numDegenerate > 0)
  {
    fprintf(stdout,"  Zero area triangles: %lld\n", (long long)numDegenerate);
  }

  for (size_t i = 0; i < worst.size(); i++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    surface->GetCellPoints(worst[i].second, npts, pts);

    double center[3] = {0.0, 0.0, 0.0};
    double corner[3][3];
    for (int j = 0; j < 3; j++)
    {
      surface->GetPoint(pts[j], corner[j]);
      for (int k = 0; k < 3; k++)
      {
        center[k] += corner[j][k]/3.0;
      }
    }

    double edges[3];
    for (int j = 0; j < 3; j++)
    {
      const double *a = corner[j];
      const double *b = corner[(j+1)%3];
      edges[j] = std::sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) +
          (b[2]-a[2])*(b[2]-a[2]));
    }

    fprintf(stdout,"  aspect ratio %.3f at (%.6g, %.6g, %.6g), edges %.6g %.6g %.6g\n",
        worst[i].first, center[0], center[1], center[2], edges[0], edges[1], edges[2]);
  }

  return SV_OK;
}
