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

#include "simvascular_tetgen.h"

#include "sv_polydatasolid_utils.h"
#include "sv_misc_utils.h"
#include "sv_vtk_utils.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
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
