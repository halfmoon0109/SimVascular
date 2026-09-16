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
#include <chrono>
#include <cmath>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
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

  // Every array read below has to be there. TetGen fills adjtetlist only
  // when it was run with neighout above one, so a caller that left that
  // unset would be read through a null pointer here.
  if (outmesh->pointlist == nullptr || (numPolys > 0 && outmesh->tetrahedronlist == nullptr) ||
      (numFaces > 0 && (outmesh->trifacelist == nullptr || outmesh->adjtetlist == nullptr)))
  {
    fprintf(stderr,"TetGen output is missing an array the conversion needs (points %s, tetrahedra %s, faces %s, face neighbours %s); the mesher has to be run with neighout = 2\n",
        outmesh->pointlist ? "yes" : "no", outmesh->tetrahedronlist ? "yes" : "no",
        outmesh->trifacelist ? "yes" : "no", outmesh->adjtetlist ? "yes" : "no");
    return SV_ERROR;
  }

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

  // The name that was asked for. This wrote to a fixed 'out.vtp' whatever it
  // was given, so every caller shared one file and the last write of a run was
  // the only one that survived - the two wall diagnostics landed on top of each
  // other, and the file named in the log never existed.
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(PData);
#else
  writer->SetInputData(PData);
#endif
  //writer->SetDataModeToAscii();

  if (!writer->Write())
  {
    fprintf(stderr,"Could not write the polydata to %s\n", filename);
    return SV_ERROR;
  }

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
        double p[3];
        surface->GetPoint(a, p);
        fprintf(stderr,"Point %lld at (%.5g, %.5g, %.5g) starts more than one boundary edge, so the surface boundary is not a set of simple loops and cannot be capped\n",
            (long long)a, p[0], p[1], p[2]);
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
        // The chain ran into a point that ends a boundary edge without starting
        // one. That is a dangling boundary, and what leaves one is a triangle
        // meeting the rest of the surface at a point rather than along an edge:
        // the edges out of the pinch are shared by two cells and so read as
        // interior. The coordinates are what locate it; the length says whether
        // the chain was a rim that nearly closed or a stray sliver.
        double s[3], e[3];
        surface->GetPoint(start, s);
        surface->GetPoint(current, e);
        fprintf(stderr,"The boundary edge chain from point %lld at (%.5g, %.5g, %.5g) ends after %zu points at point %lld at (%.5g, %.5g, %.5g) instead of closing\n",
            (long long)start, s[0], s[1], s[2], loop.size(),
            (long long)current, e[0], e[1], e[2]);
        return SV_ERROR;
      }
      current = step->second;
    }

    if (current != start)
    {
      double p[3];
      surface->GetPoint(current, p);
      fprintf(stderr,"A boundary edge chain closed onto point %lld at (%.5g, %.5g, %.5g) rather than onto its start %lld\n",
          (long long)current, p[0], p[1], p[2], (long long)start);
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
// RemoveBoundaryFlaps
// -------------------------------------
/**
 * @brief Removes triangles that hang from an edge already shared by the rest
 * of a surface, until none are left.
 * @note Two boundary edges alone do not make a triangle a flap. A valid open
 * surface can have an ear joined to the rest by one edge; that edge has one
 * neighbour and the two boundary edges remain part of a closed loop. The
 * dangling case has two boundary edges and a third edge with at least two
 * neighbours. The third edge was already interior before the extra triangle
 * was attached, so its two boundary edges form a chain that cannot close
 * through the underlying surface.
 *
 * Removing that extra triangle is safe in a way that removing a boundary ear
 * is not: the neighbours across its third edge are still there and that edge
 * remains interior.
 *
 * Removing a flap can expose another behind it, so this repeats until the
 * surface is clean, and reports how much it took. A handful is the contour
 * leaving a near-degenerate triangle where the grid could not resolve the level
 * set - one such triangle, of zero area and a whole cell long, is what sent the
 * rim walk off the end of a three point chain. Thousands would mean the surface
 * is not a surface, which the count is there to make visible.
 * @param surface The surface, repaired in place; its cells must be triangles.
 * @param numRemoved Set to the number of triangles removed.
 * @param firstRemoved Set to the centre of the first triangle removed, for the
 * log; untouched if none were.
 * @return SV_OK if the surface is left with triangles on it.
 */

static int RemoveBoundaryFlaps(vtkPolyData *surface, int &numRemoved, double firstRemoved[3])
{
  numRemoved = 0;

  // Every non-empty pass removes at least one triangle, so this must settle in
  // no more passes than there were triangles on entry.
  auto edgeNeighbors = vtkSmartPointer<vtkIdList>::New();

  while (true)
  {
    surface->BuildLinks();

    std::vector<bool> drop((size_t)surface->GetNumberOfCells(), false);
    int numDropped = 0;

    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId, npts, pts);
      if (npts != 3)
      {
        continue;
      }

      int numBoundaryEdges = 0;
      int numOverusedEdges = 0;
      for (vtkIdType j = 0; j < npts; j++)
      {
        surface->GetCellEdgeNeighbors(cellId, pts[j], pts[(j+1)%npts], edgeNeighbors);
        if (edgeNeighbors->GetNumberOfIds() == 0)
        {
          numBoundaryEdges++;
        }
        else if (edgeNeighbors->GetNumberOfIds() > 1)
        {
          numOverusedEdges++;
        }
      }

      if (numBoundaryEdges == 2 && numOverusedEdges == 1)
      {
        if (numRemoved + numDropped == 0)
        {
          for (int k = 0; k < 3; k++)
          {
            firstRemoved[k] = 0.0;
          }
          for (vtkIdType j = 0; j < npts; j++)
          {
            double p[3];
            surface->GetPoint(pts[j], p);
            for (int k = 0; k < 3; k++)
            {
              firstRemoved[k] += p[k]/3.0;
            }
          }
        }
        drop[(size_t)cellId] = true;
        numDropped++;
      }
    }

    if (numDropped == 0)
    {
      break;
    }

    auto keptCells = vtkSmartPointer<vtkCellArray>::New();
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
    {
      if (drop[(size_t)cellId])
      {
        continue;
      }
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId, npts, pts);
      if (npts != 3)
      {
        continue;
      }
      vtkIdType triangle[3] = {pts[0], pts[1], pts[2]};
      keptCells->InsertNextCell(3, triangle);
    }

    auto repaired = vtkSmartPointer<vtkPolyData>::New();
    repaired->SetPoints(surface->GetPoints());
    repaired->SetPolys(keptCells);
    surface->DeepCopy(repaired);

    numRemoved += numDropped;

    if (surface->GetNumberOfCells() == 0)
    {
      fprintf(stderr,"Removing %d triangles hanging off the boundary, the first at (%.5g, %.5g, %.5g), left no surface at all\n",
          numRemoved, firstRemoved[0], firstRemoved[1], firstRemoved[2]);
      return SV_ERROR;
    }
  }

  if (numRemoved > 0)
  {
    // The removed triangles took the last use of some points with them, and an
    // unused point is one the volume mesher renumbers its input around.
    std::vector<vtkIdType> newPointId((size_t)surface->GetNumberOfPoints(), -1);
    auto keptPoints = vtkSmartPointer<vtkPoints>::New();
    auto keptCells = vtkSmartPointer<vtkCellArray>::New();

    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId, npts, pts);
      if (npts != 3)
      {
        continue;
      }

      vtkIdType triangle[3];
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (newPointId[(size_t)pts[j]] < 0)
        {
          double p[3];
          surface->GetPoint(pts[j], p);
          newPointId[(size_t)pts[j]] = keptPoints->InsertNextPoint(p);
        }
        triangle[j] = newPointId[(size_t)pts[j]];
      }
      keptCells->InsertNextCell(3, triangle);
    }

    auto compacted = vtkSmartPointer<vtkPolyData>::New();
    compacted->SetPoints(keptPoints);
    compacted->SetPolys(keptCells);
    surface->DeepCopy(compacted);
  }

  return SV_OK;
}

// -------------------------------------
// LabelConnectedShells
// -------------------------------------
/**
 * @brief Labels each triangle of a surface with the connected shell it belongs
 * to, two triangles being connected when they share a point.
 * @note This is what vtkPolyDataConnectivityFilter computes, and reading the
 * labels off that filter is how it was done first. The filter writes them into
 * a 'RegionId' array, and the array it was read from - the one on the cells -
 * was not there in the VTK this is built against, so the pass that needed the
 * labels stopped instead of running. Which arrays that filter attaches, under
 * which extraction mode, is a detail of a library this cannot be tested
 * against; the labelling itself is a breadth-first walk over the cells and
 * costs less to write than the question costs to answer.
 * @param surface The surface; its cells must be triangles.
 * @param cellShell Set to the shell index of each cell, or -1 for a cell that
 * is not a triangle.
 * @return The number of shells.
 */

static int LabelConnectedShells(vtkPolyData *surface, std::vector<int> &cellShell)
{
  vtkIdType numCells = surface->GetNumberOfCells();
  cellShell.assign((size_t)numCells, -1);

  surface->BuildLinks();

  auto pointCells = vtkSmartPointer<vtkIdList>::New();
  std::deque<vtkIdType> queue;
  int numShells = 0;

  for (vtkIdType seed = 0; seed < numCells; seed++)
  {
    if (cellShell[(size_t)seed] >= 0)
    {
      continue;
    }

    int shell = numShells++;
    cellShell[(size_t)seed] = shell;
    queue.clear();
    queue.push_back(seed);

    while (!queue.empty())
    {
      vtkIdType cellId = queue.front();
      queue.pop_front();

      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId, npts, pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        surface->GetPointCells(pts[j], pointCells);
        for (vtkIdType n = 0; n < pointCells->GetNumberOfIds(); n++)
        {
          vtkIdType neighbor = pointCells->GetId(n);
          if (cellShell[(size_t)neighbor] < 0)
          {
            cellShell[(size_t)neighbor] = shell;
            queue.push_back(neighbor);
          }
        }
      }
    }
  }

  return numShells;
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
 * @param maxFieldBytes The memory the distance field may occupy; the spacing is
 * coarsened until the grid fits in it and the outcome is reported. It bounds
 * the field alone, not the contour or the locators built alongside it.
 * @param maxThicknessSlope The bound the caller has already imposed on how fast
 * the thickness may change along the surface. The band is sized from it, so a
 * caller that does not enforce it will see the band's own check fire.
 * @param outer Set to the contoured offset surface.
 * @param gridSpacing Set to the spacing the grid was built at, which may be
 * coarser than the one asked for. It bounds how accurately the surface locates
 * the level set, so it is what any later pass has to compare its own tolerance
 * against rather than assuming the surface is exact.
 * @return SV_OK if the offset surface is built.
 */

int TGenUtils_BuildOffsetOuterSurface(vtkPolyData *surface, vtkDoubleArray *array,
    double targetEdgeSize, double maxFieldBytes, double maxThicknessSlope, vtkPolyData *outer,
    double &gridSpacing)
{
  gridSpacing = 0.0;

  // A voxel costs a float of distance and a byte of bookkeeping, so how many
  // voxels a byte budget buys is a property of this function rather than of its
  // caller. The two arrays are declared below; working the count out anywhere
  // else leaves arithmetic elsewhere to be corrected by hand the day one of
  // them changes type, and nothing would report that it had not been.
  const size_t bytesPerVoxel = sizeof(float) + sizeof(unsigned char);
  vtkIdType maxVoxels = (vtkIdType)(maxFieldBytes/(double)bytesPerVoxel);
  if (maxVoxels < 8)
  {
    fprintf(stderr,"A field budget of %.0f bytes leaves room for %lld voxels, which is not a grid\n",
        maxFieldBytes, (long long)maxVoxels);
    return SV_ERROR;
  }

  if (maxThicknessSlope < 0.0 || maxThicknessSlope >= 1.0)
  {
    fprintf(stderr,"The wall thickness slope bound is %.5g; the band around the distance field is only finite for a bound below one\n",
        maxThicknessSlope);
    return SV_ERROR;
  }

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

  // The band a triangle marks reaches past its own edges, so the grid has to
  // clear the model by more than the band or a row would start inside it.
  double maxEdgeLength = 0.0;
  for (vtkIdType cellId = 0; cellId < closed->GetNumberOfCells(); cellId++)
  {
    vtkIdType npts;
    const vtkIdType *pts;
    closed->GetCellPoints(cellId, npts, pts);
    if (npts != 3)
    {
      continue;
    }
    for (vtkIdType j = 0; j < npts; j++)
    {
      double a[3], b[3];
      closedPoints->GetPoint(pts[j], a);
      closedPoints->GetPoint(pts[(j+1)%npts], b);
      maxEdgeLength = std::max(maxEdgeLength, std::sqrt(vtkMath::Distance2BetweenPoints(a, b)));
    }
  }

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
    // Two cells wider than the widest band a triangle can mark, so a row of the
    // grid always starts outside the band and the sign it starts with is known.
    margin = (thicknessMax + maxThicknessSlope*maxEdgeLength + 2.0*spacing)/(1.0 - maxThicknessSlope)
        + 2.0*spacing;
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
  gridSpacing = spacing;

  fprintf(stdout,"Wall outer surface by distance field offset:\n");
  fprintf(stdout,"  thickness %.5g to %.5g, %zu cap rims closed with a fan\n",
      thicknessMin, thicknessMax, loops.size());
  fprintf(stdout,"  grid %d x %d x %d = %zu voxels at spacing %.5g, which is %.2f cells per smallest wall thickness\n",
      dims[0], dims[1], dims[2], numVoxels, spacing, thicknessMin/spacing);
  if (spacing > requestedSpacing*1.001)
  {
    fprintf(stdout,"  the spacing was coarsened from %.5g to fit the %lld voxel budget (%.0f MB of field)\n",
        requestedSpacing, (long long)maxVoxels, (double)maxVoxels*bytesPerVoxel/1.0e6);
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
  const unsigned char kBand = 1;

  // This pass and the distance evaluation after it are the two that scale with
  // the grid, and how long they take decides whether the offset is affordable
  // on a larger model than the one it was written against. Timing them costs
  // nothing and answers that from the log rather than from a guess.
  auto markStart = std::chrono::steady_clock::now();

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
    double corners[3][3];
    for (vtkIdType j = 0; j < npts; j++)
    {
      double p[3];
      closedPoints->GetPoint(pts[j], p);
      for (int k = 0; k < 3; k++)
      {
        if (j == 0 || p[k] < lo[k]) { lo[k] = p[k]; }
        if (j == 0 || p[k] > hi[k]) { hi[k] = p[k]; }
        corners[j][k] = p[k];
      }
      radius = std::max(radius, closedThickness[(size_t)pts[j]]);
    }

    double longestEdge = 0.0;
    for (int j = 0; j < 3; j++)
    {
      longestEdge = std::max(longestEdge,
          std::sqrt(vtkMath::Distance2BetweenPoints(corners[j], corners[(j+1)%3])));
    }

    // How far the band has to reach is not the thickness on this triangle. A
    // voxel takes its thickness from whichever surface point is nearest, which
    // need not be one of these three, and the level set sits at that thickness.
    //
    // The gradation limiter has already bounded how fast the thickness can
    // change along the surface, at slope s. A point the band can reach is
    // within r + e of this triangle, so the thickness there is at most
    // t + s(r + e), and the band has to hold that plus a cell on each side:
    //
    //   r >= t + s(r + e) + 2h   ->   r >= (t + s*e + 2h) / (1 - s)
    //
    // Sizing the band that way makes the escape below impossible for a
    // thickness field that varies only along the surface. It can still escape
    // where two parts of the surface pass close to each other carrying very
    // different thicknesses, since the limiter bounds the slope along the
    // surface and not across the gap, which is why the check stays.
    //
    // The slope is taken from the caller rather than written in here, because
    // the caller is what enforces it and a band derived from a different number
    // than the one enforced is a band that is quietly too thin.
    radius = (radius + maxThicknessSlope*longestEdge + 2.0*spacing)/(1.0 - maxThicknessSlope);

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

  double markSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - markStart).count();

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
  auto evaluateStart = std::chrono::steady_clock::now();

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
  // a sheet between them. Which one a voxel is follows from the last band voxel
  // the row passed through: the band is closed around the model and at least
  // two cells thick, so a row cannot get from one side to the other without
  // crossing it, and the sign it crossed on is the sign it carries until it
  // crosses again. Each row starts outside, because the margin clears the model
  // by more than the band.
  //
  // Sweeping the rows this way rather than flooding the grid from a corner
  // costs one sequential pass instead of a breadth-first walk over a hundred
  // million voxels, and it removes the question of whether a seed corner was
  // available to start from.
  size_t numInside = 0;
  for (int k = 0; k < dims[2]; k++)
  {
    for (int j = 0; j < dims[1]; j++)
    {
      size_t index = (size_t)dims[0]*((size_t)j + (size_t)dims[1]*(size_t)k);
      bool inside = false;
      for (int i = 0; i < dims[0]; i++, index++)
      {
        if (state[index] == kBand)
        {
          inside = (values->GetValue((vtkIdType)index) < 0.0);
          continue;
        }
        if (inside)
        {
          values->SetValue((vtkIdType)index, (float)(-farValue));
          numInside++;
        }
        else
        {
          values->SetValue((vtkIdType)index, (float)farValue);
        }
      }
    }
  }

  double evaluateSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - evaluateStart).count();

  fprintf(stdout,"  distance evaluated at %zu band voxels; %zu enclosed and %zu outside voxels filled by sign\n",
      numBandVoxels, numInside, numVoxels - numBandVoxels - numInside);
  fprintf(stdout,"  %.1f s marking the band, %.1f s evaluating and filling it\n",
      markSeconds, evaluateSeconds);

  // The signs the sweep wrote are only right if the level set stayed inside the
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
          if ((value < 0.0) != (values->GetValue((vtkIdType)neighbor) < 0.0))
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
    fprintf(stderr,"The offset level set leaves the evaluated band at %zu voxels, so the band around the surface is too thin to hold it. The band is sized so that a thickness varying along the surface within the gradation limit cannot do this, which leaves the case it cannot cover: two parts of the surface passing within a wall of each other while carrying very different thicknesses. Evening out the local wall thickness between those faces is what fixes it.\n",
        numStraddling);
    return SV_ERROR;
  }

  auto image = vtkSmartPointer<vtkImageData>::New();
  image->SetDimensions(dims[0], dims[1], dims[2]);
  image->SetOrigin(origin[0], origin[1], origin[2]);
  image->SetSpacing(spacing, spacing, spacing);
  image->GetPointData()->SetScalars(values);

  // The bookkeeping has done its work by here, and holding it through the
  // contour would put its fifth of the budget alongside the contour's own
  // output at the one moment this function is at its largest.
  std::vector<unsigned char>().swap(state);

  auto contourStart = std::chrono::steady_clock::now();

  auto contour = vtkSmartPointer<vtkFlyingEdges3D>::New();
  contour->SetInputData(image);
  contour->SetValue(0, 0.0);
  contour->ComputeNormalsOff();
  contour->ComputeGradientsOff();
  contour->ComputeScalarsOff();
  contour->Update();

  fprintf(stdout,"  %.1f s contouring\n",
      std::chrono::duration<double>(std::chrono::steady_clock::now() - contourStart).count());

  if (contour->GetOutput()->GetNumberOfCells() == 0)
  {
    fprintf(stderr,"The offset level set is empty, so no outer wall surface was produced\n");
    return SV_ERROR;
  }

  // Marching cubes puts each vertex where the level set crosses a grid edge, so
  // wherever the level set passes close to a grid corner the vertices coming
  // off the edges that meet there arrive on top of each other and the triangle
  // between them has almost no area. On a surface of this size that is not a
  // possibility but a certainty, and the triangles are the kind a remesher
  // rejects rather than repairs. Collapsing anything closer together than a
  // hundredth of a cell removes them; the surface does not move, because the
  // points being merged were already at the same place to that tolerance, and
  // no hole is left behind, because a triangle whose vertices merge is an edge
  // collapse rather than a deletion.
  auto degenerateCleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  degenerateCleaner->SetInputData(contour->GetOutput());
  degenerateCleaner->ToleranceIsAbsoluteOn();
  degenerateCleaner->SetAbsoluteTolerance(0.01*spacing);
  degenerateCleaner->Update();

  vtkIdType numContourCells = contour->GetOutput()->GetNumberOfCells();

  outer->Initialize();
  outer->SetPoints(degenerateCleaner->GetOutput()->GetPoints());
  outer->SetPolys(degenerateCleaner->GetOutput()->GetPolys());

  if (outer->GetNumberOfCells() == 0)
  {
    fprintf(stderr,"Removing the degenerate contour triangles left no outer wall surface\n");
    return SV_ERROR;
  }

  fprintf(stdout,"  %lld of the %lld contour triangles were degenerate and collapsed\n",
      (long long)(numContourCells - outer->GetNumberOfCells()), (long long)numContourCells);

  // The outer wall is one shell per connected piece of the inner surface,
  // because dilating a connected set leaves it connected. So anything past that
  // count is debris the grid shed where it could not hold the level set, and
  // the count itself is worth reporting either way.
  //
  // Debris does not stay harmless. Each fragment is a closed shell until the
  // cap trim cuts it, and then it is a rim belonging to no vessel end - which
  // is how it was first noticed, as a boundary walk that would not close. It is
  // dropped here, and how much of the contour goes with it is the measure of
  // how badly the grid resolved the level set.
  //
  // Each shell is checked against the sign of the distance as well. The level
  // set is d = t with t positive, and d is negative inside the model, so every
  // point of every shell has to lie outside it; a shell that does not means the
  // field or the sign its unevaluated voxels were filled with is wrong, which
  // is worth seeing in the log rather than meshing.
  std::vector<int> innerShell;
  int numInnerRegions = LabelConnectedShells(surface, innerShell);

  std::vector<int> cellShell;
  int numShells = LabelConnectedShells(outer, cellShell);

  fprintf(stdout,"  offset surface has %lld points and %lld triangles in %d connected shells,"
      " against %d connected regions of the inner surface\n",
      (long long)outer->GetNumberOfPoints(), (long long)outer->GetNumberOfCells(),
      numShells, numInnerRegions);

  if (numShells > numInnerRegions)
  {
    // One shell is entirely on one side of the model, so a handful of its
    // points settle which side that is; they are taken from the cells the shell
    // happens to be labelled on first, since there is nothing to choose between
    // them.
    const int maxSamples = 8;
    std::vector<vtkIdType> shellCells((size_t)numShells, 0);
    std::vector<double> shellCentroid((size_t)numShells*3, 0.0);
    std::vector<std::vector<vtkIdType> > shellSamples((size_t)numShells);

    for (vtkIdType cellId = 0; cellId < outer->GetNumberOfCells(); cellId++)
    {
      int shell = cellShell[(size_t)cellId];
      if (shell < 0 || shell >= numShells)
      {
        continue;
      }
      vtkIdType npts;
      const vtkIdType *pts;
      outer->GetCellPoints(cellId, npts, pts);
      if (npts != 3)
      {
        continue;
      }
      shellCells[(size_t)shell]++;
      double p[3];
      outer->GetPoint(pts[0], p);
      for (int k = 0; k < 3; k++)
      {
        shellCentroid[(size_t)shell*3 + k] += p[k];
      }
      if ((int)shellSamples[(size_t)shell].size() < maxSamples)
      {
        shellSamples[(size_t)shell].push_back(pts[0]);
      }
    }

    std::vector<int> outsideShells;
    std::vector<int> insideShells;
    for (int shell = 0; shell < numShells; shell++)
    {
      if (shellCells[(size_t)shell] == 0)
      {
        continue;
      }
      for (int k = 0; k < 3; k++)
      {
        shellCentroid[(size_t)shell*3 + k] /= (double)shellCells[(size_t)shell];
      }

      int numOutsideSamples = 0;
      for (size_t s = 0; s < shellSamples[(size_t)shell].size(); s++)
      {
        double p[3];
        outer->GetPoint(shellSamples[(size_t)shell][s], p);
        if (outwardSign*implicit->EvaluateFunction(p) > 0.0)
        {
          numOutsideSamples++;
        }
      }
      if (2*numOutsideSamples > (int)shellSamples[(size_t)shell].size())
      {
        outsideShells.push_back(shell);
      }
      else
      {
        insideShells.push_back(shell);
      }
    }

    // Largest first, so that keeping one shell per region of the inner surface
    // is keeping the ones that are the wall rather than the ones that are dust.
    std::sort(outsideShells.begin(), outsideShells.end(),
        [&shellCells](int a, int b) { return shellCells[(size_t)a] > shellCells[(size_t)b]; });

    if ((int)outsideShells.size() < numInnerRegions)
    {
      fprintf(stderr,"The offset surface has %zu shells outside the model but the inner surface has %d regions, so at least one region has no outer wall at all\n",
          outsideShells.size(), numInnerRegions);
      return SV_ERROR;
    }

    std::vector<int> kept(outsideShells.begin(), outsideShells.begin() + numInnerRegions);
    std::vector<int> dropped(outsideShells.begin() + numInnerRegions, outsideShells.end());
    dropped.insert(dropped.end(), insideShells.begin(), insideShells.end());

    vtkIdType numDroppedCells = 0;
    for (size_t d = 0; d < dropped.size(); d++)
    {
      numDroppedCells += shellCells[(size_t)dropped[d]];
    }

    fprintf(stdout,"  keeping %zu shells (%lld triangles) as the outer wall and dropping %zu (%lld triangles, %.2f%% of the contour):\n",
        kept.size(), (long long)(outer->GetNumberOfCells() - numDroppedCells),
        dropped.size(), (long long)numDroppedCells,
        100.0*(double)numDroppedCells/(double)outer->GetNumberOfCells());

    if (!insideShells.empty())
    {
      fprintf(stdout,"  WARNING: %zu of them lie inside the model, which the level set d = t cannot produce. The distance field or the sign its unevaluated voxels were filled with is wrong, and the surface being kept is only as trustworthy as that.\n",
          insideShells.size());
    }

    const size_t maxListed = 12;
    for (size_t d = 0; d < dropped.size() && d < maxListed; d++)
    {
      int shell = dropped[d];
      bool isInside = std::find(insideShells.begin(), insideShells.end(), shell) != insideShells.end();
      fprintf(stdout,"    dropped %lld triangles %s the model at (%.5g, %.5g, %.5g)\n",
          (long long)shellCells[(size_t)shell], isInside ? "inside" : "outside",
          shellCentroid[(size_t)shell*3], shellCentroid[(size_t)shell*3+1],
          shellCentroid[(size_t)shell*3+2]);
    }
    if (dropped.size() > maxListed)
    {
      fprintf(stdout,"    ... %zu further dropped shells\n", dropped.size() - maxListed);
    }

    // Carrying the kept cells over by hand, for the same reason the labelling
    // is done by hand: the labels are this function's own, so there is no
    // filter to hand them back to.
    std::vector<bool> keepShell((size_t)numShells, false);
    for (size_t s = 0; s < kept.size(); s++)
    {
      keepShell[(size_t)kept[s]] = true;
    }

    // The points of the dropped shells have to go with them. An unused point is
    // one the volume mesher renumbers its input around, and the wall tagging
    // downstream depends on that numbering - so the points are renumbered here,
    // where the map from old to new is still in hand.
    std::vector<vtkIdType> newPointId((size_t)outer->GetNumberOfPoints(), -1);
    auto keptPoints = vtkSmartPointer<vtkPoints>::New();
    auto keptCells = vtkSmartPointer<vtkCellArray>::New();

    for (vtkIdType cellId = 0; cellId < outer->GetNumberOfCells(); cellId++)
    {
      int shell = cellShell[(size_t)cellId];
      if (shell < 0 || !keepShell[(size_t)shell])
      {
        continue;
      }
      vtkIdType npts;
      const vtkIdType *pts;
      outer->GetCellPoints(cellId, npts, pts);
      if (npts != 3)
      {
        continue;
      }

      vtkIdType triangle[3];
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (newPointId[(size_t)pts[j]] < 0)
        {
          double p[3];
          outer->GetPoint(pts[j], p);
          newPointId[(size_t)pts[j]] = keptPoints->InsertNextPoint(p);
        }
        triangle[j] = newPointId[(size_t)pts[j]];
      }
      keptCells->InsertNextCell(3, triangle);
    }

    auto keptSurface = vtkSmartPointer<vtkPolyData>::New();
    keptSurface->SetPoints(keptPoints);
    keptSurface->SetPolys(keptCells);
    outer->DeepCopy(keptSurface);

    if (outer->GetNumberOfCells() == 0)
    {
      fprintf(stderr,"Dropping the shells that are not the outer wall left no offset surface\n");
      return SV_ERROR;
    }

    fprintf(stdout,"  the outer wall is %lld points and %lld triangles\n",
        (long long)outer->GetNumberOfPoints(), (long long)outer->GetNumberOfCells());
  }

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
 * @param thickness The wall thickness per point of the inner surface, which is
 * what the dome over each cap was offset by. May be null, in which case every
 * cap falls back to maxThickness.
 * @param maxThickness The largest wall thickness, used for a cap whose own
 * thickness is not available.
 * @param caps Set to one entry per vessel end, holding both rims and the plane.
 * @return SV_OK if every cap was trimmed and its two rims paired.
 */

int TGenUtils_TrimOffsetSurfaceAtCaps(vtkPolyData *surface, vtkPolyData *outer,
    vtkDoubleArray *thickness, double maxThickness,
    std::vector<TGenUtilsCapRim> &caps)
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
  std::vector<double> capThickness(innerLoops.size(), maxThickness);
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

    // The dome over this cap stands as far off the end as the wall is thick
    // there, and the thickness at the rim is what the fan across the opening
    // was offset by. The largest thickness anywhere in the model says nothing
    // about this end and, on a model whose walls range over an order of
    // magnitude, sizes the cut for the thickest vessel at every cap.
    if (thickness != nullptr)
    {
      double local = 0.0;
      for (size_t m = 0; m < loop.size(); m++)
      {
        if (loop[m] >= 0 && loop[m] < thickness->GetNumberOfTuples())
        {
          local = std::max(local, thickness->GetValue(loop[m]));
        }
      }
      if (local > 0.0)
      {
        capThickness[c] = local;
      }
    }
  }

  // Clip once per cap. The scalar is positive on everything that is kept: below
  // the plane, or far enough from this rim that the plane has no business
  // reaching it. The dome is the only place both are negative.
  //
  // One clip against the smallest of the per-cap scalars would produce the same
  // set of kept points and save thirteen passes over the surface, but not the
  // same cut. The clip places a cut point by interpolating the scalar along an
  // edge linearly, and that is exact for one cap - 'below' is a linear function
  // of position, so a cut governed by it lands on the cap plane itself. The
  // smallest of several is concave rather than linear, and interpolating it
  // linearly underestimates it, so on any edge where the cap that gives the
  // smallest value changes between its two ends the cut lands short by up to an
  // edge length. That is a rim pulled off its cap plane, and the pass below
  // pairs rims to caps by how well they lie on one.
  //
  // The dome of a rim of radius R under a wall of thickness t meets the cap
  // plane at R + t, so the window has to hold that and no more than it needs
  // to. Sizing it off the thickness rather than off R keeps it tight where the
  // vessel is wide, and keeps it valid where the wall is thick relative to the
  // vessel - a window of a fixed multiple of R would fall inside the rim it is
  // meant to cut once t approached R. The t is this cap's own: a window built
  // from the model's thickest wall reaches the same distance past every rim,
  // and past a thin vessel's cap that is most of the way to its neighbour.
  for (size_t c = 0; c < caps.size(); c++)
  {
    const TGenUtilsCapRim &cap = caps[c];
    double window = capRadius[c] + 2.5*capThickness[c];

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

  // The rims are about to be walked, and a triangle joined to the surface by
  // one edge would send the walk out along it and leave it there.
  {
    int numFlaps = 0;
    double firstFlap[3] = {0.0, 0.0, 0.0};
    if (RemoveBoundaryFlaps(outer, numFlaps, firstFlap) != SV_OK)
    {
      fprintf(stderr,"Problem removing the triangles hanging off the trimmed offset surface\n");
      return SV_ERROR;
    }
    if (numFlaps > 0)
    {
      fprintf(stdout,"  removed %d triangles hanging off the trimmed surface by a single edge, the first at (%.5g, %.5g, %.5g)\n",
          numFlaps, firstFlap[0], firstFlap[1], firstFlap[2]);
    }
  }

  std::vector<std::vector<vtkIdType> > outerLoops;
  if (TGenUtils_ExtractBoundaryLoops(outer, outerLoops) != SV_OK)
  {
    // The walk reports where it broke, but a point on a surface of millions is
    // not something that can be looked at from a coordinate alone. Writing the
    // surface out is what makes the failure visible, and it is written from
    // here rather than from the walk because this is the pass that produced it.
    char trimmedFile[] = "wall_offset_trimmed.vtp";
    if (TGenUtils_WriteVTP(trimmedFile, outer) == SV_OK)
    {
      fprintf(stderr,"Problem extracting the trimmed rims of the offset surface; it has been written to %s as it stood when the walk failed\n",
          trimmedFile);
    }
    else
    {
      fprintf(stderr,"Problem extracting the trimmed rims of the offset surface, and it could not be written out to be looked at\n");
    }
    return SV_ERROR;
  }

  if (outerLoops.size() != caps.size())
  {
    fprintf(stderr,"The trimmed offset surface has %zu rims but the wall has %zu cap openings. A cut has taken more than the dome off its end, which happens when a vessel curves back within a rim radius of another vessel's cap.\n",
        outerLoops.size(), caps.size());

    // Which rim is the spare decides what to do about it, and the two cases
    // read differently here. A spare that lies on a cap plane is a second loop
    // left by that cap's own cut, so the dome it cut was joined to something -
    // the offset of two vessels fused before the cut reached them. A spare far
    // off every plane is a cut that landed somewhere it does not belong, and
    // then the window that bounds it is what needs tightening.
    for (size_t l = 0; l < outerLoops.size(); l++)
    {
      double centre[3] = {0.0, 0.0, 0.0};
      for (size_t m = 0; m < outerLoops[l].size(); m++)
      {
        double x[3];
        outer->GetPoint(outerLoops[l][m], x);
        for (int k = 0; k < 3; k++)
        {
          centre[k] += x[k]/(double)outerLoops[l].size();
        }
      }

      double radius = 0.0;
      for (size_t m = 0; m < outerLoops[l].size(); m++)
      {
        double x[3];
        outer->GetPoint(outerLoops[l][m], x);
        radius = std::max(radius, std::sqrt(vtkMath::Distance2BetweenPoints(x, centre)));
      }

      size_t nearest = 0;
      double nearestDeviation = 0.0;
      for (size_t c = 0; c < caps.size(); c++)
      {
        double deviation = 0.0;
        for (size_t m = 0; m < outerLoops[l].size(); m++)
        {
          double x[3];
          outer->GetPoint(outerLoops[l][m], x);
          double offset[3];
          vtkMath::Subtract(x, caps[c].origin, offset);
          deviation = std::max(deviation, std::abs(vtkMath::Dot(offset, caps[c].outward)));
        }
        if (c == 0 || deviation < nearestDeviation)
        {
          nearest = c;
          nearestDeviation = deviation;
        }
      }

      fprintf(stderr,"  rim %zu: %zu points, centre (%.5g, %.5g, %.5g), radius %.5g; flattest against the cap at (%.5g, %.5g, %.5g), off its plane by %.5g against that cap's rim radius %.5g%s\n",
          l, outerLoops[l].size(), centre[0], centre[1], centre[2], radius,
          caps[nearest].origin[0], caps[nearest].origin[1], caps[nearest].origin[2],
          nearestDeviation, capRadius[nearest],
          nearestDeviation <= 0.05*capRadius[nearest] ? " -- on that plane" : "");
    }

    char trimmedFile[] = "wall_offset_trimmed.vtp";
    if (TGenUtils_WriteVTP(trimmedFile, outer) == SV_OK)
    {
      fprintf(stderr,"  the trimmed surface has been written to %s\n", trimmedFile);
    }

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
    fprintf(stdout,"    cap at (%.5g, %.5g, %.5g): inner rim %zu points, trimmed rim %zu points, cut within %.5g of the rim centre (radius %.5g, wall %.5g)\n",
        caps[c].origin[0], caps[c].origin[1], caps[c].origin[2],
        caps[c].innerLoop.size(), caps[c].outerLoop.size(),
        capRadius[c] + 2.5*capThickness[c], capRadius[c], capThickness[c]);
  }

  return SV_OK;
}

// -------------------------------------
// ClosestPointOnTriangle
// -------------------------------------
/**
 * @brief The point of a triangle nearest to a point, with its barycentric
 * weights.
 * @note Region by region over the triangle's Voronoi diagram (Ericson,
 * Real-Time Collision Detection, 5.1.5), so it never divides by the area of a
 * sliver.
 * @return The squared distance.
 */

static double ClosestPointOnTriangle(const double p[3], const double a[3],
    const double b[3], const double c[3], double closest[3], double weights[3])
{
  double ab[3], ac[3], ap[3], bp[3], cp[3];
  vtkMath::Subtract(b, a, ab);
  vtkMath::Subtract(c, a, ac);
  vtkMath::Subtract(p, a, ap);
  double d1 = vtkMath::Dot(ab, ap);
  double d2 = vtkMath::Dot(ac, ap);
  double u = 0.0, v = 0.0, w = 0.0;
  if (d1 <= 0.0 && d2 <= 0.0)
  {
    u = 1.0;
  }
  else
  {
    vtkMath::Subtract(p, b, bp);
    double d3 = vtkMath::Dot(ab, bp);
    double d4 = vtkMath::Dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
      v = 1.0;
    }
    else
    {
      double vc = d1*d4 - d3*d2;
      if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
      {
        v = (d1 - d3 != 0.0) ? d1/(d1 - d3) : 0.0;
        u = 1.0 - v;
      }
      else
      {
        vtkMath::Subtract(p, c, cp);
        double d5 = vtkMath::Dot(ab, cp);
        double d6 = vtkMath::Dot(ac, cp);
        if (d6 >= 0.0 && d5 <= d6)
        {
          w = 1.0;
        }
        else
        {
          double vb = d5*d2 - d1*d6;
          if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
          {
            w = (d2 - d6 != 0.0) ? d2/(d2 - d6) : 0.0;
            u = 1.0 - w;
          }
          else
          {
            double va = d3*d6 - d5*d4;
            if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
            {
              double denom = (d4 - d3) + (d5 - d6);
              w = (denom != 0.0) ? (d4 - d3)/denom : 0.0;
              v = 1.0 - w;
            }
            else
            {
              double denom = 1.0/(va + vb + vc);
              v = vb*denom;
              w = vc*denom;
              u = 1.0 - v - w;
            }
          }
        }
      }
    }
  }
  for (int k = 0; k < 3; k++)
  {
    closest[k] = u*a[k] + v*b[k] + w*c[k];
  }
  weights[0] = u;
  weights[1] = v;
  weights[2] = w;
  return vtkMath::Distance2BetweenPoints(p, closest);
}

// -------------------------------------
// SegmentCrossesTriangle
// -------------------------------------
/**
 * @brief Whether an open segment passes through the interior of a triangle.
 * @note Moller-Trumbore. Touching at an end of the segment or on the edge of
 * the triangle does not count, so triangles that share a corner or an edge
 * are not reported against each other by this alone; the caller leaves those
 * out anyway.
 */

static bool SegmentCrossesTriangle(const double p0[3], const double p1[3],
    const double a[3], const double b[3], const double c[3])
{
  double d[3], e1[3], e2[3], pvec[3], tvec[3], qvec[3];
  vtkMath::Subtract(p1, p0, d);
  vtkMath::Subtract(b, a, e1);
  vtkMath::Subtract(c, a, e2);
  vtkMath::Cross(d, e2, pvec);
  double det = vtkMath::Dot(e1, pvec);
  if (std::abs(det) < 1.0e-14)
  {
    return false;
  }
  double inv = 1.0/det;
  vtkMath::Subtract(p0, a, tvec);
  double u = vtkMath::Dot(tvec, pvec)*inv;
  if (u <= 1.0e-9 || u >= 1.0 - 1.0e-9)
  {
    return false;
  }
  vtkMath::Cross(tvec, e1, qvec);
  double v = vtkMath::Dot(d, qvec)*inv;
  if (v <= 1.0e-9 || u + v >= 1.0 - 1.0e-9)
  {
    return false;
  }
  double t = vtkMath::Dot(e2, qvec)*inv;
  return t > 1.0e-9 && t < 1.0 - 1.0e-9;
}

// -------------------------------------
// TrianglesCross
// -------------------------------------
/**
 * @brief Whether two triangles pass through each other.
 * @note Two triangles that share no corner cross when an edge of either
 * passes through the other. Two that share one corner can still cross, with
 * the edge opposite the shared corner of one passing through the other - a
 * fan folded over at its apex does this - and that edge is the only one
 * tested for them, because the two edges at the shared corner touch the other
 * triangle there by construction. Two that share an edge are not tested:
 * they can only overlap, coplanar, which the wall has no reason to produce.
 * @param f The corners of one triangle, nine values.
 * @param fp Its three point ids.
 * @param g The corners of the other, nine values.
 * @param gp Its three point ids.
 */

static bool TrianglesCross(const double *f, const vtkIdType *fp,
    const double *g, const vtkIdType *gp)
{
  int sharedF = -1, sharedG = -1, numShared = 0;
  for (int m = 0; m < 3; m++)
  {
    for (int n = 0; n < 3; n++)
    {
      if (fp[m] == gp[n])
      {
        sharedF = m;
        sharedG = n;
        numShared++;
      }
    }
  }
  if (numShared >= 2)
  {
    return false;
  }
  if (numShared == 1)
  {
    int f1 = (sharedF+1)%3, f2 = (sharedF+2)%3;
    int g1 = (sharedG+1)%3, g2 = (sharedG+2)%3;
    return SegmentCrossesTriangle(&f[3*f1], &f[3*f2], &g[0], &g[3], &g[6]) ||
        SegmentCrossesTriangle(&g[3*g1], &g[3*g2], &f[0], &f[3], &f[6]);
  }
  for (int k = 0; k < 3; k++)
  {
    if (SegmentCrossesTriangle(&f[3*k], &f[3*((k+1)%3)], &g[0], &g[3], &g[6]) ||
        SegmentCrossesTriangle(&g[3*k], &g[3*((k+1)%3)], &f[0], &f[3], &f[6]))
    {
      return true;
    }
  }
  return false;
}

// -------------------------------------
// TriangulateLoopByLeastArea
// -------------------------------------
/**
 * @brief Fills a closed loop of points with the triangulation of least total
 * area.
 * @note The loop is the edge of the hole left where a fold was cut out of the
 * outer wall, and the hole is thin: the two sides of the loop are the cut
 * edges of the two sheets that met at the crease, and they run alongside each
 * other a few percent of a wall apart. The least-area triangulation of a loop
 * like that is the one that steps back and forth across the gap, which is the
 * zipper wanted, and it needs no notion of which side is which - which matters,
 * because at a junction the two sides are one loop joined at the ends of the
 * crease. It is the textbook cubic dynamic programme over the loop order, so
 * the loop has to be short enough for that; a crease is.
 *
 * The triangles are wound against the loop's own order, which is the order
 * the surface's boundary walk gives, so they face the same way as the surface
 * they close.
 * @param points The points the loop indexes into.
 * @param loop The loop, in boundary walk order.
 * @param triangles The triangles are appended here, three ids each.
 * @param quiet Whether to say nothing when the loop cannot be filled, for a
 * caller trying this as one of several closures.
 * @return SV_OK if the loop was filled.
 */

static int TriangulateLoopByLeastArea(vtkPoints *points, const std::vector<vtkIdType> &loop,
    std::vector<vtkIdType> &triangles, bool quiet)
{
  size_t n = loop.size();
  if (n < 3)
  {
    if (!quiet)
    {
      fprintf(stderr,"Cannot fill a loop of %zu points\n", n);
    }
    return SV_ERROR;
  }
  // The cost is cubic in the loop, and this is the loop a crease reaches
  // before the cost is seconds.
  const size_t maxLoop = 1500;
  if (n > maxLoop)
  {
    if (!quiet)
    {
      double p[3];
      points->GetPoint(loop[0], p);
      fprintf(stderr,"A hole in the outer wall at (%.5g, %.5g, %.5g) has a %zu point edge, more than the %zu the least-area fill is bounded to\n",
          p[0], p[1], p[2], n, maxLoop);
    }
    return SV_ERROR;
  }
  std::vector<double> xyz(3*n);
  double loopBounds[6];
  for (size_t i = 0; i < n; i++)
  {
    points->GetPoint(loop[i], &xyz[3*i]);
    for (int c = 0; c < 3; c++)
    {
      loopBounds[2*c] = (i == 0) ? xyz[c] : std::min(loopBounds[2*c], xyz[3*i+c]);
      loopBounds[2*c+1] = (i == 0) ? xyz[c] : std::max(loopBounds[2*c+1], xyz[3*i+c]);
    }
  }
  // A triangle with no area is charged the whole loop's extent squared on
  // top, so that the fill takes any triangulation without one before one
  // with. The least-area choice alone does not avoid them: every
  // triangulation of a flat loop has the same area, and a triangle of three
  // points in a line, whose long edge passes through the middle point, then
  // costs nothing.
  double extent = 0.0;
  for (int c = 0; c < 3; c++)
  {
    extent = std::max(extent, loopBounds[2*c+1] - loopBounds[2*c]);
  }
  const double degeneratePenalty = extent*extent;
  auto area = [&](size_t i, size_t k, size_t j)
  {
    double e1[3], e2[3], cross[3];
    for (int c = 0; c < 3; c++)
    {
      e1[c] = xyz[3*k+c] - xyz[3*i+c];
      e2[c] = xyz[3*j+c] - xyz[3*i+c];
    }
    vtkMath::Cross(e1, e2, cross);
    double a = 0.5*vtkMath::Norm(cross);
    // Flat to within the length of its longest edge times a hair.
    double longest2 = std::max(vtkMath::Dot(e1, e1), vtkMath::Dot(e2, e2));
    if (a <= 1.0e-6*longest2)
    {
      a += degeneratePenalty;
    }
    return a;
  };
  std::vector<double> cost(n*n, 0.0);
  std::vector<int> split(n*n, -1);
  for (size_t len = 2; len < n; len++)
  {
    for (size_t i = 0; i + len < n; i++)
    {
      size_t j = i + len;
      double best = std::numeric_limits<double>::max();
      int bestK = -1;
      for (size_t k = i+1; k < j; k++)
      {
        double c = cost[i*n+k] + cost[k*n+j] + area(i, k, j);
        if (c < best)
        {
          best = c;
          bestK = (int)k;
        }
      }
      cost[i*n+j] = best;
      split[i*n+j] = bestK;
    }
  }
  std::vector<std::pair<size_t,size_t> > pending;
  pending.push_back(std::make_pair((size_t)0, n-1));
  while (!pending.empty())
  {
    size_t i = pending.back().first;
    size_t j = pending.back().second;
    pending.pop_back();
    if (j - i < 2)
    {
      continue;
    }
    size_t k = (size_t)split[i*n+j];
    triangles.push_back(loop[j]);
    triangles.push_back(loop[k]);
    triangles.push_back(loop[i]);
    pending.push_back(std::make_pair(i, k));
    pending.push_back(std::make_pair(k, j));
  }
  return SV_OK;
}

// -------------------------------------
// BuildSeamBands
// -------------------------------------
/**
 * @brief Builds the two bands of triangles that join two closed loops running
 * alongside each other, one for each way round the second loop, the way two
 * rims are joined.
 * @note Where two vessels run closer than twice the wall, or a branch leaves a
 * vessel whose wall is thicker than the branch is wide, the wall of one
 * crosses the wall of the other along a closed curve, and cutting both back
 * leaves a loop on each sheet with nothing between them. Those two loops have
 * to be joined to each other. Filling either one on its own would seal its
 * sheet shut along the seam, which is a surface inside the wall and not its
 * boundary.
 *
 * The loops are walked together from their closest pair of points, advancing
 * whichever loop leaves the shorter diagonal, which is the merge that suits two
 * curves a small fraction of a wall apart. Which way round the second loop has
 * to be walked cannot be told from the neighbours of the closest pair alone,
 * so both bands are built and handed back with their areas; the caller
 * chooses. Each triangle is wound against the loop edge it uses, so the band
 * faces the way the sheets do.
 * @param points The points both loops index into.
 * @param first One loop, in boundary walk order.
 * @param second The other loop, in boundary walk order.
 * @param bands Set to the two bands, three ids per triangle: the second loop
 * walked forwards in the first, backwards in the second.
 * @param areas Set to the total area of each band.
 * @return SV_OK if the bands were built.
 */

static int BuildSeamBands(vtkPoints *points, const std::vector<vtkIdType> &first,
    const std::vector<vtkIdType> &second, std::vector<vtkIdType> bands[2], double areas[2])
{
  size_t n = first.size();
  size_t m = second.size();
  if (n < 3 || m < 3)
  {
    fprintf(stderr,"Cannot join loops of %zu and %zu points\n", n, m);
    return SV_ERROR;
  }
  auto distance2 = [&](vtkIdType a, vtkIdType b)
  {
    double p[3], q[3];
    points->GetPoint(a, p);
    points->GetPoint(b, q);
    return vtkMath::Distance2BetweenPoints(p, q);
  };
  size_t i0 = 0, j0 = 0;
  double closest = std::numeric_limits<double>::max();
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < m; j++)
    {
      double d2 = distance2(first[i], second[j]);
      if (d2 < closest)
      {
        closest = d2;
        i0 = i;
        j0 = j;
      }
    }
  }
  auto F = [&](size_t step) { return first[(i0 + step) % n]; };
  // Deciding the direction from the neighbours of the closest pair alone is
  // wrong when that pair sits at the end of a seam, where the loop turns
  // round and both neighbours are equally close; the band then joins one side
  // of the seam to the other and crosses everything between. So both are
  // built.
  auto band = [&](int dir, std::vector<vtkIdType> &triangles)
  {
    triangles.clear();
    auto S = [&](size_t step)
    {
      size_t offset = step % m;
      return (dir > 0) ? second[(j0 + offset) % m] : second[(j0 + m - offset) % m];
    };
    double area = 0.0;
    size_t i = 0, j = 0;
    while (i < n || j < m)
    {
      bool advanceFirst;
      if (i >= n)
      {
        advanceFirst = false;
      }
      else if (j >= m)
      {
        advanceFirst = true;
      }
      else
      {
        advanceFirst = distance2(F(i+1), S(j)) <= distance2(F(i), S(j+1));
      }
      vtkIdType triangle[3];
      if (advanceFirst)
      {
        triangle[0] = F(i+1);
        triangle[1] = F(i);
        triangle[2] = S(j);
        i++;
      }
      else if (dir < 0)
      {
        triangle[0] = S(j);
        triangle[1] = S(j+1);
        triangle[2] = F(i);
        j++;
      }
      else
      {
        triangle[0] = S(j+1);
        triangle[1] = S(j);
        triangle[2] = F(i);
        j++;
      }
      if (triangle[0] == triangle[1] || triangle[1] == triangle[2] || triangle[0] == triangle[2])
      {
        continue;
      }
      double p0[3], p1[3], p2[3], e1[3], e2[3], cross[3];
      points->GetPoint(triangle[0], p0);
      points->GetPoint(triangle[1], p1);
      points->GetPoint(triangle[2], p2);
      vtkMath::Subtract(p1, p0, e1);
      vtkMath::Subtract(p2, p0, e2);
      vtkMath::Cross(e1, e2, cross);
      area += 0.5*vtkMath::Norm(cross);
      for (int k = 0; k < 3; k++)
      {
        triangles.push_back(triangle[k]);
      }
    }
    return area;
  };
  areas[0] = band(1, bands[0]);
  areas[1] = band(-1, bands[1]);
  return SV_OK;
}

// -------------------------------------
// TGenUtils_BuildTrimmedExtrudedOuterSurface
// -------------------------------------
/**
 * @brief Builds the outer wall surface by extruding the inner surface along its
 * normals and cutting out every part of the result that lies inside the wall.
 * @note Dilating the solid by the wall thickness rounds its convex features and
 * creases its concave ones. Pushing every surface point out along its normal
 * gives the same surface everywhere except at the creases, where the extruded
 * sheets of the two sides run on through each other, and the parts that run on
 * are not on the boundary of the dilated solid: every one of them lies within
 * the wall of some other part of the inner surface. That is a test that can be
 * made at each extruded point alone, against the triangles of the inner
 * surface around it, so this is the offset surface without a distance field
 * over a grid. A grid has to resolve the thinnest wall in the model over the
 * whole model, and on a model whose walls span an order of magnitude it cannot
 * afford to; this resolves everything at the resolution of the surface.
 *
 * An extruded point is measured against every inner triangle within the
 * largest wall of it, other than its own sheet - the triangles at its own
 * point, and those facing within a few degrees of its own extrusion direction,
 * which a sheet cannot have folded through at that angle. Each such triangle
 * carries its own wall, and the point is cut if it stands closer to any of
 * them than that wall and a small clearance. It is measured against each triangle's own
 * wall rather than against the wall at the nearest point because the nearest
 * point is the wrong one where a thin vessel leaves a thick one: the thin
 * wall is nearest, and the thick wall is the one the point is inside. A point
 * inside a lumen, or on a triangle the extrusion turned inside out, is cut
 * whatever the distances say.
 *
 * The cut runs through the triangles between kept and cut points, with the
 * cut point on each edge put where the edge crosses the clearance, so the
 * surface is left with holes whose edges stand a small clearance off the
 * sheet they ran into, and the two sheets stop short of each other instead of
 * crossing. The holes are then closed: a hole on its
 * own is the gap at a crease and is zipped shut across it, and two holes that
 * run alongside each other are the two sides of a seam where one vessel's wall
 * crosses another's, and are joined to each other. Every closure is tested
 * for passing through the surface, which is what the volume mesher refuses,
 * and a pair whose band does so is closed the other way round, or each on its
 * own, whichever crosses least.
 *
 * The extruded rim of each cap is kept as the outer rim of that vessel end, so
 * the two rims are one to one. A cap whose rim is cut into by a neighbouring
 * wall is refused, because the end face of that vessel would no longer be an
 * annulus.
 *
 * The inner surface is never touched. Its points are the fluid/wall interface.
 * @param surface The inner surface with its 'Normals' point data.
 * @param array The wall thickness per point of the inner surface.
 * @param clearance An extruded point standing less than this fraction of
 * another sheet's wall from that sheet's surface is cut, and the cut points
 * on the edges are put where the edge crosses this fraction. Above one, so
 * the sheets that met end short of each other.
 * @param outer Set to the trimmed and closed outer surface.
 * @param caps Set to one entry per vessel end, pairing its inner rim with its
 * extruded rim.
 * @param numUnresolved Set to the number of triangles the surface is left
 * with that the volume mesher will refuse: fragments still turned over when
 * the cut stopped, and triangles passing through the surface. The surface is
 * returned whatever this is, so that it can be measured and looked at; the
 * caller decides whether to go on.
 * @return SV_OK if the outer surface was built and every hole closed.
 */

int TGenUtils_BuildTrimmedExtrudedOuterSurface(vtkPolyData *surface, vtkDoubleArray *array,
    double clearance, vtkPolyData *outer, std::vector<TGenUtilsCapRim> &caps, int &numUnresolved)
{
  caps.clear();
  numUnresolved = 0;

  if (surface == nullptr || array == nullptr || outer == nullptr)
  {
    fprintf(stderr,"Cannot build the extruded outer surface without a surface, a thickness array and an output\n");
    return SV_ERROR;
  }
  if (!(clearance > 1.0))
  {
    fprintf(stderr,"The extrusion trim needs a clearance above one, not %.5g\n", clearance);
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
    fprintf(stderr,"The surface has no 'Normals' point data to extrude along\n");
    return SV_ERROR;
  }
  auto start = std::chrono::steady_clock::now();

  double largestThickness = 0.0;
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    largestThickness = std::max(largestThickness, array->GetValue(ptId));
  }
  if (largestThickness <= 0.0)
  {
    fprintf(stderr,"Every wall thickness is zero or negative, so there is no wall to extrude\n");
    return SV_ERROR;
  }

  // Close the cap openings with a fan each so the distance can be signed. The
  // fans are not walls: nothing is measured against them.
  std::vector<std::vector<vtkIdType> > rims;
  if (TGenUtils_ExtractBoundaryLoops(surface, rims) != SV_OK)
  {
    fprintf(stderr,"Problem extracting the cap rims of the wall surface\n");
    return SV_ERROR;
  }
  auto closedPoints = vtkSmartPointer<vtkPoints>::New();
  closedPoints->DeepCopy(surface->GetPoints());
  auto closedCells = vtkSmartPointer<vtkCellArray>::New();
  std::vector<vtkIdType> wallCellPts;
  std::vector<double> wallCellNormal;
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
    double p0[3], p1[3], p2[3], e1[3], e2[3], n[3];
    surface->GetPoint(pts[0], p0);
    surface->GetPoint(pts[1], p1);
    surface->GetPoint(pts[2], p2);
    vtkMath::Subtract(p1, p0, e1);
    vtkMath::Subtract(p2, p0, e2);
    vtkMath::Cross(e1, e2, n);
    vtkMath::Normalize(n);
    for (int j = 0; j < 3; j++)
    {
      wallCellPts.push_back(pts[j]);
      wallCellNormal.push_back(n[j]);
    }
  }
  vtkIdType numWallCells = (vtkIdType)(wallCellPts.size()/3);
  std::vector<int> rimOfPoint((size_t)numPts, -1);
  for (size_t r = 0; r < rims.size(); r++)
  {
    const std::vector<vtkIdType> &rim = rims[r];
    double centroid[3] = {0.0, 0.0, 0.0};
    for (size_t m = 0; m < rim.size(); m++)
    {
      double p[3];
      closedPoints->GetPoint(rim[m], p);
      for (int k = 0; k < 3; k++)
      {
        centroid[k] += p[k];
      }
      rimOfPoint[(size_t)rim[m]] = (int)r;
    }
    for (int k = 0; k < 3; k++)
    {
      centroid[k] /= (double)rim.size();
    }
    vtkIdType centroidId = closedPoints->InsertNextPoint(centroid);
    for (size_t m = 0; m < rim.size(); m++)
    {
      vtkIdType a = rim[m];
      vtkIdType b = rim[(m+1)%rim.size()];
      vtkIdType fan[3] = {b, a, centroidId};
      closedCells->InsertNextCell(3, fan);
    }
  }
  auto closed = vtkSmartPointer<vtkPolyData>::New();
  closed->SetPoints(closedPoints);
  closed->SetPolys(closedCells);
  closed->BuildLinks();

  auto implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  implicit->SetInput(closed);
  // Which side is negative is the filter's convention; a corner well outside
  // the model says which, and its magnitude says the reading is a distance.
  double bounds[6];
  closed->GetBounds(bounds);
  double margin = 0.0;
  for (int k = 0; k < 3; k++)
  {
    margin = std::max(margin, bounds[2*k+1] - bounds[2*k]);
  }
  margin = 0.5*margin + 1.0;
  double corner[3] = {bounds[0] - margin, bounds[2] - margin, bounds[4] - margin};
  double cornerDistance = implicit->EvaluateFunction(corner);
  if (std::abs(cornerDistance) < 0.9*margin)
  {
    fprintf(stderr,"The signed distance at a corner %.5g outside the model reads %.5g, so the distance is not measuring what the trim needs\n",
        margin, cornerDistance);
    return SV_ERROR;
  }
  double outwardSign = (cornerDistance > 0.0) ? 1.0 : -1.0;

  // Extrude.
  std::vector<double> extruded((size_t)3*numPts, 0.0);
  std::vector<double> direction((size_t)3*numPts, 0.0);
  int numNoThickness = 0;
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double p[3], n[3];
    surface->GetPoint(ptId, p);
    normals->GetTuple(ptId, n);
    vtkMath::Normalize(n);
    double t = array->GetValue(ptId);
    if (t <= 0.0)
    {
      numNoThickness++;
      t = 0.0;
    }
    for (int k = 0; k < 3; k++)
    {
      extruded[(size_t)3*ptId + k] = p[k] + t*n[k];
      direction[(size_t)3*ptId + k] = n[k];
    }
  }

  // A triangle the extrusion turned over is a fold whatever the distances at
  // its corners say, so its corners are cut.
  std::vector<bool> inverted((size_t)numPts, false);
  int numInvertedCells = 0;
  for (vtkIdType cellId = 0; cellId < numWallCells; cellId++)
  {
    const vtkIdType *pts = &wallCellPts[(size_t)3*cellId];
    const double *q0 = &extruded[(size_t)3*pts[0]];
    const double *q1 = &extruded[(size_t)3*pts[1]];
    const double *q2 = &extruded[(size_t)3*pts[2]];
    double e1[3], e2[3], outerNormal[3];
    for (int k = 0; k < 3; k++)
    {
      e1[k] = q1[k] - q0[k];
      e2[k] = q2[k] - q0[k];
    }
    vtkMath::Cross(e1, e2, outerNormal);
    if (vtkMath::Dot(&wallCellNormal[(size_t)3*cellId], outerNormal) <= 0.0)
    {
      numInvertedCells++;
      for (int j = 0; j < 3; j++)
      {
        inverted[(size_t)pts[j]] = true;
      }
    }
  }

  // The surface's own structure: the triangles at each point, and the
  // neighbour of each triangle across each of its edges. Both are what "the
  // same sheet" is measured on.
  std::vector<int> pointCellStart((size_t)numPts + 1, 0);
  for (vtkIdType cellId = 0; cellId < numWallCells; cellId++)
  {
    for (int j = 0; j < 3; j++)
    {
      pointCellStart[(size_t)wallCellPts[(size_t)3*cellId + j] + 1]++;
    }
  }
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    pointCellStart[(size_t)ptId + 1] += pointCellStart[(size_t)ptId];
  }
  std::vector<vtkIdType> pointCells((size_t)pointCellStart[(size_t)numPts], -1);
  {
    std::vector<int> cursor(pointCellStart.begin(), pointCellStart.end() - 1);
    for (vtkIdType cellId = 0; cellId < numWallCells; cellId++)
    {
      for (int j = 0; j < 3; j++)
      {
        vtkIdType ptId = wallCellPts[(size_t)3*cellId + j];
        pointCells[(size_t)cursor[(size_t)ptId]++] = cellId;
      }
    }
  }
  std::vector<vtkIdType> cellNeighbor((size_t)3*numWallCells, -1);
  {
    struct EdgeSide
    {
      vtkIdType a, b, cell;
      int side;
    };
    std::vector<EdgeSide> edges;
    edges.reserve((size_t)3*numWallCells);
    for (vtkIdType cellId = 0; cellId < numWallCells; cellId++)
    {
      for (int j = 0; j < 3; j++)
      {
        vtkIdType a = wallCellPts[(size_t)3*cellId + j];
        vtkIdType b = wallCellPts[(size_t)3*cellId + (j+1)%3];
        EdgeSide e = {std::min(a, b), std::max(a, b), cellId, j};
        edges.push_back(e);
      }
    }
    std::sort(edges.begin(), edges.end(), [](const EdgeSide &l, const EdgeSide &r)
    {
      return (l.a != r.a) ? (l.a < r.a) : (l.b < r.b);
    });
    for (size_t e = 0; e + 1 < edges.size(); e++)
    {
      if (edges[e].a == edges[e+1].a && edges[e].b == edges[e+1].b)
      {
        cellNeighbor[(size_t)3*edges[e].cell + edges[e].side] = edges[e+1].cell;
        cellNeighbor[(size_t)3*edges[e+1].cell + edges[e+1].side] = edges[e].cell;
        e++;
      }
    }
  }

  // Bin the triangles on a grid at the reach, so the triangles within reach
  // of a point are the ones in the bins around it. Each triangle also carries
  // how far it can matter - the clearance times the largest wall at its
  // corners - so the thin-walled triangles that fill most of a model are
  // dismissed on their box before any distance is worked out.
  const double reach = clearance*largestThickness;
  double gridOrigin[3], gridDims[3];
  int gridSize[3];
  // The bins are the reach across, so that the triangles within reach of a
  // point are in the bins around it, but never more of them than this: the
  // grid covers the bounding box, empty space included, and a thin wall on a
  // large model would otherwise ask for bounding box over reach cubed of
  // them, which is memory in the tens of gigabytes. A bin wider than the
  // reach just holds more triangles per bin, and the box test on each
  // triangle still dismisses them cheaply.
  const double maxBins = 4.0e6;
  double binSize = reach;
  {
    double volume = 1.0;
    for (int k = 0; k < 3; k++)
    {
      volume *= bounds[2*k+1] - bounds[2*k] + 2.0*reach;
    }
    binSize = std::max(reach, std::cbrt(volume/maxBins));
  }
  for (int k = 0; k < 3; k++)
  {
    gridOrigin[k] = bounds[2*k] - reach;
    gridDims[k] = bounds[2*k+1] + reach - gridOrigin[k];
    gridSize[k] = std::max(1, (int)std::ceil(gridDims[k]/binSize));
  }
  if ((double)gridSize[0]*gridSize[1]*gridSize[2] > 2.0*maxBins)
  {
    fprintf(stderr,"The trim grid would be %d x %d x %d bins over a %.5g x %.5g x %.5g box, more than it is bounded to\n",
        gridSize[0], gridSize[1], gridSize[2], gridDims[0], gridDims[1], gridDims[2]);
    return SV_ERROR;
  }
  auto binOf = [&](double value, int axis)
  {
    int bin = (int)std::floor((value - gridOrigin[axis])/binSize);
    return std::min(std::max(bin, 0), gridSize[axis] - 1);
  };
  auto binIndex = [&](int i, int j, int k)
  {
    return ((size_t)k*gridSize[1] + (size_t)j)*gridSize[0] + (size_t)i;
  };
  size_t numBins = (size_t)gridSize[0]*gridSize[1]*gridSize[2];
  std::vector<double> cellBox((size_t)6*numWallCells, 0.0);
  std::vector<double> cellReach2((size_t)numWallCells, 0.0);
  std::vector<int> binStart(numBins + 1, 0);
  for (vtkIdType cellId = 0; cellId < numWallCells; cellId++)
  {
    double *box = &cellBox[(size_t)6*cellId];
    box[0] = box[2] = box[4] = std::numeric_limits<double>::max();
    box[1] = box[3] = box[5] = -std::numeric_limits<double>::max();
    double thickest = 0.0;
    for (int j = 0; j < 3; j++)
    {
      vtkIdType ptId = wallCellPts[(size_t)3*cellId + j];
      double p[3];
      surface->GetPoint(ptId, p);
      for (int k = 0; k < 3; k++)
      {
        box[2*k] = std::min(box[2*k], p[k]);
        box[2*k+1] = std::max(box[2*k+1], p[k]);
      }
      thickest = std::max(thickest, array->GetValue(ptId));
    }
    double r = clearance*thickest;
    cellReach2[(size_t)cellId] = r*r;
    for (int k = binOf(box[4], 2); k <= binOf(box[5], 2); k++)
    {
      for (int j = binOf(box[2], 1); j <= binOf(box[3], 1); j++)
      {
        for (int i = binOf(box[0], 0); i <= binOf(box[1], 0); i++)
        {
          binStart[binIndex(i, j, k) + 1]++;
        }
      }
    }
  }
  for (size_t b = 0; b < numBins; b++)
  {
    binStart[b + 1] += binStart[b];
  }
  std::vector<vtkIdType> binCells((size_t)binStart[numBins], -1);
  {
    std::vector<int> cursor(binStart.begin(), binStart.end() - 1);
    for (vtkIdType cellId = 0; cellId < numWallCells; cellId++)
    {
      const double *box = &cellBox[(size_t)6*cellId];
      for (int k = binOf(box[4], 2); k <= binOf(box[5], 2); k++)
      {
        for (int j = binOf(box[2], 1); j <= binOf(box[3], 1); j++)
        {
          for (int i = binOf(box[0], 0); i <= binOf(box[1], 0); i++)
          {
            binCells[(size_t)cursor[binIndex(i, j, k)]++] = cellId;
          }
        }
      }
    }
  }

  // How far a point stands inside the wall of another sheet, as the smallest
  // ratio of its distance to an inner triangle over that triangle's own wall,
  // taken over every triangle within reach of it that is not its own sheet.
  //
  // Its own sheet is the part of the surface it was extruded from, and that is
  // a matter of connection, not of distance: the triangles reached from its
  // own point by walking across edges without leaving a ball around it. A
  // wall standing behind the sheet is not reached that way however close it
  // is, which is what cuts a sheet that has run through the vessel behind it.
  // But connection alone is too generous where the walk crosses a junction,
  // because at a junction the other vessel's wall is connected too, so the
  // walk's triangles count as own only while they face within an angle of the
  // extrusion direction; the other side of a crotch faces away by the angle
  // of the crotch. The angle is what a sheet cannot have folded through: its
  // own facets lean away from a vertex normal by the warp smoothing's ten
  // degrees and a little curvature, and two vessels meet at more than this.
  const double cosOwnSheet = std::cos(vtkMath::RadiansFromDegrees(30.0));
  int stamp = 0;
  std::vector<int> ownStamp((size_t)numWallCells, 0);
  std::vector<int> seenStamp((size_t)numWallCells, 0);
  std::vector<vtkIdType> walk;
  long long numQueries = 0;
  auto markOwn = [&](vtkIdType seedPt, double radius)
  {
    double seed[3];
    surface->GetPoint(seedPt, seed);
    double radius2 = radius*radius;
    walk.clear();
    for (int c = pointCellStart[(size_t)seedPt]; c < pointCellStart[(size_t)seedPt + 1]; c++)
    {
      vtkIdType cellId = pointCells[(size_t)c];
      if (ownStamp[(size_t)cellId] != stamp)
      {
        ownStamp[(size_t)cellId] = stamp;
        walk.push_back(cellId);
      }
    }
    for (size_t w = 0; w < walk.size(); w++)
    {
      vtkIdType cellId = walk[w];
      for (int j = 0; j < 3; j++)
      {
        vtkIdType next = cellNeighbor[(size_t)3*cellId + j];
        if (next < 0 || ownStamp[(size_t)next] == stamp)
        {
          continue;
        }
        bool inside = false;
        for (int m = 0; m < 3 && !inside; m++)
        {
          double p[3];
          surface->GetPoint(wallCellPts[(size_t)3*next + m], p);
          inside = vtkMath::Distance2BetweenPoints(p, seed) <= radius2;
        }
        if (inside)
        {
          ownStamp[(size_t)next] = stamp;
          walk.push_back(next);
        }
      }
    }
  };
  auto standing = [&](const double x[3], vtkIdType ownA, vtkIdType ownB,
      const double along[3], double closest[3], double &wall)
  {
    numQueries++;
    stamp++;
    // A triangle within reach of the point is within the point's extrusion
    // plus the reach of its own point, plus the edge a cut point sits on.
    double pa[3], pb[3];
    surface->GetPoint(ownA, pa);
    surface->GetPoint(ownB, pb);
    double ownRadius = reach + clearance*std::max(array->GetValue(ownA), array->GetValue(ownB)) +
        std::sqrt(vtkMath::Distance2BetweenPoints(pa, pb));
    markOwn(ownA, ownRadius);
    if (ownB != ownA)
    {
      markOwn(ownB, ownRadius);
    }
    double best = std::numeric_limits<double>::max();
    wall = 0.0;
    for (int k = binOf(x[2] - reach, 2); k <= binOf(x[2] + reach, 2); k++)
    {
      for (int j = binOf(x[1] - reach, 1); j <= binOf(x[1] + reach, 1); j++)
      {
        for (int i = binOf(x[0] - reach, 0); i <= binOf(x[0] + reach, 0); i++)
        {
          size_t bin = binIndex(i, j, k);
          for (int c = binStart[bin]; c < binStart[bin + 1]; c++)
          {
            vtkIdType cellId = binCells[(size_t)c];
            if (seenStamp[(size_t)cellId] == stamp)
            {
              continue;
            }
            seenStamp[(size_t)cellId] = stamp;
            const vtkIdType *pts = &wallCellPts[(size_t)3*cellId];
            if (pts[0] == ownA || pts[1] == ownA || pts[2] == ownA ||
                pts[0] == ownB || pts[1] == ownB || pts[2] == ownB)
            {
              continue;
            }
            if (ownStamp[(size_t)cellId] == stamp &&
                vtkMath::Dot(&wallCellNormal[(size_t)3*cellId], along) >= cosOwnSheet)
            {
              continue;
            }
            const double *box = &cellBox[(size_t)6*cellId];
            double boxDistance2 = 0.0;
            for (int m = 0; m < 3; m++)
            {
              double d = (x[m] < box[2*m]) ? box[2*m] - x[m] : ((x[m] > box[2*m+1]) ? x[m] - box[2*m+1] : 0.0);
              boxDistance2 += d*d;
            }
            if (boxDistance2 > cellReach2[(size_t)cellId])
            {
              continue;
            }
            double a[3], b[3], cc[3], foot[3], weights[3];
            closedPoints->GetPoint(pts[0], a);
            closedPoints->GetPoint(pts[1], b);
            closedPoints->GetPoint(pts[2], cc);
            double d2 = ClosestPointOnTriangle(x, a, b, cc, foot, weights);
            double t = weights[0]*array->GetValue(pts[0]) + weights[1]*array->GetValue(pts[1]) +
                weights[2]*array->GetValue(pts[2]);
            if (t <= 0.0)
            {
              continue;
            }
            double ratio = std::sqrt(d2)/t;
            if (ratio < best)
            {
              best = ratio;
              wall = t;
              for (int m = 0; m < 3; m++)
              {
                closest[m] = foot[m];
              }
            }
          }
        }
      }
    }
    return best;
  };

  std::vector<double> standingOf((size_t)numPts, 1.0);
  std::vector<bool> removed((size_t)numPts, false);
  int numInLumen = 0, numUnderWall = 0, numInvertedPts = 0, numNoThicknessOnly = 0;
  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    double x[3] = {extruded[(size_t)3*ptId], extruded[(size_t)3*ptId+1], extruded[(size_t)3*ptId+2]};
    double closest[3], wall = 0.0;
    double f = standing(x, ptId, ptId, &direction[(size_t)3*ptId], closest, wall);
    // A point can be deep inside a lumen whose wall is thin, well over a wall
    // from that wall, so the sign is a separate question from the ratio.
    if (outwardSign*implicit->EvaluateFunction(x) < 0.0)
    {
      f = -1.0;
    }
    standingOf[(size_t)ptId] = f;
    if (inverted[(size_t)ptId])
    {
      removed[(size_t)ptId] = true;
      numInvertedPts++;
    }
    else if (f < 0.0)
    {
      removed[(size_t)ptId] = true;
      numInLumen++;
    }
    else if (f < clearance)
    {
      removed[(size_t)ptId] = true;
      numUnderWall++;
    }
    else if (array->GetValue(ptId) <= 0.0)
    {
      removed[(size_t)ptId] = true;
      numNoThicknessOnly++;
    }
  }
  double firstPassSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  // Cut, and look at what that made. Nothing is moved off the extruded
  // sheets: a cut point is put where its edge crosses the clearance, found by
  // bisection on the standing measure along the edge, so every outer point
  // is on a sheet or on an edge of one and stands at least the clearance from
  // every other sheet. Moving points off the sheets instead, measured, fanned
  // slivers over each other at a point whose neighbours were all cut, and
  // where the other wall stood beyond the sheet rather than across it the
  // move away from that wall went into the sheet's own solid.
  //
  // A kept point on no whole triangle is left holding a fan of slivers
  // between cut points, and a kept point whose edges cross the clearance
  // within a twentieth of their length would leave slivers of that size; both
  // are cut in turn and the cut made again, until the cut is stable. A cap
  // rim point is never cut this way: the rim has to survive whole for the
  // vessel end to be closed, and cutting one would fail the build outright
  // where holding the cut point off it leaves a sliver at worst.
  auto outerPoints = vtkSmartPointer<vtkPoints>::New();
  auto outerCells = vtkSmartPointer<vtkCellArray>::New();
  std::vector<vtkIdType> newId((size_t)numPts, -1);
  std::vector<double> scale;
  std::vector<vtkIdType> cutFrom, cutTo, sourceCell, outerTris;
  vtkIdType numKeptOriginal = 0, numOuterPts = 0, numCutPts = 0;
  int numCutCells = 0, numDroppedCells = 0;
  int numPeninsula = 0, numDemoted = 0, numFolded = 0, numRounds = 0, lastFolded = 0;
  int numRimHeld = 0;
  bool converged = false;
  const int maxRounds = 8;
  const double edgeFloor = 0.05;
  while (true)
  {
    numRounds++;
    numRimHeld = 0;
    // Peninsulas: kept points on no whole triangle.
    while (true)
    {
      int found = 0;
      for (vtkIdType ptId = 0; ptId < numPts; ptId++)
      {
        if (removed[(size_t)ptId])
        {
          continue;
        }
        bool onWhole = false;
        for (int c = pointCellStart[(size_t)ptId]; c < pointCellStart[(size_t)ptId + 1] && !onWhole; c++)
        {
          const vtkIdType *pts = &wallCellPts[(size_t)3*pointCells[(size_t)c]];
          onWhole = !removed[(size_t)pts[0]] && !removed[(size_t)pts[1]] && !removed[(size_t)pts[2]];
        }
        if (!onWhole)
        {
          removed[(size_t)ptId] = true;
          found++;
        }
      }
      numPeninsula += found;
      if (found == 0)
      {
        break;
      }
    }

    // Every cap rim has to survive whole: its extruded rim is the outer rim of
    // the vessel end, and the end face is the annulus between the two.
    for (size_t r = 0; r < rims.size(); r++)
    {
      int numCut = 0;
      double where[3] = {0.0, 0.0, 0.0};
      for (size_t m = 0; m < rims[r].size(); m++)
      {
        if (removed[(size_t)rims[r][m]])
        {
          if (numCut == 0)
          {
            surface->GetPoint(rims[r][m], where);
          }
          numCut++;
        }
      }
      if (numCut > 0)
      {
        fprintf(stderr,"The wall of a neighbouring vessel cuts into the rim of the cap at (%.5g, %.5g, %.5g): %d of its %zu extruded rim points lie inside another wall, so the end face of that vessel is not an annulus and the wall cannot be closed there\n",
            where[0], where[1], where[2], numCut, rims[r].size());
        return SV_ERROR;
      }
    }

    // Cut the triangles between kept and cut points, keeping the kept side.
    // The cut point on an edge is shared by the two triangles on that edge,
    // which is what keeps the cut edges a closed chain.
    outerPoints = vtkSmartPointer<vtkPoints>::New();
    outerCells = vtkSmartPointer<vtkCellArray>::New();
    std::fill(newId.begin(), newId.end(), -1);
    scale.clear();
    cutFrom.clear();
    cutTo.clear();
    sourceCell.clear();
    outerTris.clear();
    auto addOuterTriangle = [&](const vtkIdType triangle[3], vtkIdType source)
    {
      outerCells->InsertNextCell(3, triangle);
      for (int j = 0; j < 3; j++)
      {
        outerTris.push_back(triangle[j]);
      }
      sourceCell.push_back(source);
    };
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      if (removed[(size_t)ptId])
      {
        continue;
      }
      newId[(size_t)ptId] = outerPoints->InsertNextPoint(&extruded[(size_t)3*ptId]);
      scale.push_back(array->GetValue(ptId));
      cutFrom.push_back(ptId);
      cutTo.push_back(ptId);
    }
    numKeptOriginal = outerPoints->GetNumberOfPoints();
    std::vector<bool> demote((size_t)numPts, false);
    int demotedThisRound = 0;
    std::map<std::pair<vtkIdType,vtkIdType>, vtkIdType> cutPoints;
    auto cutPointOn = [&](vtkIdType kept, vtkIdType cut)
    {
      std::pair<vtkIdType,vtkIdType> key(std::min(kept, cut), std::max(kept, cut));
      std::map<std::pair<vtkIdType,vtkIdType>, vtkIdType>::iterator found = cutPoints.find(key);
      if (found != cutPoints.end())
      {
        return found->second;
      }
      const double *qk = &extruded[(size_t)3*kept];
      const double *qc = &extruded[(size_t)3*cut];
      const double *along = &direction[(size_t)3*kept];
      double u = 0.5;
      // A point cut for its distance has the clearance crossing somewhere on
      // the edge; one cut for another reason (turned over, or left on no
      // whole triangle) need not, and the middle of the edge does for it.
      if (standingOf[(size_t)cut] >= 0.0 && standingOf[(size_t)cut] < clearance)
      {
        double lo = 0.0, hi = 1.0;
        for (int step = 0; step < 8; step++)
        {
          double mid = 0.5*(lo + hi);
          double x[3], closest[3], wall = 0.0;
          for (int k = 0; k < 3; k++)
          {
            x[k] = (1.0 - mid)*qk[k] + mid*qc[k];
          }
          if (standing(x, kept, cut, along, closest, wall) >= clearance)
          {
            lo = mid;
          }
          else
          {
            hi = mid;
          }
        }
        u = lo;
      }
      else if (standingOf[(size_t)cut] < 0.0)
      {
        // Inside a lumen: the wall between is crossed somewhere on the edge,
        // and the clearance beyond it. Bisect on the sign first, then stand
        // clear of that wall.
        double lo = 0.0, hi = 1.0;
        for (int step = 0; step < 8; step++)
        {
          double mid = 0.5*(lo + hi);
          double x[3], closest[3], wall = 0.0;
          for (int k = 0; k < 3; k++)
          {
            x[k] = (1.0 - mid)*qk[k] + mid*qc[k];
          }
          bool clear = outwardSign*implicit->EvaluateFunction(x) > 0.0 &&
              standing(x, kept, cut, along, closest, wall) >= clearance;
          if (clear)
          {
            lo = mid;
          }
          else
          {
            hi = mid;
          }
        }
        u = lo;
      }
      if (u < edgeFloor)
      {
        // The crossing is at the kept end: the fragment there would be a
        // sliver a twentieth of an edge across. The kept point goes next
        // round, unless it is on a cap rim; for this one the cut point is
        // held off it either way.
        if (rimOfPoint[(size_t)kept] >= 0)
        {
          numRimHeld++;
        }
        else if (!demote[(size_t)kept])
        {
          demote[(size_t)kept] = true;
          demotedThisRound++;
        }
        u = edgeFloor;
      }
      u = std::min(u, 1.0 - edgeFloor);
      double x[3];
      for (int k = 0; k < 3; k++)
      {
        x[k] = (1.0 - u)*qk[k] + u*qc[k];
      }
      vtkIdType id = outerPoints->InsertNextPoint(x);
      scale.push_back(0.5*(array->GetValue(kept) + array->GetValue(cut)));
      cutFrom.push_back(kept);
      cutTo.push_back(cut);
      cutPoints[key] = id;
      return id;
    };
    numCutCells = 0;
    numDroppedCells = 0;
    for (vtkIdType cellId = 0; cellId < numWallCells; cellId++)
    {
      const vtkIdType *pts = &wallCellPts[(size_t)3*cellId];
      int numKept = 0;
      for (int j = 0; j < 3; j++)
      {
        if (!removed[(size_t)pts[j]])
        {
          numKept++;
        }
      }
      if (numKept == 3)
      {
        vtkIdType triangle[3] = {newId[(size_t)pts[0]], newId[(size_t)pts[1]], newId[(size_t)pts[2]]};
        addOuterTriangle(triangle, cellId);
        continue;
      }
      if (numKept == 0)
      {
        numDroppedCells++;
        continue;
      }
      numCutCells++;
      // Rotate so the odd one out is first: the one kept corner, or the one
      // cut corner. Rotation keeps the winding.
      vtkIdType a = pts[0], b = pts[1], c = pts[2];
      for (int j = 0; j < 3; j++)
      {
        bool odd = (numKept == 1) ? !removed[(size_t)pts[j]] : removed[(size_t)pts[j]];
        if (odd)
        {
          a = pts[j];
          b = pts[(j+1)%3];
          c = pts[(j+2)%3];
          break;
        }
      }
      if (numKept == 1)
      {
        vtkIdType ab = cutPointOn(a, b);
        vtkIdType ac = cutPointOn(a, c);
        vtkIdType triangle[3] = {newId[(size_t)a], ab, ac};
        addOuterTriangle(triangle, cellId);
      }
      else
      {
        vtkIdType ba = cutPointOn(b, a);
        vtkIdType ca = cutPointOn(c, a);
        vtkIdType first[3] = {ba, newId[(size_t)b], newId[(size_t)c]};
        vtkIdType second[3] = {ba, newId[(size_t)c], ca};
        addOuterTriangle(first, cellId);
        addOuterTriangle(second, cellId);
      }
    }
    numOuterPts = outerPoints->GetNumberOfPoints();
    numCutPts = numOuterPts - numKeptOriginal;

    // Anything that came out turned over against its source triangle is cut
    // at its kept corners. With nothing moved off the sheets this should not
    // happen, and it is counted so that the log says whether it did.
    lastFolded = 0;
    for (size_t outerCellId = 0; outerCellId < sourceCell.size(); outerCellId++)
    {
      const vtkIdType *pts = &outerTris[3*outerCellId];
      double q0[3], q1[3], q2[3], e1[3], e2[3], normal[3];
      outerPoints->GetPoint(pts[0], q0);
      outerPoints->GetPoint(pts[1], q1);
      outerPoints->GetPoint(pts[2], q2);
      vtkMath::Subtract(q1, q0, e1);
      vtkMath::Subtract(q2, q0, e2);
      vtkMath::Cross(e1, e2, normal);
      if (vtkMath::Dot(normal, &wallCellNormal[(size_t)3*sourceCell[outerCellId]]) > 0.0)
      {
        continue;
      }
      for (int j = 0; j < 3; j++)
      {
        if (pts[j] >= numKeptOriginal)
        {
          continue;
        }
        vtkIdType corner = cutFrom[(size_t)pts[j]];
        if (rimOfPoint[(size_t)corner] >= 0)
        {
          numRimHeld++;
        }
        else if (!demote[(size_t)corner])
        {
          demote[(size_t)corner] = true;
          lastFolded++;
        }
      }
    }
    numFolded += lastFolded;
    numDemoted += demotedThisRound;
    converged = (lastFolded == 0 && demotedThisRound == 0);
    if (converged || numRounds >= maxRounds)
    {
      break;
    }
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      if (demote[(size_t)ptId])
      {
        removed[(size_t)ptId] = true;
      }
    }
  }
  int numRemoved = numInvertedPts + numInLumen + numUnderWall + numPeninsula + numDemoted + numFolded +
      numNoThicknessOnly;

  // The boundary of what is left: the extruded cap rims, and the holes. The
  // walk builds links, so it is given a surface that is not added to after.
  auto clipped = vtkSmartPointer<vtkPolyData>::New();
  clipped->SetPoints(outerPoints);
  clipped->SetPolys(outerCells);
  std::vector<std::vector<vtkIdType> > loops;
  if (TGenUtils_ExtractBoundaryLoops(clipped, loops) != SV_OK)
  {
    fprintf(stderr,"Problem walking the edges of the trimmed outer surface\n");
    return SV_ERROR;
  }
  std::vector<int> rimOfOuter((size_t)numOuterPts, -1);
  for (vtkIdType ptId = 0; ptId < numKeptOriginal; ptId++)
  {
    rimOfOuter[(size_t)ptId] = rimOfPoint[(size_t)cutFrom[(size_t)ptId]];
  }
  caps.resize(rims.size());
  std::vector<bool> capFound(rims.size(), false);
  std::vector<std::vector<vtkIdType> > holes;
  for (size_t l = 0; l < loops.size(); l++)
  {
    const std::vector<vtkIdType> &loop = loops[l];
    int rim = -1;
    bool mixed = false;
    for (size_t m = 0; m < loop.size(); m++)
    {
      int r = rimOfOuter[(size_t)loop[m]];
      if (m == 0)
      {
        rim = r;
      }
      else if (r != rim)
      {
        mixed = true;
      }
    }
    if (mixed || (rim >= 0 && loop.size() != rims[(size_t)rim].size()))
    {
      double p[3];
      outerPoints->GetPoint(loop[0], p);
      fprintf(stderr,"A boundary loop of the trimmed outer surface at (%.5g, %.5g, %.5g) runs through a cap rim and a cut, so a neighbouring wall reaches a vessel end and the wall cannot be closed there\n",
          p[0], p[1], p[2]);
      return SV_ERROR;
    }
    if (rim >= 0)
    {
      caps[(size_t)rim].outerLoop = loop;
      capFound[(size_t)rim] = true;
    }
    else
    {
      holes.push_back(loop);
    }
  }
  for (size_t r = 0; r < rims.size(); r++)
  {
    if (!capFound[r])
    {
      fprintf(stderr,"The extruded rim of a cap did not come out as a boundary loop of the outer surface\n");
      return SV_ERROR;
    }
    TGenUtilsCapRim &cap = caps[r];
    cap.innerLoop = rims[r];
    for (int k = 0; k < 3; k++)
    {
      cap.origin[k] = 0.0;
    }
    for (size_t m = 0; m < rims[r].size(); m++)
    {
      double p[3];
      surface->GetPoint(rims[r][m], p);
      for (int k = 0; k < 3; k++)
      {
        cap.origin[k] += p[k];
      }
    }
    for (int k = 0; k < 3; k++)
    {
      cap.origin[k] /= (double)rims[r].size();
    }
    // Newell's normal is the rim's own normal for the order it is stored in;
    // the outward direction is its reverse.
    double normal[3] = {0.0, 0.0, 0.0};
    for (size_t m = 0; m < rims[r].size(); m++)
    {
      double p[3], q[3];
      surface->GetPoint(rims[r][m], p);
      surface->GetPoint(rims[r][(m+1)%rims[r].size()], q);
      normal[0] += (p[1]-q[1])*(p[2]+q[2]);
      normal[1] += (p[2]-q[2])*(p[0]+q[0]);
      normal[2] += (p[0]-q[0])*(p[1]+q[1]);
    }
    if (vtkMath::Normalize(normal) <= 0.0)
    {
      fprintf(stderr,"A cap rim of %zu points at (%.5g, %.5g, %.5g) encloses no area, so its plane cannot be found\n",
          rims[r].size(), cap.origin[0], cap.origin[1], cap.origin[2]);
      return SV_ERROR;
    }
    for (int k = 0; k < 3; k++)
    {
      cap.outward[k] = -normal[k];
    }
  }

  // Pair the holes that run alongside each other. Two loops are a pair when
  // most of the points of each have a point of the other within a wall and a
  // half - the two edges of a seam are a few percent of a wall apart across
  // the seam, more where the sheets cross at a shallow angle, and nothing else
  // comes that close along most of its length.
  size_t numHoles = holes.size();
  std::vector<double> holeBounds(6*numHoles, 0.0);
  std::vector<double> holeReach(numHoles, 0.0);
  for (size_t h = 0; h < numHoles; h++)
  {
    double *hb = &holeBounds[6*h];
    hb[0] = hb[2] = hb[4] = std::numeric_limits<double>::max();
    hb[1] = hb[3] = hb[5] = -std::numeric_limits<double>::max();
    for (size_t m = 0; m < holes[h].size(); m++)
    {
      double p[3];
      outerPoints->GetPoint(holes[h][m], p);
      for (int k = 0; k < 3; k++)
      {
        hb[2*k] = std::min(hb[2*k], p[k]);
        hb[2*k+1] = std::max(hb[2*k+1], p[k]);
      }
      holeReach[h] = std::max(holeReach[h], 1.5*scale[(size_t)holes[h][m]]);
    }
  }
  auto nearFraction = [&](size_t a, size_t b)
  {
    int numNear = 0;
    for (size_t i = 0; i < holes[a].size(); i++)
    {
      double p[3];
      outerPoints->GetPoint(holes[a][i], p);
      double within = 1.5*scale[(size_t)holes[a][i]];
      double within2 = within*within;
      for (size_t j = 0; j < holes[b].size(); j++)
      {
        double q[3];
        outerPoints->GetPoint(holes[b][j], q);
        if (vtkMath::Distance2BetweenPoints(p, q) <= within2)
        {
          numNear++;
          break;
        }
      }
    }
    return (double)numNear/(double)holes[a].size();
  };
  std::vector<int> partner(numHoles, -1);
  std::vector<double> partnerScore(numHoles, 0.0);
  for (size_t a = 0; a < numHoles; a++)
  {
    for (size_t b = a+1; b < numHoles; b++)
    {
      bool apart = false;
      double gap = std::max(holeReach[a], holeReach[b]);
      for (int k = 0; k < 3 && !apart; k++)
      {
        apart = holeBounds[6*a+2*k] > holeBounds[6*b+2*k+1] + gap ||
                holeBounds[6*b+2*k] > holeBounds[6*a+2*k+1] + gap;
      }
      if (apart)
      {
        continue;
      }
      double score = std::min(nearFraction(a, b), nearFraction(b, a));
      if (score < 0.5)
      {
        continue;
      }
      if (score > partnerScore[a] && score > partnerScore[b])
      {
        if (partner[a] >= 0)
        {
          partnerScore[(size_t)partner[a]] = 0.0;
          partner[(size_t)partner[a]] = -1;
        }
        if (partner[b] >= 0)
        {
          partnerScore[(size_t)partner[b]] = 0.0;
          partner[(size_t)partner[b]] = -1;
        }
        partner[a] = (int)b;
        partner[b] = (int)a;
        partnerScore[a] = partnerScore[b] = score;
      }
    }
  }

  // Close the holes. Each outer triangle is tagged with what it is - whole,
  // fragment, crease zip or seam band - and which hole it closes, for the
  // crossing check and for the surface written out after it.
  const int roleWhole = 0, roleFragment = 1, roleZip = 2, roleSeam = 3;
  std::vector<int> cellRole, cellHole;
  for (size_t outerCellId = 0; outerCellId < sourceCell.size(); outerCellId++)
  {
    const vtkIdType *pts = &outerTris[3*outerCellId];
    bool whole = pts[0] < numKeptOriginal && pts[1] < numKeptOriginal && pts[2] < numKeptOriginal;
    cellRole.push_back(whole ? roleWhole : roleFragment);
    cellHole.push_back(-1);
  }

  // Whether a triangle passes through the surface as it stands, which is what
  // the volume mesher refuses. The whole and fragment triangles are binned on
  // the same grid as the inner surface; the closing triangles accepted so far
  // are kept in a sparse map over the same bins, since they are few and are
  // added one hole at a time.
  vtkIdType numBaseCells = outerCells->GetNumberOfCells();
  std::vector<double> tri;
  std::vector<vtkIdType> triPts;
  std::vector<int> testedStamp;
  tri.reserve((size_t)9*numBaseCells);
  triPts.reserve((size_t)3*numBaseCells);
  auto appendTri = [&](const vtkIdType pts[3])
  {
    for (int j = 0; j < 3; j++)
    {
      triPts.push_back(pts[j]);
      double p[3];
      outerPoints->GetPoint(pts[j], p);
      tri.insert(tri.end(), p, p + 3);
    }
    testedStamp.push_back(0);
    return (vtkIdType)(triPts.size()/3) - 1;
  };
  for (vtkIdType cellId = 0; cellId < numBaseCells; cellId++)
  {
    appendTri(&outerTris[(size_t)3*cellId]);
  }
  auto binRange = [&](const double *t, int lo[3], int hi[3])
  {
    for (int k = 0; k < 3; k++)
    {
      double low = std::min(t[k], std::min(t[3+k], t[6+k]));
      double high = std::max(t[k], std::max(t[3+k], t[6+k]));
      lo[k] = binOf(low, k);
      hi[k] = binOf(high, k);
    }
  };
  std::vector<int> outerBinStart(numBins + 1, 0);
  for (vtkIdType cellId = 0; cellId < numBaseCells; cellId++)
  {
    int lo[3], hi[3];
    binRange(&tri[(size_t)9*cellId], lo, hi);
    for (int k = lo[2]; k <= hi[2]; k++)
    {
      for (int j = lo[1]; j <= hi[1]; j++)
      {
        for (int i = lo[0]; i <= hi[0]; i++)
        {
          outerBinStart[binIndex(i, j, k) + 1]++;
        }
      }
    }
  }
  for (size_t b = 0; b < numBins; b++)
  {
    outerBinStart[b + 1] += outerBinStart[b];
  }
  std::vector<vtkIdType> outerBinCells((size_t)outerBinStart[numBins], -1);
  {
    std::vector<int> cursor(outerBinStart.begin(), outerBinStart.end() - 1);
    for (vtkIdType cellId = 0; cellId < numBaseCells; cellId++)
    {
      int lo[3], hi[3];
      binRange(&tri[(size_t)9*cellId], lo, hi);
      for (int k = lo[2]; k <= hi[2]; k++)
      {
        for (int j = lo[1]; j <= hi[1]; j++)
        {
          for (int i = lo[0]; i <= hi[0]; i++)
          {
            outerBinCells[(size_t)cursor[binIndex(i, j, k)]++] = cellId;
          }
        }
      }
    }
  }
  std::map<size_t, std::vector<vtkIdType> > extraBins;
  auto insertExtra = [&](vtkIdType cellId)
  {
    int lo[3], hi[3];
    binRange(&tri[(size_t)9*cellId], lo, hi);
    for (int k = lo[2]; k <= hi[2]; k++)
    {
      for (int j = lo[1]; j <= hi[1]; j++)
      {
        for (int i = lo[0]; i <= hi[0]; i++)
        {
          extraBins[binIndex(i, j, k)].push_back(cellId);
        }
      }
    }
  };
  int testStamp = 0;
  // How many triangles of the surface a triangle passes through; the ids of
  // those are appended to hits when it is given. A triangle already in the
  // surface is skipped by its own index.
  auto crossings = [&](const double *f, const vtkIdType *fp, vtkIdType self,
      std::vector<vtkIdType> *hits)
  {
    testStamp++;
    int count = 0;
    int lo[3], hi[3];
    binRange(f, lo, hi);
    auto test = [&](vtkIdType other)
    {
      if (other == self || testedStamp[(size_t)other] == testStamp)
      {
        return;
      }
      testedStamp[(size_t)other] = testStamp;
      if (TrianglesCross(f, fp, &tri[(size_t)9*other], &triPts[(size_t)3*other]))
      {
        count++;
        if (hits != nullptr)
        {
          hits->push_back(other);
        }
      }
    };
    for (int k = lo[2]; k <= hi[2]; k++)
    {
      for (int j = lo[1]; j <= hi[1]; j++)
      {
        for (int i = lo[0]; i <= hi[0]; i++)
        {
          size_t bin = binIndex(i, j, k);
          for (int c = outerBinStart[bin]; c < outerBinStart[bin + 1]; c++)
          {
            test(outerBinCells[(size_t)c]);
          }
          std::map<size_t, std::vector<vtkIdType> >::const_iterator found = extraBins.find(bin);
          if (found != extraBins.end())
          {
            for (size_t c = 0; c < found->second.size(); c++)
            {
              test(found->second[c]);
            }
          }
        }
      }
    }
    return count;
  };
  // How many crossings a candidate closure would make: its triangles against
  // the surface, and against each other, since a fill can fold over itself.
  auto closureCrossings = [&](const std::vector<vtkIdType> &triangles)
  {
    size_t numTriangles = triangles.size()/3;
    std::vector<double> xyz(9*numTriangles);
    for (size_t t = 0; t < 3*numTriangles; t++)
    {
      outerPoints->GetPoint(triangles[t], &xyz[3*t]);
    }
    int total = 0;
    for (size_t t = 0; t < numTriangles; t++)
    {
      total += crossings(&xyz[9*t], &triangles[3*t], -1, nullptr);
    }
    for (size_t t = 0; t < numTriangles; t++)
    {
      for (size_t u = t+1; u < numTriangles; u++)
      {
        if (TrianglesCross(&xyz[9*t], &triangles[3*t], &xyz[9*u], &triangles[3*u]))
        {
          total++;
        }
      }
    }
    return total;
  };

  // Each hole is closed by whichever of its closures passes through the
  // surface least. A hole on its own is the gap at a crease and has one
  // closure, the zip. A pair of holes are the two sides of a seam, and have
  // three: the band with the second loop walked forwards, the band with it
  // walked backwards, and - when the pairing was wrong and they are two
  // creases that happen to run close - each zipped on its own. Preference
  // among equals is the band of less area, then the other band, then the
  // zips.
  struct Closure
  {
    std::vector<vtkIdType> triangles;
    std::vector<int> holeOf;
    int role;
    const char *how;
    int crossings;
  };
  int numZipped = 0, numJoined = 0, numFellBack = 0;
  size_t largestHole = 0;
  std::vector<std::pair<double,vtkIdType> > holeSeeds;
  std::vector<std::string> closureLines;
  for (size_t h = 0; h < numHoles; h++)
  {
    largestHole = std::max(largestHole, holes[h].size());
    holeSeeds.push_back(std::make_pair(-(double)holes[h].size(), holes[h][0]));
    if (partner[h] >= 0 && (size_t)partner[h] < h)
    {
      continue;
    }
    std::vector<Closure> options;
    if (partner[h] < 0)
    {
      Closure zip;
      if (TriangulateLoopByLeastArea(outerPoints, holes[h], zip.triangles, false) != SV_OK)
      {
        fprintf(stderr,"Problem zipping a crease in the outer wall\n");
        return SV_ERROR;
      }
      zip.holeOf.assign(zip.triangles.size()/3, (int)h);
      zip.role = roleZip;
      zip.how = "zipped";
      options.push_back(zip);
    }
    else
    {
      size_t g = (size_t)partner[h];
      std::vector<vtkIdType> bands[2];
      double areas[2];
      if (BuildSeamBands(outerPoints, holes[h], holes[g], bands, areas) != SV_OK)
      {
        fprintf(stderr,"Problem joining the two sides of a seam in the outer wall\n");
        return SV_ERROR;
      }
      int lesser = (areas[0] <= areas[1]) ? 0 : 1;
      for (int which = 0; which < 2; which++)
      {
        Closure band;
        int b = (which == 0) ? lesser : 1 - lesser;
        band.triangles = bands[b];
        band.holeOf.assign(band.triangles.size()/3, (int)h);
        band.role = roleSeam;
        band.how = (which == 0) ? "seam band" : "seam band the other way round";
        options.push_back(band);
      }
      Closure zips;
      std::vector<vtkIdType> second;
      if (TriangulateLoopByLeastArea(outerPoints, holes[h], zips.triangles, true) == SV_OK &&
          TriangulateLoopByLeastArea(outerPoints, holes[g], second, true) == SV_OK)
      {
        zips.holeOf.assign(zips.triangles.size()/3, (int)h);
        zips.holeOf.insert(zips.holeOf.end(), second.size()/3, (int)g);
        zips.triangles.insert(zips.triangles.end(), second.begin(), second.end());
        zips.role = roleZip;
        zips.how = "zipped each on its own";
        options.push_back(zips);
      }
    }
    size_t chosen = 0;
    for (size_t o = 0; o < options.size(); o++)
    {
      options[o].crossings = closureCrossings(options[o].triangles);
      if (options[o].crossings < options[chosen].crossings)
      {
        chosen = o;
      }
    }
    const Closure &closure = options[chosen];
    if (chosen > 0)
    {
      numFellBack++;
    }
    if (chosen > 0 || closure.crossings > 0)
    {
      double p[3];
      outerPoints->GetPoint(holes[h][0], p);
      std::string line;
      char head[256];
      snprintf(head, sizeof(head), "    hole %zu (%zu points) at (%.5g, %.5g, %.5g): %s, %d triangles crossing the surface",
          h, holes[h].size(), p[0], p[1], p[2], closure.how, closure.crossings);
      line = head;
      int numOthers = 0;
      for (size_t o = 0; o < options.size(); o++)
      {
        if (o == chosen)
        {
          continue;
        }
        char other[128];
        snprintf(other, sizeof(other), "%s %s would cross %d", (numOthers == 0) ? ";" : ",",
            options[o].how, options[o].crossings);
        line += other;
        numOthers++;
      }
      closureLines.push_back(line);
    }
    if (closure.role == roleSeam)
    {
      numJoined++;
    }
    else
    {
      numZipped += (partner[h] < 0) ? 1 : 2;
    }
    for (size_t t = 0; t < closure.triangles.size()/3; t++)
    {
      const vtkIdType *pts = &closure.triangles[3*t];
      outerCells->InsertNextCell(3, pts);
      cellRole.push_back(closure.role);
      cellHole.push_back(closure.holeOf[t]);
      insertExtra(appendTri(pts));
    }
  }
  outer->Initialize();
  outer->SetPoints(outerPoints);
  outer->SetPolys(outerCells);

  // What is left crossing once every hole is closed: every triangle is
  // tested against the whole surface, and what it crosses is marked with it.
  // The closing triangles are counted per hole, so that the log says which
  // holes the volume mesher will refuse and how they were closed; the whole
  // and fragment triangles are counted on their own, because a sheet crossing
  // another sheet is a fold the cut did not reach - the cut only looks at
  // where the points stand, and a triangle whose three points all stand
  // clear can still run through a wall between them.
  int numCrossing = 0, numSheetCrossing = 0;
  vtkIdType firstSheetCrossing = -1;
  std::vector<int> cellCrossing(cellRole.size(), 0);
  std::vector<int> holeCrossing(numHoles, 0);
  std::vector<vtkIdType> holeFirstCrossing(numHoles, -1);
  {
    vtkIdType numOuterCells = (vtkIdType)cellRole.size();
    std::vector<vtkIdType> hits;
    for (vtkIdType cellId = 0; cellId < numOuterCells; cellId++)
    {
      hits.clear();
      int count = crossings(&tri[(size_t)9*cellId], &triPts[(size_t)3*cellId], cellId, &hits);
      if (count == 0)
      {
        continue;
      }
      cellCrossing[(size_t)cellId] = 1;
      for (size_t c = 0; c < hits.size(); c++)
      {
        cellCrossing[(size_t)hits[c]] = 1;
      }
      if (cellId < numBaseCells)
      {
        numSheetCrossing++;
        if (firstSheetCrossing < 0)
        {
          firstSheetCrossing = cellId;
        }
        continue;
      }
      numCrossing++;
      int hole = cellHole[(size_t)cellId];
      if (hole >= 0)
      {
        holeCrossing[(size_t)hole]++;
        if (holeFirstCrossing[(size_t)hole] < 0)
        {
          holeFirstCrossing[(size_t)hole] = cellId;
        }
      }
    }
  }

  // The surface with its tags, for looking at what the log can only list.
  {
    auto roleArray = vtkSmartPointer<vtkIntArray>::New();
    roleArray->SetName("TrimRole");
    auto holeArray = vtkSmartPointer<vtkIntArray>::New();
    holeArray->SetName("TrimHole");
    auto crossingArray = vtkSmartPointer<vtkIntArray>::New();
    crossingArray->SetName("TrimCrossing");
    for (size_t cellId = 0; cellId < cellRole.size(); cellId++)
    {
      roleArray->InsertNextValue(cellRole[cellId]);
      holeArray->InsertNextValue(cellHole[cellId]);
      crossingArray->InsertNextValue(cellCrossing[cellId]);
    }
    outer->GetCellData()->AddArray(roleArray);
    outer->GetCellData()->AddArray(holeArray);
    outer->GetCellData()->AddArray(crossingArray);
    char trimmedFile[] = "wall_outer_trimmed.vtp";
    TGenUtils_WriteVTP(trimmedFile, outer);
  }

  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  fprintf(stdout,"Wall outer surface by extrusion, trimmed where it runs inside the wall:\n");
  fprintf(stdout,"  %lld points extruded; %d cut: %d on a triangle the extrusion turned over (%d triangles), %d inside a lumen, %d within %.3g of another sheet's wall, %d left on no whole triangle, %d whose edges crossed the clearance within %.2g of their length, %d on a fragment that came out turned over, %d for having no thickness alone; %d had no thickness\n",
      (long long)numPts, numRemoved, numInvertedPts, numInvertedCells, numInLumen, numUnderWall,
      clearance, numPeninsula, numDemoted, edgeFloor, numFolded, numNoThicknessOnly, numNoThickness);
  fprintf(stdout,"  %d rounds of cutting%s; the last left %d fragments turned over%s, and held %d cap rim corners that would otherwise have been cut\n",
      numRounds, converged ? ", converged" : " - the bound, not convergence: the last round's cuts were not made",
      lastFolded, (lastFolded > 0) ? " - the volume mesher will refuse them" : "", numRimHeld);
  fprintf(stdout,"  %d triangles cut through, %d dropped whole, %lld cut points added on the clearance crossing of their edges\n",
      numCutCells, numDroppedCells, (long long)numCutPts);
  fprintf(stdout,"  %lld standing queries over a %d x %d x %d grid of %.4g bins (reach %.4g), %.1f s for the first pass over every point\n",
      numQueries, gridSize[0], gridSize[1], gridSize[2], binSize, reach, firstPassSeconds);
  fprintf(stdout,"  %zu cap rims kept whole; %zu holes: %d zipped across a crease, %d pairs joined as the two sides of a seam, the largest %zu points around; %d closed by a fallback because the preferred closure crossed the surface more\n",
      rims.size(), numHoles, numZipped, numJoined, largestHole, numFellBack);
  for (size_t l = 0; l < closureLines.size(); l++)
  {
    fprintf(stdout,"%s\n", closureLines[l].c_str());
  }
  fprintf(stdout,"  %d closing triangles and %d sheet triangles pass through the surface%s; the surface is written to wall_outer_trimmed.vtp with TrimRole (0 whole, 1 fragment, 2 zip, 3 seam band), TrimHole and TrimCrossing on its cells\n",
      numCrossing, numSheetCrossing, (numCrossing + numSheetCrossing > 0) ? " - the volume mesher will refuse them" : "");
  if (numSheetCrossing > 0)
  {
    const double *f = &tri[(size_t)9*firstSheetCrossing];
    fprintf(stdout,"    a sheet triangle crossing the surface is a fold the cut did not reach, the first at (%.5g, %.5g, %.5g)\n",
        (f[0]+f[3]+f[6])/3.0, (f[1]+f[4]+f[7])/3.0, (f[2]+f[5]+f[8])/3.0);
  }
  for (size_t h = 0; h < numHoles; h++)
  {
    if (holeCrossing[h] == 0)
    {
      continue;
    }
    const double *f = &tri[(size_t)9*holeFirstCrossing[h]];
    fprintf(stdout,"    hole %zu (%s, %zu points): %d of its closing triangles cross the surface, the first at (%.5g, %.5g, %.5g)\n",
        h, (cellRole[(size_t)holeFirstCrossing[h]] == roleZip) ? "zipped" : "seam band", holes[h].size(),
        holeCrossing[h], (f[0]+f[3]+f[6])/3.0, (f[1]+f[4]+f[7])/3.0, (f[2]+f[5]+f[8])/3.0);
  }
  if (!holeSeeds.empty())
  {
    std::sort(holeSeeds.begin(), holeSeeds.end());
    size_t shown = std::min(holeSeeds.size(), (size_t)8);
    for (size_t h = 0; h < shown; h++)
    {
      double p[3];
      outerPoints->GetPoint(holeSeeds[h].second, p);
      fprintf(stdout,"    hole of %d points at (%.5g, %.5g, %.5g)\n",
          (int)(-holeSeeds[h].first), p[0], p[1], p[2]);
    }
  }
  fprintf(stdout,"  the outer wall is %lld points and %lld triangles, %.1f s\n",
      (long long)outer->GetNumberOfPoints(), (long long)outer->GetNumberOfCells(), seconds);
  numUnresolved = lastFolded + numCrossing + numSheetCrossing;
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
    fprintf(stdout,"    this is the construction, not the shape: every outer point was put a wall from the interface, so a shortfall is the construction giving it back - a vertex normal leaning against its facets reads a percent or two short, and a point left inside a fold reads far shorter\n");

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
