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

#include "SimVascular.h"

#include "sv_TetGenMeshObject.h"
#include "sv_SolidModel.h"
#include "sv_misc_utils.h"
#include "sv_polydatasolid_utils.h"

#include "sv_tetgenmesh_utils.h"

#include "sv_sys_geom.h"
#include "vtkGeometryFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkThreshold.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkCenterOfMass.h"

#ifdef SV_USE_VMTK
  #include "sv_vmtk_utils.h"
  #include "vtkvmtkPolyDataToUnstructuredGridFilter.h"
#endif

#ifdef SV_USE_MMG
  #include "sv_mmg_mesh_utils.h"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <utility>
#include <vector>
#include <math.h>

// -----------
// cvTetGenMeshObject for python
// -----------
cvTetGenMeshObject::cvTetGenMeshObject() : cvMeshObject()
{
  inmesh_ = nullptr;
  outmesh_ = nullptr;
  polydatasolid_ = nullptr;
  inputug_ = nullptr;
  meshloaded_ = 0;
  loadedVolumeMesh_ = 0;
  //nodemap_ = nullptr;
  pts_ = nullptr;

  originalpolydata_ = nullptr;
  surfacemesh_ = nullptr;
  volumemesh_ = nullptr;
  boundarylayermesh_ = nullptr;
  innerblmesh_ = nullptr;
  wallmesh_ = nullptr;
  holelist_ = nullptr;
  regionlist_ = nullptr;
  regionsizelist_ = nullptr;

  meshFileName_[0] = '\0';
  solidFileName_[0] = '\0';

  solidmodeling_kernel_ = SM_KT_POLYDATA;
  numModelRegions_ = 0;
  numBoundaryRegions_ = 0;

  //All the different mesh options. Originally set to zero. Changed by calls to object from GUI
  //Specifically ->SetMeshOptions()
  meshoptions_.surfacemeshflag=0;
  meshoptions_.volumemeshflag=0;
  meshoptions_.nomerge=0;
  meshoptions_.quiet=0;
  meshoptions_.docheck=0;
  meshoptions_.verbose=0;
  meshoptions_.diagnose=0;
  meshoptions_.nobisect=0;
  meshoptions_.optlevel=0;
  meshoptions_.maxedgesize=0;
  meshoptions_.epsilon=0;
  meshoptions_.minratio=0;
  meshoptions_.mindihedral=0;
  meshoptions_.coarsenpercent=0;
  meshoptions_.boundarylayermeshflag=0;
  meshoptions_.numsublayers=0;
  meshoptions_.blthicknessfactor=0;
  meshoptions_.sublayerratio=0;
  meshoptions_.useconstantblthickness=0;
  meshoptions_.newregionboundarylayer=0;
  meshoptions_.boundarylayerdirection=1;
  meshoptions_.wallmeshflag=0;
  meshoptions_.walltetgenshell=0;
  meshoptions_.wallthickness=0.0;
  meshoptions_.numwallsublayers=2;
  meshoptions_.wallthicknesssmoothingiterations=5;
  meshoptions_.wallthicknesscurvaturefactor=0.8;
  meshoptions_.refinement=0;
  meshoptions_.refinedsize=0;
  meshoptions_.sphereradius=0;
  meshoptions_.cylinderradius=0;
  meshoptions_.cylinderlength=0;
  meshoptions_.functionbasedmeshing=0;
  meshoptions_.secondarrayfunction=0;
  meshoptions_.meshwallfirst=0;
  meshoptions_.startwithvolume=0;
  meshoptions_.refinecount=0;
  meshoptions_.numberofholes=0;
  meshoptions_.numberofregions=0;
  meshoptions_.allowMultipleRegions=false;
#ifdef SV_USE_MMG
  meshoptions_.usemmg=1;
#else
  meshoptions_.usemmg=0;
#endif
  meshoptions_.hausd=0;
  for (int i=0;i<3;i++)
  {
    meshoptions_.spherecenter[i] = 0;
    meshoptions_.cylindercenter[i] = 0;
    meshoptions_.cylindernormal[i] = 0;
  }
}

// -----------
// cvTetGenMeshObject
// -----------
//
cvTetGenMeshObject::cvTetGenMeshObject( const cvTetGenMeshObject& sm ) : cvMeshObject()
{

  // Copy( sm );   // relying on automatic upcast to cvMeshObject*

}

// ------------
// ~cvTetGenMeshObject
// ------------
/**
 * @brief Destructor for cvTetGenMeshObject
 */

cvTetGenMeshObject::~cvTetGenMeshObject()
{
  if (inmesh_ != nullptr)
    delete inmesh_;

  if (outmesh_ != nullptr)
    delete outmesh_;

  if (surfacemesh_ != nullptr)
    surfacemesh_->Delete();

  if (volumemesh_ != nullptr)
    volumemesh_->Delete();

  if (polydatasolid_ != nullptr)
    polydatasolid_->Delete();

  if (inputug_ != nullptr)
    inputug_->Delete();

  if (originalpolydata_ != nullptr)
    originalpolydata_->Delete();

  if (boundarylayermesh_ != nullptr)
    boundarylayermesh_->Delete();

  if (innerblmesh_ != nullptr)
    innerblmesh_->Delete();

  if (wallmesh_ != nullptr)
    wallmesh_->Delete();

  if (holelist_ != nullptr)
    holelist_->Delete();

  if (regionlist_ != nullptr)
    regionlist_->Delete();

  if (regionsizelist_ != nullptr)
    regionsizelist_->Delete();
}

//-----------------------------
// EnableSizeFunctionBasedMesh
//-----------------------------
// Enable meshing based on a size function.
//
// This is only called by the Python API.
//
void cvTetGenMeshObject::EnableSizeFunctionBasedMesh()
{
  meshoptions_.functionbasedmeshing = 1;
}

//------------------------------
// DisableSizeFunctionBasedMesh  
//------------------------------
// Disable meshing based on a size function.
//
// This is only called by the Python API.
//
void cvTetGenMeshObject::DisableSizeFunctionBasedMesh()
{
  meshoptions_.functionbasedmeshing = 0;
}

//--------------------------------
// SizeFunctionBasedMeshIsEnabled
//--------------------------------
//
bool cvTetGenMeshObject::SizeFunctionBasedMeshIsEnabled()
{
  return (meshoptions_.functionbasedmeshing == 1);
}

int cvTetGenMeshObject::SetMeshFileName( const char* meshFileName )
{
  if (meshFileName != nullptr)
    sprintf(meshFileName_, "%s", meshFileName);
  else
    meshFileName_[0] = '\0';

  return SV_OK;
}

void cvTetGenMeshObject::SetAllowMultipleRegions(bool value) 
{
  meshoptions_.allowMultipleRegions = value;
}

int cvTetGenMeshObject::SetSolidFileName( const char* solidFileName )
{
  if (solidFileName != nullptr)
    sprintf(solidFileName_, "%s", solidFileName);
  else
    solidFileName_[0] = '\0';

  return SV_OK;
}

int cvTetGenMeshObject::Print()
{
  int num_nodes = 0;
  int nMeshFaces = 0;
  int num_elems = 0;
  int nMeshEdges = 0;
  int nRegion = 1;
  int nFace = 0;
  int nEdge = 0;
  int nVertex = 0;

  if (polydatasolid_ == nullptr)
  {
    fprintf(stderr,"Solid has not been loaded\n");
    return SV_ERROR;
  }
  if(meshoptions_.surfacemeshflag && !meshoptions_.volumemeshflag)
  {
    num_nodes = polydatasolid_->GetNumberOfPoints();
    nMeshFaces = polydatasolid_->GetNumberOfCells();
  }
  else if(meshoptions_.boundarylayermeshflag)
  {
    if (volumemesh_ == nullptr)
    {
      fprintf(stderr,"Mesh has not been created\n");
      return SV_ERROR;
    }
    num_nodes = volumemesh_->GetNumberOfPoints();
    num_elems = volumemesh_->GetNumberOfCells();
    nMeshFaces = 4.0*(volumemesh_->GetNumberOfCells());
    nMeshEdges = 3.0*(nMeshFaces);
  }
  else
  {
    if (outmesh_ == nullptr)
    {
      fprintf(stderr,"Mesh has not been created\n");
      return SV_ERROR;
    }
    num_nodes = outmesh_->numberofpoints;
    num_elems = outmesh_->numberoftetrahedra;
    nMeshFaces = outmesh_->numberoftrifaces;
    nMeshEdges = outmesh_->numberofedges;
    nFace = outmesh_->numberoftrifaces;
    nEdge = outmesh_->numberofregions;
    nVertex = outmesh_->numberofpoints;
  }
   /* output the statistics */

  nFace = originalpolydata_->GetNumberOfPolys();
  nVertex = originalpolydata_->GetNumberOfPoints();
  nEdge = 3*(originalpolydata_->GetNumberOfPolys());

  fprintf(stdout,"\nMESH STATISTICS:\n");
  fprintf(stdout,"  elements         = %i\n",num_elems);
  fprintf(stdout,"  nodes            = %i\n",num_nodes);
  fprintf(stdout,"  mesh edges       = %i\n",nMeshEdges);
  fprintf(stdout,"  mesh faces       = %i\n",nMeshFaces);
  fprintf(stdout,"\nMODEL STATISTICS:\n");
  fprintf(stdout,"  material regions = %i\n",nRegion);
  fprintf(stdout,"  edges            = %i\n",nEdge);
  fprintf(stdout,"  vertices         = %i\n\n",nVertex);

  //Reset all the meshoptions to be ready for next mesh
  meshoptions_.volumemeshflag = 0;
  meshoptions_.surfacemeshflag = 0;
  meshoptions_.boundarylayermeshflag = 0;

  return SV_OK;
}

// ----
// Copy
// ----

cvMeshObject *cvTetGenMeshObject::Copy() const
{
  cvTetGenMeshObject *result = new cvTetGenMeshObject( *this );
  return result;
}


// -----------
// Update
// -----------
/**
 * @brief Function that updates in between routines to make sure solid
 * is loaded
 * @return SV_OK if executed correctly
 */

int cvTetGenMeshObject::Update() {

  //I'm not sure why this is used. Would really mesh with some of the
  //meshing functions to Update the mesh or solid periodically.
  //However, must return SV_OK, otherwise running of anything breaks.

  return SV_OK;

}

// -----------------------
// GetSolid
// -----------------------
//
/**
 * @brief Function that returns the polydatasolid member data loaded
 * @return SV_OK if executed correctly
 */
cvPolyData* cvTetGenMeshObject::GetSolid() {

  if (polydatasolid_ == nullptr)
  {
    //Mesh must be created first
    return SV_ERROR;
  }

  cvPolyData* result =nullptr;

  result = new cvPolyData(polydatasolid_);

  return result;

}

//-----------
// HasSolid
//-----------
// Check if the mesh has a solid model defined for it.
//
bool cvTetGenMeshObject::HasSolid() {

  if (polydatasolid_ == nullptr) {
    return false;
  }

  return true;
}


// -----------------------
// GetPolyData
// -----------------------
//
/**
 * @brief Function that returns the polydatasolid member data loaded
 * @return SV_OK if executed correctly
 */
cvPolyData* cvTetGenMeshObject::GetPolyData() {

  if (surfacemesh_ == nullptr)
  {
    //Mesh must be created first
    return SV_ERROR;
  }

  cvPolyData* result =nullptr;

  result = new cvPolyData(surfacemesh_);

  return result;

}

//---------------
// HasSurfaceMesh
//---------------
// Check if a surface mesh has been generated.
//
bool cvTetGenMeshObject::HasSurfaceMesh()
{
  return (surfacemesh_ != nullptr);
}

//---------------
// HasVolumeMesh
//---------------
// Check if a volume mesh has been generated.
//
bool cvTetGenMeshObject::HasVolumeMesh()
{
  return (volumemesh_ != nullptr);
}

// --------------------
//  GetUnstructuredGrid
// --------------------
/**
 * @brief Function that returns the volume mesh if the mesh has been
 * calculated
 * @return the mesh if executed correctly
 */

cvUnstructuredGrid* cvTetGenMeshObject::GetUnstructuredGrid() {

  // recall the node numbers start at 1 in the P_id,
  // but in vtkPolyData file they start at 0.
  if (volumemesh_ == nullptr)
  {
    //Mesh must be created first
    return nullptr;
  }

  cvUnstructuredGrid *result = nullptr;

  result = new cvUnstructuredGrid(volumemesh_);

  return result;

}
/**
 * @brief Function that writes the adjacency between tetrahedral elements
 * @param *filename char holding the name of the file to be written to
 * @return SV_OK if executed correctly
 * @note This function is not necessarily needed anymore. The functionality
 * to extract the adjacency from the vtu is added in the presolver
 */

int cvTetGenMeshObject::WriteMetisAdjacency(char *filename) {

  //No longer need to write adjacency
  return SV_OK;

  if (filename == nullptr) {
        return SV_ERROR;
  }

  if (meshoptions_.volumemeshflag)
  {
    // open the output file
    if (openOutputFile(filename) != SV_OK) return SV_ERROR;

    //Different way to write Adjacency, used for comparison
    int i;
    int numCells;
    int *xadj;
    int *adjacency;
    vtkIdType cellId;
    vtkIdType meshCellId;
    vtkIdType p1,p2,p3;
    vtkIdType npts = 0;
    const vtkIdType *pts;
    vtkSmartPointer<vtkCellArray> volCells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkIntArray> globalIds = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    volumemesh_->BuildLinks();

    if (VtkUtils_UGCheckArrayName(volumemesh_,1,"GlobalElementID") != SV_OK)
    {
      fprintf(stderr,"Array name 'GlobalElementID' does not exist in volume mesh. \
		      Something wrong with ids on mesh");
      return SV_ERROR;
    }
    globalIds = vtkIntArray::SafeDownCast(volumemesh_->GetCellData()->
      GetScalars("GlobalElementID"));
    numCells = volumemesh_->GetNumberOfCells();
    volCells = volumemesh_->GetCells();

    xadj = new int[numCells+1];
    adjacency = new int[4*numCells];
    int adj = 0;
    int xcheck = 0;
    xadj[xcheck] = 0;

    ptIds->SetNumberOfIds(3);
    for (cellId = 0;cellId<numCells;cellId++)
    {
      meshCellId = globalIds->LookupValue(cellId+1);
      volumemesh_->GetCellPoints(meshCellId,npts,pts);
      for (i=0;i < npts; i++)
      {
	p1 = pts[i];
	p2 = pts[(i+1)%(npts)];
	p3 = pts[(i+2)%(npts)];

	ptIds->InsertId(0,p1);
	ptIds->InsertId(1,p2);
	ptIds->InsertId(2,p3);

	volumemesh_->GetCellNeighbors(meshCellId,ptIds,cellIds);

	//If it is zero, it is a face on the exterior. Otherwise, it has
	//neighbors
	if (cellIds->GetNumberOfIds() != 0)
	{
	  adjacency[adj++] = (int) globalIds->GetValue(cellIds->GetId(0))-1;
	}

      }
      xadj[++xcheck] = adj;
    }

    gzprintf(fp_,"xadj: %i\n",numCells+1);
    gzprintf(fp_,"adjncy: %i\n",adj);

    for (i=0;i < numCells+1; i++)
    {
	gzprintf(fp_,"%i\n",xadj[i]);
    }
    for (i=0;i < adj; i++)
    {
	gzprintf(fp_,"%i\n",adjacency[i]);
    }

    delete [] xadj;
    delete [] adjacency;

    return closeOutputFile();
  }
  else
  {
  //  fprintf(stdout,"No volumemesh, not writing adjacency file\n");
    return SV_OK;
  }
}

int cvTetGenMeshObject::GetNodeCoords(int node)
{
  if (volumemesh_ == 0)
  {
    fprintf(stderr,"Mesh needs to be computed before node coords can be retrieved\n");
    return SV_ERROR;
  }
  nodeID_ = node;
  nodeX_ = 0; nodeY_ = 0; nodeZ_ = 0;
  double pts[3];

  volumemesh_->GetPoint(node,pts);
  nodeX_ = pts[0];
  nodeY_ = pts[1];
  nodeZ_ = pts[2];

  return SV_OK;
}

// --------------------
//  LoadModel
// --------------------
/**
 * @brief Function that loads a model from solid and store it in the
 * member data polydatasolid
 * @param *filename char holding the name of the file to read
 * @return SV_OK if executed correctly
 * @note uses the PolyDataUtils to read the solid: see cv_polydatasolid_utils
 */

int cvTetGenMeshObject::LoadModel(char *filename) {

  if (filename == nullptr) {
    return SV_ERROR;
  }

  // must load model before mesh!
  if (inmesh_ != nullptr) {
    return SV_ERROR;
  }

  //polydatasolid cannot already exist
  if (polydatasolid_ != nullptr)
  {
    polydatasolid_->Delete();
  }
  if (originalpolydata_ != nullptr) {
    originalpolydata_->Delete();
  }

  polydatasolid_ = vtkPolyData::New();
  originalpolydata_ = vtkPolyData::New();
  if (PlyDtaUtils_ReadNative(filename,polydatasolid_) != SV_OK) {
    return SV_ERROR;
  }

  originalpolydata_->DeepCopy(polydatasolid_);
  return SV_OK;

}

int cvTetGenMeshObject::LoadModel(vtkPolyData *pd) {

  if (pd == nullptr) {
    return SV_ERROR;
  }

  // must load model before mesh!
  if (inmesh_ != nullptr) {
    return SV_ERROR;
  }

  //polydatasolid cannot already exist
  if (polydatasolid_ != nullptr)
  {
    polydatasolid_->Delete();
  }
  if (originalpolydata_ != nullptr)
  {
    originalpolydata_->Delete();
  }

  originalpolydata_=vtkPolyData::New();
  originalpolydata_->DeepCopy(pd);

  polydatasolid_ = vtkPolyData::New();
  polydatasolid_->DeepCopy(pd);

  return SV_OK;

}

// --------------------
//  GetBoundaryFaces
// --------------------
/**
 * @brief Function to extract the boundaries for the member vtkPolyData
 * @param angle double that specifies the extraction angle. Any faces
 * with a difference between face normals larger than this angle will be
 * considered a separate face
 * @return *result: SV_ERROR is member data hasn't been loaded, or if the
 * GetBoundaryFaces function does not work properly. SV_OK is executed
 * properly
 */

int cvTetGenMeshObject::GetBoundaryFaces(double angle)
{
  if (polydatasolid_ == nullptr) {
    return SV_ERROR;
  }

  if(PlyDtaUtils_GetBoundaryFaces(polydatasolid_,angle,
	&numBoundaryRegions_) != SV_OK) {
    return SV_ERROR;
  }

  return SV_OK;
}

// --------------------
//  LoadMesh
// --------------------
/**
 * @brief Function to load the vtkUnstructuredGrid of a Mesh
 * @return SV_OK if function executes properly.
 * @note If a current volume
 * mesh exists, the current volumemesh is deleted and the new one is loaded
 */
int cvTetGenMeshObject::LoadMesh(char *filename,char *surfilename) {

  if (filename == nullptr) {
    return SV_ERROR;
  }
  if (volumemesh_ != nullptr)
    volumemesh_->Delete();
  volumemesh_ = vtkUnstructuredGrid::New();
  if (TGenUtils_LoadMesh(filename,volumemesh_) != SV_OK)
    return SV_ERROR;
  if (surfilename != 0)
  {
    if (surfacemesh_ != nullptr)
      surfacemesh_->Delete();

    surfacemesh_ = vtkPolyData::New();
    if (PlyDtaUtils_ReadNative(surfilename,surfacemesh_) != SV_OK)
      return SV_ERROR;
  }
  return SV_OK;

}

// --------------------
//  NewMesh
// --------------------
/**
 * @brief Function to prepare a mesh from the polydatasolid. Takes input
 * mesh and converts to tetgen structures in preparation.
 * @return *result: SV_ERROR if solid hasn't been loaded. If a current mesh
 * exists, it is deleted so that a new one can be made.
 */

int cvTetGenMeshObject::NewMesh() {

  // cant overwrite mesh
  if (inmesh_ != nullptr) {
    delete inmesh_;
  }
  if (outmesh_ != nullptr)
  {
    delete outmesh_;
  }

  // need solid to convert to new tetgen mesh
  if (polydatasolid_ == nullptr) {
    return SV_ERROR;
  }

  //Create new tetgen mesh objects and set first number of output mesh to 0
  inmesh_ = new tetgenio;
  inmesh_->firstnumber = 0;
  outmesh_ = new tetgenio;
  outmesh_->firstnumber = 0;

  //MarkerListName
  std::string markerListName;
  int useSizingFunction;
  int useBoundary;

  if (meshoptions_.boundarylayermeshflag)
  {
    useSizingFunction = 1;
    useBoundary = 0;
    markerListName = "CellEntityIds";
  }
  else if (meshoptions_.functionbasedmeshing || meshoptions_.refinement)
  {
    fprintf(stdout,"Using size function\n");
    useSizingFunction = 1;
    useBoundary = 1;
    markerListName = "ModelFaceID";
  }
  else
  {
    useSizingFunction = 0;
    useBoundary = 1;
    markerListName = "ModelFaceID";
  }

  vtkSmartPointer<vtkCleanPolyData> cleaner =
   vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(polydatasolid_);
  cleaner->Update();
  polydatasolid_->DeepCopy(cleaner->GetOutput());
  fprintf(stderr,"Converting to TetGen...\n");
  //Convert the polydata to tetgen for meshing with given option
  if (TGenUtils_ConvertSurfaceToTetGen(inmesh_,polydatasolid_) != SV_OK)
  {
    fprintf(stderr,"Error converting surface to tetgen object\n");
    return SV_ERROR;
  }

  // Add mesh sizing function
  if (useSizingFunction)
  {
    if (TGenUtils_AddPointSizingFunction(inmesh_,polydatasolid_,
          "MeshSizingFunction", meshoptions_.maxedgesize) != SV_OK)
    {
      fprintf(stderr,"Could not add mesh sizing function to mesh\n");
      return SV_ERROR;
    }
  }

  // Add facet markers
  if (useBoundary)
  {
    if(TGenUtils_AddFacetMarkers(inmesh_,polydatasolid_,
      markerListName) != SV_OK)
    {
      fprintf(stderr,"Could not add facet markers to mesh\n");
      return SV_ERROR;
    }
  }

  // Add holes
  if (meshoptions_.numberofholes > 0)
  {
    if (TGenUtils_AddHoles(inmesh_, holelist_) != SV_OK)
    {
      fprintf(stderr,"Could not add hole to mesh\n");
      return SV_ERROR;
    }
  }

  if (meshoptions_.numberofregions > 0)
  {
    if (TGenUtils_AddRegions(inmesh_, regionlist_, regionsizelist_) != SV_OK)
    {
      fprintf(stderr,"Could not add region to mesh\n");
      return SV_ERROR;
    }
  }

  //The mesh is now loaded, and TetGen is ready to be called
  meshloaded_ = 1;

  return SV_OK;
}

//------------------------------
// GenerateLocalSizeSizingArray
//------------------------------
//
int
cvTetGenMeshObject::GenerateLocalSizeSizingArray(int faceID, double edgeSize)
{
  meshoptions_.functionbasedmeshing = 1;

  if (TGenUtils_SetLocalMeshSize(polydatasolid_, faceID, edgeSize) != SV_OK) {
      return SV_ERROR;
  }   
  meshoptions_.secondarrayfunction = 1;

  return SV_OK;
}

// --------------------
//  SetMeshOptions
// --------------------
/**
 * @brief Function to set the options for tetgen. Store temporarily in
 * meshoptions_ object until the mesh is run
 * @param *flag char containing the flag to set
 * @param value if the flag requires a value, this double contains that
 * value to be set
 * @return *result: SV_ERROR if the mesh doesn't exist. New Mesh must be
 * called before the options can be set
 */

int cvTetGenMeshObject::SetMeshOptions(char *flags,int numValues,double *values) 
{
  if(!strncmp(flags,"GlobalEdgeSize",14)) {            //Global edge size
     if (numValues < 1)
       return SV_ERROR;

    meshoptions_.maxedgesize=values[0];
  }
  else if(!strncmp(flags,"LocalEdgeSize",13)) {

    if (numValues < 2)
    {
      fprintf(stderr,"Must give face id and local edge size\n");
      return SV_ERROR;
    }
    meshoptions_.functionbasedmeshing = 1;
    //Create a new mesh sizing function and call TGenUtils to compute function.
    //Store in the member data vtkDouble Array meshsizingfunction
    if (TGenUtils_SetLocalMeshSize(polydatasolid_,values[0],values[1]) != SV_OK)
      return SV_ERROR;
    meshoptions_.secondarrayfunction = 1;
  }
  else if(!strncmp(flags,"SurfaceMeshFlag",15)) {
#ifdef SV_USE_VMTK
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.surfacemeshflag = values[0];
#else
      fprintf(stderr,"Plugin VMTK is not being used!\
	  In order to use surface meshing, plugin VMTK must be available!\n");
      return SV_ERROR;
#endif
  }
  else if(!strncmp(flags,"VolumeMeshFlag",14)) {
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.volumemeshflag = values[0];
  }
  else if(!strncmp(flags,"QualityRatio",12)) {//q
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.minratio=values[0];
  }
  else if(!strncmp(flags,"MinDihedral",11)) {//q
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.mindihedral=values[0];
  }
  else if(!strncmp(flags,"Optimization",12)) {//O
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.optlevel=(int)values[0];
  }
  else if(!strncmp(flags,"Epsilon",7)) {//T
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.epsilon=values[0];
  }
  else if(!strncmp(flags,"CoarsenPercent",14)) {//R
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.coarsenpercent=values[0]/100;
  }
  else if(!strncmp(flags,"AddHole",7)) {
    if (numValues < 3)
    {
      fprintf(stderr,"Must provide x,y,z coordinate of hole\n");
      return SV_ERROR;
    }
    meshoptions_.numberofholes++;
    if (holelist_ == nullptr)
      holelist_ = vtkPoints::New();
    holelist_->InsertNextPoint(values[0], values[1], values[2]);
  }
  else if(!strncmp(flags,"AddSubDomain",12)) {
    if (numValues < 4)
    {
      fprintf(stderr,"Must provide x,y,z, size of region, and coordinate of region\n");
      return SV_ERROR;
    }
    meshoptions_.numberofregions++;
    if (regionsizelist_ == nullptr)
      regionsizelist_ = vtkDoubleArray::New();
    regionsizelist_->InsertNextTuple1(values[0]);
    if (regionlist_ == nullptr)
      regionlist_ = vtkPoints::New();
    regionlist_->InsertNextPoint(values[1], values[2], values[3]);

  }
  else if(!strncmp(flags,"Verbose",7)) {//V
    meshoptions_.verbose=1;
  }
  else if(!strncmp(flags,"NoMerge",7)) {//M
    meshoptions_.nomerge=1;
  }
  else if(!strncmp(flags,"Check",5)) {//C
    meshoptions_.docheck=1;
  }
  else if(!strncmp(flags,"NoBisect",8)) {//Y
    meshoptions_.nobisect=1;
  }
  else if(!strncmp(flags,"Quiet",5)) {//Q
    meshoptions_.quiet=1;
  }
  else if(!strncmp(flags,"Diagnose",8)) {//d
    meshoptions_.diagnose=1;
  }
  else if(!strncmp(flags,"MeshWallFirst",13)) {//k
    meshoptions_.meshwallfirst=1;
  }
  else if(!strncmp(flags,"StartWithVolume",15)) {//r
    meshoptions_.startwithvolume=1;
  }
  else if(!strncmp(flags,"Hausd",5)) {
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.hausd=values[0];
  }
  else if (!strncmp(flags,"UseMMG",6)){
      if (numValues < 1)
        return SV_ERROR;
      meshoptions_.usemmg=values[0];
  }
  else if (!strncmp(flags,"NewRegionBoundaryLayer",22)) {
    meshoptions_.newregionboundarylayer=1;
  }
  else if (!strncmp(flags,"BoundaryLayerDirection",22)) {
    if (numValues < 1)
      return SV_ERROR;
    meshoptions_.boundarylayerdirection=values[0];
  }
  else if (!strncmp(flags,"GenerateWallMesh",16)) {
    // With no value the option acts as a flag enabling wall meshing.
    if (numValues < 1)
    {
      meshoptions_.wallmeshflag=1;
    }
    else
    {
      meshoptions_.wallmeshflag=(values[0] != 0.0);
    }
  }
  else if (!strncmp(flags,"WallMeshTetGenShell",19)) {
    // With no value the option acts as a flag selecting the TetGen shell fill.
    if (numValues < 1)
    {
      meshoptions_.walltetgenshell=1;
    }
    else
    {
      meshoptions_.walltetgenshell=(values[0] != 0.0);
    }
  }
  else if (!strncmp(flags,"LocalWallThickness",18)) {
    if (numValues < 2)
    {
      fprintf(stderr,"Must give face id and local wall thickness\n");
      return SV_ERROR;
    }
    // The GUI, Python API and .msh file paths all set options through here,
    // so finiteness (a NaN passes a '<= 0' check) and the integer face id
    // are enforced once in this core path.
    if (!std::isfinite(values[0]) || values[0] != std::floor(values[0]))
    {
      fprintf(stderr,"The local wall thickness face id must be an integer\n");
      return SV_ERROR;
    }
    if (!std::isfinite(values[1]) || values[1] <= 0.0)
    {
      fprintf(stderr,"Local wall thickness must be a finite value greater than zero\n");
      return SV_ERROR;
    }
    localWallThickness_[(int)values[0]]=values[1];
  }
  // Must be checked before the shorter 'WallThickness' prefix.
  else if (!strncmp(flags,"WallThicknessSmoothingIterations",32)) {
    if (numValues < 1)
      return SV_ERROR;
    // The GUI, Python API and .msh file paths all set options through here,
    // so the valid range (integer 0-50, matching the GUI spin box) is
    // enforced once in this core path.
    if (!std::isfinite(values[0]) || values[0] < 0.0 || values[0] > 50.0 ||
        values[0] != std::floor(values[0]))
    {
      fprintf(stderr,"The number of wall thickness smoothing iterations must be an integer between 0 and 50\n");
      return SV_ERROR;
    }
    meshoptions_.wallthicknesssmoothingiterations=(int)values[0];
  }
  // Must be checked before the shorter 'WallThickness' prefix.
  else if (!strncmp(flags,"WallThicknessCurvatureFactor",28)) {
    if (numValues < 1)
      return SV_ERROR;
    // The GUI, Python API and .msh file paths all set options through here,
    // so the valid range (0.0-1.0, matching the GUI spin box) is enforced
    // once in this core path. A value of zero disables the clamp.
    if (!std::isfinite(values[0]) || values[0] < 0.0 || values[0] > 1.0)
    {
      fprintf(stderr,"The wall thickness curvature factor must be between 0 and 1\n");
      return SV_ERROR;
    }
    meshoptions_.wallthicknesscurvaturefactor=values[0];
  }
  else if (!strncmp(flags,"WallThickness",13)) {
    if (numValues < 1)
      return SV_ERROR;
    // A value of zero means the wall thickness is not set; whether a
    // positive wall thickness has been set is checked in GenerateMesh().
    // A NaN passes a '< 0' check, so finiteness is checked explicitly.
    if (!std::isfinite(values[0]) || values[0] < 0.0)
    {
      fprintf(stderr,"Wall thickness must be a finite value greater than zero\n");
      return SV_ERROR;
    }
    meshoptions_.wallthickness=values[0];
  }
  else if (!strncmp(flags,"NumberOfWallLayers",18)) {
    if (numValues < 1)
      return SV_ERROR;
    // Converting a NaN to int is undefined and a fractional value would
    // silently be truncated, so both are rejected here.
    if (!std::isfinite(values[0]) || values[0] < 1.0 ||
        values[0] != std::floor(values[0]))
    {
      fprintf(stderr,"The number of wall layers must be an integer of at least one\n");
      return SV_ERROR;
    }
    meshoptions_.numwallsublayers=(int)values[0];
  }
  else if (!strncmp(flags,"AllowMultipleRegions",20)) {
      meshoptions_.allowMultipleRegions = (int(values[0]) == 1);
  }
  else {
      fprintf(stderr,"%s: flag is not recognized\n",flags);
  }
  return SV_OK;
}

// --------------------------
//  ClearLocalWallThickness
// --------------------------
/**
 * @brief Removes all per-face wall thickness values set with the
 * 'LocalWallThickness' option. The values otherwise accumulate in the
 * mesher, so a caller applying a new set of options (e.g. the Python API
 * reusing a mesher) must clear them first or values from a previous mesh
 * generation would still be applied.
 */
void cvTetGenMeshObject::ClearLocalWallThickness()
{
  localWallThickness_.clear();
}

// --------------------
//  SetBoundaryLayer
// --------------------
/**
 * @brief Function to set the boundary layer of the mesh. In this case
 * @brief the only things set are the meshoptions_ parameters.
 * @note The polydatasolid_ object is prepared for meshing by separating the
 * @note surface that is to have the boundary layer away from the rest of
 * @note the mesh.
 * @note Must have VMTK to use Boundary Layer!
 * @param type UNUSED
 * @param id This is the value of the region on which we want the BL.
 * @param side UNUSED
 * @param nL This is the number of layers for the BL mesh.
 * @param H This is sublayer ratio or how much one layer decreases from last.
 * @return SV_OK if surface exists and is extracted with Threshold correctly
 */
int cvTetGenMeshObject::SetBoundaryLayer(int type, int id, int side,
    int nL, double* H)
{
#ifdef SV_USE_VMTK
  meshoptions_.boundarylayermeshflag = 1;
  meshoptions_.numsublayers = nL;
  meshoptions_.blthicknessfactor = *H;
  meshoptions_.sublayerratio = *(H+1);
  meshoptions_.useconstantblthickness = *(H+2);
  meshoptions_.meshwallfirst = 1;
#else
  fprintf(stderr,"Plugin VMTK is not being used! \
      In order to use boundary layer meshing, \
      plugin VMTK must be available!\n");
  return SV_ERROR;
#endif

  return SV_OK;
}

// --------------------
//  SetWalls
// --------------------
/**
 * @brief Function to set the walls of the mesh with an integer array
 * @param numWalls number of walls being set with value
 * @param walls, integer list for ModelFaceIds in the wall.
 * @return SV_OK if surface exists and array is set correctly
 */
int cvTetGenMeshObject::SetWalls(int numWalls, int *walls)
{
  int max=0;
  double range[2];
  vtkIntArray *wallArray = vtkIntArray::New();
  meshoptions_.meshwallfirst = 1;

  if (VtkUtils_PDCheckArrayName(polydatasolid_,1,"ModelFaceID") != SV_OK)
  {
    fprintf(stderr,"ModelFaceID array not on object, so cannot set walls\n");
    return SV_ERROR;
  }
  vtkIntArray *modelIds;
  modelIds = vtkIntArray::SafeDownCast( polydatasolid_->GetCellData()->GetArray("ModelFaceID"));
  modelIds->GetRange(range);
  max = range[1];

  int *isWall = new int[max];
  for (int i=0; i < max; i++)
    isWall[i] = 0;
  for (int i=0; i < numWalls; i++)
  {
    int wallid = *(walls+i);
    isWall[wallid-1] = 1;
    wallFaceIDs_.insert(wallid);
  }

  int numCells = polydatasolid_->GetNumberOfCells();
  for (vtkIdType cellId = 0; cellId < numCells ; cellId++)
  {
    int value = modelIds->GetValue(cellId);
    if (isWall[value - 1] == 1)
      wallArray->InsertValue(cellId,1);
    else
      wallArray->InsertValue(cellId,0);
  }

  wallArray->SetName("WallID");
  polydatasolid_->GetCellData()->AddArray(wallArray);
  wallArray->Delete();

#ifdef SV_USE_MMG
  if (meshoptions_.usemmg == 0)
  {
#endif
    auto thresholder = vtkSmartPointer<vtkThreshold>::New();
    thresholder->SetInputData(polydatasolid_);
    thresholder->SetInputArrayToProcess(0,0,0,1,"WallID");
    thresholder->SetLowerThreshold(1.0);
    thresholder->SetUpperThreshold(1.0);
    thresholder->Update();

    auto surfacer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfacer->SetInputData(thresholder->GetOutput());
    surfacer->Update();

    polydatasolid_->DeepCopy(surfacer->GetOutput());
#ifdef SV_USE_MMG
  }
#endif

  delete [] isWall;
  return SV_OK;
}

// --------------------
//  SetCylinderRefinement
// --------------------
/**
 * @brief Function to set the region to refine based on input cylinder
 * @param size This is the smaller refined of the edges within cylinder region.
 * @param radius This is the radius of the refinement cylinder.
 * @param center This is the center of the refinement cylinder.
 * Halfway along the length.
 * @param normal This is the normal direction from the center that the length
 * of the cylinder follows.
 * @return SV_OK if the mesh sizing function based on the circle is computed
 * correctly
 */
int cvTetGenMeshObject::SetCylinderRefinement(double size, double radius,
    double length, double* center, double *normal)
{
  //Set meshoptions_ parameters based on input.
  int i;
  meshoptions_.refinement = 1;
  meshoptions_.refinedsize = size;
  meshoptions_.cylinderradius = radius;
  meshoptions_.cylinderlength = length;
  for (i=0;i<3;i++)
  {
    meshoptions_.cylindercenter[i] = center[i];
    meshoptions_.cylindernormal[i] = normal[i];
  }

  //Create a new mesh sizing function and call TGenUtils to compute function.
  //Store in the member data vtkDouble Array meshsizingfunction
  if (TGenUtils_SetRefinementCylinder(polydatasolid_,"MeshSizingFunction",
	size,radius,center,length,normal,meshoptions_.secondarrayfunction,
	meshoptions_.maxedgesize,"RefineID",meshoptions_.refinecount) != SV_OK)
  {
    return SV_ERROR;
  }

  meshoptions_.secondarrayfunction = 1;
  meshoptions_.refinecount += 1;
  return SV_OK;
}


// --------------------
//  SetSphereRefinement
// --------------------
/**
 * @brief Function to set the region to refine based on input sphere
 * @param size This is the smaller refined of the edges within sphere region.
 * @param radius This is the radius of the refinement sphere.
 * @param center This is the center of the refinement sphere.
 * @return SV_OK if the mesh sizing function based on the circle is computed
 * correctly
 */
int cvTetGenMeshObject::SetSphereRefinement(double size, double radius,
    double* center)
{
  //Set meshoptions_ parameters based on input.
  int i;
  meshoptions_.refinement = 1;
  meshoptions_.refinedsize = size;
  meshoptions_.sphereradius = radius;
  for (i=0;i<3;i++)
  {
    meshoptions_.spherecenter[i] = center[i];
  }

  //Create a new mesh sizing function and call TGenUtils to compute function.
  //Store in the member data vtkDouble Array meshsizingfunction
  if (TGenUtils_SetRefinementSphere(polydatasolid_,"MeshSizingFunction",
	size,radius,center,meshoptions_.secondarrayfunction,
	meshoptions_.maxedgesize,"RefineID",meshoptions_.refinecount) != SV_OK)
  {
    return SV_ERROR;
  }

  meshoptions_.secondarrayfunction = 1;
  meshoptions_.refinecount += 1;
  return SV_OK;
}

// --------------------
//  SetSizeFunctionBasedMesh
// --------------------
/**
 * @brief Function to set up the mesh edge size for VMTK and TetGen
 * based on an input size array. For radius-based, this will take
 * the array, normalize it based on the minimum value and apply the
 * metric to the surface
 * @param size This is the size to be specified on selected surface if
 * there is no value specified in the array.
 * @param sizefunctionname This is the name of the function name
 * on the surface
 * @return SV_OK if the mesh sizing function is computed correctly
 */
int cvTetGenMeshObject::SetSizeFunctionBasedMesh(double size,char *sizefunctionname)
{

  fprintf(stderr,"Setting size based function...\n");
  //Set meshoptions_ parameters based on input.
  int i;
  meshoptions_.functionbasedmeshing = 1;

  //Create a new mesh sizing function and call TGenUtils to compute function.
  //Store in the member data vtkDouble Array meshsizingfunction
  if (TGenUtils_SetSizeFunctionArray(polydatasolid_,"MeshSizingFunction",
	size,sizefunctionname,meshoptions_.secondarrayfunction) != SV_OK)
  {
    return SV_ERROR;
  }

  //originalpolydata_->DeepCopy(polydatasolid_);
  meshoptions_.secondarrayfunction = 1;
  return SV_OK;
}

/**
 * @brief Function to generate a mesh
 * @return *result: SV_ERROR if the mesh doesn't exist. New Mesh must be
 * called before a mesh can be generated
 * @note Function checks to see if any of the mesh options have been set.
 * It they have, the corresponding tetgenbehavior object values are set.
 */

int cvTetGenMeshObject::GenerateMesh() {

  // Check the options controlling the generation of a solid vessel wall mesh.
  //
  if (meshoptions_.wallmeshflag)
  {
    // The wall mesh is generated inside the surface remeshing / boundary
    // layer path, so without these flags GenerateMesh() would silently
    // finish without a wall mesh.
    if (!meshoptions_.surfacemeshflag || !meshoptions_.volumemeshflag)
    {
      fprintf(stderr,"Surface and volume meshing must be enabled to generate a vessel wall mesh\n");
      return SV_ERROR;
    }
    if (!meshoptions_.boundarylayermeshflag)
    {
      fprintf(stderr,"Boundary layer meshing must be enabled to generate a vessel wall mesh\n");
      return SV_ERROR;
    }
    if (meshoptions_.boundarylayerdirection != 1)
    {
      fprintf(stderr,"The boundary layer must extrude inward (BoundaryLayerDirection 1)\
 when generating a vessel wall mesh\n");
      return SV_ERROR;
    }
    if (meshoptions_.wallthickness <= 0.0)
    {
      fprintf(stderr,"A wall thickness greater than zero (option WallThickness) must be set\
 to generate a vessel wall mesh\n");
      return SV_ERROR;
    }
  }

  if (surfacemesh_ != nullptr)
  {
    surfacemesh_->Delete();
    surfacemesh_ = nullptr;
  }

  if (volumemesh_ != nullptr)
  {
    volumemesh_->Delete();
    volumemesh_ = nullptr;
  }

//All these complicated options exist if using VMTK. Should be stopped prior
//to this if trying to use VMTK options and don't have VMTK.
#ifdef SV_USE_VMTK
  //If doing surface remeshing!
  if (meshoptions_.surfacemeshflag)
  {
     if (GenerateSurfaceRemesh() != SV_OK)
       return SV_ERROR;

    //If we are doing a volumemesh based off the surface mesh!
    if (meshoptions_.volumemeshflag)
    {
      //If we are doing boundary layer mesh, it gets complicated!
      if (meshoptions_.boundarylayermeshflag)
      {
        if (GenerateBoundaryLayerMesh() != SV_OK)
          return SV_ERROR;

        if (GenerateAndMeshCaps() != SV_OK)
          return SV_ERROR;
      }

      if (meshoptions_.boundarylayermeshflag || meshoptions_.functionbasedmeshing
        || meshoptions_.refinement)
      {
        if (GenerateMeshSizingFunction() != SV_OK)
          return SV_ERROR;
      }
      if (NewMesh() != SV_OK)
        return SV_ERROR;
    }
    //In this case, we only are doing a surface remesh, and we are essentially
    //done
    else
    {
      surfacemesh_ = vtkPolyData::New();
      surfacemesh_->DeepCopy(polydatasolid_);
    }
  }
#endif
  //If we are doing sphere refinement and only a volumemesh, then we need
  //to re-set up the mesh based on sizing function for sphere refinement
 
  if (meshoptions_.volumemeshflag && !meshoptions_.surfacemeshflag)
  {
    int meshInfo[3];
    if (TGenUtils_CheckSurfaceMesh(polydatasolid_, meshInfo) != SV_OK)
    {
      fprintf(stderr,"Error checking surface\n");
      return SV_ERROR;
    }

    if (!(meshoptions_.numberofholes > 0 || meshoptions_.numberofregions > 0))
    {
      if (!meshoptions_.allowMultipleRegions && meshInfo[0] > 1)
      {
        fprintf(stderr,"There are too many regions here!\n");
        fprintf(stderr,"Terminating meshing!\n");
        return SV_ERROR;
      }
      if (meshInfo[1] > 0 && !meshoptions_.boundarylayermeshflag)
      {
        fprintf(stderr,"There are free edes on surface!\n");
        fprintf(stderr,"Terminating meshing!\n");
        return SV_ERROR;
      }
      if (meshInfo[2] > 0)
      {
        fprintf(stderr,"There are bad edes on surface!\n");
        fprintf(stderr,"Terminating meshing!\n");
        return SV_ERROR;
      }
    }
    if (NewMesh() != SV_OK)
      return SV_ERROR;
  }

  //Here we set all the mesh flags for TetGen!
  if (meshoptions_.volumemeshflag)
  {
    tetgenbehavior* tgb = new tetgenbehavior;
    //Default flags: plc and neihgborlist/adjtetlist output
    tgb->plc=1;
    tgb->neighout=2;

    //User defined options below
    if (meshoptions_.maxedgesize != 0)
    {
      double mES = meshoptions_.maxedgesize;
      //Volume of an element is approximately (a^3)/(6*sqrt(2))
      double maxvol = (mES*mES*mES)/(6*sqrt(2.));
      tgb->fixedvolume=1;
      tgb->maxvolume=maxvol;
    }
    if (meshoptions_.minratio != 0)
    {
      tgb->quality=1;
      tgb->minratio=meshoptions_.minratio;
    }
    if (meshoptions_.mindihedral != 0.0)
    {
      tgb->quality=1;
      tgb->mindihedral = meshoptions_.mindihedral;
    }

    if (meshoptions_.optlevel != 0)
    {
      tgb->optlevel=meshoptions_.optlevel;
    }
    if (meshoptions_.epsilon != 0)
    {
      tgb->epsilon=meshoptions_.epsilon;
    }
    if (meshoptions_.verbose)
    {
      tgb->verbose=1;
    }
    if (meshoptions_.docheck)
    {
      tgb->docheck=1;
    }
    if (meshoptions_.quiet)
    {
      tgb->quiet=1;
    }
    if (meshoptions_.nobisect)
    {
      tgb->nobisect=1;
    }
    if (meshoptions_.diagnose)
    {
      tgb->diagnose=1;
    }
    if (meshoptions_.boundarylayermeshflag)
    {
      tgb->quality = 3;
      tgb->metric = 1;
      tgb->mindihedral = 10.0;
      tgb->nobisect=1;
    }
    if (meshoptions_.refinement)
    {
      tgb->quality = 3;
      tgb->metric = 1;
      tgb->mindihedral = 10.0;
    }
    if (meshoptions_.functionbasedmeshing)
    {
      tgb->quality = 3;
      tgb->metric = 1;
      tgb->mindihedral = 10.0;
    }
    if (meshoptions_.numberofregions > 0)
    {
      tgb->fixedvolume = 0;
      tgb->varvolume = 1;
      tgb->regionattrib = 1;
    }
#if defined(TETGEN150) || defined(TETGEN151)
    if (meshoptions_.coarsenpercent != 0)
    {
      tgb->coarsen=1;
      tgb->coarsen_percent=meshoptions_.coarsenpercent;
    }
    if (meshoptions_.nomerge)
    {
      tgb->nomergefacet=1;
      tgb->nomergevertex=1;
    }
#endif
#ifdef TETGEN143
    if (meshoptions_.boundarylayermeshflag)
    {
      tgb->goodratio = 2.0;
    }
    else
    {
      tgb->goodratio = 4.0;
    }
    tgb->goodangle = 0.88;
    tgb->useshelles = 1;
#endif

    if (meshloaded_ != 1)
    {
      fprintf(stderr,"For some reason, mesh is not loaded! TetGen cannot\
	  be run.\n");
    }

    fprintf(stdout,"TetGen Meshing Started...\n");
    try
    {
//      std::freopen("mesh_stats.txt","w",stdout);
      tetrahedralize(tgb, inmesh_, outmesh_);
//      std::fclose(stdout);
    }
    catch (int r)
    {
      fprintf(stderr,"ERROR: TetGen quit and returned error code %d\n",r);
      return SV_ERROR;
    }
    fprintf(stdout,"TetGen Meshing Finished...\n");
  }

  else
  {
    fprintf(stdout,"Only surface\n");
  }

#ifdef SV_USE_VMTK
  // This is a post meshing step that needs to be done for boundary layer mesh.
  if (meshoptions_.boundarylayermeshflag)
  {
    if (AppendBoundaryLayerMesh() != SV_OK)
    {
      fprintf(stderr,"Problem appending the boundary layer mesh\n");
      return SV_ERROR;
    }
  }
#endif

  // must have created mesh
  if (meshoptions_.volumemeshflag && !meshoptions_.boundarylayermeshflag)
  {
    if (outmesh_ == nullptr) {
      return SV_ERROR;
    }
    if (surfacemesh_ != nullptr)
    {
      surfacemesh_->Delete();
    }
    if (volumemesh_ != nullptr)
    {
      volumemesh_->Delete();
    }

    surfacemesh_ = vtkPolyData::New();
    volumemesh_ = vtkUnstructuredGrid::New();

    if (TGenUtils_ConvertToVTK(outmesh_,volumemesh_,surfacemesh_,
	  &numBoundaryRegions_,1) != SV_OK)
      return SV_ERROR;
  }

  // Boundary layer meshing creates 'ModelFaceID' IDs for regions of local mesh
  // refinement (e.g. local sphere refinement). Set the 'ModelFaceID' IDs to
  // the IDs from the original input model. This is not done for boundary layer
  // meshing extruded outward.
  //
  // This is also not done when a wall mesh is generated: the wall surface lies
  // outside of the original model so its face IDs (e.g. the IDs on the wall
  // end caps) cannot be set from the original model surface.
  //
  if (meshoptions_.boundarylayermeshflag && (meshoptions_.boundarylayerdirection == 1)
    && !meshoptions_.wallmeshflag) {
    if (TGenUtils_ResetOriginalRegions(surfacemesh_ ,originalpolydata_, "ModelFaceID") != SV_OK) {
      std::cout << "Failed to reset original face IDs for boundary layer mesh." << std::endl;
    }
  }

  // Report the quality of the volume mesh. The element aspect ratios are
  // also stored in an 'AspectRatio' cell data array on the volume mesh.
  //
  if (meshoptions_.volumemeshflag && (volumemesh_ != nullptr)) {
    TGenUtils_ReportMeshQuality(volumemesh_);
  }

  return SV_OK;
}

//-----------
// WriteMesh
//-----------
// Write the volume mesh to a file.
//
// [TODO:DaveP] What is 'smsver' ?
//
int cvTetGenMeshObject::WriteMesh(char *filename, int smsver) 
{
  if (volumemesh_ == nullptr) {
    return SV_ERROR;
  }

  TGenUtils_WriteVTU(filename, volumemesh_);

  return SV_OK;
}

int cvTetGenMeshObject::WriteStats(char *filename) {
  // must have created mesh
  if (inmesh_ == nullptr) {
    return SV_ERROR;
  }
  return SV_OK;
}

/**
 * @brief Procedure gets one face based on faceid defined in PolyDataUtils
 * @param orgfaceid int which is the number of the face to extract
 * @return *result: cvPolyData containg the face vtkPolyData
 */
cvPolyData* cvTetGenMeshObject::GetFacePolyData (int orgfaceid) {

  // recall the node numbers start at 1 in the P_id,
  // but in vtkPolyData file they start at 0.
  vtkPolyData *face = vtkPolyData::New();
  cvPolyData *result = nullptr;

  if (TGenUtils_GetFacePolyData(orgfaceid,surfacemesh_,face) != SV_OK)
  {
    return SV_ERROR;
  }

  result = new cvPolyData(face);

  face->Delete();

  return result;

}

//------------------
// GetModelFaceInfo
//------------------
// Get the face information used by the mesh.
//
// Returns:
//   desc: The description of the information returned.
//   faceInfo: A list of strings containing face IDs.
//
int cvTetGenMeshObject::GetModelFaceInfo(std::map<std::string,std::vector<std::string>>& faceInfo)
{
  faceInfo.clear();

  // [TODO:DaveP] What else can the kernel be?
  //
  if (solidmodeling_kernel_ == SM_KT_POLYDATA) {
      if (VtkUtils_PDCheckArrayName(originalpolydata_,1,"ModelFaceID") != SV_OK) {
          fprintf(stderr,"ModelFaceID does not exist\n");
          return SV_ERROR;
      }

      // Get face IDs.
      int *faces;
      int numFaces = 0;
      if (PlyDtaUtils_GetFaceIds(originalpolydata_, &numFaces, &faces) != SV_OK) {
          fprintf(stderr,"Could not get face ids\n");
          return SV_ERROR;
      }

      // Store face IDs into the returned string.
      for (int i = 0; i < numFaces; i++) {
          faceInfo[ModelFaceInfo::ID].push_back(std::to_string(faces[i]));
      }

      delete [] faces;
  }

  return SV_OK;
}

int cvTetGenMeshObject::GetModelFaceIDs(std::vector<int>& faceIDs)
{
  faceIDs.clear();

  // [TODO:DaveP] What else can the kernel be?
  //
  if (solidmodeling_kernel_ == SM_KT_POLYDATA) {
      if (VtkUtils_PDCheckArrayName(originalpolydata_,1,"ModelFaceID") != SV_OK) {
          fprintf(stderr,"ModelFaceID does not exist\n");
          return SV_ERROR;
      }

      // Get face IDs.
      int *faces;
      int numFaces = 0;
      if (PlyDtaUtils_GetFaceIds(originalpolydata_, &numFaces, &faces) != SV_OK) {
          fprintf(stderr,"Could not get face ids\n");
          return SV_ERROR;
      }

      // Store face IDs into the returned string.
      for (int i = 0; i < numFaces; i++) {
          faceIDs.push_back(faces[i]);
      }

      delete [] faces;
  }

  return SV_OK;
}


/**
 * @brief Function to set the PolyData member object
 * @param *newPolyData Pointer to vtkPolyData object that you want to be
 * set as the class member data
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::SetVtkPolyDataObject(vtkPolyData *newPolyData)
{
  if (polydatasolid_ != nullptr)
  {
    polydatasolid_->Delete();
  }
  polydatasolid_ = vtkPolyData::New();
  polydatasolid_->DeepCopy(newPolyData);

  return SV_OK;
}

/**
 * @brief Function to set an input unstructured grid
 * @param *newPolyData Pointer to vtkPolyData object that you want to be
 * set as the class member data
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::SetInputUnstructuredGrid(vtkUnstructuredGrid *ug)
{
  if (inputug_ != nullptr)
  {
    inputug_->Delete();
  }
  inputug_ = vtkUnstructuredGrid::New();
  inputug_->DeepCopy(ug);

  return SV_OK;
}

/**
 * @brief Helper function to generate surface mesh
 * @note This is a helper function. It is called from GenerateMesh
 * and it calls the VMTK utils for generating surface meshes
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::GenerateSurfaceRemesh()
{
#ifdef SV_USE_VMTK
  int meshcapsonly = 0;
  int preserveedges;
  int trianglesplitfactor;
  int useSizingFunction = 0;
  double collapseanglethreshold;
  std::string markerListName;
  vtkSmartPointer<vtkDoubleArray> meshsizingfunction =
    vtkSmartPointer<vtkDoubleArray>::New();
  //If we are doing a boundary layer mesh, we do not want to retain edges
  //for our surface remeshing
  //Else, we would like to preserve the edges of the boundary layer mesh
  if (meshoptions_.meshwallfirst)
  {
    preserveedges = 0;
    trianglesplitfactor = 5.0;
    collapseanglethreshold = 0.2;
    markerListName = "CellEntityIds";
  }
  else
  {
    preserveedges = 1;
    trianglesplitfactor = NULL;
    collapseanglethreshold = NULL;
    markerListName = "ModelFaceID";
  }
  //If doing sphere refinement, we need to base surface mesh on mesh
  //sizing function
  if (meshoptions_.refinement || meshoptions_.functionbasedmeshing)
  {
    useSizingFunction = 1;
    if (VtkUtils_PDCheckArrayName(polydatasolid_,0,"MeshSizingFunction") != SV_OK)
    {
      fprintf(stderr,"Array name 'MeshSizingFunction' does not exist. \
	              Something may have gone wrong when setting up BL");
      return SV_ERROR;
    }
    fprintf(stdout,"Getting sizing function Surface\n");
    meshsizingfunction = vtkDoubleArray::SafeDownCast(polydatasolid_->\
	  GetPointData()->GetScalars("MeshSizingFunction"));

  }
  //If not doing sphere refinement or function based meshing,
  //we do not base surface mesh on sizing function
  else
  {
    useSizingFunction = 0;
    meshsizingfunction = nullptr;
  }

  // Report the triangle shape quality on both sides of the remeshing. A
  // sliver reaching the wall extrusion cannot be fixed by any wall thickness
  // (see TGenUtils_LimitThicknessToPreventFold), so it has to be traced to
  // either the model surface or the remeshing itself.
  TGenUtils_ReportSurfaceTriangleQuality(polydatasolid_, "before surface remesh");

#ifdef SV_USE_MMG
  if (meshoptions_.usemmg)
  {
    double meshFactor = 0.8;
    double meshsize = meshFactor*meshoptions_.maxedgesize;
    double mmg_maxsize = 1.5*meshsize;
    double mmg_minsize = 0.5*meshsize;
    if (meshoptions_.hausd == 0)
      meshoptions_.hausd = 10.0*meshsize;
    double hausd = meshoptions_.hausd;
    double dumAng = 45.0;
    double hgrad = 1.01;

    //Generate Surface Remeshing
    if(MMGUtils_SurfaceRemeshing(polydatasolid_, mmg_minsize,
	  mmg_maxsize, hausd, dumAng, hgrad,
	  useSizingFunction, meshsizingfunction, meshoptions_.refinecount) != SV_OK)
    {
      fprintf(stderr,"Problem with surface meshing\n");
      return SV_ERROR;
    }
  }
  else
  {
#endif
    //Generate Surface Remeshing
    if(VMTKUtils_SurfaceRemeshing(polydatasolid_,meshoptions_.maxedgesize,
          meshcapsonly,preserveedges,trianglesplitfactor,
          collapseanglethreshold,nullptr,markerListName,
          useSizingFunction,meshsizingfunction) != SV_OK)
    {
      fprintf(stderr,"Problem with surface meshing\n");
      return SV_ERROR;
    }
    ResetOriginalRegions("ModelFaceID");
#ifdef SV_USE_MMG
  }
#endif

  TGenUtils_ReportSurfaceTriangleQuality(polydatasolid_, "after surface remesh");

  int meshInfo[3];
  if (TGenUtils_CheckSurfaceMesh(polydatasolid_, meshInfo) != SV_OK)
  {
    fprintf(stderr,"Mesh surface is bad\n");
    return SV_ERROR;
  }

  if (!(meshoptions_.numberofholes > 0 ||
        meshoptions_.numberofregions > 0))
  {
    if (!meshoptions_.allowMultipleRegions && meshInfo[0] > 1)
    {
      fprintf(stderr,"There are too many regions here!\n");
      fprintf(stderr,"Terminating meshing!\n");
      return SV_ERROR;
    }
    if (meshInfo[1] > 0 && !meshoptions_.meshwallfirst)
    {
      fprintf(stderr,"There are free edes on surface!\n");
      fprintf(stderr,"Terminating meshing!\n");
      return SV_ERROR;
    }
    if (meshInfo[2] > 0)
    {
      fprintf(stderr,"There are bad edes on surface!\n");
      fprintf(stderr,"Terminating meshing!\n");
      return SV_ERROR;
    }
  }

  if (meshoptions_.meshwallfirst && !meshoptions_.boundarylayermeshflag)
  {
    if (GenerateAndMeshCaps() != SV_OK)
      return SV_ERROR;
    ResetOriginalRegions("ModelFaceID");
  }

#else
  fprintf(stderr,"Cannot do a surface remesh without using VMTK\n");
  return SV_ERROR;
#endif

  return SV_OK;
}

//-----------------------
// SetCapBoundaryNormals 
//-----------------------
// Set the normals for caps boundary points.
//
// The point normals for the mesh at the cap ends typically do not lie in the
// cap plane due to the averaging of polygon face normals. This causes the
// boundary layer mesh at the cap not to be flat.
//
// Cap boundary points normals are set for 'surface' which is used as the surface 
// for computing the boundary layer mesh. The normals are computed as the vector 
// from a boundary point to the cap boundary center.
//
// Data modified:
//   surface - 'Normals' data array. 
//
void cvTetGenMeshObject::SetCapBoundaryNormals(vtkPolyData* surface) 
{
  // Extract surface cap boundaries. 
  //
  auto feature_edges = vtkFeatureEdges::New();
  feature_edges->SetInputData(surface);
  feature_edges->BoundaryEdgesOn();
  feature_edges->ManifoldEdgesOff();
  feature_edges->NonManifoldEdgesOff();
  feature_edges->FeatureEdgesOff();
  feature_edges->Update();

  auto boundary_edges = feature_edges->GetOutput();
  auto clean_filter = vtkCleanPolyData::New();
  clean_filter->SetInputData(boundary_edges);
  clean_filter->Update();
  auto cleaned_edges = clean_filter->GetOutput();

  auto conn_filter = vtkPolyDataConnectivityFilter::New();
  conn_filter->SetInputData(cleaned_edges);
  conn_filter->SetExtractionModeToSpecifiedRegions();
  std::vector<vtkPolyData*> cap_boundaries;
  int edge_id = 0;

  while (true) {
    conn_filter->AddSpecifiedRegion(edge_id);
    conn_filter->Update();
    auto component = vtkPolyData::New();
    component->DeepCopy(conn_filter->GetOutput());
    if (component->GetNumberOfCells() <= 0) {
      break;
    }
    auto clean_filter = vtkCleanPolyData::New();
    clean_filter->SetInputData(component);
    clean_filter->Update();
    auto cleaned_edges = clean_filter->GetOutput();
    cap_boundaries.push_back(cleaned_edges);
    conn_filter->DeleteSpecifiedRegion(edge_id);
    edge_id += 1;
  }

  // Set cap boundary normals.
  //
  // For each point on cap boundaries find the corresponding
  // point on the surface, compute the normal there, and set it
  // for the surface.
  //
  auto surface_points = surface->GetPoints();
  auto surface_normals = surface->GetPointData()->GetArray("Normals");

  auto pointLocator = vtkSmartPointer<vtkPointLocator>::New();
  pointLocator->SetDataSet(surface);
  pointLocator->BuildLocator();

  for (auto& cap_boundary : cap_boundaries) { 
    auto com_filter = vtkSmartPointer<vtkCenterOfMass>::New();
    com_filter->SetInputData(cap_boundary);
    com_filter->Update();
    auto cap_center = com_filter->GetCenter();
    auto cap_points = cap_boundary->GetPoints();

    for (int i = 0; i < cap_boundary->GetNumberOfPoints(); i++) {
      auto cap_point = cap_boundary->GetPoint(i);
      // Find the cap point on the surface.
      int point_id = pointLocator->FindClosestPoint(cap_point);
      auto surface_point = surface_points->GetPoint(point_id);

      // Compute the normal.
      double normal[3];
      double mag = 0.0;
      for (int j = 0; j < 3; j++) {
        normal[j] = surface_point[j] - cap_center[j];
        mag += normal[j]*normal[j];
      }
      mag = sqrt(mag);
      for (int j = 0; j < 3; j++) {
        normal[j] /= mag; 
      }

      // Set the surface normal. 
      surface_normals->SetComponent(point_id, 0, normal[0]);
      surface_normals->SetComponent(point_id, 1, normal[1]);
      surface_normals->SetComponent(point_id, 2, normal[2]);
    }
  }
}

/**
 * @brief Helper function to generate boundary layer mesh
 * @note This is a helper function. It is called from GenerateMesh
 * and it calls the VMTK utils for generating a boundary layer mesh
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::GenerateBoundaryLayerMesh()
{
  std::cout << "[GenerateBoundaryLayerMesh] " << std::endl;
  std::cout << "[GenerateBoundaryLayerMesh] ========== cvTetGenMeshObject::GenerateBoundaryLayerMesh ==========" << std::endl;

#ifdef SV_USE_VMTK
  if (boundarylayermesh_ != nullptr)
  {
    boundarylayermesh_->Delete();
  }

  if (innerblmesh_ != nullptr)
  {
    innerblmesh_->Delete();
  }

  std::string markerListName;

  if (meshoptions_.boundarylayermeshflag)
  {
    markerListName = "CellEntityIds";
  }
  else
  {
    markerListName = "ModelFaceID";
  }
  std::cout << "[GenerateBoundaryLayerMesh] markerListName: " << markerListName << std::endl;

  // Clean and convert the polydata to a vtu.
  //
  auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
  normaler->SetInputData(polydatasolid_);
  normaler->SetConsistency(1);
  normaler->SetAutoOrientNormals(1);
  normaler->SetFlipNormals(0);
  normaler->SetComputeCellNormals(0);
  normaler->SplittingOff();
  normaler->Update();

  auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(normaler->GetOutput());
  cleaner->Update();

  auto originalsurfpd = vtkSmartPointer<vtkPolyData>::New();
  originalsurfpd->DeepCopy(cleaner->GetOutput());

  // This surface is shared: the fluid boundary layer, the TetGen fluid volume
  // and the extruded solid wall are all derived from it, so a sliver here is
  // in all three. It is reported after the cleaning because merging coincident
  // points can itself change the triangles.
  TGenUtils_ReportSurfaceTriangleQuality(originalsurfpd, "fluid/wall interface surface");

  // Create a sizing function vtk data array named 'MeshSizingFunction' for the 
  // current mesh used to define the size of the mesh at each node. 
  //
  // This modifies 'originalsurfpd'.
  //
  if (VMTKUtils_ComputeSizingFunction(originalsurfpd, NULL, "MeshSizingFunction") != SV_OK) {
    fprintf(stderr,"Problem when computing sizing function");
    return SV_ERROR;
  }

  // Set cap boundary normals.
  SetCapBoundaryNormals(originalsurfpd);

  // Convert the surface to an vtkUnstructuredGrid.
  auto converter = vtkSmartPointer<vtkvmtkPolyDataToUnstructuredGridFilter>::New();
  converter->SetInputData(originalsurfpd);
  converter->Update();

  // Copy the vtkUnstructuredGrid.
  innerblmesh_ = vtkUnstructuredGrid::New();
  innerblmesh_->DeepCopy(converter->GetOutput());
  boundarylayermesh_ = vtkUnstructuredGrid::New();
  boundarylayermesh_->DeepCopy(converter->GetOutput());


  // Compute the boundary layer mesh using vmtk. 
  //
  // Modifies: 'boundarylayermesh_' and 'innerSurface'.
  //
  int negateWarpVectors = meshoptions_.boundarylayerdirection;
  int innerSurfaceCellId = 1;
  int sidewallCellEntityId = 9999;
  int useConstantThickness = meshoptions_.useconstantblthickness;
  std::string layerThicknessArrayName = "MeshSizingFunction";
  auto innerSurface = vtkSmartPointer<vtkUnstructuredGrid>::New();
  std::cout << "[GenerateBoundaryLayerMesh] negateWarpVectors: " << negateWarpVectors << std::endl;
  std::cout << "[GenerateBoundaryLayerMesh] useConstantThickness: " << useConstantThickness << std::endl;

  if (!boundarylayermesh_->GetPointData()->GetArray("Normals")) {
    std::cout << "[GenerateBoundaryLayerMesh] **** boundarylayermesh_ does not have point normals." << std::endl;
  } else {
    std::cout << "[GenerateBoundaryLayerMesh] boundarylayermesh_ has point normals." << std::endl;
  }

  // Measure the fluid boundary layer against its own concave curvature, the
  // same way the solid wall is measured. The two are extruded from this one
  // surface in opposite directions, so a junction that is concave for the wall
  // is convex for the boundary layer and vice versa, and the ratio that
  // decides whether an extrusion can succeed (t/R) is therefore a different
  // number on each side. Reporting only the wall side leaves the claim that
  // the boundary layer is safe as an assumption; this measures it.
  {
    double blThickness = meshoptions_.maxedgesize*meshoptions_.blthicknessfactor;
    fprintf(stdout,"Fluid boundary layer options in effect: MaxEdgeSize %g, BLThicknessFactor %g, "
        "constant thickness %d (total layer thickness %g), NumSubLayers %d, SubLayerRatio %g\n",
        meshoptions_.maxedgesize, meshoptions_.blthicknessfactor, useConstantThickness,
        blThickness, meshoptions_.numsublayers, meshoptions_.sublayerratio);

    // vmtk builds the layer thickness per point as the sizing function scaled
    // by the thickness factor, or as one constant when constant thickness is
    // set; rebuild the same value here so the ratio is reported against the
    // thickness that is actually extruded.
    auto blThicknessArray = vtkSmartPointer<vtkDoubleArray>::New();
    blThicknessArray->SetName("FluidLayerThickness");
    blThicknessArray->SetNumberOfComponents(1);
    blThicknessArray->SetNumberOfTuples(originalsurfpd->GetNumberOfPoints());
    auto sizing = originalsurfpd->GetPointData()->GetArray(layerThicknessArrayName.c_str());
    for (vtkIdType ptId = 0; ptId < originalsurfpd->GetNumberOfPoints(); ptId++)
    {
      double thickness = blThickness;
      if (!useConstantThickness && sizing != nullptr)
      {
        thickness = sizing->GetComponent(ptId,0)*meshoptions_.blthicknessfactor;
      }
      blThicknessArray->SetValue(ptId, thickness);
    }

    if (TGenUtils_ReportConcaveCurvatureVsThickness(originalsurfpd, blThicknessArray, 1,
          "fluid boundary layer, inward") != SV_OK)
    {
      fprintf(stderr,"Problem reporting the fluid boundary layer concave curvature\n");
      return SV_ERROR;
    }

    originalsurfpd->GetPointData()->RemoveArray("ThicknessOverRadius");
    originalsurfpd->GetPointData()->RemoveArray("ConcaveRadiusTypical");
    originalsurfpd->GetPointData()->RemoveArray("ConcaveRadiusSmallest");
  }

  if (VMTKUtils_BoundaryLayerMesh(boundarylayermesh_, innerSurface, meshoptions_.maxedgesize, meshoptions_.blthicknessfactor,
	meshoptions_.numsublayers, meshoptions_.sublayerratio, sidewallCellEntityId, innerSurfaceCellId, negateWarpVectors,
	markerListName, useConstantThickness, layerThicknessArrayName) != SV_OK)
  {
    fprintf(stderr,"Problem with boundary layer meshing\n");
    return SV_ERROR;
  }

  // We take the inside surface of the boundary layer mesh and set the
  // member vtkPolyData polydatasolid_ to be equal to this. We will
  // use this to create the volume mesh with TetGen.
  //
  auto surfacer = vtkSmartPointer<vtkGeometryFilter>::New();
  surfacer->SetInputData(innerSurface);
  surfacer->Update();

  if (meshoptions_.boundarylayerdirection)
  {
    polydatasolid_->DeepCopy(surfacer->GetOutput());
  }
  else
  {
    // The remeshing excludes cell entity ids with value of 1. Need to
    // set all to one so that only caps are remeshed.
    //
    polydatasolid_->DeepCopy(originalsurfpd);
    polydatasolid_->GetCellData()->RemoveArray("CellEntityIds");
    auto newEntityIds = vtkSmartPointer<vtkIntArray>::New();
    newEntityIds->SetNumberOfTuples(polydatasolid_->GetNumberOfCells());
    newEntityIds->FillComponent(0, 1);
    newEntityIds->SetName("CellEntityIds");
    polydatasolid_->GetCellData()->AddArray(newEntityIds);
  }

  // Generate a solid vessel wall mesh by extruding the original model
  // surface outward by the wall thickness. The wall mesh is appended
  // to the fluid mesh as a separate region (see AppendBoundaryLayerMesh).
  //
  if (meshoptions_.wallmeshflag)
  {
    if (GenerateWallMesh(originalsurfpd, markerListName) != SV_OK)
    {
      fprintf(stderr,"Problem with vessel wall meshing\n");
      return SV_ERROR;
    }
  }

#else
  fprintf(stderr,"Cannot generate a boundary layer mesh without VMTK\n");
  return SV_ERROR;
#endif

  return SV_OK;
}

//------------------
// GenerateWallMesh
//------------------
// Generate a solid vessel wall mesh by extruding the given surface (the
// original model wall surface, no caps) outward using the vmtk boundary
// layer generator.
//
// The wall thickness at each surface node is taken from the global
// 'WallThickness' option; faces with a 'LocalWallThickness' option override
// the global value in the same way local edge sizes override the global
// edge size.
//
// The generated mesh is stored in the member data 'wallmesh_'.
//
// Arguments:
//   wallSurface - The model wall surface with outward oriented point normals.
//   markerListName - The name of the cell array used to mark the extruded
//     mesh cells (see GenerateBoundaryLayerMesh).
//
int cvTetGenMeshObject::GenerateWallMesh(vtkPolyData* wallSurface, std::string markerListName)
{
#ifdef SV_USE_VMTK
  if (wallmesh_ != nullptr)
  {
    wallmesh_->Delete();
    wallmesh_ = nullptr;
  }

  if (meshoptions_.wallthickness <= 0.0)
  {
    fprintf(stderr,"A wall thickness greater than zero must be set to generate a wall mesh\n");
    return SV_ERROR;
  }

  // Echo the wall mesh options actually in effect so the log unambiguously
  // shows which values (in particular the concave curvature factor) drove
  // this run, independent of how the GUI or command history set them.
  fprintf(stdout,"Wall mesh options in effect: WallThickness %g, CurvatureFactor %g, SmoothingIterations %d, NumberOfWallLayers %d\n",
      meshoptions_.wallthickness, meshoptions_.wallthicknesscurvaturefactor,
      meshoptions_.wallthicknesssmoothingiterations, meshoptions_.numwallsublayers);

  // Every thickness pass below except the gradation limit exists to keep a
  // one-to-one outward extrusion valid, and each one buys that validity by
  // taking thickness away - measured, the clamp, the rounding and the fold
  // limit between them left hundreds of points under half the wall asked for.
  // The shell fill does not extrude. It offsets the surface as a whole, where
  // a concave junction creases instead of folding, so there is nothing for
  // those passes to prevent and handing them the field first would only thin
  // it for a problem the fill does not have.
  const bool extrudeWedges = !meshoptions_.walltetgenshell;

  // Create a point data array giving the wall thickness at each node
  // of the surface.
  //
  auto surface = vtkSmartPointer<vtkPolyData>::New();
  surface->DeepCopy(wallSurface);

  int numPts = surface->GetNumberOfPoints();
  auto thicknessArray = vtkSmartPointer<vtkDoubleArray>::New();
  thicknessArray->SetNumberOfComponents(1);
  thicknessArray->SetNumberOfTuples(numPts);
  thicknessArray->SetName("WallThickness");

  // Build the requested per-node wall thickness: the global wall thickness,
  // or the local wall thickness of the face the node belongs to when one is
  // set for it. A point shared by faces with different thicknesses gets the
  // arithmetic mean of the unique face thicknesses, so the value at a face
  // boundary does not depend on how many triangles each face is divided into
  // or on the order the cells are visited.
  std::vector<double> requestedThickness(numPts, meshoptions_.wallthickness);
  if (!localWallThickness_.empty())
  {
    auto faceIds = vtkIntArray::SafeDownCast(surface->GetCellData()->GetArray("ModelFaceID"));
    if (faceIds == nullptr)
    {
      fprintf(stderr,"No 'ModelFaceID' array on the surface; cannot set local wall thickness\n");
      return SV_ERROR;
    }
    std::vector<std::set<int>> pointFaceIds(numPts);
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); cellId++)
    {
      int faceId = faceIds->GetValue(cellId);
      vtkIdType npts;
      const vtkIdType *pts;
      surface->GetCellPoints(cellId,npts,pts);
      for (vtkIdType j = 0; j < npts; j++)
      {
        pointFaceIds[pts[j]].insert(faceId);
      }
    }
    // Points not connected to any cell keep the global thickness.
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      if (pointFaceIds[ptId].empty())
      {
        continue;
      }
      double thicknessSum = 0.0;
      for (auto faceId : pointFaceIds[ptId])
      {
        auto localThickness = localWallThickness_.find(faceId);
        thicknessSum += (localThickness == localWallThickness_.end()) ?
          meshoptions_.wallthickness : localThickness->second;
      }
      requestedThickness[ptId] = thicknessSum / pointFaceIds[ptId].size();
    }
  }

  std::vector<double> baseThickness(requestedThickness);

  for (vtkIdType ptId = 0; ptId < numPts; ptId++)
  {
    thicknessArray->SetValue(ptId, baseThickness[ptId]);
  }

  // Smooth the extrusion warp vectors (the point normals) in concave regions
  // so the outward wall extrusion does not dip inward and twist where the
  // normals converge (such as the crotch where two vessels merge). Unlike the
  // thickness clamp and fold-prevention pass below, which only change the
  // thickness magnitude, this changes the extrusion direction, which is what
  // actually causes the inward dip and the twisted wall elements at junctions.
  // Only the normal direction changes; the wall thickness and the surface
  // points (the fluid/wall interface) never move, and the cap-boundary normals
  // set by SetCapBoundaryNormals are pinned. It runs before the clamp and the
  // fold-prevention pass so those act on the smoothed extrusion directions.
  const double warpVectorRelaxation = 0.5;
  if (TGenUtils_SmoothWarpVectorsInConcaveRegions(surface, "Normals",
        meshoptions_.wallthicknesssmoothingiterations, warpVectorRelaxation) != SV_OK)
  {
    fprintf(stderr,"Problem smoothing the wall extrusion warp vectors\n");
    return SV_ERROR;
  }

  // Diagnostic: report the requested thickness against the local concave
  // radius of curvature (t/R) before any pass reduces the thickness. A junction
  // with t/R > 1 cannot carry the requested thickness whatever the extrusion
  // does, so this separates a junction whose shape is too sharp for the
  // thickness from a thinning the passes below produce on a junction that could
  // have carried it. It runs after the warp vectors are smoothed because the
  // heights are measured along those normals, and it is independent of the
  // curvature clamp, which is skipped entirely when its factor is zero.
  // Nothing is modified.
  if (TGenUtils_ReportConcaveCurvatureVsThickness(surface, thicknessArray, 0, "solid wall, outward") != SV_OK)
  {
    fprintf(stderr,"Problem reporting the concave curvature against the wall thickness\n");
    return SV_ERROR;
  }

  // The largest change in wall thickness allowed per unit distance along the
  // surface. The outer surface then tilts at most atan(slope) away from the
  // inner one where the thickness varies, so 0.5 is a 26.6 degree taper - on
  // the same order as the extrusion tilt that the depression diagnostic
  // already treats as suspicious above 15 degrees, and well under the cliffs
  // the thickness passes produce on their own (the curvature clamp was
  // measured dropping a thickness by 2.3 times across a single edge).
  const double wallThicknessMaxSlope = 0.5;

  // Clamp the thickness values in concave regions (such as the crotch
  // where two vessels merge) where a thickness larger than the concave
  // radius of curvature would make the outward extruded outer wall fold
  // over and self-intersect. Only the thickness values change; the surface
  // points (the fluid/wall interface) never move.
  if (extrudeWedges &&
      TGenUtils_ClampThicknessToConcaveCurvature(surface, thicknessArray,
        meshoptions_.wallthicknesscurvaturefactor) != SV_OK)
  {
    fprintf(stderr,"Problem clamping the wall thickness array to the surface curvature\n");
    return SV_ERROR;
  }

  // Smooth the thickness values so the wall thickness transitions gradually
  // across face boundaries where neighboring faces have different local
  // thicknesses and around clamped concave regions; an abrupt step in the
  // thickness would otherwise show up as a step in the extruded outer wall
  // surface. Only the thickness values are smoothed; the surface points
  // (the fluid/wall interface) never move. Without local overrides, the
  // curvature clamp the thickness is uniform and smoothing is a no-op.
  if (!localWallThickness_.empty() || meshoptions_.wallthicknesscurvaturefactor > 0.0)
  {
    if (TGenUtils_SmoothPointArray(surface, thicknessArray,
          meshoptions_.wallthicknesssmoothingiterations) != SV_OK)
    {
      fprintf(stderr,"Problem smoothing the wall thickness array\n");
      return SV_ERROR;
    }

    // The smoothing pulls the clamped values back up toward their
    // unclamped neighbors, so the clamp is applied once more to restore
    // the curvature limit; the thickness stays smooth away from the
    // reclamped points.
    if (extrudeWedges &&
        TGenUtils_ClampThicknessToConcaveCurvature(surface, thicknessArray,
          meshoptions_.wallthicknesscurvaturefactor) != SV_OK)
    {
      fprintf(stderr,"Problem clamping the wall thickness array to the surface curvature\n");
      return SV_ERROR;
    }
  }

  // The reclamp above restores the curvature limit at the clamped points and
  // with it the step in the thickness the smoothing had just removed, because
  // averaging cannot both respect a ceiling and remove a step at the point the
  // ceiling applies to. Bound the gradient of the thickness field instead: a
  // point standing too far above a neighbor is brought down to within the
  // allowed slope. Because this only ever lowers, it keeps the curvature limit
  // and every local thickness as ceilings, and no further reclamp is needed.
  //
  // A step in the thickness is not cosmetic. The fold prevention pass below
  // levels a folded triangle to its smallest thickness, so a step feeds it
  // points to level and each levelling makes the next step; measured, that
  // turned 15 infeasible junction regions into 37 thinned ones. Bounding the
  // gradient removes the fuel rather than the fire.
  if (TGenUtils_LimitThicknessGradation(surface, thicknessArray,
        wallThicknessMaxSlope, "requested wall thickness") != SV_OK)
  {
    fprintf(stderr,"Problem limiting the wall thickness gradation\n");
    return SV_ERROR;
  }

  // Preserve the assigned thickness at concave junctions by rounding the outer
  // wall outward into a smooth convex fillet instead of letting the fold
  // prevention pass thin it there (which reaches less far outward and so caves
  // in as the junction "depression"). The naive outer surface (each point at
  // its thickness along the normal) self-intersects at a concave crotch even
  // though every point is at the full thickness, because the outward normals
  // converge; this moves those outer points outward and apart into a fillet
  // that keeps the thickness and does not self-intersect, the way the outer
  // side of a thick welded junction fills with material. The inner surface
  // (the fluid/wall interface) is fixed; only the outer surface moves, encoded
  // back into the normal (extrusion direction) and thickness (extrusion
  // magnitude) so the existing extrusion reproduces it. It runs after the
  // thickness is finalized and before the fold prevention pass, which stays as
  // the safety net for a degenerate input sliver that no rounding can carry.
  const double outerRoundingRelaxation = 0.5;
  const double maxFilletRatio = 3.0;
  if (extrudeWedges &&
      TGenUtils_RoundOuterWallToPreserveThickness(surface, thicknessArray,
        meshoptions_.wallthicknesssmoothingiterations, outerRoundingRelaxation, maxFilletRatio) != SV_OK)
  {
    fprintf(stderr,"Problem rounding the outer wall to preserve the junction thickness\n");
    return SV_ERROR;
  }

  // Final safety pass: the rounding above keeps the thickness at concave
  // junctions, but a degenerate input sliver cannot carry a wall in any
  // direction, so check the actual extruded outer geometry (point +
  // thickness*normal) and reduce the thickness wherever an outer triangle
  // would still invert, so the wall mesh does not self-intersect. This runs
  // unconditionally because it only reduces thickness where a fold would
  // otherwise occur. Only the thickness values change; the surface points
  // (the fluid/wall interface) never move.
  // 30 iterations of the 0.8 reduction factor reach 0.1% of the requested
  // thickness, enough to clear a fold at a sliver triangle, whose thickness
  // has to drop to the order of the sliver's altitude; each iteration only
  // walks the surface triangles once and costs nothing against the meshing.
  if (extrudeWedges &&
      TGenUtils_LimitThicknessToPreventFold(surface, thicknessArray, 30) != SV_OK)
  {
    fprintf(stderr,"Problem limiting the wall thickness array to prevent the outer wall folding over\n");
    return SV_ERROR;
  }

  // The fold prevention pass leaves the same kind of cliff behind: it levels a
  // folded triangle to its smallest thickness and steps that value down, so a
  // point pulled far down sits next to neighbors it never touched. Bound the
  // gradient once more on the final field. It can only lower the thickness, and
  // lowering it in a concave region reduces rather than creates folding, so
  // this cannot undo the pass above. What it does undo is the crater: the
  // levelling propagates one ring per iteration for up to thirty iterations,
  // and a bounded gradient makes the transition it leaves a taper instead of a
  // step that the next extrusion has to absorb.
  if (extrudeWedges &&
      TGenUtils_LimitThicknessGradation(surface, thicknessArray,
        wallThicknessMaxSlope, "final wall extrusion length") != SV_OK)
  {
    fprintf(stderr,"Problem limiting the final wall thickness gradation\n");
    return SV_ERROR;
  }

  // Diagnostic: report how far the thickness passes (curvature clamp,
  // smoothing, fold prevention) reduced the wall below the requested value.
  // A wall much thinner than requested is the thin part of a junction
  // depression, so this locates the thinning from the log alone.
  {
    int numBelow90 = 0, numBelow50 = 0, numBelow25 = 0;
    std::vector<std::pair<double,vtkIdType> > thinPts;
    for (vtkIdType ptId = 0; ptId < numPts; ptId++)
    {
      double requested = baseThickness[ptId];
      if (requested <= 0.0)
      {
        continue;
      }
      double ratio = thicknessArray->GetValue(ptId) / requested;
      if (ratio < 0.90) { numBelow90++; thinPts.push_back(std::make_pair(ratio,ptId)); }
      if (ratio < 0.50) { numBelow50++; }
      if (ratio < 0.25) { numBelow25++; }
    }
    fprintf(stdout,"Wall thickness reduction (final vs requested): points below 90%%/50%%/25%%: %d/%d/%d of %d\n",
        numBelow90, numBelow50, numBelow25, numPts);

    // Report the thinned points as spatially separated regions instead of a
    // single worst point. One severe local defect otherwise hides every other
    // junction, so the log cannot distinguish "one bad spot" from "every
    // junction is thinned" - which is the question that decides whether the
    // fix belongs in the surface or in the thickness passes. The worst point
    // has the smallest ratio, which is the order the clustering seeds in.
    // A radius of 2% of the model diagonal keeps separate junctions apart while
    // still absorbing all the points belonging to one depression.
    if (!thinPts.empty())
    {
      const int maxRegions = 8;
      const double radiusFraction = 0.02;
      std::vector<TGenUtilsPointRegion> regions;
      double radius = 0.0;
      int numOutside = 0;
      int numRegionsTotal = 0;
      if (TGenUtils_ClusterPointsIntoRegions(surface, thinPts, maxRegions, radiusFraction,
            regions, radius, numOutside, numRegionsTotal) != SV_OK)
      {
        fprintf(stderr,"Problem clustering the thinned wall thickness points into regions\n");
        return SV_ERROR;
      }

      fprintf(stdout,"  thinned regions: %d in total (separated by %.4g), worst %d shown:\n",
          numRegionsTotal, radius, (int)regions.size());
      for (size_t i = 0; i < regions.size(); i++)
      {
        vtkIdType seedId = regions[i].seedId;
        double seed[3];
        surface->GetPoint(seedId, seed);
        fprintf(stdout,"    [%d] ratio %.3f (final %.5g / requested %.5g) at (%.5g, %.5g, %.5g), %d points\n",
            (int)(i+1), thicknessArray->GetValue(seedId)/baseThickness[seedId],
            thicknessArray->GetValue(seedId), baseThickness[seedId],
            seed[0], seed[1], seed[2], regions[i].numPoints);
      }
      if (numOutside > 0)
      {
        fprintf(stdout,"    ... %d further thinned points in the remaining %d regions\n",
            numOutside, numRegionsTotal - (int)regions.size());
      }
    }
  }

  // The report above divides the surviving extrusion length by the requested
  // one, which is not the wall thickness: at a concave junction an outer point
  // moves toward the inner surface of the vessel it is merging with, so the
  // wall can be far thinner than the extrusion is long. Measure the wall that
  // is actually produced, as the clearance from each outer point to the whole
  // inner surface, so the two can be compared and the gap between them is
  // visible rather than assumed away.
  //
  // It measures the clearance from each extruded outer point, so it only says
  // anything where those points are the outer surface. The shell fill's outer
  // surface is offset as a whole and has no such point-for-point relation to
  // the inner one, so this would be measuring a surface that is never built.
  if (extrudeWedges &&
      TGenUtils_ReportAchievedWallThickness(surface, thicknessArray, baseThickness,
        "solid wall, outward") != SV_OK)
  {
    fprintf(stderr,"Problem reporting the achieved wall thickness\n");
    return SV_ERROR;
  }

  // The diagnostic fields (t/R, the concave radii and the achieved thickness
  // ratio) live on this surface, so writing it once gives the whole picture as
  // a field that can be viewed, which is the only way to answer whether every
  // junction is affected or only the few the log has room to list. The write
  // is a polydata because the surface is one; it happens before the arrays are
  // stripped for the extrusion.
  {
    char wallDiagnosticsFile[] = "wall_diagnostics.vtp";
    TGenUtils_WriteVTP(wallDiagnosticsFile, surface);
  }

  // The diagnostic fields have been written out, so drop them here: everything
  // left on this surface is carried through the extrusion into the wall mesh
  // and then into the appended volume mesh, and a report has no business
  // ending up in the mesh that is handed to the solver.
  surface->GetPointData()->RemoveArray("ThicknessOverRadius");
  surface->GetPointData()->RemoveArray("ConcaveRadiusTypical");
  surface->GetPointData()->RemoveArray("ConcaveRadiusSmallest");
  surface->GetPointData()->RemoveArray("AchievedThicknessRatio");

  // Fill the wall with unstructured tetrahedra instead of extruding wedges,
  // when asked. Both fills leave the inner surface points untouched, which is
  // the only part of the wall the solver constrains, so the choice is free.
  if (meshoptions_.walltetgenshell)
  {
    return FillWallMeshWithTetGen(surface, thicknessArray);
  }

  surface->GetPointData()->RemoveArray("WallThickness");
  surface->GetPointData()->AddArray(thicknessArray);

  // Convert the surface to a vtkUnstructuredGrid as required by the vmtk
  // boundary layer generator.
  auto converter = vtkSmartPointer<vtkvmtkPolyDataToUnstructuredGridFilter>::New();
  converter->SetInputData(surface);
  converter->Update();

  wallmesh_ = vtkUnstructuredGrid::New();
  wallmesh_->DeepCopy(converter->GetOutput());

  // Extrude the wall mesh outward (negateWarpVectors = 0) with the thickness
  // given by the 'WallThickness' point data array. A sublayer ratio of 1.0
  // gives element layers of uniform thickness through the wall.
  //
  int negateWarpVectors = 0;
  int innerSurfaceCellId = 1;
  int sidewallCellEntityId = 9999;
  int useConstantThickness = 0;
  double thicknessFactor = 1.0;
  double wallSublayerRatio = 1.0;
  auto outerSurface = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (VMTKUtils_BoundaryLayerMesh(wallmesh_, outerSurface, meshoptions_.wallthickness, thicknessFactor,
	meshoptions_.numwallsublayers, wallSublayerRatio, sidewallCellEntityId, innerSurfaceCellId,
	negateWarpVectors, markerListName, useConstantThickness, "WallThickness") != SV_OK)
  {
    fprintf(stderr,"Problem with wall mesh extrusion\n");
    return SV_ERROR;
  }

  return SV_OK;
#else
  fprintf(stderr,"Cannot generate a wall mesh without VMTK\n");
  return SV_ERROR;
#endif
}

//------------------------
// FillWallMeshWithTetGen
//------------------------
/**
 * @brief Fills the solid wall with tetrahedra generated between the inner
 * surface and the outer surface, instead of extruding wedge layers from it.
 * @note The wedge extrusion gives every outer node exactly one inner node. At
 * a concave junction the outer nodes converge on each other, and with the node
 * correspondence fixed the only way to keep the extrusion valid is to shorten
 * it, which is what thins the wall there. Filling the volume between the two
 * surfaces has no such correspondence, so the mesher is free to put whatever
 * nodes it needs where the two vessels merge.
 *
 * The inner surface is passed through unchanged. It is the fluid/wall
 * interface, and the solver matches its nodes one to one against the fluid
 * mesh with a tolerance of a few units in the last place, so nothing may move
 * them - not the fill, and not the mesher, which is run with boundary
 * splitting off so it cannot insert a node on the input facets either.
 *
 * @param surface The inner surface with its final 'Normals' and thickness.
 * @param thicknessArray The final extrusion length per point.
 * @return SV_OK if the wall is filled.
 */

int cvTetGenMeshObject::FillWallMeshWithTetGen(vtkPolyData* surface, vtkDoubleArray* thicknessArray)
{
  // The outer surface is the true offset of the inner one at the requested
  // thickness, not the inner one pushed along its normals. Dilating a solid by
  // t rounds its convex features and creases its concave ones, and the parts of
  // the naive offset that run past that crease are not on the boundary at all;
  // the points whose offset lands there have no outer point, which is why every
  // pass that insisted on giving them one had to pay in thickness.
  //
  // The ceiling on the distance field: a float per voxel and a byte of
  // bookkeeping, so 160 million is about 800 MB, paid once per mesh.
  //
  // It has to be this large because a vascular bounding box is mostly empty. A
  // model a few tens of units across and a few hundred long needs a few hundred
  // cells on its long axis before the spacing even reaches the wall thickness,
  // and almost all of those cells are air. The builder coarsens the spacing to
  // fit the budget and reports how many voxels are left per wall thickness; the
  // smooth part of the offset survives a coarse grid, but the crease at a
  // junction is rounded over about half a cell, so a count near one means the
  // grid is eating the very thing the offset is for. If that happens on a model
  // this budget cannot hold, the fix is to contour in slabs rather than to
  // spend more memory.
  const vtkIdType maxOffsetVoxels = 160000000;
  auto offsetOuter = vtkSmartPointer<vtkPolyData>::New();
  if (TGenUtils_BuildOffsetOuterSurface(surface, thicknessArray,
        meshoptions_.maxedgesize, maxOffsetVoxels, offsetOuter) != SV_OK)
  {
    fprintf(stderr,"Problem building the offset outer wall surface\n");
    return SV_ERROR;
  }

  // The smallest wall the offset was asked for. It bounds how far the remesh
  // below may move the surface, because that is the quantity the whole offset
  // exists to keep.
  double smallestThickness = 0.0, largestThickness = 0.0;
  for (vtkIdType ptId = 0; ptId < thicknessArray->GetNumberOfTuples(); ptId++)
  {
    double t = thicknessArray->GetValue(ptId);
    if (t <= 0.0)
    {
      continue;
    }
    if (smallestThickness == 0.0 || t < smallestThickness)
    {
      smallestThickness = t;
    }
    if (t > largestThickness)
    {
      largestThickness = t;
    }
  }
  if (smallestThickness <= 0.0)
  {
    fprintf(stderr,"Every wall thickness is zero or negative, so there is no wall to fill\n");
    return SV_ERROR;
  }

  // The contour is at grid resolution and full of the slivers marching cubes
  // leaves behind, which the volume mesher would inherit. Remesh it to the
  // mesh edge size first.
  //
  // Two things about this remesher decide whether the offset survives it.
  //
  // It reads 'ModelFaceID' and returns without building anything if the array
  // is absent, so the contour has to be given one or the remesh silently
  // produces nothing. One face for the whole surface is the right answer here:
  // the offset has no correspondence to the model's faces, and the wall
  // boundary is tagged from the shell afterwards.
  //
  // Ridges - the edges MMG will not move - come only from the boundaries
  // between different 'ModelFaceID' values; angle detection is switched off in
  // the wrapper, which is why the angle argument is named 'dumAng' at the call
  // site that uses it. With one face there are no ridges, so nothing pins the
  // crease at the junctions except the Hausdorff distance. The wrapper's usual
  // value is ten times the mesh size, which here would let MMG move the surface
  // by many times the wall it is carrying; a tenth of the smallest wall keeps
  // it. That is not a refinement constraint at these radii - a chord of the
  // vessel deviates that far only at edges far longer than hmax - so it costs
  // no elements.
#ifdef SV_USE_MMG
  {
    double meshFactor = 0.8;
    double meshsize = meshFactor*meshoptions_.maxedgesize;
    double mmg_maxsize = 1.5*meshsize;
    double mmg_minsize = 0.5*meshsize;
    double offsetHausd = 0.1*smallestThickness;
    double dumAng = 45.0;
    double hgrad = 1.01;
    int useSizingFunction = 0;
    auto meshsizingfunction = vtkSmartPointer<vtkDoubleArray>::New();

    auto offsetFaceIds = vtkSmartPointer<vtkIntArray>::New();
    offsetFaceIds->SetName("ModelFaceID");
    offsetFaceIds->SetNumberOfComponents(1);
    offsetFaceIds->SetNumberOfTuples(offsetOuter->GetNumberOfCells());
    offsetFaceIds->FillComponent(0, 1);
    offsetOuter->GetCellData()->AddArray(offsetFaceIds);

    fprintf(stdout,"  remeshing the offset surface to edge sizes %.5g..%.5g, holding it within %.5g of where the level set put it\n",
        mmg_minsize, mmg_maxsize, offsetHausd);

    TGenUtils_ReportSurfaceTriangleQuality(offsetOuter, "offset outer wall, before remesh");
    if (MMGUtils_SurfaceRemeshing(offsetOuter, mmg_minsize, mmg_maxsize,
          offsetHausd, dumAng, hgrad, useSizingFunction, meshsizingfunction, 0) != SV_OK)
    {
      fprintf(stderr,"Problem remeshing the offset outer wall surface\n");
      return SV_ERROR;
    }
    if (offsetOuter->GetNumberOfCells() == 0)
    {
      fprintf(stderr,"The remesh of the offset outer wall surface produced nothing\n");
      return SV_ERROR;
    }
    TGenUtils_ReportSurfaceTriangleQuality(offsetOuter, "offset outer wall, after remesh");
  }
#else
  fprintf(stdout,"  built without MMG, so the offset surface is filled at its grid resolution; expect many more wall elements than the mesh edge size asks for\n");
#endif

  // Trim the dome the capping left over each vessel end, and pair the rim it
  // leaves with the inner cap rim so the wall can be closed between them.
  std::vector<TGenUtilsCapRim> caps;
  if (TGenUtils_TrimOffsetSurfaceAtCaps(surface, offsetOuter, largestThickness, caps) != SV_OK)
  {
    fprintf(stderr,"Problem trimming the offset outer wall surface at the caps\n");
    return SV_ERROR;
  }

  // Measure the wall the offset actually makes, now that it has been remeshed
  // and trimmed and is the surface the fill will use. Both directions are
  // reported because they answer different questions: outward checks that the
  // construction kept the thickness through the grid and the remesh, inward is
  // the wall standing over the interface.
  if (TGenUtils_ReportOffsetWallThickness(surface, thicknessArray, offsetOuter,
        "solid wall, offset") != SV_OK)
  {
    fprintf(stderr,"Problem reporting the offset wall thickness\n");
    return SV_ERROR;
  }

  // The ratio is a field over the interface, and the log has room for eight
  // regions of it. Writing it out is the only way to see whether a shortfall is
  // one junction or all of them.
  {
    char offsetDiagnosticsFile[] = "wall_offset_diagnostics.vtp";
    TGenUtils_WriteVTP(offsetDiagnosticsFile, surface);
  }
  surface->GetPointData()->RemoveArray("OffsetThicknessRatio");

  auto shell = vtkSmartPointer<vtkPolyData>::New();
  int numDegenerate = 0;
  if (TGenUtils_BuildWallShellSurface(surface, offsetOuter, caps, shell, numDegenerate) != SV_OK)
  {
    fprintf(stderr,"Problem building the wall shell surface\n");
    return SV_ERROR;
  }

  fprintf(stdout,"Filling the wall with TetGen tetrahedra: shell surface has %lld points and %lld triangles closing %zu vessel ends\n",
      (long long)shell->GetNumberOfPoints(), (long long)shell->GetNumberOfCells(), caps.size());
  if (numDegenerate > 0)
  {
    fprintf(stdout,"  %d of the vessel end triangles have no area; they are kept because dropping one would leave a hole in the wall, but the mesher may refuse them\n",
        numDegenerate);
  }

  auto shellInMesh = new tetgenio;
  auto shellOutMesh = new tetgenio;

  if (TGenUtils_ConvertSurfaceToTetGen(shellInMesh, shell) != SV_OK)
  {
    fprintf(stderr,"Problem converting the wall shell surface to TetGen\n");
    delete shellInMesh;
    delete shellOutMesh;
    return SV_ERROR;
  }

  // A closed inner surface encloses the lumen as well as the wall, so the
  // lumen has to be marked as a hole or it would be filled with wall elements.
  // An inner surface left open at the caps has a rim pair at each end and is
  // closed by the annulus between them, so it encloses the wall alone and there
  // is nothing to exclude.
  if (caps.empty())
  {
    double holePoint[3];
    if (TGenUtils_FindLumenHolePoint(surface, holePoint) != SV_OK)
    {
      fprintf(stderr,"Problem finding a point inside the lumen to exclude it from the wall\n");
      delete shellInMesh;
      delete shellOutMesh;
      return SV_ERROR;
    }
    fprintf(stdout,"  the inner surface is closed, so the lumen is marked as a hole at (%.5g, %.5g, %.5g)\n",
        holePoint[0], holePoint[1], holePoint[2]);

    auto holeList = vtkSmartPointer<vtkPoints>::New();
    holeList->InsertNextPoint(holePoint);
    if (TGenUtils_AddHoles(shellInMesh, holeList) != SV_OK)
    {
      fprintf(stderr,"Problem marking the lumen as a hole in the wall shell\n");
      delete shellInMesh;
      delete shellOutMesh;
      return SV_ERROR;
    }
  }

  auto shellBehavior = new tetgenbehavior;
  shellBehavior->plc = 1;
  // The input facets carry the fluid/wall interface nodes, so no node on them
  // may be added or moved; without this the interface stops matching the fluid
  // mesh and the solver refuses the case.
  shellBehavior->nobisect = 1;
  shellBehavior->quality = 1;
  shellBehavior->minratio = 1.414;
  shellBehavior->mindihedral = 10.0;

  fprintf(stdout,"  TetGen wall fill started...\n");
  try
  {
    tetrahedralize(shellBehavior, shellInMesh, shellOutMesh);
  }
  catch (int r)
  {
    fprintf(stderr,"ERROR: TetGen quit with error code %d while filling the wall. The shell it was\
 given is the inner surface, the offset surface, and the annulus closing the two at each vessel end;\
 TetGen prints the coordinates of the intersection above, so look them up to see which of the three it\
 is in. It is not the junction problem the extrusion had: the offset surface is a level set and cannot\
 fold onto itself, and where two walls would run into each other it merges them rather than crossing\
 them. Check the annulus first. It is the only part not built by the offset, it is triangulated by\
 angle about the cap axis, and that ordering only exists for a rim that winds once about its own\
 centre; a rim that got past that check but is close to failing it shows up as the zero-area end\
 triangles counted above. Otherwise the inner surface is self-intersecting and the wall has inherited\
 it - see the interface triangle quality report. Note that two vessels closer together than twice the\
 wall thickness do not fail here at all, they fuse into one solid, which reads in the offset thickness\
 report as an interface region carrying far more wall than was asked for\n", r);
    delete shellBehavior;
    delete shellInMesh;
    delete shellOutMesh;
    return SV_ERROR;
  }
  fprintf(stdout,"  TetGen wall fill finished\n");

  if (wallmesh_ != nullptr)
  {
    wallmesh_->Delete();
  }
  wallmesh_ = vtkUnstructuredGrid::New();

  auto wallSurfaceMesh = vtkSmartPointer<vtkPolyData>::New();
  int totRegions = 0;
  if (TGenUtils_ConvertToVTK(shellOutMesh, wallmesh_, wallSurfaceMesh, &totRegions, 0) != SV_OK)
  {
    fprintf(stderr,"Problem converting the filled wall mesh from TetGen\n");
    delete shellBehavior;
    delete shellInMesh;
    delete shellOutMesh;
    return SV_ERROR;
  }

  fprintf(stdout,"  wall filled with %lld tetrahedra on %lld nodes\n",
      (long long)wallmesh_->GetNumberOfCells(), (long long)wallmesh_->GetNumberOfPoints());

  // The wedge extrusion hands downstream a mesh holding both the volume cells
  // and the surface cells that bound them, tagged with 'CellEntityIds' and
  // 'ModelFaceID'; VMTKUtils_CreateBoundaryLayerSurfaceAndCaps splits the two
  // on the cell type and then reads 'CellEntityIds' without checking that it
  // exists. A mesh of tetrahedra alone therefore has to be given the same
  // shape, or that split produces an empty surface and the read dereferences
  // null. The boundary triangles are exactly the shell surface triangles,
  // because the mesher was run with boundary splitting off and so preserved
  // the input facets, and its output keeps the input points at their input
  // indices, so the shell cells can be inserted unchanged.
  {
    vtkIdType numInner = surface->GetNumberOfPoints();
    vtkIdType numTets = wallmesh_->GetNumberOfCells();

    auto cellEntityIds = vtkSmartPointer<vtkIntArray>::New();
    cellEntityIds->SetName("CellEntityIds");
    cellEntityIds->SetNumberOfComponents(1);
    auto modelFaceIds = vtkSmartPointer<vtkIntArray>::New();
    modelFaceIds->SetName("ModelFaceID");
    modelFaceIds->SetNumberOfComponents(1);

    for (vtkIdType cellId = 0; cellId < numTets; cellId++)
    {
      cellEntityIds->InsertNextValue(0);
      modelFaceIds->InsertNextValue(0);
    }

    // A triangle whose points are all inner points is the fluid/wall
    // interface, all outer points is the free outer wall, and a mix is a side
    // wall closing the two at a cap. 'CellEntityIds' follows the wedge
    // extrusion, which is all the downstream split reads it for: it only tests
    // for the side wall value.
    const int innerSurfaceCellId = 1;
    const int outerSurfaceCellId = 2;
    const int sidewallCellEntityId = 9999;

    // 'ModelFaceID' is a different matter, because it is the face the solver
    // and the mesh-complete output name their boundaries by, and the model
    // already uses small integers for its own faces. Writing the entity id into
    // it would give the interface the model's face 1 and the outer wall its
    // face 2. So the interface keeps the face it came from - the wedge
    // extrusion carries the same array through for the same reason - and the
    // outer wall gets one new face past the end of the model's range. A side
    // wall is left at 9999, which the downstream pass replaces with the id of
    // the cap it closes against.
    auto surfaceFaceIds = vtkIntArray::SafeDownCast(surface->GetCellData()->GetArray("ModelFaceID"));
    int outerWallFaceId = outerSurfaceCellId;
    if (surfaceFaceIds != nullptr)
    {
      double faceIdRange[2];
      surfaceFaceIds->GetRange(faceIdRange, 0);
      outerWallFaceId = (int)faceIdRange[1] + 1;
    }
    else
    {
      fprintf(stdout,"  the wall surface carries no 'ModelFaceID', so the wall boundary is tagged by its role alone\n");
    }

    auto faceLocator = vtkSmartPointer<vtkCellLocator>::New();
    faceLocator->SetDataSet(surface);
    faceLocator->BuildLocator();
    auto faceCell = vtkSmartPointer<vtkGenericCell>::New();

    int numInnerCells = 0, numOuterCells = 0, numSideCells = 0;

    for (vtkIdType cellId = 0; cellId < shell->GetNumberOfCells(); cellId++)
    {
      vtkIdType npts;
      const vtkIdType *pts;
      shell->GetCellPoints(cellId, npts, pts);
      if (npts != 3)
      {
        continue;
      }

      int numInnerPts = 0;
      for (vtkIdType j = 0; j < npts; j++)
      {
        if (pts[j] < numInner)
        {
          numInnerPts++;
        }
      }

      int entityId = sidewallCellEntityId;
      int faceId = sidewallCellEntityId;
      if (numInnerPts == 3)
      {
        entityId = innerSurfaceCellId;
        numInnerCells++;

        // The face of the interface triangle is the face of the wall surface
        // under it. It is looked up by position rather than by cell index
        // because the shell only holds the triangles of that surface, and a
        // surface that held anything else would put the two out of step.
        faceId = outerWallFaceId;
        if (surfaceFaceIds != nullptr)
        {
          double centroid[3] = {0.0, 0.0, 0.0};
          for (vtkIdType j = 0; j < npts; j++)
          {
            double p[3];
            shell->GetPoint(pts[j], p);
            for (int k = 0; k < 3; k++)
            {
              centroid[k] += p[k]/3.0;
            }
          }
          double closest[3];
          vtkIdType closestCell = -1;
          int subId = 0;
          double distanceSquared = 0.0;
          faceLocator->FindClosestPoint(centroid, closest, faceCell, closestCell, subId, distanceSquared);
          if (closestCell >= 0)
          {
            faceId = surfaceFaceIds->GetValue(closestCell);
          }
        }
      }
      else if (numInnerPts == 0)
      {
        entityId = outerSurfaceCellId;
        faceId = outerWallFaceId;
        numOuterCells++;
      }
      else
      {
        numSideCells++;
      }

      wallmesh_->InsertNextCell(VTK_TRIANGLE, npts, pts);
      cellEntityIds->InsertNextValue(entityId);
      modelFaceIds->InsertNextValue(faceId);
    }

    wallmesh_->GetCellData()->AddArray(cellEntityIds);
    wallmesh_->GetCellData()->AddArray(modelFaceIds);

    fprintf(stdout,"  wall boundary tagged: %d interface, %d outer, %d side wall triangles; the outer wall is ModelFaceID %d\n",
        numInnerCells, numOuterCells, numSideCells, outerWallFaceId);
  }

  delete shellBehavior;
  delete shellInMesh;
  delete shellOutMesh;

  return SV_OK;
}

/**
 * @brief Helper function to cap the surface and remesh them
 * @note This is a helper function. It is called from GenerateMesh
 * and it calls the VMTK utils for generating caps and remeshing them
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::GenerateAndMeshCaps()
{
#ifdef SV_USE_VMTK
  vtkSmartPointer<vtkIdList> excluded =
    vtkSmartPointer<vtkIdList>::New();
  int marker;
  int meshcapsonly = 1;
  int preserveedges = 0;
  int captype = 0;
  int trioutput = 1;
  int cellOffset = 1;
  int useSizeFunction = 0;
  int trianglesplitfactor;
  double collapseanglethreshold;
  vtkSmartPointer<vtkDoubleArray> meshsizingfunction =
    vtkSmartPointer<vtkDoubleArray>::New();
  std::string markerListName = "ModelFaceID";
  if (meshoptions_.meshwallfirst)
  {
    collapseanglethreshold = 0.2;
    trianglesplitfactor = 5.0;
  }
  else
  {
    collapseanglethreshold = NULL;
    trianglesplitfactor = NULL;
  }
  if (meshoptions_.functionbasedmeshing || meshoptions_.refinement)
  {
    useSizeFunction = 1;
    if (VtkUtils_PDCheckArrayName(polydatasolid_,0,"MeshSizingFunction") != SV_OK)
    {
      fprintf(stderr,"Array name 'MeshSizingFunctionID' does not exist. \
	              Something may have gone wrong when setting up BL");
      return SV_ERROR;
    }
    fprintf(stderr,"Getting sizing function caps\n");
    meshsizingfunction = vtkDoubleArray::SafeDownCast(polydatasolid_->\
	  GetPointData()->GetScalars("MeshSizingFunction"));
    fprintf(stderr,"Got function caps\n");
  }

  vtkDataArray *currentMarkers = polydatasolid_->GetCellData()->GetArray(markerListName.c_str());
  if (currentMarkers == nullptr)
  {
    excluded->SetNumberOfIds(1);
    excluded->InsertId(0,1);
    cellOffset = 1;
  }
  else
  {
    cellOffset = -1;
    for (int i=0; i<polydatasolid_->GetNumberOfCells(); i++)
    {
      marker = currentMarkers->GetTuple1(i);
      excluded->InsertUniqueId(marker);
      if (marker > cellOffset)
      {
        cellOffset = marker;
      }
    }
  }

  //We cap the inner surface of the boundary layer mesh
  if (VMTKUtils_Capper(polydatasolid_,captype,trioutput,cellOffset,
	markerListName) != SV_OK)
  {
    fprintf(stderr,"Problem with capping\n");
    return SV_ERROR;
  }

  //Set caps equal to caps on original
  if (TGenUtils_ResetOriginalRegions(polydatasolid_, originalpolydata_,
        markerListName, excluded) != SV_OK)
  {
    fprintf(stderr,"Problem with name resetting on cap remeshing\n");
    return SV_ERROR;
  }

  //We use VMTK surface remeshing to remesh just the caps to be
  //the same mesh size as the inner surface of the BL mesh
  if (VMTKUtils_SurfaceRemeshing(polydatasolid_,
	meshoptions_.maxedgesize,meshcapsonly,preserveedges,
	trianglesplitfactor,collapseanglethreshold,excluded,
	markerListName,useSizeFunction,meshsizingfunction) != SV_OK)
  {
    fprintf(stderr,"Problem with cap remeshing\n");
    return SV_ERROR;
  }
#else
  fprintf(stderr,"Cannot generate and mesh caps without VMTK\n");
  return SV_ERROR;
#endif

  // Add wall id back on to the surface
  vtkSmartPointer<vtkIntArray> wallIds = vtkSmartPointer<vtkIntArray>::New();
  wallIds->SetNumberOfTuples(polydatasolid_->GetNumberOfCells());
  wallIds->SetName("WallID");
  wallIds->FillComponent(0, 0);

  currentMarkers = polydatasolid_->GetCellData()->GetArray(markerListName.c_str());
  for (int i=0; i<polydatasolid_->GetNumberOfCells(); i++)
  {
    marker = currentMarkers->GetTuple1(i);
    if (excluded->IsId(marker) != -1)
    {
      wallIds->SetTuple1(i, 1);
    }
  }
  polydatasolid_->GetCellData()->AddArray(wallIds);

  return SV_OK;
}

/**
 * @brief Helper function to generate a mesh sizing function
 * @note This is a helper function. It is called from GenerateMesh
 * and it calls the VMTK utils to generate a mesh sizing function
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::GenerateMeshSizingFunction()
{
#ifdef SV_USE_VMTK
  vtkSmartPointer<vtkCleanPolyData> cleaner =
    vtkSmartPointer<vtkCleanPolyData>::New();
  //Compute a mesh sizing function to send to TetGen and create a
  //volume mesh based on. Must do this otherwise when appending
  //volume mesh and BL mesh, they won't match up!
  if (VtkUtils_PDCheckArrayName(polydatasolid_,0,"MeshSizingFunction") == 1)
  {
    fprintf(stderr,"MeshSizingFunction Name exists. Delete!\n");
    polydatasolid_->GetPointData()->RemoveArray("MeshSizingFunction");
  }
  if (VMTKUtils_ComputeSizingFunction(polydatasolid_,NULL,
	"MeshSizingFunction") != SV_OK)
  {
    fprintf(stderr,"Problem when computing sizing function");
    return SV_ERROR;
  }

  //Clean the output and make it a vtkUnstructuredGrid
  cleaner->SetInputData(polydatasolid_);
  cleaner->Update();

  polydatasolid_->DeepCopy(cleaner->GetOutput());
#else
  fprintf(stderr,"Cannot apply mesh sizing function without using VMTK\n");
  return SV_ERROR;
#endif

  return SV_OK;
}

/**
 * @brief Helper function to append meshes together for a full bl mesh
 * @note This is a helper function. It is called from GenerateMesh
 * and it calls the VMTK utils for appending meshes
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::AppendBoundaryLayerMesh()
{
#ifdef SV_USE_VMTK
  if (surfacemesh_ != nullptr)
  {
    surfacemesh_->Delete();
  }
  if (volumemesh_ != nullptr)
  {
    volumemesh_->Delete();
  }
  if (boundarylayermesh_ == nullptr)
  {
    fprintf(stderr,"Cannot append mesh without a Boundary Layer Mesh\n");
    return SV_ERROR;
  }
  if (innerblmesh_ == nullptr)
  {
    fprintf(stderr,"Cannot append mesh without inner surface from boundary layer\n");
    return SV_ERROR;
  }
  surfacemesh_ = vtkPolyData::New();
  volumemesh_ = vtkUnstructuredGrid::New();

  auto surfacetomesh = vtkSmartPointer<vtkvmtkPolyDataToUnstructuredGridFilter>::New();
  auto newnormaler = vtkSmartPointer<vtkPolyDataNormals>::New();

  if (TGenUtils_ConvertToVTK(outmesh_,volumemesh_,surfacemesh_,
    &numBoundaryRegions_,0) != SV_OK)
  {
  return SV_ERROR;
  }

  surfacetomesh->SetInputData(polydatasolid_);
  surfacetomesh->Update();

  // We append the volume mesh from tetgen, the inner surface, and the
  // boundary layer mesh all together.
  //
  auto newVolumeMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto newSurfaceMesh = vtkSmartPointer<vtkPolyData>::New();
  int giveblnewregion = meshoptions_.newregionboundarylayer;

  // When a wall mesh is generated it becomes the new region and the
  // boundary layer remains part of the fluid region.
  vtkUnstructuredGrid* wallmesh = nullptr;
  if (meshoptions_.wallmeshflag)
  {
    if (wallmesh_ == nullptr)
    {
      fprintf(stderr,"Cannot append mesh without a wall mesh\n");
      return SV_ERROR;
    }
    wallmesh = wallmesh_;
    if (giveblnewregion)
    {
      fprintf(stdout,"Note: the NewRegionBoundaryLayer option is ignored when generating a wall mesh\n");
      giveblnewregion = 0;
    }
  }

  fprintf(stdout,"Appending Boundary Layer and Volume Mesh\n");
  if (VMTKUtils_AppendData(volumemesh_,boundarylayermesh_,
    surfacetomesh->GetOutput(), newVolumeMesh, newSurfaceMesh, giveblnewregion, wallmesh) != SV_OK)
  {
    return SV_ERROR;
  }
  volumemesh_->DeepCopy(newVolumeMesh);

  //We generate normals for region detection and save in surfacemesh_
  newnormaler->SetInputData(newSurfaceMesh);
  newnormaler->SplittingOff();
  newnormaler->ComputeCellNormalsOn();
  newnormaler->ComputePointNormalsOff();
  newnormaler->AutoOrientNormalsOn();
  newnormaler->Update();

  surfacemesh_->DeepCopy(newnormaler->GetOutput());
  polydatasolid_->DeepCopy(newnormaler->GetOutput());
#else
  fprintf(stderr,"Cannot append the BL mesh without using VMTK\n");
  return SV_ERROR;
#endif

  return SV_OK;
}

/**
 * @brief Helper function to reset the original region ids
 * @note This is a helper function. It is called from GenerateMesh
 * and it calls the VMTK utils in order to set the regions back
 * to the regions originally identified on the first surface
 * @return SV_OK if executed correctly
 */
int cvTetGenMeshObject::ResetOriginalRegions(std::string regionName)
{
  if (polydatasolid_ == nullptr || originalpolydata_ == nullptr)
  {
    fprintf(stderr,"Cannot reset original regions without orignal surface \
		    or without new surface\n");
    return SV_ERROR;
  }

  if (TGenUtils_ResetOriginalRegions(polydatasolid_,originalpolydata_, regionName)
      != SV_OK)
  {
    fprintf(stderr,"Error while resetting the original region values\n");
    return SV_ERROR;
  }


 return SV_OK;
}

// --------------------
//  Adapt
// --------------------
/**
 * @brief Function to Adapt Mesh based on input adaption features etc.
 * @return SV_OK if adaptions performs correctly
 */
int cvTetGenMeshObject::Adapt()
{
  cout<<"Starting Adaptive Mesh..."<<endl;
  tetgenbehavior* newtgb = new tetgenbehavior;

  newtgb->refine=1;
  newtgb->metric=1;
  newtgb->quality = 3;
  newtgb->mindihedral = 10.0;
  newtgb->minratio = 1.2;
  newtgb->neighout=2;
  newtgb->verbose=1;
  //newtgb->coarsen=1;
  //newtgb->coarsen_param=8;
  //newtgb->coarsenpercent=1;
#if USE_TETGEN143
  newtgb->goodratio = 4.0;
  newtgb->goodangle = 0.88;
  newtgb->useshelles = 1;
#endif

  try
  {
    tetrahedralize(newtgb, inmesh_, outmesh_);
  }
  catch (int r)
  {
    fprintf(stderr,"ERROR: TetGen quit and returned error code %d\n",r);
    return SV_ERROR;
  }

  cout<<"Done with Adaptive Mesh..."<<endl;

  if (outmesh_ == nullptr) {
    return SV_ERROR;
  }
  if (surfacemesh_ != nullptr)
    surfacemesh_->Delete();

  if (volumemesh_ != nullptr)
    volumemesh_->Delete();

  surfacemesh_ = vtkPolyData::New();
  volumemesh_ = vtkUnstructuredGrid::New();
  if (TGenUtils_ConvertToVTK(outmesh_,volumemesh_,surfacemesh_,
	&numBoundaryRegions_,1) != SV_OK)
    return SV_ERROR;

  if (TGenUtils_ResetOriginalRegions(surfacemesh_,originalpolydata_,
	"ModelFaceID")
      != SV_OK)
  {
    fprintf(stderr,"Error while resetting the original region values\n");
    return SV_ERROR;
  }
  return SV_OK;
}

// --------------------
//  SetMetricOnMesh
// --------------------
int cvTetGenMeshObject::SetMetricOnMesh(double *error_indicator,int lstep,double factor, double hmax, double hmin,int strategy)
{
  // cant overwrite mesh
  if (inmesh_ != nullptr) {
    delete inmesh_;
  }
  if (outmesh_ != nullptr)
  {
    delete outmesh_;
  }

  if (inputug_ == nullptr)
    return SV_ERROR;

  if (polydatasolid_ == nullptr)
    return SV_ERROR;

  //Create new tetgen mesh objects and set first number of output mesh to 0
  inmesh_ = new tetgenio;
  inmesh_->firstnumber = 0;
  outmesh_ = new tetgenio;
  outmesh_->firstnumber = 0;

  if (TGenUtils_ConvertVolumeToTetGen(inputug_,polydatasolid_,inmesh_) != SV_OK)
  {
    fprintf(stderr,"Conversion from volume to TetGen failed\n");
    return SV_ERROR;
  }

  return SV_OK;
}

// --------------------
//  GetAdaptedMesh
// --------------------
int cvTetGenMeshObject::GetAdaptedMesh(vtkUnstructuredGrid *ug, vtkPolyData *pd)
{
  if (outmesh_ == nullptr) {
    return SV_ERROR;
  }
  if (volumemesh_ == nullptr) {
    return SV_ERROR;
  }
  if (surfacemesh_ == nullptr) {
    return SV_ERROR;
  }

  ug->DeepCopy(volumemesh_);
  pd->DeepCopy(surfacemesh_);

  return SV_OK;
}
