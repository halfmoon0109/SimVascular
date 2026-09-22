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

#include "sv4gui_MeshLegacyIO.h"
#include "sv4gui_MitkMeshIO.h"

#include "sv_polydatasolid_utils.h"
#include "sv_vtk_utils.h"

#include <QDir>

#include <set>

#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkConnectivityFilter.h>
#include <vtkErrorCode.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkThreshold.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGrid.h>
#include <string>
#include <vector>

//-----------------------
// ComputeVolumeMeshMaps
//-----------------------
//
static void ComputeVolumeMeshMaps(vtkSmartPointer<vtkUnstructuredGrid> volumeMesh, std::map<int,int>& node_map, std::map<int,int>& elem_map)
{
  auto nodeIDs = vtkIntArray::SafeDownCast(volumeMesh->GetPointData()->GetArray("GlobalNodeID"));
  for (int i = 0; i < nodeIDs->GetNumberOfTuples(); i++) {
    auto nid = nodeIDs->GetValue(i);
    if (node_map.count(nid) != 0) { 
      std::cout << "[ComputeVolumeMeshMaps] Duplicate node ID " << nid << std::endl;
    }
    node_map.insert({nid, i});
  }

  auto elemIDs = vtkIntArray::SafeDownCast(volumeMesh->GetCellData()->GetScalars("GlobalElementID"));
  for (int i = 0; i < elemIDs->GetNumberOfTuples(); i++) {
    auto eid = elemIDs->GetValue(i);
    if (elem_map.count(eid) != 0) { 
      std::cout << "[ComputeVolumeMeshMaps] Duplicate element ID " << eid << std::endl;
    }
    elem_map.insert({eid, i});
  }
}

//---------------------
// ResetFaceSurfaceIds 
//---------------------
// Reset the node and element IDs for surface faces so they 
// correctly index into their parent volume mesh.
//
// A face cell whose element id is not in this domain's volume mesh (measured
// 2026-09-22 on the solid wall's free outer face: one cell carried element
// 888626, a fluid element) is given back its element by its nodes: the one
// cell of the volume mesh on all three of them. Faces this cannot mend keep
// the id they came with, and the count is reported, since a face file whose
// element ids do not index the domain's mesh is what the solver reads.
//
static void ResetFaceSurfaceIds(vtkPolyData* surface, const std::map<int,int>& node_map, const std::map<int,int>& elem_map,
    vtkUnstructuredGrid* volumeMesh, const std::string& faceName)
{
  // Reset surface node IDs.
  //
  int num_nodes = surface->GetNumberOfPoints();
  auto node_ids = vtkIntArray::SafeDownCast(surface->GetPointData()->GetArray("GlobalNodeID"));

  auto node_ids_data = vtkSmartPointer<vtkIntArray>::New();
  node_ids_data->SetNumberOfValues(num_nodes);
  node_ids_data->SetName("GlobalNodeID");
  std::set<int> node_id_set;

  for (int i = 0; i < num_nodes; i++) { 
    auto nid = node_ids->GetValue(i);
    if (node_id_set.count(nid) != 0) {
      std::cout << "[ResetFaceSurfaceIds] Duplicate node ID " << nid << std::endl;
    }
    node_id_set.insert(nid);
    try {
      int index = node_map.at(nid);
      node_ids_data->SetValue(i, index+1);
    } catch (...) {
      std::cout << "[ResetFaceSurfaceIds] Can't find node " << nid << std::endl;
      return;
    }
  }
  surface->GetPointData()->RemoveArray("GlobalNodeID");
  surface->GetPointData()->AddArray(node_ids_data);

  // Reset surface element IDs.
  //
  int num_elems = surface->GetNumberOfCells();
  auto elem_ids = vtkIntArray::SafeDownCast(surface->GetCellData()->GetScalars("GlobalElementID"));

  auto elem_ids_data = vtkSmartPointer<vtkIntArray>::New();
  elem_ids_data->SetNumberOfValues(num_elems);
  elem_ids_data->SetName("GlobalElementID");

  int num_missing = 0, num_recovered = 0, first_missing = 0;
  auto cell_point_ids = vtkSmartPointer<vtkIdList>::New();
  auto candidate_cells = vtkSmartPointer<vtkIdList>::New();
  auto candidate_points = vtkSmartPointer<vtkIdList>::New();

  for (int i = 0; i < num_elems; i++) { 
    auto eid = elem_ids->GetValue(i);
    auto found = elem_map.find(eid);
    if (found != elem_map.end()) {
      elem_ids_data->SetValue(i, found->second+1);
      continue;
    }
    if (num_missing == 0) {
      first_missing = eid;
    }
    num_missing++;
    elem_ids_data->SetValue(i, eid);

    // The volume cell on every node of this face cell.
    if (volumeMesh == nullptr) {
      continue;
    }
    surface->GetCellPoints(i, cell_point_ids);
    std::vector<int> volume_points;
    bool have_points = true;
    for (vtkIdType k = 0; k < cell_point_ids->GetNumberOfIds(); k++) {
      auto nid = node_ids->GetValue(cell_point_ids->GetId(k));
      auto node_found = node_map.find(nid);
      if (node_found == node_map.end()) {
        have_points = false;
        break;
      }
      volume_points.push_back(node_found->second);
    }
    if (!have_points || volume_points.empty()) {
      continue;
    }
    volumeMesh->GetPointCells(volume_points[0], candidate_cells);
    for (vtkIdType c = 0; c < candidate_cells->GetNumberOfIds(); c++) {
      auto cell_id = candidate_cells->GetId(c);
      volumeMesh->GetCellPoints(cell_id, candidate_points);
      bool has_all = true;
      for (size_t k = 1; k < volume_points.size() && has_all; k++) {
        has_all = candidate_points->IsId(volume_points[k]) >= 0;
      }
      if (has_all) {
        elem_ids_data->SetValue(i, (int)cell_id+1);
        num_recovered++;
        break;
      }
    }
  }
  if (num_missing > 0) {
    std::cout << "[ResetFaceSurfaceIds] face '" << faceName << "': " << num_missing << " of " << num_elems
        << " cells carry an element id that is not in this domain's volume mesh (the first is " << first_missing
        << "); " << num_recovered << " were given their element back by their nodes"
        << (num_missing > num_recovered ? ", the rest keep the id they came with" : "") << std::endl;
  }
  surface->GetCellData()->RemoveArray("GlobalElementID");
  surface->GetCellData()->AddArray(elem_ids_data);
}

//------------
// WriteFiles
//------------
// Write mesh complete files.
//
bool sv4guiMeshLegacyIO::WriteFiles(mitk::DataNode::Pointer meshNode, sv4guiModelElement* modelElement, QString meshDir)
{
    if (meshNode.IsNull()) {
        return false;
    }

    auto mitkMesh = dynamic_cast<sv4guiMitkMesh*>(meshNode->GetData());
    if(!mitkMesh) {
        return false;
    }

    auto mesh = mitkMesh->GetMesh();
    if(!mesh) {
        return false;
    }

    std::string path="";
    meshNode->GetStringProperty("path",path);

    // Get the surface mesh from a data node or from a file.
    std::string surfaceFileName = path+"/"+meshNode->GetName()+".vtp";
    auto surfaceMesh = mesh->GetSurfaceMesh();
    if (surfaceMesh == nullptr && path != "") {
        surfaceMesh = mesh->CreateSurfaceMeshFromFile(surfaceFileName);
    }

    // Get the volume mesh from a data node or from a file.
    std::string volumeFileName = path+"/"+meshNode->GetName()+".vtu";
    auto volumeMesh = mesh->GetVolumeMesh();
    if (volumeMesh == nullptr && path != "") {
        volumeMesh = mesh->CreateVolumeMeshFromFile(volumeFileName);
    }

    // Check to see if multi domain
    double minmax[2] = {0.0, 0.0}; 
    if (surfaceMesh->GetCellData()->GetArray("ModelRegionID")) {
      surfaceMesh->GetCellData()->GetArray("ModelRegionID")->GetRange(minmax);
    }

    // Write a single-domain mesh. 
    //
    if (minmax[0] == minmax[1]) {
      QDir().mkpath(meshDir);
      return WriteFiles(surfaceMesh, volumeMesh, modelElement, meshDir);
    }

    // Extract meshes based on ModelRegionID.
    //
    for (int i = (int)minmax[0]; i <= (int) minmax[1]; i++) {

      // Extract surface mesh.
      //
      auto threshold_surface = VtkUtils_ThresholdSurface(i, i, "ModelRegionID", surfaceMesh);
      if (threshold_surface->GetNumberOfCells() == 0) { 
        continue;
      }

      // Extract volume mesh.
      auto threshold_volume = VtkUtils_ThresholdUgrid(i, i, "ModelRegionID", volumeMesh);
      if (threshold_volume->GetNumberOfCells() == 0) {
        continue;
      }

      // Write meshes.
      //
      fprintf(stdout, "[sv4guiMeshLegacyIO::WriteFiles] Writing domain %d\n", i);
      QString newDir = meshDir+"_domain-" + QString::number(i);
      QDir().mkpath(newDir);
      if (WriteFiles(threshold_surface, threshold_volume, modelElement, newDir) == false) {
        return false;
      }
    }

    return true;
}

//------------
// WriteFiles
//------------
// Write mesh complete files for the given domain (region).
//
// These mesh files can be used for FSI simulations. 
//
bool sv4guiMeshLegacyIO::WriteFiles(vtkSmartPointer<vtkPolyData> surfaceMesh, vtkSmartPointer<vtkUnstructuredGrid> volumeMesh, 
    sv4guiModelElement* modelElement, QString meshDir)
{
    if (!surfaceMesh || !volumeMesh || !modelElement) {
        return false;
    }

    QString vtuFilePath = meshDir + "/mesh-complete.mesh.vtu";
    vtuFilePath = QDir::toNativeSeparators(vtuFilePath);
    auto vtuWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    vtuWriter->SetCompressorTypeToZLib();
    vtuWriter->EncodeAppendedDataOff();
    vtuWriter->SetInputData(volumeMesh);
    vtuWriter->SetFileName(vtuFilePath.toStdString().c_str());
    vtuWriter->Write();

    // Get mappings between node and element IDs to an 
    // index into the data array location.
    std::map<int,int> node_map;
    std::map<int,int> elem_map;
    ComputeVolumeMeshMaps(volumeMesh, node_map, elem_map);
    volumeMesh->BuildLinks();

    QString vtpFilePath = meshDir + "/mesh-complete.exterior.vtp";
    vtpFilePath = QDir::toNativeSeparators(vtpFilePath);
    auto vtpWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    vtpWriter->SetCompressorTypeToZLib();
    vtpWriter->EncodeAppendedDataOff();
    vtpWriter->SetInputData(surfaceMesh);
    vtpWriter->SetFileName(vtpFilePath.toStdString().c_str());
    vtpWriter->Write();

    // Write face meshes used for bcs to the mesh-surfaces directory.
    //
    bool wallFound = false;
    auto wallAppender = vtkSmartPointer<vtkAppendPolyData>::New();
    wallAppender->UserManagedInputsOff();

    QDir mDir(meshDir);
    mDir.mkdir("mesh-surfaces");
    auto faces = modelElement->GetFaces();
    std::set<int> modelFaceIdents;
    std::set<std::string> writtenNames;

    for (int i = 0; i < faces.size(); i++) {
      auto face = faces[i];
      if (face == nullptr) {
        continue;
      }

      auto facepd = vtkSmartPointer<vtkPolyData>::New();
      int ident = modelElement->GetFaceIdentifierFromInnerSolid(face->id);
      modelFaceIdents.insert(ident);
      writtenNames.insert(face->name);
      PlyDtaUtils_GetFacePolyData(surfaceMesh.GetPointer(), &ident, facepd);

      ResetFaceSurfaceIds(facepd, node_map, elem_map, volumeMesh.GetPointer(), face->name);

      vtpFilePath = meshDir + "/mesh-surfaces/" + QString::fromStdString(face->name) + ".vtp";
      vtpFilePath = QDir::toNativeSeparators(vtpFilePath);
      vtpWriter->SetInputData(facepd);
      vtpWriter->SetFileName(vtpFilePath.toStdString().c_str());
      vtpWriter->Write();

      if (face->type == "wall") {
        wallAppender->AddInputData(facepd);
        wallFound = true;
      }
    }

    // A face id on the mesh that is no face of the model is one the mesher
    // made: the solid wall mesh tags its free outer surface with an id past
    // the model's, because the model has no face there. It is written under
    // its own name so the solid domain has a boundary file for it too. The
    // name is wall_outer, with the id appended when there are several or
    // when a model face already took the name, so nothing written above is
    // overwritten.
    //
    {
      std::set<int> extraFaceIdents;
      auto faceIds = vtkIntArray::SafeDownCast(surfaceMesh->GetCellData()->GetArray("ModelFaceID"));
      if (faceIds != nullptr) {
        for (vtkIdType cellId = 0; cellId < faceIds->GetNumberOfTuples(); cellId++) {
          int ident = faceIds->GetValue(cellId);
          if (modelFaceIdents.count(ident) == 0) {
            extraFaceIdents.insert(ident);
          }
        }
      }

      for (auto ident : extraFaceIdents) {
        auto facepd = vtkSmartPointer<vtkPolyData>::New();
        PlyDtaUtils_GetFacePolyData(surfaceMesh.GetPointer(), &ident, facepd);
        ResetFaceSurfaceIds(facepd, node_map, elem_map, volumeMesh.GetPointer(), std::string("face id ") + std::to_string(ident));

        QString name = (extraFaceIdents.size() == 1) ? QString("wall_outer") :
            QString("wall_outer_") + QString::number(ident);
        if (writtenNames.count(name.toStdString()) > 0) {
          name = QString("wall_outer_") + QString::number(ident);
          while (writtenNames.count(name.toStdString()) > 0) {
            name += "_";
          }
        }
        writtenNames.insert(name.toStdString());
        vtpFilePath = meshDir + "/mesh-surfaces/" + name + ".vtp";
        vtpFilePath = QDir::toNativeSeparators(vtpFilePath);
        vtpWriter->SetInputData(facepd);
        vtpWriter->SetFileName(vtpFilePath.toStdString().c_str());
        vtpWriter->Write();
        fprintf(stdout, "[sv4guiMeshLegacyIO::WriteFiles] Face id %d is not a face of the model; written as mesh-surfaces/%s.vtp\n",
            ident, name.toStdString().c_str());
      }
    }

    // If there are wall faces then extract separate faces from the
    // mesh surface based on RegionId. 
    //
    if (wallFound) {
      wallAppender->Update();
      auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
      cleaner->PointMergingOn();
      cleaner->PieceInvariantOff();
      cleaner->SetInputData(wallAppender->GetOutput());
      cleaner->Update();

      // Determine if wall surfaces are connected.
      auto connectFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
      connectFilter->SetInputData(cleaner->GetOutput());
      connectFilter->SetExtractionModeToAllRegions();
      connectFilter->ColorRegionsOn();
      connectFilter->Update();

      // If wall surfaces are not connected then for each regiond extract 
      // and write them to separate files.
      //
      // Note that these surfaces don't need to have their IDs remapped.
      //
      if (connectFilter->GetNumberOfExtractedRegions() > 1) {

        for (int j = 0; j < connectFilter->GetNumberOfExtractedRegions(); j++) {
          auto region_surface = VtkUtils_ThresholdSurface(j, j, "RegionId", connectFilter->GetOutput());
          region_surface->GetCellData()->RemoveArray("RegionId");
          region_surface->GetPointData()->RemoveArray("RegionId");
          /* dp
          auto thresholder = vtkSmartPointer<vtkThreshold>::New();
          thresholder->SetInputData(connectFilter->GetOutput());
          thresholder->SetInputArrayToProcess(0,0,0,1,"RegionId");
          thresholder->ThresholdBetween(j, j);
          thresholder->Update();

          auto surfacer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
          surfacer->SetInputData(thresholder->GetOutput());
          surfacer->Update();
          surfacer->GetOutput()->GetCellData()->RemoveArray("RegionId");
          surfacer->GetOutput()->GetPointData()->RemoveArray("RegionId");
          auto region_surface = surfacer->GetOutput(); 
          */

          vtpFilePath = meshDir + "/walls_combined_connected_region_" + QString::number(j) + ".vtp";
          vtpFilePath = QDir::toNativeSeparators(vtpFilePath);
          vtpWriter->SetInputData(region_surface);
          vtpWriter->SetFileName(vtpFilePath.toStdString().c_str());
          vtpWriter->Write();
        }

      // Write walls to a single file.
      //
      } else {
        // The faces appended above already had their node and element ids
        // reset to this domain's indices, so the combined surface is not
        // reset again: a second pass would read those indices as global ids
        // (measured 2026-09-22 on the solid wall domain, whose element ids
        // do not start at one: all 356293 cells of walls_combined came back
        // "not in this domain", and eleven were even handed a wrong element
        // by their likewise doubly mapped nodes). The fluid domain, numbered
        // from one, never showed it because there the two are the same.
        auto cleaned_surface = vtkSmartPointer<vtkPolyData>::New();
        cleaned_surface = cleaner->GetOutput(); 

        vtpFilePath = meshDir + "/walls_combined.vtp";
        vtpFilePath = QDir::toNativeSeparators(vtpFilePath);
        vtpWriter->SetInputData(cleaned_surface);
        vtpWriter->SetFileName(vtpFilePath.toStdString().c_str());
        vtpWriter->Write();
      }
    }

    return true;
}


