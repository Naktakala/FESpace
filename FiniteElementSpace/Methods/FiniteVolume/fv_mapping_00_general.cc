#include "fv_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/chi_mesh_utils.h"

using namespace chi_math::finite_element;

//###################################################################
/**Constructs a Finite Volume mapping of a cell.*/
FiniteVolume::
  FiniteVolume(const chi_mesh::Cell& cell,
               const chi_mesh::MeshContinuum& grid,
               std::vector<NodeInfo>& node_list) :
               FiniteElementMapping(cell, grid)
{
  const size_t num_nodes = 1;

  SetNumNodesAndLocalRegister(num_nodes, node_list.size());

  //======================================== Add the nodes to node_list
  node_list.emplace_back(NodeType::INTERNAL,
                         IdentifyingInfo({{},{m_cell.global_id}}),
                         m_cell.centroid);

  //======================================== Compute mandatories
  auto cell_info = chi_mesh::ComputeCellVolumeAndFaceAreas(cell, grid);
  m_volume = cell_info.volume;
  m_face_areas = std::move(cell_info.face_areas);
}

size_t FiniteVolume::FaceNumNodes(const chi_mesh::Cell& cell, const size_t f) const
{
  return 0;
}

std::vector<chi_mesh::Vector3> FiniteVolume::
  CellNodeLocations(const chi_mesh::Cell& cell) const
{
  return {cell.centroid};
}

size_t FiniteVolume::MapFaceNodeToCellNode(const chi_mesh::Cell& cell,
                                           size_t face_index,
                                           size_t face_node_index) const
{
  return 0;
}


