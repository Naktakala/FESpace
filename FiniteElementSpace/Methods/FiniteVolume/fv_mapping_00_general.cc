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

size_t FiniteVolume::FaceNumNodes(const size_t face_index) const
{
  return 0;
}

size_t FiniteVolume::MapFaceNodeToCellNode(const size_t face_index,
                                           const size_t face_node_index) const
{
  return 0;
}


