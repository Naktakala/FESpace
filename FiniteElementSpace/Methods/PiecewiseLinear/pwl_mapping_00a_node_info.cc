#include "pwl_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

using namespace chi_math::finite_element;

//std::vector<NodeInfo> PiecewiseLinear::GetNodeInfo() const
//{
//  std::vector<NodeInfo> cell_node_info;
//  cell_node_info.reserve(NumNodes());
//
//  for (uint64_t cvid : m_cell.vertex_ids)
//    cell_node_info.emplace_back(NodeType::CORNER,
//                                IdentifyingInfo({{cvid}}),
//                                m_grid.vertices[cvid]);
//
//  return cell_node_info;
//}