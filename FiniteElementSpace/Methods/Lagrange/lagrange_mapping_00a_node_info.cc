#include "lagrange_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

using namespace chi_math::finite_element;

//std::vector<NodeInfo> LagrangeQ2::GetNodeInfo() const
//{
//  std::vector<NodeInfo> cell_node_info;
//  cell_node_info.reserve(NumNodes());
//
//  //================================================== Add cell vertices
//  // All of the cell vertices participate in Q2 vertices
//  // so we can add all them here abstractly.
//  for (uint64_t cvid : m_cell.vertex_ids)
//    cell_node_info.emplace_back(NodeType::CORNER,
//                                IdentifyingInfo({{cvid}}),
//                                m_grid.vertices[cvid]);
//
//  //================================================== Cell-subtype specifics
//  // For a slab we only have to add the cell-centroid.
//  if (m_cell.SubType() == chi_mesh::CellType::SLAB)
//    cell_node_info.emplace_back(NodeType::INTERNAL,
//                                IdentifyingInfo({{},{m_cell.global_id}}),
//                                m_cell.centroid);
//  // For a triangle we need to add side centroids.
//  else if (m_cell.SubType() == chi_mesh::CellType::TRIANGLE)
//  {
//    for (auto& face : m_cell.faces)
//      cell_node_info.emplace_back(NodeType::FACE,
//                                  IdentifyingInfo(
//                                  {std::set<uint64_t>(face.vertex_ids.begin(),
//                                                      face.vertex_ids.end())}),
//                                  face.centroid);
//  }
//  // For a quadrilateral we add the edge centers and cell centroid.
//  else if (m_cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
//  {
//    for (auto& face : m_cell.faces)
//      cell_node_info.emplace_back(NodeType::FACE,
//                                  IdentifyingInfo(
//                                  {std::set<uint64_t>(face.vertex_ids.begin(),
//                                                      face.vertex_ids.end())}),
//                                  face.centroid);
//
//    cell_node_info.emplace_back(NodeType::INTERNAL,
//                                IdentifyingInfo({{},{m_cell.global_id}}),
//                                m_cell.centroid);
//  }
//  else if (m_cell.SubType() == chi_mesh::CellType::TETRAHEDRON or
//           m_cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
//  {
//    // First we build the edges
//    std::set<std::pair<uint64_t, uint64_t>> edges;
//
//    for (auto& face : m_cell.faces)
//    {
//      const size_t num_face_verts = face.vertex_ids.size();
//      for (size_t fv=0; fv<num_face_verts; ++fv)
//      {
//        size_t fvp1 = (fv <(num_face_verts - 1))? fv+1 : 0;
//
//        uint64_t vid0 = face.vertex_ids[fv  ];
//        uint64_t vid1 = face.vertex_ids[fvp1];
//
//        edges.insert(std::make_pair(std::min(vid0, vid1),
//                                    std::max(vid0, vid1)));
//      }
//    }//for face
//
//    // Now we insert the edge centers.
//    for (const auto& edge : edges)
//    {
//      const auto& v0 = m_grid.vertices[edge.first];
//      const auto& v1 = m_grid.vertices[edge.second];
//
//      auto edge_centroid = 0.5*(v0+v1);
//
//      cell_node_info.emplace_back(NodeType::EDGE,
//                                  IdentifyingInfo({{edge.first,edge.second}}),
//                                  edge_centroid);
//    }
//
//    for (auto& face : m_cell.faces)
//      cell_node_info.emplace_back(NodeType::FACE,
//                                  IdentifyingInfo(
//                                  {std::set<uint64_t>(face.vertex_ids.begin(),
//                                                      face.vertex_ids.end())}),
//                                  face.centroid);
//
//    cell_node_info.emplace_back(NodeType::INTERNAL,
//                                IdentifyingInfo({{},{m_cell.global_id}}),
//                                m_cell.centroid);
//  }
//
//  return cell_node_info;
//}