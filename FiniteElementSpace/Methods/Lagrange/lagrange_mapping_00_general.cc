#include "lagrange_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

using namespace chi_math::finite_element;

LagrangeQ2::LagrangeQ2(const chi_mesh::Cell &cell,
                       const chi_mesh::MeshContinuum &grid,
                       size_t& node_register_size) :
                       FiniteElementMapping(cell, grid)
{
  if (cell.SubType() == chi_mesh::CellType::SLAB)
  {
    node_register_size = SetNumNodesAndLocalRegister(3, node_register_size);
  }
  else if (cell.SubType() == chi_mesh::CellType::TRIANGLE)
  {
    node_register_size = SetNumNodesAndLocalRegister(6, node_register_size);
  }
  else if (cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
  {
    node_register_size = SetNumNodesAndLocalRegister(9, node_register_size);
  }
  else if (cell.SubType() == chi_mesh::CellType::TETRAHEDRON)
  {
    node_register_size = SetNumNodesAndLocalRegister(10, node_register_size);
  }
  else if (cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
  {
    node_register_size = SetNumNodesAndLocalRegister(27, node_register_size);
  }
  else
    throw std::logic_error(std::string(__FUNCTION__) +
                           "Unsupported cell-type encountered.");

}

LagrangeQ2::LagrangeQ2(const chi_mesh::Cell &cell,
                       const chi_mesh::MeshContinuum &grid,
                       std::vector<NodeInfo>& node_list) :
                       FiniteElementMapping(cell, grid)
{
  size_t reserve_size;
  if (cell.SubType() == chi_mesh::CellType::SLAB)
  {
    reserve_size = SetNumNodesAndLocalRegister(3, node_list.size());
  }
  else if (cell.SubType() == chi_mesh::CellType::TRIANGLE)
  {
    reserve_size = SetNumNodesAndLocalRegister(6, node_list.size());
  }
  else if (cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
  {
    reserve_size = SetNumNodesAndLocalRegister(9, node_list.size());
  }
  else if (cell.SubType() == chi_mesh::CellType::TETRAHEDRON)
  {
    reserve_size = SetNumNodesAndLocalRegister(10, node_list.size());
  }
  else if (cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
  {
    reserve_size = SetNumNodesAndLocalRegister(27, node_list.size());
  }
  else
    throw std::logic_error(std::string(__FUNCTION__) +
                           "Unsupported cell-type encountered.");

  //======================================== Add the nodes to node_list
//  node_list.reserve(reserve_size);
  // Add cell vertices
  // All of the cell vertices participate in Q2 vertices
  // so we can add all them here abstractly.
  for (uint64_t cvid : m_cell.vertex_ids)
    node_list.emplace_back(NodeType::CORNER,
                           IdentifyingInfo({{cvid}}),
                           m_grid.vertices[cvid]);

  // Cell-subtype specifics
  // For a slab we only have to add the cell-centroid.
  if (m_cell.SubType() == chi_mesh::CellType::SLAB)
    node_list.emplace_back(NodeType::INTERNAL,
                           IdentifyingInfo({{},{m_cell.global_id}}),
                           m_cell.centroid);
  // For a triangle we need to add side centroids.
  else if (m_cell.SubType() == chi_mesh::CellType::TRIANGLE)
  {
    for (auto& face : m_cell.faces)
      node_list.emplace_back(NodeType::FACE,
                             IdentifyingInfo(
                               {std::set<uint64_t>(face.vertex_ids.begin(),
                                                   face.vertex_ids.end())}),
                             face.centroid);
  }
  // For a quadrilateral we add the edge centers and cell centroid.
  else if (m_cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
  {
    for (auto& face : m_cell.faces)
      node_list.emplace_back(NodeType::FACE,
                             IdentifyingInfo(
                             {std::set<uint64_t>(face.vertex_ids.begin(),
                                                 face.vertex_ids.end())}),
                             face.centroid);

    node_list.emplace_back(NodeType::INTERNAL,
                           IdentifyingInfo({{},{m_cell.global_id}}),
                           m_cell.centroid);
  }
  else if (m_cell.SubType() == chi_mesh::CellType::TETRAHEDRON or
           m_cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
  {
    // First we build the edges
    std::set<std::pair<uint64_t, uint64_t>> edges;

    for (auto& face : m_cell.faces)
    {
      const size_t num_face_verts = face.vertex_ids.size();
      for (size_t fv=0; fv<num_face_verts; ++fv)
      {
        size_t fvp1 = (fv <(num_face_verts - 1))? fv+1 : 0;

        uint64_t vid0 = face.vertex_ids[fv  ];
        uint64_t vid1 = face.vertex_ids[fvp1];

        edges.insert(std::make_pair(std::min(vid0, vid1),
                                    std::max(vid0, vid1)));
      }
    }//for face

    // Now we insert the edge centers.
    for (const auto& edge : edges)
    {
      const auto& v0 = m_grid.vertices[edge.first];
      const auto& v1 = m_grid.vertices[edge.second];

      auto edge_centroid = 0.5*(v0+v1);

      node_list.emplace_back(NodeType::EDGE,
                             IdentifyingInfo({{edge.first,edge.second}}),
                             edge_centroid);
    }

    // Now we insert the face centers
    for (auto& face : m_cell.faces)
      node_list.emplace_back(NodeType::FACE,
                             IdentifyingInfo(
                             {std::set<uint64_t>(face.vertex_ids.begin(),
                                                 face.vertex_ids.end())}),
                             face.centroid);

    // Now we insert the cell center
    node_list.emplace_back(NodeType::INTERNAL,
                           IdentifyingInfo({{},{m_cell.global_id}}),
                           m_cell.centroid);
  }
}

size_t LagrangeQ2::
  FaceNumNodes(const chi_mesh::Cell& cell, const size_t f) const
{
  return cell.faces.at(f).vertex_ids.size();
}

std::vector<chi_mesh::Vector3> LagrangeQ2::
  CellNodeLocations(const chi_mesh::Cell& cell) const
{
  std::vector<chi_mesh::Vector3> node_locations;
  node_locations.reserve(NumNodes());

//  const auto node_info = GetNodeInfo();
//  for (auto& node : node_info)
//    node_locations.push_back(node.location);

  return node_locations;
}

size_t LagrangeQ2::MapFaceNodeToCellNode(const chi_mesh::Cell& cell,
                                              size_t face_index,
                                              size_t face_node_index) const
{
  return 0; //TODO: Fix
}

