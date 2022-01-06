#include "lagrange_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_quadrilateral.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"
#include "ChiMath/Quadratures/quadrature_hexahedron.h"

using namespace chi_math::finite_element;

LagrangeQ2::LagrangeQ2(const chi_mesh::Cell &cell,
                       const chi_mesh::MeshContinuum &grid,
                       std::vector<NodeInfo>& node_list) :
                       FiniteElementMapping(cell, grid)
{
  const std::string fname = __FUNCTION__;

  if (cell.SubType() == chi_mesh::CellType::SLAB)
  {
    SetNumNodesAndLocalRegister(3, node_list.size());
    m_face_num_nodes.assign(cell.faces.size(), 1);
  }
  else if (cell.SubType() == chi_mesh::CellType::TRIANGLE)
  {
    SetNumNodesAndLocalRegister(6, node_list.size());
    m_face_num_nodes.assign(cell.faces.size(), 3);
  }
  else if (cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
  {
    SetNumNodesAndLocalRegister(9, node_list.size());
    m_face_num_nodes.assign(cell.faces.size(), 3);
  }
  else if (cell.SubType() == chi_mesh::CellType::TETRAHEDRON)
  {
    SetNumNodesAndLocalRegister(10, node_list.size());
    m_face_num_nodes.assign(cell.faces.size(), 6);
  }
  else if (cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
  {
    SetNumNodesAndLocalRegister(27, node_list.size());
    m_face_num_nodes.assign(cell.faces.size(), 9);
  }
  else
    throw std::logic_error(fname + "Unsupported cell-subtype encountered.");

  const size_t nls_before = node_list.size();

  //======================================== Add the nodes to node_list
  //==================== Add cell vertices
  // All of the cell vertices participate in Q2 vertices
  // so we can add all them here abstractly.
  for (uint64_t cvid : m_cell.vertex_ids)
    node_list.emplace_back(NodeType::CORNER,
                           IdentifyingInfo({{cvid}}),
                           m_grid.vertices[cvid]);

  //==================== Cell-subtype specifics
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

  const size_t nls_after = node_list.size();

  if ((nls_after - nls_before) != NumNodes())
    throw std::logic_error(fname + ": Cell node accounting error.");

  //=================================== Initialize face_2_cell_map
  m_face_2_cell_map.assign(cell.faces.size(), std::vector<uint>());
  for (size_t f=0; f<cell.faces.size(); ++f)
    m_face_2_cell_map[f].assign(m_face_num_nodes[f], 0);

  //=================================== Map face node
  // The Slab is so simple we can just hard code
  if (cell.SubType() == chi_mesh::CellType::SLAB)
  {
    m_face_2_cell_map[0][0] = 0;
    m_face_2_cell_map[1][0] = 1;
  }
  // The faces of a triangle and a quadrilateral are both just
  // edges with an additional node at the center. We can use the
  // same logic for both of them
  else if (cell.SubType() == chi_mesh::CellType::TRIANGLE or
           cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
  {
    size_t f=0;
    for (auto& face : m_cell.faces)
    {
      std::vector<NodeInfo> fnod_list;
      // Add face corners
      for (uint64_t fvid : face.vertex_ids)
        fnod_list.emplace_back(NodeType::CORNER,
                               IdentifyingInfo({{fvid}}),
                               m_grid.vertices[fvid]);

      // Add face centroid
      fnod_list.emplace_back(NodeType::FACE,
                             IdentifyingInfo(
                               {std::set<uint64_t>(face.vertex_ids.begin(),
                                                   face.vertex_ids.end())}),
                             face.centroid);

      if (fnod_list.size() != m_face_num_nodes[f])
        throw std::logic_error(fname + "Face node mapping error.");

      for (size_t fn=0; fn<fnod_list.size(); ++fn)
        for (size_t cn=0; cn<NumNodes(); ++cn)
          if (node_list[nls_before + cn] == fnod_list[fn])
          {
            m_face_2_cell_map[f][fn] = cn;
            break;
          }

      ++f;
    }//for face
  }
  // The faces of a tet and a hex have their edges split.
  // The faces of a hex have an additional node at the centroid.
  else if (cell.SubType() == chi_mesh::CellType::TETRAHEDRON or
           cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
  {
    size_t f=0;
    for (auto& face : m_cell.faces)
    {
      std::vector<NodeInfo> fnod_list;
      // Add face corners
      for (uint64_t fvid : face.vertex_ids)
        fnod_list.emplace_back(NodeType::CORNER,
                               IdentifyingInfo({{fvid}}),
                               m_grid.vertices[fvid]);

      // Add edge centers
      const size_t num_face_verts = face.vertex_ids.size();
      for (size_t fv=0; fv<num_face_verts; ++fv)
      {
        size_t fvp1 = (fv <(num_face_verts - 1))? fv+1 : 0;

        uint64_t vid0 = face.vertex_ids[fv  ];
        uint64_t vid1 = face.vertex_ids[fvp1];

        std::pair<uint64_t, uint64_t> edge(std::min(vid0, vid1),
                                           std::max(vid0, vid1));

        const auto& v0 = m_grid.vertices[edge.first];
        const auto& v1 = m_grid.vertices[edge.second];

        auto edge_centroid = 0.5*(v0+v1);

        fnod_list.emplace_back(NodeType::EDGE,
                               IdentifyingInfo({{edge.first,edge.second}}),
                               edge_centroid);
      }

      // Add face centroid only for hex
      if (cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
        fnod_list.emplace_back(NodeType::FACE,
                               IdentifyingInfo(
                                 {std::set<uint64_t>(face.vertex_ids.begin(),
                                                     face.vertex_ids.end())}),
                               face.centroid);

      if (fnod_list.size() != m_face_num_nodes[f])
        throw std::logic_error(fname + "Face node mapping error.");

      for (size_t fn=0; fn<fnod_list.size(); ++fn)
        for (size_t cn=0; cn<NumNodes(); ++cn)
          if (node_list[nls_before + cn] == fnod_list[fn])
          {
            m_face_2_cell_map[f][fn] = cn;
            break;
          }

      ++f;
    }//for face
  }

}

void LagrangeQ2::
  AddRequiredQuadratures(
    const chi_math::QuadratureOrder q_order,
    std::map<QuadratureKey, QuadraturePtr> &quadrature_stack) const
{
  const std::string fname = __FUNCTION__;

  const auto cell_sub_type = m_cell.SubType();

  const QuadratureKey volm_quad_key = {cell_sub_type, q_order};

  //=================================== Check if stack already has required
  //                                    quadrature
  if (quadrature_stack.count(volm_quad_key) > 0)
    return;

  //=================================== Load volumetric quadrature
  QuadraturePtr volm_quad;
  QuadraturePtr surf_quad;
  QuadratureKey surf_quad_key;

  using namespace chi_mesh; using namespace chi_math; using namespace std;

  if      (cell_sub_type == CellType::SLAB)
  {
    volm_quad = make_unique<QuadratureGaussLegendre>(q_order);
    surf_quad = make_unique<QuadratureGaussLegendre>(QuadratureOrder::CONSTANT);
    surf_quad_key = {CellType::SLAB, QuadratureOrder::CONSTANT};
  }
  else if (cell_sub_type == CellType::TRIANGLE)
  {
    volm_quad = make_unique<QuadratureTriangle     >(q_order);
    surf_quad = make_unique<QuadratureGaussLegendre>(q_order);
    surf_quad_key = {CellType::SLAB, q_order};
  }
  else if (cell_sub_type == CellType::QUADRILATERAL)
  {
    volm_quad = make_unique<QuadratureQuadrilateral>(q_order);
    surf_quad = make_unique<QuadratureGaussLegendre>(q_order);
    surf_quad_key = {CellType::SLAB, q_order};
  }
  else if (cell_sub_type == CellType::TETRAHEDRON)
  {
    volm_quad = make_unique<QuadratureTetrahedron  >(q_order);
    surf_quad = make_unique<QuadratureTriangle     >(q_order);
    surf_quad_key = {CellType::TRIANGLE, q_order};
  }
  else if (cell_sub_type == CellType::HEXAHEDRON)
  {
    volm_quad = make_unique<QuadratureHexahedron   >(q_order);
    surf_quad = make_unique<QuadratureQuadrilateral>(q_order);
    surf_quad_key = {CellType::QUADRILATERAL, q_order};
  }
  else
    throw std::logic_error(fname + "Unsupported cell-subtype encountered.");

  quadrature_stack.insert( make_pair(volm_quad_key, move(volm_quad)) );
  quadrature_stack.insert( make_pair(surf_quad_key, move(surf_quad)) );

}//AddRequiredQuadratures

size_t LagrangeQ2::
  FaceNumNodes(const size_t face_index) const
{
  return m_face_2_cell_map.at(face_index).size();
}

size_t LagrangeQ2::MapFaceNodeToCellNode(const size_t face_index,
                                         const size_t face_node_index) const
{
  return m_face_2_cell_map.at(face_index).at(face_node_index);
}

