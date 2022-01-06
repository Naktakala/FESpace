#include "pwl_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/chi_mesh_utils.h"

#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"

using namespace chi_math::finite_element;

//###################################################################
/**Constructs a Piecewise linear mapping of a cell.*/
PiecewiseLinear::
  PiecewiseLinear(const chi_mesh::Cell& cell,
                  const chi_mesh::MeshContinuum& grid,
                  std::vector<NodeInfo>& node_list) :
                  FiniteElementMapping(cell, grid)
{
  const size_t num_nodes = cell.vertex_ids.size();

  SetNumNodesAndLocalRegister(num_nodes, node_list.size());

  //======================================== Add the nodes to node_list
  for (uint64_t cvid : m_cell.vertex_ids)
    node_list.emplace_back(NodeType::CORNER,
                           IdentifyingInfo({{cvid}}),
                           m_grid.vertices[cvid]);

  //======================================== Compute mandatories
  // Volume and Area
  auto cell_info = chi_mesh::ComputeCellVolumeAndFaceAreas(cell, grid);
  m_volume = cell_info.volume;
  m_face_areas = std::move(cell_info.face_areas);

  // m_face_2_cell node map
  {
    // Making vertex-id to cell-node map
    // This is done to make the face_2_cell
    // mapping more efficient
    std::map<uint64_t, uint> cell_node_id_2_node_index_map;
    uint node_index=0;
    for (uint64_t node_id : cell.vertex_ids)
      cell_node_id_2_node_index_map[node_id] = node_index++;

    //============================== face_2_cell node map
    m_face_2_cell_map.reserve(cell.faces.size());
    for (auto& face : cell.faces)
    {
      std::vector<uint> face_node_map;
      face_node_map.reserve(face.vertex_ids.size());
      for (uint64_t fnode_id : face.vertex_ids)
        face_node_map.push_back(cell_node_id_2_node_index_map.at(fnode_id));

      m_face_2_cell_map.push_back(std::move(face_node_map));
    }//for face
  }
}


void PiecewiseLinear::
  AddRequiredQuadratures(
    const QuadratureOrder q_order,
    std::map<QuadratureKey, QuadraturePtr> &quadrature_stack) const
{
  const std::string fname = __FUNCTION__;

  const auto cell_type = m_cell.Type();

  const QuadratureKey volm_quad_key = {cell_type, q_order};

  //=================================== Check if stack already has required
  //                                    quadrature
  if (quadrature_stack.count(volm_quad_key) > 0)
    return;

  //=================================== Load volumetric quadrature
  QuadraturePtr volm_quad;
  QuadraturePtr surf_quad;
  QuadratureKey surf_quad_key;

  using namespace chi_mesh; using namespace chi_math; using namespace std;

  if      (cell_type == CellType::SLAB)
  {
    volm_quad = make_unique<QuadratureGaussLegendre>(q_order);
    surf_quad = make_unique<QuadratureGaussLegendre>(QuadratureOrder::CONSTANT);
    surf_quad_key = {CellType::SLAB, QuadratureOrder::CONSTANT};
  }
  else if (cell_type == CellType::POLYGON)
  {
    volm_quad = make_unique<QuadratureTriangle     >(q_order);
    surf_quad = make_unique<QuadratureGaussLegendre>(q_order);
    surf_quad_key = {CellType::SLAB, q_order};
  }
  else if (cell_type == CellType::POLYHEDRON)
  {
    volm_quad = make_unique<QuadratureTetrahedron  >(q_order);
    surf_quad = make_unique<QuadratureTriangle     >(q_order);
    surf_quad_key = {CellType::TRIANGLE, q_order};
  }
  else
    throw std::logic_error(fname + "Unsupported cell-type encountered.");

  quadrature_stack.insert( make_pair(volm_quad_key, move(volm_quad)) );
  quadrature_stack.insert( make_pair(surf_quad_key, move(surf_quad)) );

}//AddRequiredQuadratures


size_t PiecewiseLinear::
  FaceNumNodes(const size_t face_index) const
{
  return m_cell.faces.at(face_index).vertex_ids.size();
}

size_t PiecewiseLinear::MapFaceNodeToCellNode(const size_t face_index,
                                              const size_t face_node_index) const
{
  return m_face_2_cell_map.at(face_index).size();
}

