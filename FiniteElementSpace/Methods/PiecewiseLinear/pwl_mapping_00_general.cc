#include "pwl_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/chi_mesh_utils.h"

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


size_t PiecewiseLinear::
  FaceNumNodes(const chi_mesh::Cell& cell, const size_t f) const
{
  return cell.faces.at(f).vertex_ids.size();
}

std::vector<chi_mesh::Vector3> PiecewiseLinear::
  CellNodeLocations(const chi_mesh::Cell& cell) const
{
  std::vector<chi_mesh::Vector3> node_locations;
  node_locations.reserve(NumNodes());

  for (uint64_t vid : cell.vertex_ids)
    node_locations.push_back(m_grid.vertices[vid]);

  return node_locations;
}

size_t PiecewiseLinear::MapFaceNodeToCellNode(const chi_mesh::Cell& cell,
                                              size_t face_index,
                                              size_t face_node_index) const
{
  return m_face_2_cell_map.at(face_index).at(face_node_index);
}

VolumeQPData PiecewiseLinear::
  BuildVolumetricQPData(const chi_mesh::Cell &cell,
                        chi_math::QuadratureOrder order) const
{
  if (cell.Type() == chi_mesh::CellType::SLAB)
    return BuildVolumetricQPDataSlab(cell, order);
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
    return BuildVolumetricQPDataPolygon(cell, order);
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    return BuildVolumetricQPDataPolyhedron(cell,order);
  else
    throw std::logic_error(std::string(__FUNCTION__) + "Unsupported cell-type"
                           " encountered.");
}