#include "pwld.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

using namespace chi_math::finite_element;

//###################################################################
size_t PiecewiseLinearDiscontinuous::
  GetCellNumNodes(const chi_mesh::Cell& cell) const
{
  return cell.vertex_ids.size();
}

std::vector<chi_mesh::Vector3> PiecewiseLinearDiscontinuous::
  GetCellNodeLocations(const chi_mesh::Cell& cell,
                       const chi_mesh::MeshContinuum& grid) const
{
  size_t num_nodes = GetCellNumNodes(cell);
  std::vector<chi_mesh::Vector3> node_locations;
  node_locations.reserve(num_nodes);

  for (uint64_t vid : cell.vertex_ids)
    node_locations.push_back(grid.vertices[vid]);

  return node_locations;
}

std::map<chi_mesh::CellType, FiniteElementMethod::QuadraturePtr>
  PiecewiseLinearDiscontinuous::
    MapRequiredQuadratures(const chi_mesh::MeshContinuum& grid,
                           const chi_math::QuadratureOrder& q_order) const
{
  using LineQuad = chi_math::QuadratureLine;
  using TriQuad = chi_math::QuadratureTriangle;
  using TetQuad = chi_math::QuadratureTetrahedron;

  std::map<chi_mesh::CellType, QuadraturePtr> qmap;

  for (const auto& cell : grid.local_cells)
  {
    //Check if cell-type already mapped
    if (qmap.find(cell.Type()) != qmap.end()) continue;

    if (cell.Type() == chi_mesh::CellType::SLAB)
      qmap[cell.Type()] = std::make_unique<LineQuad>(q_order);
    else if (cell.Type() == chi_mesh::CellType::POLYGON)
      qmap[cell.Type()] = std::make_unique<TriQuad>(q_order);
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      qmap[cell.Type()] = std::make_unique<TetQuad>(q_order);
    else
      throw std::logic_error("Unsupported cell-type encountered when"
                             " building quadratures for a PiecewiseLinear"
                             " finite element space.");
  }//for cell

  return qmap;
}