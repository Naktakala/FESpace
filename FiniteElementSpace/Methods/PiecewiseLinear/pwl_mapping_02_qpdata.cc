#include "pwl_mapping.h"

using namespace chi_math::finite_element;

VolumeQPData PiecewiseLinear::
  BuildVolumetricQPData(
    chi_math::QuadratureOrder order,
    std::map<QuadratureKey, QuadraturePtr>& quadrature_stack) const
{
  if      (m_cell.Type() == chi_mesh::CellType::SLAB)
    return BuildVolumetricQPDataSlab(order);
  else if (m_cell.Type() == chi_mesh::CellType::POLYGON)
    return BuildVolumetricQPDataPolygon(order);
  else if (m_cell.Type() == chi_mesh::CellType::POLYHEDRON)
    return BuildVolumetricQPDataPolyhedron(order);
  else
    throw std::logic_error(std::string(__FUNCTION__) + "Unsupported cell-type"
                           " encountered.");
}