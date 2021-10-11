#ifndef FINITE_ELEMENT_METHOD_PIECEWISE_LINEAR_H
#define FINITE_ELEMENT_METHOD_PIECEWISE_LINEAR_H

#include "../FEMethodBase.h"
#include "pwl_mapping.h"

#include "ChiMath/Quadratures/quadrature_line.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_quadrilateral.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"
#include "ChiMath/Quadratures/quadrature_hexahedron.h"

namespace chi_math::finite_element
{
  /**Object for mapping cells onto a piecewise-linear finite
   * elemnt space.*/
  class PiecewiseLinearDiscontinuous : public FiniteElementMethod
  {
  public:

    FEMappingPtr
      MakeCellMapping(const chi_mesh::Cell& cell,
                      const chi_mesh::MeshContinuum& grid) const override
    {
      auto map = std::make_unique<PiecewiseLinearMapping>(cell, grid);
      return std::unique_ptr<FiniteElementMapping>(map.release());
    }

    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override;

    std::vector<chi_mesh::Vector3>
      GetCellNodeLocations(const chi_mesh::Cell& cell,
                          const chi_mesh::MeshContinuum& grid) const override;

    std::map<chi_mesh::CellType, QuadraturePtr>
      MapRequiredQuadratures(
        const chi_mesh::MeshContinuum& grid,
        const chi_math::QuadratureOrder& q_order) const override;
  };

}//namespace chi_math::finite_element_mappings

#endif//FINITE_ELEMENT_METHOD_PIECEWISE_LINEAR_H