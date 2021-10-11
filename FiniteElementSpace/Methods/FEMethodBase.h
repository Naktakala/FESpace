#ifndef FINITE_ELEMENT_METHOD_BASE_H
#define FINITE_ELEMENT_METHOD_BASE_H

#include "FEMappingBase.h"

#include "ChiMath/Quadratures/quadrature.h"

#include "ChiMesh/Cell/cell.h"

#include <vector>
#include <map>
#include <memory>

namespace chi_math::finite_element
{
  /**Base class of any finite element method.*/
  class FiniteElementMethod
  {
  private:
  protected:
    std::vector<std::unique_ptr<FiniteElementMapping>> cell_mapping_info;

  public:
    FiniteElementMethod() = default;

    virtual FEMappingPtr
      MakeCellMapping(const chi_mesh::Cell& cell,
                      const chi_mesh::MeshContinuum& grid) const = 0;

    /**Returns the number of nodes associated with a given cell.*/
    virtual
    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const = 0;

    /**Returns the locations of the nodes associated with a given cell.*/
    virtual
    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell,
                         const chi_mesh::MeshContinuum& grid) const = 0;

    using QuadraturePtr = std::unique_ptr<chi_math::Quadrature>;
    /**Returns a map of cell-types, to quadrature-formulas required
     * to integrate both over the cell-volume and the cell-surface.*/
    virtual std::map<chi_mesh::CellType, QuadraturePtr>
      MapRequiredQuadratures(const chi_mesh::MeshContinuum& grid,
                             const chi_math::QuadratureOrder& q_order) const = 0;

    virtual ~FiniteElementMethod() {};
  };

}//namespace chi_math

#endif //FINITE_ELEMENT_METHOD_BASE_H