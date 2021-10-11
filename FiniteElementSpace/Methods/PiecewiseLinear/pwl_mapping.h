#ifndef FINITE_ELEMENT_METHOD_PWL_MAPPING_H
#define FINITE_ELEMENT_METHOD_PWL_MAPPING_H

#include "Methods/FEMappingBase.h"

namespace chi_math::finite_element
{

  //################################################################# Class def
  /**Piecewise linear Finite Element mapping class*/
  class PiecewiseLinear : public FiniteElementMapping
  {
  public:
    explicit PiecewiseLinear(const chi_mesh::Cell& cell,
                             const chi_mesh::MeshContinuum& grid);

    size_t CellNumNodes(const chi_mesh::Cell& cell) const override;

    std::vector<chi_mesh::Vector3>
    CellNodeLocations(const chi_mesh::Cell& cell) const override;
  };

}

#endif //FINITE_ELEMENT_METHOD_PWL_MAPPING_H