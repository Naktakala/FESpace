#ifndef FINITE_ELEMENT_METHOD_FINITE_VOLUME_H
#define FINITE_ELEMENT_METHOD_FINITE_VOLUME_H

#include "../FEMethodBase.h"
#include "fv_mapping.h"

#include <map>

namespace chi_math::finite_element
{
  //################################################################# Class def
  /**Object for mapping cells to a finite-volume space.*/
  class FiniteVolume : public FiniteElementMethod
  {
  public:

    FEMappingPtr
      MakeCellMapping(const chi_mesh::Cell& cell,
                      const chi_mesh::MeshContinuum& grid) const override
    {
      auto map = std::make_unique<FiniteVolumeMapping>(cell, grid);
      return std::unique_ptr<FiniteElementMapping>(map.release());
    }

    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override
    {return 1;}

    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell,
                         const chi_mesh::MeshContinuum& grid) const override
    {return {cell.centroid};}

    std::map<chi_mesh::CellType, QuadraturePtr>
    MapRequiredQuadratures(const chi_mesh::MeshContinuum& grid,
                           const chi_math::QuadratureOrder& q_order) const override
    {return {};}
  };

}//namespace chi_math::finite_element_mappings

#endif//FINITE_ELEMENT_METHOD_FINITE_VOLUME_H