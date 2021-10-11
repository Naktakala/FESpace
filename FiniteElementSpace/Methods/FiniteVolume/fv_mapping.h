#ifndef FINITE_ELEMENT_METHOD_FINITE_VOLUME_MAPPING_H
#define FINITE_ELEMENT_METHOD_FINITE_VOLUME_MAPPING_H

#include "Methods/FEMappingBase.h"

namespace chi_math::finite_element
{

//################################################################# Class def
/**Finite volume mapping class*/
class FiniteVolumeMapping : public FiniteElementMapping
{
public:
  explicit FiniteVolumeMapping(const chi_mesh::Cell& cell,
                               const chi_mesh::MeshContinuum& grid);

  size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override;

  std::vector<chi_mesh::Vector3>
  GetCellNodeLocations(const chi_mesh::Cell& cell,
                       const chi_mesh::MeshContinuum& grid) const override;
};

}

#endif //FINITE_ELEMENT_METHOD_FINITE_VOLUME_MAPPING_H