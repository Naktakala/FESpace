#ifndef FINITE_ELEMENT_METHOD_FINITE_VOLUME_MAPPING_H
#define FINITE_ELEMENT_METHOD_FINITE_VOLUME_MAPPING_H

#include "Methods/FEMappingBase.h"

namespace chi_math::finite_element
{

  //################################################################# Class def
  /**Finite volume mapping class*/
  class FiniteVolume : public FiniteElementMapping
  {
  public:
    FiniteVolume(const chi_mesh::Cell& cell,
                 const chi_mesh::MeshContinuum& grid,
                 std::vector<NodeInfo>& node_list);

  public:
    std::vector<chi_mesh::Vector3>
    CellNodeLocations(const chi_mesh::Cell& cell) const override;

    size_t FaceNumNodes(const chi_mesh::Cell& cell, size_t f) const override;

    size_t MapFaceNodeToCellNode(const chi_mesh::Cell& cell, size_t face_index,
                                 size_t face_node_index) const override;

    VolumeQPData BuildVolumetricQPData(
      const chi_mesh::Cell& cell,
      chi_math::QuadratureOrder order) const override
    {
      return VolumeQPData({}, {}, {}, {}, {});
    }

  };

}//namespace chi_math::finite_element

#endif //FINITE_ELEMENT_METHOD_FINITE_VOLUME_MAPPING_H