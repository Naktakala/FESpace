#ifndef FINITE_ELEMENT_METHOD_LAGRANGE_MAPPING_H
#define FINITE_ELEMENT_METHOD_LAGRANGE_MAPPING_H

#include "Methods/FEMappingBase.h"

namespace chi_math::finite_element
{

  //################################################################# Class def
  /**Lagrange mapping class.*/
  class LagrangeQ2 : public FiniteElementMapping
  {
  protected:
    typedef unsigned int uint;
    std::vector<uint>              m_face_num_nodes;
    std::vector<std::vector<uint>> m_face_2_cell_map;
  public:
    LagrangeQ2(const chi_mesh::Cell& cell,
               const chi_mesh::MeshContinuum& grid,
               std::vector<NodeInfo>& node_list);

  public: //Overriding functions
    QuadratureOrder GetMinimumQuadratureOrder() const override
    {
      return QuadratureOrder::FOURTH;
    }
    void AddRequiredQuadratures(
      chi_math::QuadratureOrder order,
      std::map<QuadratureKey, QuadraturePtr>& quadrature_stack) const override;

    size_t FaceNumNodes(size_t face_index) const override;

    size_t MapFaceNodeToCellNode(size_t face_index,
                                 size_t face_node_index) const override;

    VolumeQPData BuildVolumetricQPData(
      chi_math::QuadratureOrder order,
      std::map<QuadratureKey, QuadraturePtr>& quadrature_stack) const override;
  };

}//namespace chi_math::finite_element

#endif //FINITE_ELEMENT_METHOD_LAGRANGE_MAPPING_H