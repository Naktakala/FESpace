#ifndef FINITE_ELEMENT_METHOD_PWL_MAPPING_H
#define FINITE_ELEMENT_METHOD_PWL_MAPPING_H

#include "Methods/FEMappingBase.h"

namespace chi_math::finite_element
{

  //################################################################# Class def
  /**Piecewise linear Finite Element mapping class*/
  class PiecewiseLinear : public FiniteElementMapping
  {
  protected:
    typedef unsigned int uint;
    std::vector<std::vector<uint>> m_face_2_cell_map;
  public:
    PiecewiseLinear(const chi_mesh::Cell& cell,
                    const chi_mesh::MeshContinuum& grid,
                    std::vector<NodeInfo>& node_list);

  public: //Overriding functions
    std::vector<chi_mesh::Vector3>
    CellNodeLocations(const chi_mesh::Cell& cell) const override;

    size_t FaceNumNodes(const chi_mesh::Cell& cell, size_t f) const override;

    size_t MapFaceNodeToCellNode(const chi_mesh::Cell& cell, size_t face_index,
                                 size_t face_node_index) const override;

    VolumeQPData BuildVolumetricQPData(
      const chi_mesh::Cell& cell,
      chi_math::QuadratureOrder order) const override;

  private:
    VolumeQPData BuildVolumetricQPDataSlab(
      const chi_mesh::Cell& cell,
      chi_math::QuadratureOrder order) const;

    VolumeQPData BuildVolumetricQPDataPolygon(
      const chi_mesh::Cell& cell,
      chi_math::QuadratureOrder order) const;

    VolumeQPData BuildVolumetricQPDataPolyhedron(
      const chi_mesh::Cell& cell,
      chi_math::QuadratureOrder order) const;

  private:
    static double SlabShape(int i, const chi_mesh::Vector3& qpoint);
    static chi_mesh::Vector3 SlabGradShape(int i);
    static double TriangleShape(int i, const chi_mesh::Vector3& qpoint);
    static chi_mesh::Vector3 TriangleGradShape(int i);
    static double TetrahedronShape(int i, const chi_mesh::Vector3& qpoint);
    static chi_mesh::Vector3 TetrahedronGradShape(int i);
  };

}//namespace chi_math::finite_element

#endif //FINITE_ELEMENT_METHOD_PWL_MAPPING_H