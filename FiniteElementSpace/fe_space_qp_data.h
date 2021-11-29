#ifndef FESPACETEST_FE_SPACE_QP_DATA_H
#define FESPACETEST_FE_SPACE_QP_DATA_H

#include <vector>
#include "ChiMesh/chi_mesh.h"

namespace chi_math::finite_element
{
  typedef std::vector<double> VecDbl;
  typedef std::vector<chi_mesh::Vector3> VecVec3;

  //#############################################
  /**Stored relevant quadrature point information
   * for volumetric integrals.*/
  class VolumeQPData
  {
  protected:
    std::vector<unsigned int>     m_quadrature_point_indices; ///< qp index only
    VecVec3                       m_qpoints_xyz             ; ///< qp index only
    std::vector<VecDbl>           m_shape_value             ; ///< Node i, then qp
    std::vector<VecVec3>          m_shape_grad              ; ///< Node i, then qp
    VecDbl                        m_JxW                     ; ///< qp index only

  public:
    VolumeQPData(std::vector<unsigned int> quadrature_point_indices,
                   VecVec3                   qpoints_xyz,
                   std::vector<VecDbl>       shape_value,
                   std::vector<VecVec3>      shape_grad,
                   VecDbl                    JxW);
    const std::vector<unsigned int>& QuadraturePointIndices() const;
    chi_mesh::Vector3 QPointXYZ(unsigned int qp) const;
    double ShapeValue(unsigned int i, unsigned int qp) const;
    chi_mesh::Vector3 ShapeGrad(unsigned int i, unsigned int qp) const;
    double JxW(unsigned int qp) const;
  };
}//namespace chi_math::finite_element

#endif //FESPACETEST_FE_SPACE_QP_DATA_H
