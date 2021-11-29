#include "fe_space_qp_data.h"

namespace chi_math::finite_element
{
  VolumeQPData::VolumeQPData(std::vector<unsigned int> quadrature_point_indices,
                             VecVec3                   qpoints_xyz,
                             std::vector<VecDbl>       shape_value,
                             std::vector<VecVec3>      shape_grad,
                             VecDbl                    JxW) :
    m_quadrature_point_indices(std::move(quadrature_point_indices)),
    m_qpoints_xyz             (std::move(qpoints_xyz             )),
    m_shape_value             (std::move(shape_value             )),
    m_shape_grad              (std::move(shape_grad              )),
    m_JxW                     (std::move(JxW                     ))
  {}

  const std::vector<unsigned int>& VolumeQPData::QuadraturePointIndices() const
  {
    return m_quadrature_point_indices;
  }
  chi_mesh::Vector3 VolumeQPData::QPointXYZ(unsigned int qp) const
  {
    return m_qpoints_xyz.at(qp);
  }
  double VolumeQPData::ShapeValue(unsigned int i, unsigned int qp) const
  {
    auto& qp_data = m_shape_value.at(i);
    return qp_data.at(qp);
  }
  chi_mesh::Vector3 VolumeQPData::ShapeGrad(unsigned int i, unsigned int qp) const
  {
    auto& qp_data = m_shape_grad.at(i);
    return qp_data.at(qp);
  }
  double VolumeQPData::JxW(unsigned int qp) const
  {
    return m_JxW.at(qp);
  }
}//chi_math::finite_element
