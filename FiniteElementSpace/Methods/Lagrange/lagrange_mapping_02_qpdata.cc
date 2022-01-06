#include "lagrange_mapping.h"

using namespace chi_math::finite_element;

VolumeQPData LagrangeQ2::BuildVolumetricQPData(
      chi_math::QuadratureOrder q_order,
      std::map<QuadratureKey, QuadraturePtr>& quadrature_stack) const
{
  const std::string fname = __FUNCTION__;
  const auto cell_sub_type = m_cell.SubType();
  const QuadratureKey quadrature_key = {cell_sub_type, q_order};

  const auto& quadrature = *quadrature_stack.at(quadrature_key);

  std::vector<unsigned int>     quadrature_point_indices; ///< qp index only
  VecVec3                       qpoints_xyz             ; ///< qp index only
  std::vector<VecDbl>           shape_value             ; ///< Node i, then qp
  std::vector<VecVec3>          shape_grad              ; ///< Node i, then qp
  VecDbl                        JxW                     ; ///< qp index only

  const size_t num_q_points = quadrature.qpoints.size();
  const size_t num_nodes = NumNodes();

  //=================================== Init volumetric quadrature
  quadrature_point_indices.reserve(num_q_points);
  for (unsigned int qp=0; qp<num_q_points; ++qp)
    quadrature_point_indices.push_back(qp);

  shape_value.reserve(num_nodes);
  shape_grad.reserve(num_nodes);
  for (size_t i=0; i < num_nodes; i++)
  {
    VecDbl  node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(num_q_points);
    node_shape_grad.reserve(num_q_points);

    for (const auto& qpoint : quadrature.qpoints)
    {
      double shape_i = 0.0;
//        TriangleShape(imapped,qpoint) + beta * TriangleShape(2, qpoint); //TODO:

      auto gradshape_i = chi_mesh::Vector3(0,0,0);
//        TriangleGradShape(imapped) + beta * TriangleGradShape(2); //TODO:

      chi_mesh::Matrix3x3 J_Tinv; //TODO:

      node_shape_value.push_back(shape_i);
      node_shape_grad.push_back(J_Tinv * gradshape_i);
    }//for qp

    shape_value.push_back(node_shape_value);
    shape_grad.push_back(node_shape_grad);
  }//for i

  const chi_mesh::Vector3 v0 = m_cell.centroid;

  JxW.reserve(num_q_points);
  qpoints_xyz.reserve(num_q_points);
  for (size_t qp=0; qp<num_q_points; ++qp)
  {
    const auto w = quadrature.weights[qp];

    chi_mesh::Matrix3x3 J; //TODO:

    JxW.push_back(J.Det() * w);

    const auto& qp_xyz_tilde = quadrature.qpoints[qp];
    qpoints_xyz.push_back(v0 + J * qp_xyz_tilde);
  }//for qp

  return VolumeQPData(quadrature_point_indices,
                      qpoints_xyz,
                      shape_value,
                      shape_grad,
                      JxW);
}

