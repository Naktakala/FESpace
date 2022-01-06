#include "pwl_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/Quadratures/quadrature_line.h"

using namespace chi_math::finite_element;

VolumeQPData PiecewiseLinear::
  BuildVolumetricQPDataSlab(chi_math::QuadratureOrder order) const
{
  std::vector<unsigned int>     quadrature_point_indices; ///< qp index only
  VecVec3                       qpoints_xyz             ; ///< qp index only
  std::vector<VecDbl>           shape_value             ; ///< Node i, then qp
  std::vector<VecVec3>          shape_grad              ; ///< Node i, then qp
  VecDbl                        JxW                     ; ///< qp index only

  chi_mesh::Vector3 v0 = m_grid.vertices[m_cell.vertex_ids[0]];
  chi_math::QuadratureLine qdata(order);

  //=================================== Determine number of internal qpoints
  size_t ttl_num_vol_qpoints = qdata.qpoints.size();
  size_t num_nodes = m_cell.vertex_ids.size();

  //=================================== Init volumetric quadrature
  quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp=0; qp<ttl_num_vol_qpoints; ++qp)
    quadrature_point_indices.push_back(qp);

  const double J = m_volume;
  shape_value.reserve(num_nodes);
  shape_grad.reserve(num_nodes);
  for (int i=0; i < num_nodes; i++)
  {
    VecDbl  node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_vol_qpoints);
    node_shape_grad.reserve(ttl_num_vol_qpoints);

    for (const auto& qpoint : qdata.qpoints)
    {
      node_shape_value.push_back(SlabShape(i,qpoint));
      node_shape_grad.push_back((1.0/J) * SlabGradShape(i));
    }//for qp

    shape_value.push_back(node_shape_value);
    shape_grad.push_back(node_shape_grad);
  }//for i

  JxW.reserve(ttl_num_vol_qpoints);
  qpoints_xyz.reserve(ttl_num_vol_qpoints);
  for (size_t qp=0; qp<ttl_num_vol_qpoints; ++qp)
  {
    const double w = qdata.weights[qp];
    JxW.push_back(J * w);

    const double qp_xyz_tilde = qdata.qpoints[qp][0];
    qpoints_xyz.push_back(v0 + J * chi_mesh::Vector3(0.0,0.0,qp_xyz_tilde));
  }//for qp

  return VolumeQPData(quadrature_point_indices,
                      qpoints_xyz,
                      shape_value,
                      shape_grad,
                      JxW);
}