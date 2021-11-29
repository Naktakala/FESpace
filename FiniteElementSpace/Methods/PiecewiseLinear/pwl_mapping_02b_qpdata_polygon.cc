#include "pwl_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/Quadratures/quadrature_triangle.h"

using namespace chi_math::finite_element;

VolumeQPData PiecewiseLinear::
  BuildVolumetricQPDataPolygon(const chi_mesh::Cell &cell,
                            chi_math::QuadratureOrder order) const
{
  std::vector<unsigned int>     quadrature_point_indices; ///< qp index only
  VecVec3                       qpoints_xyz             ; ///< qp index only
  std::vector<VecDbl>           shape_value             ; ///< Node i, then qp
  std::vector<VecVec3>          shape_grad              ; ///< Node i, then qp
  VecDbl                        JxW                     ; ///< qp index only

  chi_math::QuadratureTriangle qdata(order);

  //=================================== Determine number of internal qpoints
  size_t num_tris = cell.vertex_ids.size();
  size_t num_vol_qpoints = qdata.qpoints.size();
  size_t ttl_num_vol_qpoints = num_tris * num_vol_qpoints;
  size_t num_nodes = cell.vertex_ids.size();

  //=================================== Determine triangle data
  struct TriangleData
  {
    chi_mesh::Vector3   v0;
    chi_mesh::Matrix3x3 J;
    chi_mesh::Matrix3x3 J_Tinv;
    double detJ=0.0;

    std::vector<int> imap;
  };
  std::vector<TriangleData> triangle_data(num_tris);
  const double beta = 1.0/ static_cast<double>(num_tris);
  for (size_t s=0; s<num_tris; ++s)
  {
    size_t sp1 = (s < (num_tris-1))? s+1 : 0;
    const auto& v0 = m_grid.vertices[cell.vertex_ids[  s]];
    const auto& v1 = m_grid.vertices[cell.vertex_ids[sp1]];
    const auto& v2 = cell.centroid;

    const auto v01 = v1-v0;
    const auto v02 = v2-v0;

    chi_mesh::Matrix3x3 J;
    J.SetColJVec(0,v01);
    J.SetColJVec(1,v02);

    const auto J_Tinv = J.Transpose().Inverse();

    std::vector<int> imap(num_nodes, -1);
    imap[  s] = 0;
    imap[sp1] = 1;

    TriangleData data;
    data.v0 = v0;
    data.J = J;
    data.J_Tinv = J_Tinv;
    data.detJ = J.Det();
    data.imap = imap;

    triangle_data.push_back(std::move(data));
  }

  //=================================== Init volumetric quadrature
  quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp=0; qp<ttl_num_vol_qpoints; ++qp)
    quadrature_point_indices.push_back(qp);

  shape_value.reserve(num_nodes);
  shape_grad.reserve(num_nodes);
  for (size_t i=0; i < num_nodes; i++)
  {
    VecDbl  node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_vol_qpoints);
    node_shape_grad.reserve(ttl_num_vol_qpoints);

    for (size_t s=0; s < num_tris; s++)
    {
      int imapped = triangle_data[s].imap[i];
      for (const auto& qpoint : qdata.qpoints)
      {
        double shape_i =
          TriangleShape(imapped,qpoint) + beta * TriangleShape(2, qpoint);

        auto gradshape_i =
          TriangleGradShape(imapped) + beta * TriangleGradShape(2);

        node_shape_value.push_back(shape_i);
        node_shape_grad.push_back(triangle_data[s].J_Tinv * gradshape_i);
      }//for qp
    } //for side

    shape_value.push_back(node_shape_value);
    shape_grad.push_back(node_shape_grad);
  }//for i

  JxW.reserve(ttl_num_vol_qpoints);
  qpoints_xyz.reserve(ttl_num_vol_qpoints);
  for (const auto& side : triangle_data)
  {
    for (size_t qp=0; qp<num_vol_qpoints; ++qp)
    {
      const auto w = qdata.weights[qp];
      JxW.push_back(side.detJ * w);

      const auto& qp_xyz_tilde = qdata.qpoints[qp];
      qpoints_xyz.push_back(side.v0 + side.J * qp_xyz_tilde);
    }//for qp
  } //for side

  return VolumeQPData(quadrature_point_indices,
                      qpoints_xyz,
                      shape_value,
                      shape_grad,
                      JxW);
}