#include "pwl_mapping.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"

using namespace chi_math::finite_element;

VolumeQPData PiecewiseLinear::
  BuildVolumetricQPDataPolyhedron(chi_math::QuadratureOrder order) const
{
  std::vector<unsigned int>     quadrature_point_indices; ///< qp index only
  VecVec3                       qpoints_xyz             ; ///< qp index only
  std::vector<VecDbl>           shape_value             ; ///< Node i, then qp
  std::vector<VecVec3>          shape_grad              ; ///< Node i, then qp
  VecDbl                        JxW                     ; ///< qp index only

  typedef std::vector<int> VecInt;
  chi_math::QuadratureTetrahedron qdata(order);

  //=================================== Determine number of internal qpoints
  size_t num_tets=0;
  for (auto& face : m_cell.faces)
    for (auto& side : face.vertex_ids)
      ++num_tets;

  size_t num_vol_qpoints = qdata.qpoints.size();
  size_t ttl_num_vol_qpoints = num_tets * num_vol_qpoints;
  size_t num_nodes = m_cell.vertex_ids.size();
  size_t num_faces = m_cell.faces.size();

  //=================================== Determine alpha and betaa
  double alpha = 1.0/static_cast<double>(m_cell.vertex_ids.size());
  std::vector<double> beta;
  beta.reserve(num_faces);
  for (auto& face : m_cell.faces)
    beta.push_back(1.0/static_cast<double>(face.vertex_ids.size()));

  //=================================== Determine node to face influence map
  std::vector<VecInt> ifmap(num_nodes, VecInt(num_faces, -1));
  for (size_t i=0; i<num_nodes; ++i)
    for (size_t f=0; f<num_faces; ++f)
      for (uint64_t fvid : m_cell.faces[f].vertex_ids)
        if (m_cell.vertex_ids[i] == fvid)
          ifmap[i][f] = 2;

  //=================================== Determine tetrahedron data
  struct TetrahedronData
  {
    chi_mesh::Vector3 v0;
    chi_mesh::Matrix3x3 J;
    chi_mesh::Matrix3x3 J_Tinv;
    double detJ=0.0;

    std::vector<int> imap;
  };
  std::vector<std::vector<TetrahedronData>> tet_face_side_data;

  for (auto& face : m_cell.faces)
  {
    const size_t num_tris = face.vertex_ids.size();
    std::vector<TetrahedronData> tet_side_data;
    tet_side_data.reserve(num_tris);
    for (size_t s=0; s<num_tris; ++s)
    {
      size_t sp1 = (s < (num_tris-1))? s+1 : 0;
      const auto& v0 = m_grid.vertices[face.vertex_ids[  s]];
      const auto& v1 = m_grid.vertices[face.vertex_ids[sp1]];
      const auto& v2 = face.centroid;
      const auto& v3 = m_cell.centroid;

      const auto v01 = v1-v0;
      const auto v02 = v2-v0;
      const auto v03 = v3-v0;

      chi_mesh::Matrix3x3 J;
      J.SetColJVec(0, v01);
      J.SetColJVec(0, v02);
      J.SetColJVec(0, v03);

      const auto J_Tinv = J.Transpose().Inverse();

      std::vector<int> imap(num_nodes, -1);
      imap[  s] = 0;
      imap[sp1] = 1;

      TetrahedronData data;
      data.v0 = v0;
      data.J = J;
      data.J_Tinv = J_Tinv;
      data.detJ = J.Det();
      data.imap = imap;

      tet_side_data.push_back(std::move(data));
    }//for side

    tet_face_side_data.push_back(std::move(tet_side_data));
  }//for face

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

    for (size_t f=0; f < num_faces; f++)
    {
      auto& face = m_cell.faces[f];
      const size_t num_sides = face.vertex_ids.size();
      const double beta_f = beta[f];
      for (size_t s=0; s < num_sides; s++)
      {
        const auto& tet_data = tet_face_side_data[f][s];
        const int imapped  = tet_data.imap[i];
        const int ifmapped = ifmap[i][f];

        for (const auto& qpoint : qdata.qpoints)
        {
          double shape_i =
            alpha  * TetrahedronShape(3, qpoint) +
            beta_f * TetrahedronShape(ifmapped, qpoint) +
                      TetrahedronShape(imapped, qpoint);

          auto gradshape_i =
            alpha  * TetrahedronGradShape(3) +
            beta_f * TetrahedronGradShape(ifmapped) +
                     TetrahedronGradShape(imapped);

          node_shape_value.push_back(shape_i);
          node_shape_grad.emplace_back(tet_data.J_Tinv * gradshape_i);  //z
        }//for qp
      } //for side
    } //for face

    shape_value.push_back(node_shape_value);
    shape_grad.push_back(node_shape_grad);
  }//for i

  JxW.reserve(ttl_num_vol_qpoints);
  qpoints_xyz.reserve(ttl_num_vol_qpoints);
  for (const auto& tet_face_data : tet_face_side_data)
  {
    for (const auto& side_tet : tet_face_data)
    {
      for (size_t qp=0; qp<num_vol_qpoints; ++qp)
      {
        const auto w = qdata.weights[qp];
        JxW.push_back(side_tet.detJ * w);

        const auto& qp_xyz_tilde = qdata.qpoints[qp];
        qpoints_xyz.push_back(side_tet.v0 + side_tet.J * qp_xyz_tilde);
      }//for qp
    } //for side
  } //for face

  return VolumeQPData(quadrature_point_indices,
                      qpoints_xyz,
                      shape_value,
                      shape_grad,
                      JxW);
}