#include "pwl_mapping.h"

using namespace chi_math::finite_element;

double PiecewiseLinear::SlabShape(int i, const chi_mesh::Vector3& qpoint)
{
  if (i < 0) return 0.0;

  if (qpoint.x < 0.0) return 0.0;
  if (qpoint.x > 1.0) return 0.0;

  double value = 0.0;
  if (i == 0) value = 1.0 - qpoint.x;
  if (i == 1) value = qpoint.x;

  return value;
}

chi_mesh::Vector3 PiecewiseLinear::SlabGradShape(int i)
{
  if (i < 0) return chi_mesh::Vector3();

  double value = 0.0;
  if (i == 0) value = -1.0;
  if (i == 1) value =  1.0;

  return chi_mesh::Vector3(0.0,0.0,value);
}

double PiecewiseLinear::TriangleShape(int i, const chi_mesh::Vector3& qpoint)
{
  if (i < 0) return 0.0;

  double value = 0.0;
  if (i == 0) value = 1 - qpoint.x - qpoint.y;
  if (i == 1) value = qpoint.x;
  if (i == 2) value = qpoint.y;

  return value;
}

chi_mesh::Vector3 PiecewiseLinear::TriangleGradShape(int i)
{
  if (i < 0) return chi_mesh::Vector3();

  chi_mesh::Vector3 grad_shape;

  if (i == 0) grad_shape = chi_mesh::Vector3(-1.0,-1.0,0.0);
  if (i == 1) grad_shape = chi_mesh::Vector3(-1.0, 0.0,0.0);
  if (i == 2) grad_shape = chi_mesh::Vector3( 0.0,-1.0,0.0);

  return grad_shape;
}

double PiecewiseLinear::
  TetrahedronShape(int i, const chi_mesh::Vector3 &qpoint)
{
  if (i < 0) return 0.0;

  double value = 0.0;
  if (i == 0) value = 1.0 - qpoint.x - qpoint.y - qpoint.z;
  if (i == 1) value = qpoint.x;
  if (i == 2) value = qpoint.y;
  if (i == 3) value = qpoint.z;

  return value;
}

chi_mesh::Vector3 PiecewiseLinear::TetrahedronGradShape(int i)
{
  if (i < 0) return chi_mesh::Vector3();

  if (i == 0) return chi_mesh::Vector3(-1.0,-1.0,-1.0);
  if (i == 1) return chi_mesh::Vector3( 1.0, 0.0, 0.0);
  if (i == 2) return chi_mesh::Vector3( 0.0, 1.0, 0.0);
  if (i == 3) return chi_mesh::Vector3( 0.0, 0.0, 1.0);

  return chi_mesh::Vector3(0.0,0.0,0.0);
}


