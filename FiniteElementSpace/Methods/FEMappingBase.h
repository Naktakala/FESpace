#ifndef FINITE_ELEMENT_MAPPING_BASE_H
#define FINITE_ELEMENT_MAPPING_BASE_H

#include <vector>

#include "ChiMesh/chi_mesh.h"

namespace chi_math::finite_element
{

/**Base class for any finite element mapping*/
class FiniteElementMapping
{
protected:
  size_t              m_num_nodes=0;
  double              m_volume=0.0;
  std::vector<double> m_face_areas;

  FiniteElementMapping() {}

public:
  /**Returns the number of nodes associated with a given cell.*/
  virtual
  size_t GetCellNumNodes(const chi_mesh::Cell& cell) const = 0;

  /**Returns the locations of the nodes associated with a given cell.*/
  virtual
  std::vector<chi_mesh::Vector3>
  GetCellNodeLocations(const chi_mesh::Cell& cell,
                       const chi_mesh::MeshContinuum& grid) const = 0;

  /**Returns the volume of the cell associated with this mapping.*/
  double Volume() const {return m_volume;}

  /**Returns the area of a specific face of the cell associated
   * with this mapping.*/
  double FaceArea(size_t face_index) const {return m_face_areas.at(face_index);}

  /** Virtual function evaluation of the shape function. */
  virtual double ShapeValue(const int i, const chi_mesh::Vector3& xyz) const
  {
    return 0.0;
  }

  /** Virtual function returning the all the shape function evaluations
   * at the point.*/
  virtual void ShapeValues(const chi_mesh::Vector3& xyz,
                           std::vector<double>& shape_values) const
  {
    shape_values.resize(m_num_nodes, 0.0);
  }

  /** Virtual function evaluation of the grad-shape function. */
  virtual chi_mesh::Vector3 GradShapeValue(const int i,
                                           const chi_mesh::Vector3& xyz) const
  {
    return chi_mesh::Vector3(0.0, 0.0, 0.0);
  }

  /** Virtual function evaluation of the grad-shape function. */
  virtual void
    GradShapeValues(const chi_mesh::Vector3& xyz,
                    std::vector<chi_mesh::Vector3>& gradshape_values) const
  {
    gradshape_values.resize(m_num_nodes, chi_mesh::Vector3());
  }

  /** Destructor. */
  virtual ~FiniteElementMapping() = default;
};

typedef std::unique_ptr<FiniteElementMapping> FEMappingPtr;



}//namespace chi_math

#endif //FINITE_ELEMENT_MAPPING_BASE_H