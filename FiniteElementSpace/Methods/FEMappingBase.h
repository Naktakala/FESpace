#ifndef FINITE_ELEMENT_MAPPING_BASE_H
#define FINITE_ELEMENT_MAPPING_BASE_H

#include "FENodeInfo.h"

#include <map>

#include "fe_space_qp_data.h"

#include "ChiMath/Quadratures/quadrature.h"

#include "ChiMesh/Cell/cell.h"

namespace chi_math::finite_element
{

/**Base class for any finite element mapping*/
class FiniteElementMapping
{
protected:
  const chi_mesh::Cell&          m_cell;
  const chi_mesh::MeshContinuum& m_grid;

private: //Mandatory via mutator
  size_t                         m_num_nodes=0;
  std::vector<uint64_t>          m_local_node_register;

protected: //Mandatory to compute
  double                         m_volume=0.0;
  std::vector<double>            m_face_areas;

protected:
  explicit FiniteElementMapping(const chi_mesh::Cell& cell,
                                const chi_mesh::MeshContinuum& grid) :
                                m_cell(cell),
                                m_grid(grid)
                                {}

  size_t SetNumNodesAndLocalRegister(const size_t num_nodes,
                                     const size_t node_register_size)
  {
    m_num_nodes = num_nodes;
    m_local_node_register.assign(m_num_nodes, 0);

    for (uint64_t i=0; i<m_num_nodes; ++i)
      m_local_node_register[i] = node_register_size + i;

    return node_register_size + num_nodes;
  }

public: //Virtual functions - needing definitions
  virtual QuadratureOrder GetMinimumQuadratureOrder() const
  {
    return QuadratureOrder::CONSTANT;
  }

  typedef std::pair<chi_mesh::CellType,QuadratureOrder> QuadratureKey;
  typedef std::unique_ptr<Quadrature> QuadraturePtr;
  /**Adds the appropriate quadratures for volumetric and surface integrations,
   * for the referenced cell, to the quadrature map. The required quadrature
   * is only added to the map if not already there.*/
  virtual void AddRequiredQuadratures(
    chi_math::QuadratureOrder order,
    std::map<QuadratureKey, QuadraturePtr>& quadrature_stack) const
  {/*Does nothing by default*/};

  /**Returns the number of nodes associated with a given cell-face pair.*/
  virtual
  size_t FaceNumNodes(size_t face_index) const = 0;

  /**Maps a face node to a cell node.*/
  virtual
  size_t MapFaceNodeToCellNode(size_t face_index,
                               size_t face_node_index) const = 0;

  virtual
  VolumeQPData BuildVolumetricQPData(
    chi_math::QuadratureOrder order,
    std::map<QuadratureKey, QuadraturePtr>& quadrature_stack) const = 0;

public: //Virtual functions - not needing definitions

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

public: //Non-virtual functions
  /**Returns the number of nodes associated with a given cell.*/
  size_t NumNodes() const {return m_num_nodes;}

  /**Returns the volume of the cell associated with this mapping.*/
  double Volume() const {return m_volume;}

  /**Returns the area of a specific face of the cell associated
   * with this mapping.*/
  double FaceArea(size_t face_index) const {return m_face_areas.at(face_index);}

  /** Destructor. */
  virtual ~FiniteElementMapping() = default;

  const chi_mesh::Cell& ReferenceCell() const {return m_cell;}

public:
  uint64_t MapNodeRegister(size_t node_index) const
  {
    return m_local_node_register.at(node_index);
  }
};

typedef std::unique_ptr<FiniteElementMapping> FEMappingPtr;



}//namespace chi_math

#endif //FINITE_ELEMENT_MAPPING_BASE_H