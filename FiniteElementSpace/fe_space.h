#ifndef FINITE_ELEMENT_SPACE_H
#define FINITE_ELEMENT_SPACE_H

#include "ChiMath/chi_math.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "Methods/FEMethodBase.h"

namespace chi_math
{
  void FESpaceTest();
}//namespace chi_math

namespace chi_math::finite_element
{
  enum class FESpaceFlags : unsigned int
  {
    NO_FLAGS_SET                          = 0,
    PREBUILD_CELL_MAPPINGS                = (1 << 0),
    PRECOMPUTE_IntV_gradShapeI_gradShapeJ = (1 << 1),
    PRECOMPUTE_IntV_shapeI_gradShapeJ     = (1 << 2),
    PRECOMPUTE_IntV_shapeI_shapeJ         = (1 << 3),
    PRECOMPUTE_IntV_shapeI                = (1 << 4),
    PRECOMPUTE_IntV_gradShapeI            = (1 << 5),

    PRECOMPUTE_IntS_shapeI_shapeJ         = (1 << 6),
    PRECOMPUTE_IntS_shapeI                = (1 << 7),
    PRECOMPUTE_IntS_shapeI_gradshapeJ     = (1 << 8),

    PRECOMPUTE_FACE_NODE_MAPPINGS         = (1 << 9)
  };

  inline FESpaceFlags
  operator|(const FESpaceFlags f1, const FESpaceFlags f2)
  {
    return static_cast<FESpaceFlags>(static_cast<unsigned int>(f1) |
                                     static_cast<unsigned int>(f2));
  }
}//namespace chi_math

namespace chi_math::finite_element
{
  template<class ArbMethod>
  class FiniteElementSpace
  {
  private:
    std::string                      m_text_name;
  protected:
    const chi_mesh::MeshContinuumPtr m_grid;
    FiniteElementMethod&             m_fe_method;
    ArbMethod                        m_arb_fe_method;

    using QuadraturePtr = std::unique_ptr<chi_math::Quadrature>;
    std::map<chi_mesh::CellType, QuadraturePtr> m_celltype_to_quadrature_map;

  public:
    /**Constructor.*/
    explicit
    FiniteElementSpace(chi_mesh::MeshContinuumPtr& in_grid) :
      m_grid(in_grid),
      m_fe_method(m_arb_fe_method)
    {
      m_celltype_to_quadrature_map =
        m_fe_method.MapRequiredQuadratures(*m_grid, QuadratureOrder::SECOND);
    }

    FEMappingPtr GetCellMapping(const chi_mesh::Cell& cell)
    {
      return m_fe_method.MakeCellMapping(cell, *m_grid);
    }

    /**Returns the grid associated with this space.*/
    const chi_mesh::MeshContinuum& Grid() const {return *m_grid;}

    /**Gets the number of nodes (not DOFs) for the
     * requested cell.*/
    size_t GetCellNumNodes(const chi_mesh::Cell& cell)
    {return m_fe_method.GetCellNumNodes(cell);}

    /**Returns the locations of all the nodes of the cell
     * consistent with the supplied mapping.*/
    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell)
    {return m_fe_method.GetCellNodeLocations(cell, *m_grid);}

    /**Returns the text name associated with this space.*/
    std::string TextName() const {return m_text_name;}

    /**Sets the text name to be associated with this space.*/
    void SetTextName(const std::string& text_name) {m_text_name = text_name;}
  };

}//namespace chi_math

#endif //FINITE_ELEMENT_SPACE_H