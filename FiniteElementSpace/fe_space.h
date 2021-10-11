#ifndef FINITE_ELEMENT_SPACE_H
#define FINITE_ELEMENT_SPACE_H

#include "ChiMath/chi_math.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "fe_space_flags.h"
#include "Methods/FEMappingBase.h"

namespace chi_math::finite_element
{
  template<class ArbMapping>
  class FiniteElementSpace
  {
  private:
    std::string                      m_text_name;
  protected:
    const chi_mesh::MeshContinuumPtr m_grid;
    std::vector<ArbMapping> m_cell_mappings;

    using QuadraturePtr = std::unique_ptr<chi_math::Quadrature>;
    std::map<chi_mesh::CellType, QuadraturePtr> m_celltype_to_quadrature_map;

  public:
    /**Constructor.*/
    explicit
    FiniteElementSpace(chi_mesh::MeshContinuumPtr& in_grid) :
      m_grid(in_grid)
    {
      m_cell_mappings.reserve(m_grid->local_cells.size());
      for (const auto& cell : m_grid->local_cells)
        m_cell_mappings.push_back(ArbMapping(cell, *m_grid));
    }

    FEMappingPtr GetCellMapping(const chi_mesh::Cell& cell)
    {
      auto arb_mapping = std::make_unique<ArbMapping>(cell, *m_grid);
      return std::static_pointer_cast<FiniteElementMapping>(arb_mapping);
    }

    /**Returns the grid associated with this space.*/
    const chi_mesh::MeshContinuum& Grid() const {return *m_grid;}

    /**Gets the number of nodes (not DOFs) for the
     * requested cell.*/
    size_t GetCellNumNodes(const chi_mesh::Cell& cell)
    {return m_cell_mappings[cell.local_id].CellNumNodes(cell);}

    /**Returns the locations of all the nodes of the cell
     * consistent with the supplied mapping.*/
    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell)
    {return m_cell_mappings[cell.local_id].CellNodeLocations(cell);}

    /**Returns the text name associated with this space.*/
    std::string TextName() const {return m_text_name;}

    /**Sets the text name to be associated with this space.*/
    void SetTextName(const std::string& text_name) {m_text_name = text_name;}
  };

}//namespace chi_math

#endif //FINITE_ELEMENT_SPACE_H