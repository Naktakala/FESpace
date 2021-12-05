#ifndef FINITE_ELEMENT_SPACE_H
#define FINITE_ELEMENT_SPACE_H

#include "ChiMath/chi_math.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "fe_space_structures.h"
#include "fe_space_utils.h"
#include "Methods/FEMappingBase.h"

#include <unistd.h>

namespace chi_math::finite_element
{
  typedef std::pair<NodeInfo,uint64_t> NodeInfoIDPair;
  typedef std::vector<NodeInfoIDPair> VecNodeInfoIDPair;
  typedef std::map<uint64_t, std::vector<NodeInfo>> LocINodeInfoMap;

  //################################################################# Class def
  /**Base class for spatial discretizations.*/
  class SpatialDiscretization
  {
  private:
    std::string                      m_text_name;
  protected:
    const chi_mesh::MeshContinuumPtr m_grid;
    std::vector<FEMappingPtr>        m_local_cell_mappings;
    std::map<uint64_t, FEMappingPtr> m_ghost_cell_mappings;

    std::vector<chi_mesh::Vector3>   m_LNLL; ///< Local Node Location List
    std::vector<size_t>              m_LNR_id_2_LNLL_id_map;

    std::vector<chi_mesh::Vector3>   m_GNLL; ///< Ghost Node Location List
    std::vector<size_t>              m_GNR_id_2_GNLL_id_map;

    std::vector<int64_t>             m_LNR_local_ids;
    std::vector<int64_t>             m_LNR_global_ids;

    std::vector<int64_t>             m_GNR_global_ids;

  protected:
    explicit
    SpatialDiscretization(chi_mesh::MeshContinuumPtr& in_grid) :
      m_grid(in_grid) {}

    void AssembleNodes(
      const std::vector<NodeInfo>& LNR,
      const VecNodeInfoIDPair& GNR);

  private:
   static VecNodeInfoIDPair FilterGhostNodeRegister(
      const VecNodeInfoIDPair& GNR);

   static std::vector<NodeInfo> SubtractGNRFromLNR(
     NodeInfoListManager& LNR_manager,
     const VecNodeInfoIDPair& GNR);

   static VecNodeInfoIDPair IntersectGNRWithLNR(
     NodeInfoListManager& LNR_manager,
     const VecNodeInfoIDPair& GNR);

   static LocINodeInfoMap
    ConsolidateGNR(const VecNodeInfoIDPair& GNR);

   static LocINodeInfoMap
    SerializeAndCommunicateGNR(const LocINodeInfoMap& consolidated_GNR);

   static std::vector<int64_t> CreateGlobalIDMapForLNR(
     const std::vector<NodeInfo>& LNR);
  public:
    /**Returns the grid associated with this space.*/
    const chi_mesh::MeshContinuum& Grid() const {return *m_grid;}

    /**Gets the volume of the requested cell.*/
    double GetCellVolume(const chi_mesh::Cell& cell)
    {
      return m_local_cell_mappings[cell.local_id]->Volume();
    }

    double GetCellFaceArea(const chi_mesh::Cell& cell, size_t face_index)
    {
      return m_local_cell_mappings[cell.local_id]->FaceArea(face_index);
    }

    /**Gets the number of nodes (not DOFs) for the
     * requested cell.*/
    size_t GetCellNumNodes(const chi_mesh::Cell& cell)
    {
      return m_local_cell_mappings[cell.local_id]->NumNodes();
    }

    /**Gets the number of nodes (not DOFs) for the
     * requested cell-face pair.*/
    size_t GetFaceNumNodes(const chi_mesh::Cell& cell, size_t face)
    {
      return m_local_cell_mappings[cell.local_id]->FaceNumNodes(cell, face);
    }

    /**Returns the locations of all the nodes of the cell
     * consistent with the supplied mapping.*/
    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell)
    {
      return m_local_cell_mappings[cell.local_id]->CellNodeLocations(cell);
    }

    /**Returns the text name associated with this space.*/
    std::string TextName() const {return m_text_name;}

    /**Sets the text name to be associated with this space.*/
    void SetTextName(const std::string& text_name) {m_text_name = text_name;}

    //02
    int64_t MapNodeLocal(const chi_mesh::Cell& cell, size_t node_index) const;
    int64_t MapNodeGlobal(const chi_mesh::Cell& cell, size_t node_index) const;

    std::vector<size_t> CellNodeIndices(const chi_mesh::Cell& cell) const;

    virtual ~SpatialDiscretization() = default;
  };

  typedef std::unique_ptr<SpatialDiscretization> SpatialDiscretizationPtr;
}

namespace chi_math::finite_element
{
  //################################################################# Class def
  /**Mapping based spatial discretizations.*/
  template<class ArbMapping>
  class FiniteElementSpace : public SpatialDiscretization
  {
  public:
    /**Constructor.*/
    explicit
    FiniteElementSpace(chi_mesh::MeshContinuumPtr& in_grid) :
      SpatialDiscretization(in_grid)
    {
      std::vector<NodeInfo> local_node_register;
      m_local_cell_mappings.reserve(m_grid->local_cells.size());
      for (const auto& cell : m_grid->local_cells)
        m_local_cell_mappings.push_back(
          std::make_unique<ArbMapping>(cell, *m_grid, local_node_register));

      std::vector<std::pair<NodeInfo,uint64_t>> ghost_node_register;
      {
        const std::vector<uint64_t> ghost_cell_global_ids =
          m_grid->cells.GetGhostGlobalIDs();

        for (uint64_t ghost_global_id : ghost_cell_global_ids)
        {
          const auto& cell = m_grid->cells[ghost_global_id];
          std::vector<NodeInfo> ghost_node_list;
          m_ghost_cell_mappings[ghost_global_id] =
            std::make_unique<ArbMapping>(cell, *m_grid,ghost_node_list);
          for (auto& node : ghost_node_list)
            ghost_node_register.emplace_back(node, cell.partition_id);
        }
      }

      AssembleNodes(local_node_register, ghost_node_register);
    }

    static SpatialDiscretizationPtr New(chi_mesh::MeshContinuumPtr& in_grid)
    {
      auto new_space = std::make_unique<FiniteElementSpace<ArbMapping>>(in_grid);
      return std::move(new_space);
    }
  };

}//namespace chi_math

#endif //FINITE_ELEMENT_SPACE_H