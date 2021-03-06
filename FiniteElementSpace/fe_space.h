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

    std::vector<int64_t>             m_GNR_local_ids;
    std::vector<int64_t>             m_GNR_global_ids;

    std::vector<int64_t>             m_local_nodes_global_ids;
    std::vector<int64_t>             m_ghost_nodes_global_ids;

    int64_t                          m_num_local_nodes = 0;
    int64_t                          m_num_global_nodes = 0;

    typedef std::pair<chi_mesh::CellType, QuadratureOrder> QuadratureKey;
    typedef std::unique_ptr<Quadrature> QuadraturePtr;
    std::map<QuadratureKey, QuadraturePtr> m_quadrature_stack;

  protected:
    explicit
    SpatialDiscretization(chi_mesh::MeshContinuumPtr& in_grid) :
      m_grid(in_grid) {}

    const FiniteElementMapping& GetCellMapping(const chi_mesh::Cell& cell) const;

    const std::vector<int64_t>&
      GetCellRelevantLocalIDRegister(const chi_mesh::Cell& cell) const;

      const std::vector<int64_t>&
      GetCellRelevantGlobalIDRegister(const chi_mesh::Cell& cell) const;

  public:
    /**Returns the text name associated with this space.*/
    std::string TextName() const {return m_text_name;}

    /**Sets the text name to be associated with this space.*/
    void SetTextName(const std::string& text_name) {m_text_name = text_name;}

  protected:
  //01
    void AssembleNodes(
      const std::vector<NodeInfo>& LNR,
      const VecNodeInfoIDPair& GNR);

  //01a
  private:
   static VecNodeInfoIDPair FilterGhostNodeRegister(
      const VecNodeInfoIDPair& GNR);

   static std::vector<NodeInfo> SubtractGNRFromLNR(
     const NodeInfoListManager& LNR_manager,
     const VecNodeInfoIDPair& GNR);

   static VecNodeInfoIDPair IntersectGNRWithLNR(
     const NodeInfoListManager& LNR_manager,
     const VecNodeInfoIDPair& GNR);

   static LocINodeInfoMap
    ConsolidateGNR(const VecNodeInfoIDPair& GNR);

   static LocINodeInfoMap
    SerializeAndCommunicateGNR(const LocINodeInfoMap& consolidated_GNR);

   static std::vector<int64_t> CreateGlobalIDMapForLNR(
     const std::vector<NodeInfo>& LNR,
     int64_t& num_global_nodes);

   static VecNodeInfoIDPair SubtractLNRFromGNR(
     const NodeListFindManager& LNR_manager,
     const VecNodeInfoIDPair& GNR);

  public:
    /**Returns the grid associated with this space.*/
    const chi_mesh::MeshContinuum& Grid() const {return *m_grid;}

    //02 utils
    /**Gets the volume of the requested cell.*/
    double GetCellVolume(const chi_mesh::Cell& cell);

    double GetCellFaceArea(const chi_mesh::Cell& cell, size_t face_index);

    /**Gets the number of nodes (not DOFs) for the
     * requested cell.*/
    size_t GetCellNumNodes(const chi_mesh::Cell& cell);

    /**Gets the number of nodes (not DOFs) for the
     * requested cell-face pair.*/
    size_t GetFaceNumNodes(const chi_mesh::Cell& cell, size_t face);

    /**Returns the locations of all the nodes of the cell
     * consistent with the supplied mapping.*/
    std::vector<chi_mesh::Vector3>
      GetCellNodeLocations(const chi_mesh::Cell& cell);

    std::vector<size_t> CellNodeIndices(const chi_mesh::Cell& cell) const;

    //03
    std::pair<std::vector<int64_t>, std::vector<int64_t>>
      BuildNodalCFEMSparsityPattern() const;

    std::pair<std::vector<int64_t>, std::vector<int64_t>>
      BuildCFEMSparsityPattern(
        const chi_math::UnknownManager& unknown_manager) const;

//    std::pair<std::vector<int64_t>, std::vector<int64_t>>
//      BuildNodalDFEMSparsityPattern() const;
//
//    std::pair<std::vector<int64_t>, std::vector<int64_t>>
//      BuildDFEMSparsityPattern(
//        const chi_math::UnknownManager& unknown_manager) const;
//
//    std::pair<std::vector<int64_t>, std::vector<int64_t>>
//      BuildNodalFVSparsityPattern() const;

    //04
    int64_t NumLocalNodes() const;
    int64_t NumGlobalNodes() const;

    int64_t NumLocalDOFs(const chi_math::UnknownManager& unknown_manager) const;
    int64_t NumGlobalDOFs(const chi_math::UnknownManager& unknown_manager) const;

    int64_t MapNodeLocal(const chi_mesh::Cell& cell, size_t node_index) const;
    int64_t MapNodeGlobal(const chi_mesh::Cell& cell, size_t node_index) const;

    std::vector<int64_t> MapNodesLocal(const chi_mesh::Cell& cell) const;
    std::vector<int64_t> MapNodesGlobal(const chi_mesh::Cell& cell) const;

    int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                        size_t node_index,
                        const chi_math::UnknownManager& unknown_manager,
                        unsigned int unknown_id,
                        unsigned int component) const;

    int64_t MapDOFGlobal(const chi_mesh::Cell& cell,
                         size_t node_index,
                         const chi_math::UnknownManager& unknown_manager,
                         unsigned int unknown_id,
                         unsigned int component) const;

    std::vector<int64_t>
      MapDOFsLocal(const chi_mesh::Cell& cell,
                   const chi_math::UnknownManager& unknown_manager,
                   unsigned int unknown_id,
                   unsigned int component) const;
    std::vector<int64_t>
      MapDOFsGlobal(const chi_mesh::Cell& cell,
                    const chi_math::UnknownManager& unknown_manager,
                    unsigned int unknown_id,
                    unsigned int component) const;

    std::vector<int64_t>
      GetGhostNodesGlobalIDs() const {return m_ghost_nodes_global_ids;}

    std::vector<int64_t>
      GetGhostDOFsGlobalIDs(
        const chi_math::UnknownManager& unknown_manager) const;

    virtual ~SpatialDiscretization() = default;
  };

  typedef std::shared_ptr<SpatialDiscretization> SpatialDiscretizationPtr;
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
      //=============================== Make mappings for local cells
      std::vector<NodeInfo> local_node_register;
      m_local_cell_mappings.reserve(m_grid->local_cells.size());
      for (const auto& cell : m_grid->local_cells)
        m_local_cell_mappings.push_back(
          std::make_unique<ArbMapping>(cell, *m_grid, local_node_register));

      //=============================== Make mappings for ghost cells
      std::vector<NodeInfo> ghost_node_list;
      std::vector<std::pair<NodeInfo,uint64_t>> ghost_node_register;
      {
        const std::vector<uint64_t> ghost_cell_global_ids =
          m_grid->cells.GetGhostGlobalIDs();

        for (uint64_t ghost_global_id : ghost_cell_global_ids)
        {
          const auto& cell = m_grid->cells[ghost_global_id];

          const size_t first_node_location = ghost_node_list.size();
          m_ghost_cell_mappings[ghost_global_id] =
            std::make_unique<ArbMapping>(cell, *m_grid,ghost_node_list);
          const size_t last_node_location = ghost_node_list.size()-1;

          for (size_t n=first_node_location; n<=last_node_location; ++n)
            ghost_node_register.emplace_back(ghost_node_list[n],
                                             cell.partition_id);
        }
      }

      //=============================== Assemble nodes
      AssembleNodes(local_node_register, ghost_node_register);

      //=============================== Populate quadratures
      for (const auto& cell_mapping : m_local_cell_mappings)
        cell_mapping->AddRequiredQuadratures(
          cell_mapping->GetMinimumQuadratureOrder(),
          m_quadrature_stack);

      for (const auto& [cell_globl_id, cell_mapping] : m_ghost_cell_mappings)
        cell_mapping->AddRequiredQuadratures(
          cell_mapping->GetMinimumQuadratureOrder(),
          m_quadrature_stack);
    }

    static SpatialDiscretizationPtr New(chi_mesh::MeshContinuumPtr& in_grid)
    {
      auto new_space = std::make_unique<FiniteElementSpace<ArbMapping>>(in_grid);
      return std::move(new_space);
    }
  };

}//namespace chi_math

#endif //FINITE_ELEMENT_SPACE_H