#ifndef FINITE_ELEMENT_SPACE_H
#define FINITE_ELEMENT_SPACE_H

#include "ChiMath/chi_math.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "fe_space_structures.h"
#include "Methods/FEMappingBase.h"

#include <unistd.h>

namespace chi_math::finite_element
{
  //################################################################# Class def
  /**Base class for spatial discretizations.*/
  class SpatialDiscretization
  {
  private:
    std::string                      m_text_name;
  protected:
    const chi_mesh::MeshContinuumPtr m_grid;
    std::vector<FEMappingPtr>        m_local_cell_mappings;
    std::vector<FEMappingPtr>        m_ghost_cell_mappings;

    std::vector<chi_mesh::Vector3>   m_local_unique_node_locations;
    std::vector<size_t>              m_node_register_2_local_unique_map;

    std::vector<int64_t>             m_node_register_local_ids;
    std::vector<int64_t>             m_node_register_global_ids;

  protected:
    explicit
    SpatialDiscretization(chi_mesh::MeshContinuumPtr& in_grid) :
      m_grid(in_grid) {}

    void AssembleNodes(
      const std::vector<NodeInfo>& node_register_node_info,
      const std::vector<std::pair<NodeInfo,uint64_t>>& ghost_scope_nodes_pids);

  public:
    virtual FEMappingPtr GetCellMapping(const chi_mesh::Cell& cell) const = 0;

    virtual void OrderNodes()
    {
      OrderNodesContinuous();
    }
    void OrderNodesContinuous();
    void OrderNodesDisContinuous();

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
//      MPI_Barrier(MPI_COMM_WORLD);
//      std::cout << "Pause before cell mapping";
//      MPI_Barrier(MPI_COMM_WORLD);
//      usleep(10000000);

      MPI_Barrier(MPI_COMM_WORLD);
      std::cout << "Creating cell mappings\n";
      MPI_Barrier(MPI_COMM_WORLD);

      std::vector<NodeInfo> node_register_node_info;
      m_local_cell_mappings.reserve(m_grid->local_cells.size());
      for (const auto& cell : m_grid->local_cells)
        m_local_cell_mappings.push_back(
          std::make_unique<ArbMapping>(cell, *m_grid, node_register_node_info));

      MPI_Barrier(MPI_COMM_WORLD);
      std::cout << "Creating ghost cell mappings\n";
      MPI_Barrier(MPI_COMM_WORLD);

      std::vector<std::pair<NodeInfo,uint64_t>> ghost_scope_nodes_pids;
      {
        const std::vector<uint64_t> ghost_cell_global_ids =
          m_grid->cells.GetGhostGlobalIDs();

        for (uint64_t ghost_id : ghost_cell_global_ids)
        {
          const auto& cell = m_grid->cells[ghost_id];
          std::vector<NodeInfo> ghost_node_list;
          m_ghost_cell_mappings.push_back(
            std::make_unique<ArbMapping>(cell, *m_grid,ghost_node_list));
          for (auto& node : ghost_node_list)
            ghost_scope_nodes_pids.emplace_back(node, cell.partition_id);
        }
      }

//      MPI_Barrier(MPI_COMM_WORLD);
//      std::cout << "Pause before nodes assembly";
//      MPI_Barrier(MPI_COMM_WORLD);
//      usleep(10000000);

      AssembleNodes(node_register_node_info, ghost_scope_nodes_pids);

//      MPI_Barrier(MPI_COMM_WORLD);
//      std::cout << "Pause after nodes assembly";
//      MPI_Barrier(MPI_COMM_WORLD);
//      usleep(10000000);
    }

    static SpatialDiscretizationPtr New(chi_mesh::MeshContinuumPtr& in_grid)
    {
      auto new_space = std::make_unique<FiniteElementSpace<ArbMapping>>(in_grid);
      return std::move(new_space);
    }

    FEMappingPtr GetCellMapping(const chi_mesh::Cell& cell) const override
    {
      size_t temp_sizer=0;
      auto arb_mapping = std::make_unique<ArbMapping>(cell, *m_grid, temp_sizer);
      return std::move(arb_mapping);
    }
  };

}//namespace chi_math

#endif //FINITE_ELEMENT_SPACE_H