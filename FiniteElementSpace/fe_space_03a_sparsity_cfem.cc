#include "fe_space.h"

#include "chi_log.h"
extern ChiLog& chi_log;

using namespace chi_math::finite_element;

//###################################################################
/***/
std::pair<std::vector<int64_t>, std::vector<int64_t>> SpatialDiscretization::
  BuildNodalCFEMSparsityPattern() const
{
  const int64_t num_local_nodes = NumLocalNodes();

  typedef std::set<int64_t> SetInt64;
  typedef std::vector<SetInt64> VecSetInt64;

  VecSetInt64 native_connectivity(num_local_nodes);
  VecSetInt64 foreign_connectivity(num_local_nodes);

  /**Lambda for inserting a cell's connectivity into the
   * master connectivity sets*/
  auto InsertCellConnectivity = [this,
                                 &native_connectivity,
                                 &foreign_connectivity]
                                 (const chi_mesh::Cell& cell)
  {
    const auto cell_node_indices = CellNodeIndices(cell);

    for (const auto i : cell_node_indices)
    {
      const int64_t node_i_local_addr = MapNodeLocal(cell, i);

      if (node_i_local_addr >= 0) //i is in diagonal-block
      {
        for (const auto j : cell_node_indices)
        {
          const int64_t node_j_local_addr = MapNodeLocal(cell, j);
          const int64_t node_j_globl_addr = MapNodeGlobal(cell, j);

          if (node_j_local_addr >= 0) //j is in diagonal-block
            native_connectivity[node_i_local_addr].insert(node_j_globl_addr);
          else
            foreign_connectivity[node_i_local_addr].insert(node_j_globl_addr);

        }//for j
      }//if local row
    }//for i
  };

  for (const auto& cell : m_grid->local_cells)
    InsertCellConnectivity(cell);

  for (const uint64_t ghost_global_id : m_grid->cells.GetGhostGlobalIDs())
    InsertCellConnectivity(m_grid->cells[ghost_global_id]);

  std::vector<int64_t> nnz_in_diagonal(num_local_nodes, 0);
  std::vector<int64_t> nnz_off_diagonal(num_local_nodes, 0);

  for (int64_t local_node_id=0; local_node_id<num_local_nodes; ++local_node_id)
  {
    nnz_in_diagonal[local_node_id] = static_cast<int64_t>(
      native_connectivity[local_node_id].size());
    nnz_off_diagonal[local_node_id] = static_cast<int64_t>(
      foreign_connectivity[local_node_id].size());
  }

  return std::make_pair(std::move(nnz_in_diagonal), std::move(nnz_off_diagonal));
}

//###################################################################
/***/
std::pair<std::vector<int64_t>, std::vector<int64_t>> SpatialDiscretization::
  BuildCFEMSparsityPattern(const chi_math::UnknownManager& unknown_manager) const
{
  const auto storage_type = unknown_manager.dof_storage_type;
  const auto [nodal_nnz_in_diag, nodal_nnz_off_diag] =
  BuildNodalCFEMSparsityPattern();

  const size_t  num_local_nodes = nodal_nnz_in_diag.size();
  const size_t  num_nodal_dofs  = unknown_manager.GetTotalUnknownStructureSize();
  const int64_t num_local_dofs  = this->NumLocalDOFs(unknown_manager);

  std::vector<int64_t> nnz_in_diagonal(num_local_dofs, 0);
  std::vector<int64_t> nnz_off_diagonal(num_local_dofs, 0);

  if (storage_type == UnknownStorageType::BLOCK)
  {
    int64_t local_dof_id = 0;
    for (int64_t nodal_dof_id=0; nodal_dof_id<num_nodal_dofs; ++nodal_dof_id)
      for (int64_t node_id=0; node_id<num_local_nodes; ++node_id)
      {
        nnz_in_diagonal[local_dof_id] = nodal_nnz_in_diag[node_id];
        nnz_off_diagonal[local_dof_id] = nodal_nnz_off_diag[node_id];
        ++local_dof_id;
      }
  }
  else if (storage_type == UnknownStorageType::NODAL)
  {
    int64_t local_dof_id = 0;
    for (int64_t node_id=0; node_id<num_local_nodes; ++node_id)
      for (int64_t nodal_dof_id=0; nodal_dof_id<num_nodal_dofs; ++nodal_dof_id)
      {
        nnz_in_diagonal[local_dof_id] = nodal_nnz_in_diag[node_id];
        nnz_off_diagonal[local_dof_id] = nodal_nnz_off_diag[node_id];
        ++local_dof_id;
      }
  }

  return std::make_pair(std::move(nnz_in_diagonal), std::move(nnz_off_diagonal));
}