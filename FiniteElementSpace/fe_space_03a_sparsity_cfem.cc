#include "fe_space.h"

#include "chi_log.h"
extern ChiLog& chi_log;

using namespace chi_math::finite_element;

//###################################################################
/***/
std::pair<std::vector<int64_t>, std::vector<int64_t>> SpatialDiscretization::
  BuildNodalContinuousFESparsityPattern() const
{
  const int64_t num_local_nodes = NumLocalNodes();

  typedef std::set<int64_t> SetInt64;
  typedef std::vector<SetInt64> VecSetInt64;

  VecSetInt64 native_connectivity(num_local_nodes);
  VecSetInt64 foreign_connectivity(num_local_nodes);

  for (const auto& cell : m_grid->local_cells)
  {
    const auto cell_node_indices = CellNodeIndices(cell);

    for (const auto i : cell_node_indices)
    {
      const int64_t i_map = MapNodeLocal(cell, i);

      if (i_map >= 0) //i is in diagonal-block
      {
        for (const auto j : cell_node_indices)
        {
          const int64_t j_map_local = MapNodeLocal(cell, j);
          const int64_t j_map_globl = MapNodeGlobal(cell, j);

          if (j_map_local >= 0) //j is in diagonal-block
            native_connectivity[i_map].insert(j_map_globl);
          else
            foreign_connectivity[i_map].insert(j_map_globl);

        }
      }//if local row
    }//for i
  }//for cell

  const auto ghost_cell_indices = m_grid->cells.GetGhostGlobalIDs();
  for (const uint64_t ghost_global_id : ghost_cell_indices)
  {
    const auto& cell = m_grid->cells[ghost_global_id];

    const auto cell_node_indices = CellNodeIndices(cell);

    for (const auto i : cell_node_indices)
    {
      const int64_t i_map = MapNodeLocal(cell, i);

      if (i_map >= 0) //i is in diagonal-block
      {
        for (const auto j : cell_node_indices)
        {
          const int64_t j_map_local = MapNodeLocal(cell, j);
          const int64_t j_map_globl = MapNodeGlobal(cell, j);

          if (j_map_local >= 0) //j is in diagonal-block
            native_connectivity[i_map].insert(j_map_globl);
          else
            foreign_connectivity[i_map].insert(j_map_globl);

        }
      }//if local row
    }//for i
  }

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
  BuildCFEMSparsityPattern(
    const chi_math::UnknownManager& unknown_manager) const
{
  const auto storage_type = unknown_manager.dof_storage_type;
  const auto nodal_sparsity_pattern = BuildNodalContinuousFESparsityPattern();

  const auto& nodal_nnz_in_diag = nodal_sparsity_pattern.first;
  const auto& nodal_nnz_off_diag = nodal_sparsity_pattern.second;

  const size_t num_nodal_dofs = unknown_manager.GetTotalUnknownStructureSize();
  const size_t num_local_nodes = nodal_nnz_in_diag.size();
  const int64_t num_local_dofs = NumLocalDOFs(unknown_manager);

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