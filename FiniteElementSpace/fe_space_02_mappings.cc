#include "fe_space.h"

#include "chi_mpi.h"

using namespace chi_math::finite_element;

int64_t SpatialDiscretization::
  MapNodeLocal(const chi_mesh::Cell &cell, size_t node_index) const
{
  const auto& chi_mpi = ChiMPI::GetInstance();
  const std::string fname = __FUNCTION__;

  if (cell.partition_id != chi_mpi.location_id/*home*/)
    throw std::invalid_argument(fname + ": Function cannot be used on "
                                        "non-local cells.");

  const auto& cell_mapping = m_local_cell_mappings[cell.local_id];

  uint64_t register_index = cell_mapping->MapNodeRegister(node_index);

  return m_LNR_local_ids.at(register_index);
}

int64_t SpatialDiscretization::
  MapNodeGlobal(const chi_mesh::Cell &cell, size_t node_index) const
{
  const auto& chi_mpi = ChiMPI::GetInstance();
  const std::string fname = __FUNCTION__;

  if (cell.partition_id != chi_mpi.location_id/*home*/)
  {
    const auto& cell_mapping = m_local_cell_mappings[cell.local_id];

    uint64_t register_index = cell_mapping->MapNodeRegister(node_index);

    return m_LNR_global_ids.at(register_index);
  }
  else
  {
    try
    {
      const auto& cell_mapping = m_ghost_cell_mappings.at(cell.global_id);

      uint64_t register_index = cell_mapping->MapNodeRegister(node_index);

      return m_GNR_global_ids.at(register_index);
    }
    catch (const std::out_of_range& oor)
    {
      throw std::out_of_range(fname + ": Ghost cell mapping not found for "
                                      "given cell.");
    }
  }
}

std::vector<size_t> SpatialDiscretization::
  CellNodeIndices(const chi_mesh::Cell &cell) const
{
  const auto& chi_mpi = ChiMPI::GetInstance();

  size_t num_nodes = 0;
  if (cell.partition_id != chi_mpi.location_id/*home*/)
  {
    const auto& cell_mapping = m_local_cell_mappings[cell.local_id];

    num_nodes = cell_mapping->NumNodes();
  }
  else
  {
    const auto& cell_mapping = m_ghost_cell_mappings.at(cell.global_id);

    num_nodes = cell_mapping->NumNodes();
  }

  std::vector<size_t> node_indices(num_nodes, 0);
  for (size_t i=0; i<num_nodes; ++i)
    node_indices[i] = i;

  return node_indices;
}