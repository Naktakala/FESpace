#include "fe_space.h"

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

using namespace chi_math::finite_element;

//###################################################################
int64_t SpatialDiscretization::
  MapNodeLocal(const chi_mesh::Cell &cell, size_t node_index) const
{
  const std::string fname = __FUNCTION__;

  if (cell.partition_id != chi_mpi.location_id/*home*/)
    throw std::invalid_argument(fname + ": Function cannot be used on "
                                        "non-local cells.");

  const auto& cell_mapping = *m_local_cell_mappings[cell.local_id];

  uint64_t register_index = cell_mapping.MapNodeRegister(node_index);

  return m_LNR_local_ids.at(register_index);
}

//###################################################################
int64_t SpatialDiscretization::
  MapNodeGlobal(const chi_mesh::Cell &cell, size_t node_index) const
{
  const std::string fname = __FUNCTION__;

  const auto& cell_mapping = GetCellMapping(cell);

  const uint64_t register_index = cell_mapping.MapNodeRegister(node_index);

  if (cell.partition_id == chi_mpi.location_id/*home*/)
    return m_LNR_global_ids.at(register_index);
  else
  {
    try { return m_GNR_global_ids.at(register_index); }
    catch (const std::out_of_range& oor)
    {
      throw std::out_of_range(fname + ": Ghost cell mapping not found for "
                                      "given cell.");
    }
  }
}

//###################################################################
int64_t SpatialDiscretization::
  NumLocalDOFs(const chi_math::UnknownManager &unknown_manager) const
{
  const size_t num_comps = unknown_manager.GetTotalUnknownStructureSize();

  return m_num_local_nodes * static_cast<int64_t>(num_comps);
}

//###################################################################
int64_t SpatialDiscretization::
  NumGlobalDOFs(const chi_math::UnknownManager &unknown_manager) const
{
  const size_t num_comps = unknown_manager.GetTotalUnknownStructureSize();

  return m_num_global_nodes * static_cast<int64_t>(num_comps);
}

//###################################################################
int64_t SpatialDiscretization::
  MapDOFLocal(const chi_mesh::Cell &cell, size_t node_index,
              const chi_math::UnknownManager &unknown_manager,
              const unsigned int unknown_id,
              const unsigned int component) const
{
  const std::string fname = __FUNCTION__;

  if (cell.partition_id != chi_mpi.location_id/*home*/)
    throw std::invalid_argument(fname + ": Function cannot be used on "
                                        "non-local cells.");

  const auto& cell_mapping = GetCellMapping(cell);
  const uint64_t register_index = cell_mapping.MapNodeRegister(node_index);

  const int64_t node_address = m_LNR_local_ids[register_index];

  const auto& uk_man   = unknown_manager;
  const auto storage_type = unknown_manager.dof_storage_type;
  const size_t num_nodal_dofs_ = uk_man.GetTotalUnknownStructureSize();
  const size_t nodal_dof_id_ = uk_man.MapUnknown(unknown_id, component);

  const int64_t num_nodal_dofs = static_cast<int64_t>(num_nodal_dofs_);
  const int64_t nodal_dof_id      = static_cast<int64_t>(nodal_dof_id_);

  int64_t dof_address = 0;
  if (storage_type == UnknownStorageType::NODAL)
    dof_address = node_address * num_nodal_dofs + nodal_dof_id;
  else if (storage_type == UnknownStorageType::BLOCK)
    dof_address = nodal_dof_id * num_nodal_dofs + node_address;

  return dof_address;
}

//###################################################################
int64_t SpatialDiscretization::
  MapDOFGlobal(const chi_mesh::Cell &cell, size_t node_index,
              const chi_math::UnknownManager &unknown_manager,
              const unsigned int unknown_id,
              const unsigned int component) const
{
  const std::string fname = __FUNCTION__;

  const auto& cell_mapping = GetCellMapping(cell);
  const uint64_t register_index = cell_mapping.MapNodeRegister(node_index);

  int64_t temp_node_address = 0;
  if (cell.partition_id == chi_mpi.location_id)
    temp_node_address = m_LNR_global_ids[register_index];
  else
    temp_node_address = m_GNR_global_ids[register_index];

  const int64_t node_address = temp_node_address;

  const auto& uk_man   = unknown_manager;
  const auto storage_type = unknown_manager.dof_storage_type;
  const size_t num_nodal_dofs_ = uk_man.GetTotalUnknownStructureSize();
  const size_t nodal_dof_id_   = uk_man.MapUnknown(unknown_id, component);

  const int64_t num_nodal_dofs = static_cast<int64_t>(num_nodal_dofs_);
  const int64_t nodal_dof_id   = static_cast<int64_t>(nodal_dof_id_);

  int64_t dof_address = 0;
  if (storage_type == UnknownStorageType::NODAL)
    dof_address = node_address * num_nodal_dofs + nodal_dof_id;
  else if (storage_type == UnknownStorageType::BLOCK)
    dof_address = nodal_dof_id * num_nodal_dofs + node_address;

  return dof_address;
}

//###################################################################
std::vector<int64_t> SpatialDiscretization::
  MapDOFsLocal(const chi_mesh::Cell &cell,
               const chi_math::UnknownManager &unknown_manager,
               unsigned int unknown_id,
               unsigned int component) const
{
  const std::string fname = __FUNCTION__;

  if (cell.partition_id != chi_mpi.location_id/*home*/)
    throw std::invalid_argument(fname + ": Function cannot be used on "
                                        "non-local cells.");

  const auto& cell_mapping = GetCellMapping(cell);
  const size_t cell_num_nodes = cell_mapping.NumNodes();

  std::vector<int64_t> temp_node_addresses(cell_num_nodes, 0);
  for (size_t node_index=0; node_index < cell_num_nodes; ++node_index)
  {
    const uint64_t register_index = cell_mapping.MapNodeRegister(node_index);
    temp_node_addresses[node_index] = m_LNR_local_ids[register_index];
  }

  const auto& node_addresses = temp_node_addresses;

  const auto& uk_man   = unknown_manager;
  const auto storage_type = unknown_manager.dof_storage_type;
  const size_t num_nodal_dofs_ = uk_man.GetTotalUnknownStructureSize();
  const size_t nodal_dof_id_ = uk_man.MapUnknown(unknown_id, component);

  const int64_t num_nodal_dofs = static_cast<int64_t>(num_nodal_dofs_);
  const int64_t nodal_dof_id      = static_cast<int64_t>(nodal_dof_id_);

  std::vector<int64_t> dof_addresses(cell_num_nodes, 0);

  if (storage_type == UnknownStorageType::NODAL)
    for (size_t node_index=0; node_index < cell_num_nodes; ++node_index)
      dof_addresses[node_index] = node_addresses[node_index] * num_nodal_dofs
                                  + nodal_dof_id;
  else if (storage_type == UnknownStorageType::BLOCK)
    for (size_t node_index=0; node_index < cell_num_nodes; ++node_index)
      dof_addresses[node_index] = nodal_dof_id * num_nodal_dofs
                                  + node_addresses[node_index];

  return dof_addresses;
}

//###################################################################
std::vector<int64_t> SpatialDiscretization::
  MapDOFsGlobal(const chi_mesh::Cell &cell,
               const chi_math::UnknownManager &unknown_manager,
               unsigned int unknown_id,
               unsigned int component) const
{
  const std::string fname = __FUNCTION__;

  const auto& cell_mapping = GetCellMapping(cell);
  const size_t cell_num_nodes = cell_mapping.NumNodes();

  std::vector<int64_t> temp_node_addresses(cell_num_nodes, 0);
  for (size_t node_index=0; node_index < cell_num_nodes; ++node_index)
  {
    const uint64_t register_index = cell_mapping.MapNodeRegister(node_index);
    if (cell.partition_id == chi_mpi.location_id)
      temp_node_addresses[node_index] = m_LNR_global_ids[register_index];
    else
      temp_node_addresses[node_index] = m_GNR_global_ids.at(register_index);
  }

  const auto& node_addresses = temp_node_addresses;

  const auto& uk_man   = unknown_manager;
  const auto storage_type = unknown_manager.dof_storage_type;
  const size_t num_nodal_dofs_ = uk_man.GetTotalUnknownStructureSize();
  const size_t nodal_dof_id_   = uk_man.MapUnknown(unknown_id, component);

  const int64_t num_nodal_dofs = static_cast<int64_t>(num_nodal_dofs_);
  const int64_t nodal_dof_id   = static_cast<int64_t>(nodal_dof_id_);

  std::vector<int64_t> dof_addresses(cell_num_nodes, 0);

  if (storage_type == UnknownStorageType::NODAL)
    for (size_t node_index=0; node_index < cell_num_nodes; ++node_index)
      dof_addresses[node_index] = node_addresses[node_index] * num_nodal_dofs
                                  + nodal_dof_id;
  else if (storage_type == UnknownStorageType::BLOCK)
    for (size_t node_index=0; node_index < cell_num_nodes; ++node_index)
      dof_addresses[node_index] = nodal_dof_id * num_nodal_dofs
                                  + node_addresses[node_index];

  return dof_addresses;
}

//###################################################################
std::vector<int64_t> SpatialDiscretization::
  GetGhostDOFsGlobalIDs(const chi_math::UnknownManager &unknown_manager) const
{
  const auto storage_type = unknown_manager.dof_storage_type;
  const size_t num_nodal_dofs_ = unknown_manager.GetTotalUnknownStructureSize();
  const size_t num_ghost_nodes = m_ghost_global_ids.size();
  const size_t num_ghost_dofs = num_ghost_nodes * num_nodal_dofs_;

  const int64_t num_nodal_dofs = static_cast<int64_t>(num_nodal_dofs_);

  std::vector<int64_t> ghost_dofs_global_ids(num_ghost_dofs, -1);

  if (storage_type == UnknownStorageType::NODAL)
  {
    size_t dof_id = 0;
    for (size_t g=0; g<num_ghost_nodes; ++g)
    {
      const int64_t node_global_address = m_ghost_global_ids[g];
      for (int64_t nodal_dof_id=0; nodal_dof_id<num_nodal_dofs; ++nodal_dof_id)
      {
        ghost_dofs_global_ids[dof_id] =
          node_global_address * num_nodal_dofs + nodal_dof_id;
        ++dof_id;
      }
    }
  }
  else if (storage_type == UnknownStorageType::BLOCK)
  {
    size_t dof_id = 0;
    for (size_t g=0; g<num_ghost_nodes; ++g)
    {
      const int64_t node_global_address = m_ghost_global_ids[g];
      for (int64_t nodal_dof_id=0; nodal_dof_id<num_nodal_dofs; ++nodal_dof_id)
      {
        ghost_dofs_global_ids[dof_id] =
          nodal_dof_id * num_nodal_dofs + node_global_address;
        ++dof_id;
      }
    }
  }

  return ghost_dofs_global_ids;
}