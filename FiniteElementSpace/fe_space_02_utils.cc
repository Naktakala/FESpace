#include "fe_space.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

using namespace chi_math::finite_element;

const FiniteElementMapping& SpatialDiscretization::
  GetCellMapping(const chi_mesh::Cell& cell) const
{
  if (cell.partition_id == chi_mpi.location_id/*home*/)
    return *m_local_cell_mappings[cell.local_id];
  else
    return *m_ghost_cell_mappings.at(cell.global_id);
}

const std::vector<int64_t>& SpatialDiscretization::
  GetCellRelevantLocalIDRegister(const chi_mesh::Cell &cell) const
{
  if (cell.partition_id == chi_mpi.location_id)
    return m_LNR_local_ids;
  else
    return m_GNR_local_ids;
}

const std::vector<int64_t>& SpatialDiscretization::
  GetCellRelevantGlobalIDRegister(const chi_mesh::Cell &cell) const
{
  if (cell.partition_id == chi_mpi.location_id)
    return m_LNR_global_ids;
  else
    return m_GNR_global_ids;
}

/**Gets the volume of the requested cell.*/
double SpatialDiscretization::GetCellVolume(const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = GetCellMapping(cell);
  return cell_mapping.Volume();
}


double SpatialDiscretization::GetCellFaceArea(const chi_mesh::Cell& cell,
                                              size_t face_index)
{
  const auto& cell_mapping = GetCellMapping(cell);
  return cell_mapping.FaceArea(face_index);
}

/**Gets the number of nodes (not DOFs) for the requested cell.*/
size_t SpatialDiscretization::GetCellNumNodes(const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = GetCellMapping(cell);
  return cell_mapping.NumNodes();
}

/**Gets the number of nodes (not DOFs) for the requested cell-face pair.*/
size_t SpatialDiscretization::GetFaceNumNodes(const chi_mesh::Cell& cell, size_t face)
{
  const auto& cell_mapping = GetCellMapping(cell);
  return cell_mapping.FaceNumNodes(face);
}

/**Returns the locations of all the nodes of the cell consistent with the
 * supplied mapping.*/
std::vector<chi_mesh::Vector3> SpatialDiscretization::
  GetCellNodeLocations(const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = GetCellMapping(cell);

  const size_t num_nodes = cell_mapping.NumNodes();

  std::vector<chi_mesh::Vector3> node_locations(num_nodes);

  if (cell.partition_id == chi_mpi.location_id)
    for (size_t n=0; n<num_nodes; ++n)
    {
      const size_t register_id = cell_mapping.MapNodeRegister(n);
      const size_t register_id_map = m_LNR_id_2_LNLL_id_map[register_id];

      node_locations[n] = m_LNLL[register_id_map];
    }
  else
    for (size_t n=0; n<num_nodes; ++n)
    {
      const size_t register_id = cell_mapping.MapNodeRegister(n);
      const size_t register_id_map = m_GNR_id_2_GNLL_id_map[register_id];

      node_locations[n] = m_GNLL[register_id_map];
    }

  return node_locations;
}

//###################################################################
std::vector<size_t> SpatialDiscretization::
  CellNodeIndices(const chi_mesh::Cell &cell) const
{
  const auto& cell_mapping = GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  std::vector<size_t> node_indices(num_nodes, 0);
  for (size_t i=0; i<num_nodes; ++i)
    node_indices[i] = i;

  return node_indices;
}