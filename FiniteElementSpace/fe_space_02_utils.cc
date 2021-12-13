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

/**Gets the volume of the requested cell.*/
double SpatialDiscretization::GetCellVolume(const chi_mesh::Cell& cell)
{
  return m_local_cell_mappings[cell.local_id]->Volume();
}


double SpatialDiscretization::GetCellFaceArea(const chi_mesh::Cell& cell,
                                              size_t face_index)
{
  return m_local_cell_mappings[cell.local_id]->FaceArea(face_index);
}

/**Gets the number of nodes (not DOFs) for the requested cell.*/
size_t SpatialDiscretization::GetCellNumNodes(const chi_mesh::Cell& cell)
{
  size_t num_nodes = 0;

  return m_local_cell_mappings[cell.local_id]->NumNodes();
}

/**Gets the number of nodes (not DOFs) for the requested cell-face pair.*/
size_t SpatialDiscretization::GetFaceNumNodes(const chi_mesh::Cell& cell, size_t face)
{
  return m_local_cell_mappings[cell.local_id]->FaceNumNodes(cell, face);
}

/**Returns the locations of all the nodes of the cell consistent with the
 * supplied mapping.*/
std::vector<chi_mesh::Vector3> SpatialDiscretization::
  GetCellNodeLocations(const chi_mesh::Cell& cell)
{
  return m_local_cell_mappings[cell.local_id]->CellNodeLocations(cell);
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