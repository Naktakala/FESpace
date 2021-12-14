#include "fe_lua_utils.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "FiniteElementSpace/Methods/PiecewiseLinear/pwl_mapping.h"
#include "FiniteElementSpace/Methods/Lagrange/lagrange_mapping.h"
#include "FiniteElementSpace/Methods/FiniteVolume/fv_mapping.h"

#include "DiffusionTest/diffusion_solver.h"

#include "chi_log.h"
#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

int chi_math::lua_utils::chiFESpaceTest(lua_State *L)
{
  auto& chi_log = ChiLog::GetInstance();
  chi_log.Log() << "\nExecuting chiFESpaceTest\n\n";

  auto handler = chi_mesh::GetCurrentHandler();
  auto grid = handler->GetGrid();

  using namespace chi_math::finite_element;

  auto pwl_sdm  = FiniteElementSpace<PiecewiseLinear>::New(grid);
  auto pwld_sdm = FiniteElementSpace<PiecewiseLinear>::New(grid);
  auto q2_sdm   = FiniteElementSpace<LagrangeQ2>::New(grid);
  auto fv_sdm   = FiniteElementSpace<FiniteVolume>::New(grid);

  //=================================== Unit testing stuff
  pwl_sdm->SetTextName("PWL-Spatial Discretization");
  chi_log.Log() << pwl_sdm->TextName();
  {
    const auto& pwl_grid = pwl_sdm->Grid();
    double total_volume = 0.0;
    double total_area = 0.0;
    size_t total_face_nodes = 0;
    size_t total_num_nodes = 0;
    for (const auto& cell : pwl_grid.local_cells)
    {
      total_volume += pwl_sdm->GetCellVolume(cell);
      for (size_t f=0; f<cell.faces.size(); ++f)
      {
        total_area += pwl_sdm->GetCellFaceArea(cell, f);
        total_face_nodes += pwl_sdm->GetFaceNumNodes(cell, f);
      }
      total_num_nodes += pwl_sdm->GetCellNumNodes(cell);

      const auto node_locations = pwl_sdm->GetCellNodeLocations(cell);
    }
  }

//  pwl_sdm.ComputeUnitIntegrals();
//  pwl_sdm.InitializeQuadraturePoints();

  DiffusionSolver diffusion_solver;
  diffusion_solver.basic_options["SDM"].SetStringValue("PWL");
  diffusion_solver.SetSpatialDiscretization(q2_sdm);

  diffusion_solver.Initialize();

  chi_log.Log() << "\nDone executing chiFESpaceTest\n\n";

  auto WriteCellsInfo = [&chi_log, &grid, &q2_sdm]()
  {
    for (const auto& cell : grid->local_cells)
    {
      std::stringstream outstr;
      outstr << cell.global_id << " " << cell.centroid.PrintStr() << "\n";

      const auto node_locations = q2_sdm->GetCellNodeLocations(cell);
      size_t n=0;
      for (const auto& node : node_locations)
      {
        outstr
          << "cell " << cell.global_id
          << " node " << n << " " << node.PrintStr()
          << " " << q2_sdm->MapNodeLocal(cell, n)
          << " " << q2_sdm->MapNodeGlobal(cell, n) << "\n";
        ++n;
      }

      chi_log.Log(LOG_ALL) << outstr.str();
    }//for cell

    const auto ghost_cell_ids = grid->cells.GetGhostGlobalIDs();

    for (const uint64_t cell_global_id : ghost_cell_ids)
    {
      const auto& cell = grid->cells[cell_global_id];

      std::stringstream outstr;
      outstr << cell.global_id << " " << cell.centroid.PrintStr() << "\n";

      const auto node_locations = q2_sdm->GetCellNodeLocations(cell);
      size_t n=0;
      for (const auto& node : node_locations)
      {
        outstr
          << "cell " << cell.global_id
          << " node " << n << " " << node.PrintStr()
          << " " << q2_sdm->MapNodeLocal(cell, n)
          << " " << q2_sdm->MapNodeGlobal(cell, n) << "\n";
        ++n;
      }

      chi_log.Log(LOG_ALL) << outstr.str();
    }//for global-id
  };

//  MPI_Barrier(MPI_COMM_WORLD);
//  if (chi_mpi.location_id == 0)
//  {
//    WriteCellsInfo();
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  if (chi_mpi.location_id == 1)
//  {
//    WriteCellsInfo();
//  }
//  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}