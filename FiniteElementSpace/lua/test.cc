#include "fe_lua_utils.h"

#include "fe_space.h"
#include "Methods/FiniteVolume/fv_mapping.h"
#include "Methods/PiecewiseLinear/pwl_mapping.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "ChiLog/chi_log.h"

int chi_math::lua_utils::chiFESpaceTest(lua_State *L)
{
  auto& chi_log = ChiLog::GetInstance();
  chi_log.Log() << "Executing chiFESpaceTest\n\n";

  auto handler = chi_mesh::GetCurrentHandler();
  auto grid = handler->GetGrid();

  using namespace chi_math::finite_element;

  //============================================= FV
  {
    FiniteElementSpace<FiniteVolume> fe(grid);

    size_t c=0;
    for (const auto& cell : fe.Grid().local_cells)
    {
      size_t num_nodes = fe.GetCellNumNodes(cell);
      auto node_positions = fe.GetCellNodeLocations(cell);
      for (size_t i=0; i<num_nodes; ++i)
        chi_log.Log()
        << "cell " << cell.local_id
        << " node " << i
        << " position " << node_positions[i].PrintS();
      ++c;
      if (c>3) break;
    }

    fe.SetTextName("FV Space");

    chi_log.Log() << "Created " << fe.TextName();
  }

  //============================================= PWL
  {
    FiniteElementSpace<PiecewiseLinear> fe(grid);

    size_t c=0;
    for (const auto& cell : fe.Grid().local_cells)
    {
      size_t num_nodes = fe.GetCellNumNodes(cell);
      auto node_positions = fe.GetCellNodeLocations(cell);
      for (size_t i=0; i<num_nodes; ++i)
        chi_log.Log()
        << cell.local_id << " "
        << i << " "
        << node_positions[i].PrintS();
      ++c;
      if (c>3) break;
    }

    fe.SetTextName("FE Space");

    chi_log.Log() << "Created " << fe.TextName();
  }


  return 0;
}