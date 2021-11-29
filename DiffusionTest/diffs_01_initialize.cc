#include "diffusion_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "FiniteElementSpace/fe_space.h"
#include "FiniteElementSpace/Methods/PiecewiseLinear/pwl_mapping.h"



//###################################################################
/**Initialization routine.*/
void DiffusionSolver::Initialize()
{
  //============================================= Get grid
  auto handler = chi_mesh::GetCurrentHandler();
  m_grid_ptr = handler->GetGrid();

  //============================================= Initialize spatial
  //                                              discretization
  {
    using namespace chi_math::finite_element;
    if (basic_options("SDM").StringValue() == "PWLC")
      m_sdm = FiniteElementSpace<PiecewiseLinear>::New(m_grid_ptr);
    else
      throw std::invalid_argument("DiffusionSolver::Initialize. Unsupported"
                                  " option for SDM: " +
                                  basic_options("SDM").StringValue());
  }

  //============================================= Initialize unknown manager
  m_dof_handler.AddUnknown(chi_math::UnknownType::SCALAR);
  m_dof_handler.SetUnknownTextName(0, "Temperature");

  //============================================= Get number of DOFs
  size_t num_local_dofs = 0;
  size_t num_globl_dofs = 0;
}