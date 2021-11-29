#ifndef DIFFUSION_SOLVER_H
#define DIFFUSION_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/UnknownManager/unknown_manager.h"

#include "FiniteElementSpace/fe_space.h"

//################################################################### Class def
/**Test diffusion solver for different methods.*/
class DiffusionSolver : public chi_physics::Solver
{
protected:
  chi_math::finite_element::SpatialDiscretizationPtr m_sdm;
  chi_mesh::MeshContinuumPtr m_grid_ptr;
  chi_math::UnknownManager   m_dof_handler;

public:
  DiffusionSolver();

  void Initialize() override;

  virtual ~DiffusionSolver() = default;
};

#endif //DIFFUSION_SOLVER_H